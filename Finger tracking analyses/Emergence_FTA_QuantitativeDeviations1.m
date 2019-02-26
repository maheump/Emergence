% This script implements analyses of finger tracking that aim at showing
% quantitative deviations that subjects exhibit compared to the ideal
% inference scenario. In particular, we show that subjects are biased in
% the way they report their probabilistic beliefs in a similar manner of 
% what has been shown by economic theory: they overestimate small
% probabilities and underestimate large ones.
% 
% Copyright (c) 2018 Maxime Maheu

% Compute error over parameters grid
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Define parameters grid
GammaGrid = 0.4:0.01:1;
P0Grid = 0.2:0.01:0.6;

% Get size of grids
nGamma = numel(GammaGrid);
nP0 = numel(P0Grid);

% Prepare output variables
MSE = NaN(nGamma, nP0, nSub);

% For each subjects
for iSub = 1:nSub
    fprintf('Fitting parameter for subject #%2.0f/%2.0f.\n', iSub, nSub);
	
    % Get posterior probabilities from the ideal observer
    trueproba = cell2mat(cellfun(@(x) x.BarycCoord, IO(:,iSub), 'UniformOutput', 0));
    trueproba = reshape(trueproba(:), [1,1,numel(trueproba)]);
    
    % Apply the probability weighting
    transfproba = probaweight(trueproba, GammaGrid', P0Grid);
    
    % Get posterior probabilities from the subjects
    estproba = cell2mat(cellfun(@(x) x.BarycCoord, G(:,iSub), 'UniformOutput', 0));
    estproba = reshape(estproba(:), [1,1,numel(estproba)]);
    
    % Compute the mean squared difference between transformed probabilities
    % and reported probabilities
    error = (estproba - transfproba) .^ 2;
    MSE(:,:,iSub) = mean(error,3);
end

% Find the best parameters set
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Get best parameters set (i.e. the set that minimizes mean squared error)
% independently for each subject
[err,idx] = cellfun(@(x) min(x(:)), squeeze(mat2cell(MSE, nGamma, nP0, ones(1,1,nSub))));
[pI,gI] = ind2sub(size(MSE), idx);
PWparams = [GammaGrid(pI); P0Grid(gI)];

% Get the error level of a non-weighted scenario
errornonw = MSE(GammaGrid == 1, 1, :);

% Test whether the error corresponding to the best parameters set is
% significantly smaller than the one corresponding to gamma = 1 (no matter
% of the value of p0) which corresponds to unbiased (i.e. non weighted)
% probability matching
[~,pval,tci,stats] = ttest(err - squeeze(errornonw));
Emergence_PrintTstats(pval,tci,stats);

% Average best parameters over subjects
avgPWparams = mean(PWparams, 2);
semPWparams = sem(PWparams, 2);

% Prepare a new window
figure('Position', [1 905 200 200]);

% Define parameters' name and null values
pwparamname = {'\gamma', 'p_{0}'};
nullparamval = [1, 1/2];

% For each parameter (gamma and p0) of the probability weighting function
for iParam = 1:2
    subplot(1,2,iParam); hold('on');
    
    % Display null value of the parameter aginst which it will be compared
    plot([0,2], repmat(nullparamval(iParam),1,2), 'k--');
    
    % Display distribution of parameter value
    Emergence_PlotSubGp(PWparams(iParam,:), 'k');
    
    % Customize the axes
    xlim([0,2]);
    if     iParam == 1, ylim([0,2]);
    elseif iParam == 2, ylim([0,1]);
    end
    set(gca, 'Box', 'Off', 'XTick', [], 'XColor', 'none');
    
    % Display whether the difference is significant or not
    Emergence_DispStatTest(PWparams(iParam,:));
    
    % Add some text labels
    ylabel(['$', pwparamname{iParam}, '$'], 'Interpreter', 'LaTeX', ...
        'Rotation', 0, 'VerticalAlignment', 'Bottom');
    
    % Display the output of the statistical comparison against null value
    [~,pval,tci,stats] = ttest(PWparams(iParam,:)' - nullparamval(iParam));
    Emergence_PrintTstats(pval,tci,stats);
end

% Save the figure
save2pdf(fullfile(ftapath, 'figs', 'F_QD_GpPWFit.pdf'));

% Display the average error map over parameters grid
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Center the error matrices on that error level
ctrerrmap = MSE - errornonw;

% Average over subjects
avgerrmap = mean(ctrerrmap, 3);

% Prepare a new window
figure('Position', [202 905 180 200]);

% Diplay the average map of error over the parameters grid
imagesc(P0Grid, GammaGrid, avgerrmap); hold('on');
contour(P0Grid, GammaGrid, avgerrmap, 11, 'k-', 'LineWidth', 1/2);

% Find the minipum of that averaged error map
[~,idx] = min(avgerrmap(:));
[pI,gI] = ind2sub(size(avgerrmap), idx);
plot(P0Grid(gI), GammaGrid(pI), 'k.', 'MarkerSize', 10);

% Display the group 
plot(repmat(avgPWparams(2),1,2), avgPWparams(1)+[-1,1].*semPWparams(1), 'k-');
plot(avgPWparams(2)+[-1,1].*semPWparams(2), repmat(avgPWparams(1),1,2), 'k-');
plot(avgPWparams(2), avgPWparams(1), 'k.', 'MarkerSize', 20);

% Customize the axes
axis('xy'); caxis([-1,1] .* max(abs(caxis)));

% Add a colorbar
colorbar; colormap(flipud(cbrewer2('Spectral', 1e5)));

% Add some text labels
xlabel('p_{0} parameter');
ylabel('\gamma parameter');

% Save the figure
save2pdf(fullfile(ftapath, 'figs', 'F_QD_AvgErrMap.pdf'));

% Display binned probabilities from subjects versus the ideal observer
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Define the number of bins to use
nBin = 12;

% Create bins of probability values
binmeth = 'unif'; % type of binning method
if ~isempty(nBin) || ~isnan(nBin)
    if strcmpi(binmeth, 'unif') % bins uniformed on the probability space
        bingrid = linspace(0, 1, nBin+1);
    elseif strcmpi(binmeth, 'equil') % bins uniformed on the probability space
        iobel = cellfun(@(x) x.BarycCoord(:), IO, 'UniformOutput', 0);
        iobel = cell2mat(iobel(:));
        bingrid = prctile(iobel, linspace(0, 100, nBin+1));
    else, error('Please check the binnig method that is provided');
    end
end

% Prepare the output
binbel = NaN(nBin,2,nSub);

% For each subject
for iSub = 1:nSub
    
    % Get beliefs of the subject and those of the ideal observer in all the
    % sequences
    subbel = cellfun(@(x) x.BarycCoord(:), G(:,iSub), 'UniformOutput', 0);
    subbel = cell2mat(subbel(:));
    iobel = cellfun(@(x) x.BarycCoord(:), IO(:,iSub), 'UniformOutput', 0);
    iobel = cell2mat(iobel(:));
    
    % Average probabilities within each bin
    [~,~,idx] = histcounts(iobel, bingrid);
    for iBin = 1:nBin
        binbel(iBin,1,iSub) = mean(iobel( idx == iBin));
        binbel(iBin,2,iSub) = mean(subbel(idx == iBin));
    end
end

% Prepare window
figure('Position', [383 905 220 200]);

% Display identity line and help lines
plot([0,1], [0,1], '-', 'Color', g); hold('on');
text(0.15, 0.15, 'Identity', 'Color', g, 'VerticalAlignment', 'Top', 'Rotation', 45);

% Average over subjects
avgbinbel = mean(binbel, 3);
sembinbel = sem(binbel, 3);

% Transform posterior probabilities from the ideal observer
trueproba = avgbinbel(:,1);
trueproba = reshape(trueproba, [1 1 numel(trueproba)]);
transfproba = probaweight(trueproba, GammaGrid', P0Grid);

% Get probabilities estimated by subjects
estproba = avgbinbel(:,2);
estproba = reshape(estproba, [1 1 numel(estproba)]);

% Compute error between transformed 
error = (estproba - transfproba) .^ 2;
avgerror = mean(error, 3);

% Find the parameters set that induce the smaller error
[~,idx] = min(avgerror(:));
[pI,gI] = ind2sub(size(avgerror), idx);
P = [GammaGrid(pI); P0Grid(gI)];

% Display probability weighting function obtained with group-average
% parameters
pgrid = linspace(0, 1, 1001);
fp = probaweight(pgrid, P(1), P(2));
plot(pgrid, fp, 'k-', 'LineWidth', 2);

% Display averaged bins and related error
plot(repmat(avgbinbel(:,1), 1, 2)', avgbinbel(:,2)'+sembinbel(:,2)'.*[-1,1]', 'k-')
plot(avgbinbel(:,1), avgbinbel(:,2), 'ko', 'MarkerFaceColor', g, 'MarkerSize', 8);

% Customize the axes
axis(repmat([0,1], 1, 2)); axis('square');
set(gca, 'XTick', 0:0.2:1, 'YTick', 0:0.2:1);
set(gca, 'Box', 'Off');

% Add some text labels
xlabel('Beliefs from the ideal observer'); ylabel('Beliefs from subjects');

% Save the figure
save2pdf(fullfile(ftapath, 'figs', 'F_QD_AvgPWFit.pdf'));

% Probability weighting function
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% See: https://en.wikipedia.org/wiki/Prospect_theory
% See: Gonzalez, R., & Wu, G. (1999). On the shape of the probability
%   weighting function. Cognitive psychology, 38(1), 129-166.
function y = probaweight(x, gamma, p0)
y = logitinv(gamma .* logit(x) + (1 - gamma) .* logit(p0));
end

% Inverse of the logit transformation
function u = logitinv(v)
maxcut = -log(eps);
mincut = -log(1/realmin - 1);
u = 1 ./ (1 + exp(-max(min(v,maxcut),mincut)));
end

% Logit function
function a = logit(b)
a = log(b./(1-b));
a(real(a)~=a) = NaN;
end
