% This script implements analyses of the false alarms that are observed in
% fully-stochastic parts. We show that the strength of (false) beliefs is
% predicted by the ideal observer, meaning that subjects' false alarms
% reflect genuine local (probabilistic) regularities occuring by chance in
% actually entirely stochastic parts of the sequences.
% 
% Copyright (c) 2018 Maxime Maheu

%% COMPARE FALSE ALARMS TOWARD THE PROBABILISTIC HYPOTHESIS IN SUBJECTS AND IN THE IO
%  ==================================================================================

% Define options
%�~~~~~~~~~~~~~~

% Focus on false alarms toward 
iHyp = 3;

% Define the number of bins to use
nBin = 10;

% Define the type of binning method
binmeth = 'equil';

% Select the sequences to look at
% 1: all fully-stochastic parts
% 2: only fully-stochastic sequences
% 3: only fully-stochastic sequences that were correctly labeled
% 4: only first-part of stochastic-to-regular sequences
restopt = 3;

% Average subjects' trajectories in fully-stochastic parts according to the
% ideal observer's beliefs in the probabilistic hypothesis
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Select moment in which sequences do not entail any regularities
randidx = cellfun(@(x) find(x.Gen == 1), G, 'UniformOutput', 0);

% Restric to fully-stochastic sequences
if restopt == 2 || restopt == 3
    randidx(cell2mat(cidx(1:2)),:) = {[]};
    
    % Restric to those that were correctly labeled
    if restopt == 3
        detecmask = cellfun(@(x) isnan(x.Questions(2)), G(cidx{3},:));
        randidx(~detecmask) = {[]};
    end
    
% Restric to first part of stochastic-to-regular sequences
elseif restopt == 4
    randidx(cidx{3},:) = {[]};
end

% Get probability bins
if strcmpi(binmeth, 'unif') % bins of the same amplitude
    pgrid = linspace(0, 1, nBin+1);
elseif strcmpi(binmeth, 'equil') % bins with the same number of observations
    iobel = cellfun(@(x,y) x.BarycCoord(y,iHyp), IO, randidx, 'UniformOutput', 0);
    iobel = cell2mat(iobel(:));
    pgrid = prctile(iobel, linspace(0, 100, nBin+1));
else, error('Please check the binnig method that is provided');
end

% Prepare output variables
binsubtraj = NaN(nBin,3,nSub);
biniotraj  = NaN(nBin,3,nSub);
coef       = NaN(nSub,1);
binsubn    = NaN(nBin,nSub);

% For each subject
for iSub = 1:nSub
    
    % Get beliefs of both subject and ideal observer in the
    % fully-stochastic parts of the sequence
    subbel = cell2mat(cellfun(@(x,y) x.BarycCoord(y,:), ...
        G(:,iSub), randidx(:,iSub), 'UniformOutput', 0));
    iobel = cell2mat(cellfun(@(x,y) x.BarycCoord(y,:), ...
        IO(:,iSub), randidx(:,iSub), 'UniformOutput', 0));
    
    % Correlate subject's and ideal observer's beliefs in the probabilistic
    % hypothesis
    coef(iSub) = Emergence_Regress(subbel(:,iHyp), iobel(:,iHyp), 'CC', 'r');
    
    % Get indices of moments at which the ideal observer's beliefs in the
    % probabilistic hypothesis fell into some binned probability values
    condidx = cellfun(@(i,j) iobel(:,iHyp) >= i & iobel(:,iHyp) < j, ...
        num2cell(pgrid(1:end-1)), num2cell(pgrid(2:end)), 'UniformOutput', 0)';
    
    % Average both subject's and ideal observer's beliefs over those
    % moments for each bin
    binsubn(:,iSub) = cellfun(@(i) size(iobel(i,:), 1), condidx, 'UniformOutput', 1);
    binsubtraj(:,:,iSub) = cell2mat(cellfun(@(i) mean(subbel(i,:), 1), ...
        condidx, 'UniformOutput', 0));
    biniotraj(:,:,iSub)  = cell2mat(cellfun(@(i) mean(iobel(i,:), 1), ...
        condidx, 'UniformOutput', 0));
end

% Average over subjects
avgsubtraj = mean(binsubtraj, ndims(binsubtraj), 'OmitNaN');
avgiotraj  = mean(biniotraj,  ndims(binsubtraj), 'OmitNaN');
semsubtraj = sem( binsubtraj, ndims(biniotraj));
semiotraj  = sem( biniotraj,  ndims(biniotraj));

% Display the position of different false alarms in the triangle depending
% on the corresponding ideal observer's beliefs in that false alarm
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Prepare a new window
figure('Position', [1 805 200 300]);

% Create a colormap whose length equals the number of bins
cmap = {'Blues', 'Reds', 'Greens'};
cmap = cbrewer2(cmap{iHyp}, nBin);

% Customize the colormap and add a colorbar
cbr = colorbar('Location', 'SouthOutside');
cbr.Label.String = sprintf('p(M_%s|y) from the ideal observer', proclab{iHyp}(1));
colormap(cmap); caxis([0, 1]);

% Display the triangle
Emergence_PlotTrajOnTri; alpha(0);
Emergence_PlotGridOnTri(10, iHyp, tricol(iHyp,:));

% Convert beliefs in each hypothesis to cartesian coordinates
tcn = [0, sqrt(3)/2; 1, sqrt(3)/2; 1/2, 0];
cartavgcoord = flipud(avgsubtraj*tcn);
cartsemcoord = flipud(semsubtraj*tcn);

% Display error bars
plot(repmat(cartavgcoord(:,1), 1, 2)', cartavgcoord(:,2)' + [-1;1] .* cartsemcoord(:,2)', 'k-', 'LineWidth', 1/2);
plot(cartavgcoord(:,1)' + ([-1;1] .* cartsemcoord(:,1)'), repmat(cartavgcoord(:,2), 1, 2)', 'k-', 'LineWidth', 1/2);

% Get average size of the bin
avgbinsize = mean(binsubn, 2);
dotsize = 1 + flipud(log(avgbinsize+1).*10);

% Display binned averages
scatter(cartavgcoord(:,1), cartavgcoord(:,2), dotsize, ...
    flipud(cmap), 'filled', 'MarkerEdgeColor', 'k');

% Save the figure
save2pdf(fullfile(ftapath, 'figs', 'F_FA_Tri.pdf'));

% Display the correlation between subjects' and ideal observer's beliefs
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Prepare a new window
figure('Position', [202 905 200 200]);

% Display the identity line
plot([0,1], [0,1], '-', 'Color', g); hold('on');
text(0.15, 0.15, 'Identity', 'Color', g, 'VerticalAlignment', 'Top', 'Rotation', 45);

% Display averaged beliefs in each probability bin with its error bars
plot(avgiotraj(:,iHyp), avgsubtraj(:,iHyp), 'k-', 'LineWidth', 2);
plot(repmat(avgiotraj(:,iHyp)', [2,1]), avgsubtraj(:,iHyp)'+semsubtraj(:,iHyp)'.*[-1;1], 'k-');
plot(avgiotraj(:,iHyp)'+semiotraj(:,iHyp)'.*[-1;1], repmat(avgsubtraj(:,iHyp)', [2,1]), 'k-');
scatter(avgiotraj(:,iHyp), avgsubtraj(:,iHyp), dotsize, cmap, 'filled', 'MarkerEdgeColor', 'k');

% Customize the axes
axis('square'); set(gca, 'Box', 'Off');
axis([0,1,0,1]);
set(gca, 'XTick', get(gca, 'YTick'));

% Add some text labels
xlabel('Posterior beliefs from the IO');
ylabel('Posterior beliefs from the subjects');

% Save the figure
save2pdf(fullfile(ftapath, 'figs', 'F_FA_Corr.pdf'));

% Display the distribution of correlation coefficients over subjects
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Prepare a new window
figure('Position', [403 905 120 200]);

% Display chance level
plot([0,2], zeros(1,2), '-', 'Color', g); hold('on');

% Display distribution of correlation coefficients
Emergence_PlotSubGp(coef, tricol(iHyp,:));

% Customize the axes
set(gca, 'Box', 'Off', 'XTick', [], 'XColor', 'None');
axis([0,2,-1,1]);

% Display whether the difference is significant or not
Emergence_DispStatTest(coef);

% Add some text label
ylabel('Correlation coefficient');

% Save the figure
save2pdf(fullfile(ftapath, 'figs', 'F_FA_Gp.pdf'));

% Perform a t-test (against chance) on correlation coefficients
[~,pval,tci,stats] = ttest(coef);
disptstats(pval,tci,stats);