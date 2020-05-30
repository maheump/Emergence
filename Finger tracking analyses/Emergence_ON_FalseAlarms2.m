% This script implements analyses of the false alarms that are observed in
% fully-stochastic parts. We show that the strength of (false) beliefs is
% predicted by the ideal observer, meaning that subjects' false alarms
% reflect genuine local (probabilistic) regularities occuring by chance in
% actually entirely stochastic parts of the sequences. We show that is true
% after controling for the effect of time elapsed since the beginning od
% the sequence on the increase of false alarms.
% 
% Copyright (c) 2020 Maxime Maheu

%% COMPARE FALSE ALARMS IN SUBJECTS AND IN THE IO
%  ==============================================

% Define options
% ~~~~~~~~~~~~~~

% The likelihood of the fully-stochastic hypothesis also quantifies the
% extent of false alarms
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
restopt = 1;

% Average subjects' trajectories in fully-stochastic parts according to the
% ideal observer's beliefs in the probabilistic hypothesis
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Select observations
randidx = Emergence_SelectFullyStochSeq(G, filter, restopt);

% Prepare output variables
binsubn    = NaN(nBin,  nSub);
binsubtraj = NaN(nBin,3,nSub);
biniotraj  = NaN(nBin,3,nSub);

% For each subject
for iSub = 1:nSub
    
    % Get beliefs of both subject and ideal observer in the
    % fully-stochastic parts of the sequence
    subbel = cell2mat(cellfun(@(x,y) x.BarycCoord(y,:), ...
        G(:,iSub), randidx(:,iSub), 'uni', 0));
    iobel = cell2mat(cellfun(@(x,y) x.BarycCoord(y,:), ...
        IO(:,iSub), randidx(:,iSub), 'uni', 0));
    
    % Get probability bins
    if strcmpi(binmeth, 'unif') % bins of the same amplitude
        edges = linspace(0, 1, nBin+1);
    elseif strcmpi(binmeth, 'equil') % bins with the same number of observations
        edges = prctile(iobel(:,iHyp), linspace(0, 100, nBin+1));
    else, error('Please check the binnig method that is provided');
    end
    
    % Get in which bin falls each observation
    [binsubn(:,iSub),~,bins] = histcounts(iobel(:,iHyp), edges);
    
    % Average likelihoods in each bin
    binsubtraj(:,:,iSub) = cell2mat(arrayfun(@(i) ...
        mean(subbel(bins == i,:)), 1:nBin, 'uni', 0)');
    biniotraj(:,:,iSub)  = cell2mat(arrayfun(@(i) ...
        mean(iobel( bins == i,:)), 1:nBin, 'uni', 0)');
end

% Average over subjects
avgsubtraj = mean(binsubtraj, 3, 'OmitNaN');
avgiotraj  = mean(biniotraj,  3, 'OmitNaN');
semsubtraj = sem( binsubtraj, 3);
semiotraj  = sem( biniotraj,  3);

% Display the position of different false alarms in the triangle depending
% on the corresponding ideal observer's beliefs in that false alarm
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Prepare a new window
figure('Position', [323 805 200 300]);

% Create a colormap whose length equals the number of bins
cmap = {'Blues', 'Reds', 'Greens'};
prec = 100;
gege = Emergence_Colormap(cmap{iHyp}, prec);
grid = linspace(0, 1, prec);
xm = mean(biniotraj(:,iHyp,:), 3);
[~,colidx] = min(abs(xm-grid), [], 2);
cmap = gege(colidx,:);

% Customize the colormap and add a colorbar
cbr = colorbar('Location', 'SouthOutside');
cbr.Label.String = sprintf('p(M_%s|y) from the ideal observer', proclab{iHyp}(1));
colormap(cmap); caxis([0, 1]);

% Display the triangle
Emergence_PlotTrajOnTri; alpha(0);
Emergence_PlotGridOnTri(10, iHyp, tricol(iHyp,:));

% Convert beliefs in each hypothesis to cartesian coordinates
cartavgcoord = flipud(avgsubtraj*tricc);
cartsemcoord = flipud(semsubtraj*tricc);

% Display error bars
plot(repmat(cartavgcoord(:,1), 1, 2)', cartavgcoord(:,2)' + ...
    [-1;1] .* cartsemcoord(:,2)', 'k-', 'LineWidth', 1/2);
plot(cartavgcoord(:,1)' + ([-1;1] .* cartsemcoord(:,1)'), ...
    repmat(cartavgcoord(:,2), 1, 2)', 'k-', 'LineWidth', 1/2);

% Get average size of the bin
avgbinsize = mean(binsubn, 2);
dotsize = 1 + flipud(log(avgbinsize+1).*10);

% Display binned averages
scatter(cartavgcoord(:,1), cartavgcoord(:,2), dotsize, ...
    flipud(cmap), 'filled', 'MarkerEdgeColor', 'k');

% Zoom in the bottom part of the triangle
axis('tight');
ylim([0,(sqrt(3)/2)/2]);

% Save the figure
save2pdf(fullfile(ftapath, 'figs', 'F_FA_Tri.pdf'));

%% CONTROL FOR TIME ELAPSED SINCE THE BEGINNING OF THE SEQUENCE
%  ============================================================

% Prepare output variables
regcoef      = NaN(2,3,nSub);
resbinsubFAL = NaN(nBin,2,nSub);
resbinioFAL  = NaN(nBin,2,nSub);

% For each subject
for iSub = 1:nSub
    
    % Get beliefs of both subject and ideal observer in the
    % fully-stochastic parts of the sequence
    subbel = cell2mat(cellfun(@(x,y) x.BarycCoord(y,:), ...
        G(:,iSub), randidx(:,iSub), 'uni', 0));
    iobel = cell2mat(cellfun(@(x,y) x.BarycCoord(y,:), ...
        IO(:,iSub), randidx(:,iSub), 'uni', 0));
    
    % Get probability bins
    if strcmpi(binmeth, 'unif') % bins of the same amplitude
        edges = linspace(0, 1, nBin+1);
    elseif strcmpi(binmeth, 'equil') % bins with the same number of observations
        edges = prctile(iobel(:,iHyp), linspace(0, 100, nBin+1));
    else, error('Please check the binnig method that is provided');
    end
    
    % Get in which bin falls each observation
    [binsubn(:,iSub),~,bins] = histcounts(iobel(:,iHyp), edges);
    
    % Build the design matrix
    obsidx = cell2mat(cellfun(@find, randidx(:,iSub), 'uni', 0));
    desmat = [ones(numel(obsidx),1), iobel(:,iHyp), obsidx];
    
    % Center predictors
    desmat = desmat - [0, mean(desmat(:,2:end), 1, 'OmitNaN')];
    
    % Regress subjects' and ideal observer's beliefs in the
    % fully-stochastic hypothesis by taking, or not, into account the
    % effect of the number of observations received since the beginning of
    % the sequence (that was shown to induce false alarms)
    for iMod = 1:2
        selecpred = 1:(1+iMod);
        beta = regress(subbel(:,iHyp), desmat(:,selecpred));
        regcoef(iMod,selecpred,iSub) = beta;
        
        % Remove variance from all the predictors except the IO's false
        % alarm likelihood
        pidx = setdiff(selecpred, 2);
        pred = desmat(:,pidx) * beta(pidx); % predictions based on beta values
        ressubbel = subbel(:,iHyp) - pred; % pseudo-residuals
        
        % Average IO's false alarm likelihood and subject's residual false
        % alarms
        resbinsubFAL(:,iMod,iSub) = arrayfun(@(i) mean(ressubbel(bins == i     )), 1:nBin);
        resbinioFAL(:,iMod,iSub)  = arrayfun(@(i) mean(iobel(    bins == i,iHyp)), 1:nBin);
    end
end

% Average over subjects
avgresbinsubFAL = mean(resbinsubFAL, 3, 'OmitNaN');
avgresbinioFAL  = mean(resbinioFAL,  3, 'OmitNaN');
semresbinsubFAL = sem( resbinsubFAL, 3);
semresbinioFAL  = sem( resbinioFAL,  3);

% Display the correlation between subjects' and ideal observer's beliefs
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Prepare a new window
figure('Position', [524 905 200 200]);

% Display predictions from the first simpler model
iMod = 1;
xm = avgresbinioFAL(:,iMod);
ym = avgresbinsubFAL(:,iMod);
B = Emergence_Regress(ym, xm, 'TLS', {'beta0', 'beta1'});
xval = [min(xm), max(xm)];
plot(xval, xval.*B(2) + B(1), 'k-', 'LineWidth', 1); hold('on');

% Residuals from which regression to plot
iMod = 2;
xm = avgresbinioFAL(:,iMod);
ym = avgresbinsubFAL(:,iMod);
xs = semresbinioFAL(:,iMod);
ys = semresbinsubFAL(:,iMod);

% Display the regression line
B = Emergence_Regress(ym, xm, 'TLS', {'beta0', 'beta1'});
confint = Emergence_Regress(ym, xm, 'TLS', 'confint');
xval = Emergence_Regress(ym, xm, 'TLS', 'confintx');
fill([xval, fliplr(xval)], [confint(1,:), fliplr(confint(2,:))], 'k', ...
       'EdgeColor', 'none', 'FaceColor', g, 'FaceAlpha', 1/2); hold('on');
plot(xval, xval.*B(2) + B(1), 'k-', 'LineWidth', 3);

% Display averaged beliefs in each probability bin with its error bars
plot(xm, ym, 'k-', 'LineWidth', 1/2); 
plot(xm' + xs'.*[-1;1], repmat(ym',2,1), 'k-');
plot(repmat(xm',2,1), ym' + ys'.*[-1;1], 'k-');
scatter(xm, ym, dotsize, cmap, 'filled', 'MarkerEdgeColor', 'k')

% Customize the axes
axis([0.4,1,-0.15,0.1]);
axis('square'); set(gca, 'Box', 'Off');

% Add some text labels
xlabel('Beliefs from the IO');
ylabel('Beliefs from subjects (residuals)');

% Save the figure
save2pdf(fullfile(ftapath, 'figs', 'F_FA_Corr.pdf'));

% Display the distribution of correlation coefficients over subjects
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Which regression coefficients to analyse
coefofint = squeeze(regcoef(:,2,:))';

% Prepare a new window
figure('Position', [725 905 140 200]);

% Average regression coefficients over subjects
m = mean(coefofint);
s = sem(coefofint);

% Display averaged regression coefficients
bar(1:2, m, 'FaceColor', g, 'EdgeColor', 'k'); hold('on');
plot(repmat(1:2, 2, 1), m+[-s;s], 'k-');

% Customize the axes
set(gca, 'Box', 'Off', 'XTick', 1:2, 'XTickLabel', {'M_{1}', 'M_{2}'});
axis([0,3,0,0.7]);

% Display whether the difference is significant or not
Emergence_DispStatTest(coefofint);

% Add some text label
xlabel('Regression models');
ylabel('Regression coefficient for p(H_{S}|y)');

% Save the figure
save2pdf(fullfile(ftapath, 'figs', 'F_FA_Gp.pdf'));

% Perform a t-test (against chance) on correlation coefficients
[~,pval,tci,stats] = ttest(coefofint);
Emergence_PrintTstats(pval,tci,stats);
