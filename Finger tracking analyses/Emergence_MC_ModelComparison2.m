% This script shows the error of regressions between belief difference in
% subjects and different models.
% 
% Copyright (c) 2020 Maxime Maheu

% Define options
% ~~~~~~~~~~~~~~

% Get model simulations, i.e. choose among:
% - PredictabilityLin
% - PredictabilityMax
% - PredictabilityUshape
% - PredictabilityAshape
% - PseudoDeterministic
% - BiasedPseudoDeterministic
% - DifferentPriors
% - IndependentDiffDiscrete
% - IndependentDiffContinuous
% - IndependentRatioDiscrete
% - IndependentRatioContinuous
% - Leak
% - TreeDepth
SimuType = 'IndependentDiffContinuous';

% Define properties of the plots based on the model to simulate
lab = {'Continuous', 'Shape', 'Leak'};
dispcont = any(contains(SimuType, lab, 'IgnoreCase', true));
useXlog = false;
useYlog = false;
lab = {'Ashape', 'Independent'};
if any(contains(SimuType, lab, 'IgnoreCase', true)), useXlog = true; end
lab = {'Predictability', 'PseudoDeterministic', 'DifferentPriors'};
if any(contains(SimuType, lab, 'IgnoreCase', true)), useYlog = true; end

% Define the number of bins to use when plotting the results
nBin = 10;

% Define the type of binning method
binmeth = 'equil';

% Whether to display the fit between subjects' and model's belief
% difference
dispfit = false;

% Select the sequences to look at
% 1: all fully-stochastic parts
% 2: only fully-stochastic sequences
% 3: only fully-stochastic sequences that were correctly labeled
% 4: only first-part of stochastic-to-regular sequences
restopt = 1;

% Run simulations
% ~~~~~~~~~~~~~~~
Emergence_MC_ModelSimulations;

% Perform subject-specific regressions
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Select observations
randidx = Emergence_SelectFullyStochSeq(G, filter, restopt);

% Get subjects' P/D ratio
pHpgYpHdgYsub = cellfun(@(x,i) ...
    (x.BarycCoord(i,2) - x.BarycCoord(i,1)) ./ (x.BarycCoord(i,2) + x.BarycCoord(i,1)), ...
    G, randidx, 'uni', 0);

% Prepare output variable
MSE     = NaN(     nSub,nMod);
subavg  = NaN(nBin,nSub,nMod);
modavg  = NaN(nBin,nSub,nMod);

% For each model
for iMod = 1:nMod
    
    % Get ideal observer's P/D ratio
    pHpgYpHdgYmod = cellfun(@(x,i) (x(i,2) - x(i,1)) ./ (x(i,2) + x(i,1)), ...
        pMgY(:,:,iMod), randidx, 'uni', 0);
    
    % Get ideal observer's belief in the fully-stochastic hypothesis
    pHsgYmod = cellfun(@(x,i) x(i,3), pMgY(:,:,iMod), randidx, 'uni', 0);
    
    % For each subject
    for iSub = 1:nSub
        
        % Select variables
        subratioPD = cell2mat(pHpgYpHdgYsub(:,iSub));
        modratioPD = cell2mat(pHpgYpHdgYmod(:,iSub));
        
        % When the division is not possible (i.e. P+D = 0) set the ratio to 0
        subratioPD(isnan(subratioPD)) = 0;
        modratioPD(isnan(modratioPD)) = 0;
        
        % Create design matrix by adding additional confound regressors
        iobelinS = cell2mat(pHsgYmod(:,iSub));
        iointerac = modratioPD .* iobelinS;
        offset = ones(numel(modratioPD), 1);    
        desmat = [offset, modratioPD, iobelinS, iointerac];
        
        % Center predictors
        desmat = desmat - [0, mean(desmat(:,2:end), 1)];
        coi = 2; % regression coefficient of interest
        
        % Remove variance from all the predictors except the P/D ratio
        beta = regress(subratioPD, desmat);
        pidx = setdiff(1:numel(beta), coi);
        pred = desmat(:,pidx) * beta(pidx); % predictions based on beta values
        subratioPD2 = subratioPD - pred; % pseudo-residuals
        
        % Fit a regression model which restricts the regression coefficient
        % to be null or positive
        demeanmodratioPD = modratioPD - mean(modratioPD); % demean
        beta = regress(subratioPD2, demeanmodratioPD); % regress
        beta = max([0,beta]); % restrict the coefficient to be null or positive
        pred = demeanmodratioPD .* beta; % predicted values
        error = (subratioPD2 - pred) .^ 2; % compute SE
        MSE(iSub,iMod) = mean(error); % compute MSE
        
        % Get probability bins to later plot subjects vs. model
        % N.B. This deals with model predictions that lack variance: e.g.
        % in the case of the max version of the independent two-system
        % model, belief difference can only take 3 values: –1, 0, or 1
        if strcmpi(binmeth, 'unif') % bins of the same amplitude
            edges = linspace(-1, 1, nBin+1);
        elseif strcmpi(binmeth, 'equil') % bins with the same number of observations
            if numel(unique(modratioPD)) < nBin
                edges = movmean(unique(modratioPD), 2);
                edges = [edges(1)-1; edges(2:end); edges(end)+1];
            else, edges = prctile(modratioPD, linspace(0, 100, nBin+1));
            end
        else, error('Please check the binnig method that is provided');
        end
        prec = 0.001;
        edges = uniquetol(edges, prec);
        ncb = numel(edges) - 1;
        
        % Average IO P/D ratio and subject's residual P/D ratio
        [~,~,bins] = histcounts(modratioPD, edges);
        subavg(1:ncb,iSub,iMod) = arrayfun(@(i) mean(subratioPD2(bins == i)), 1:ncb);
        modavg(1:ncb,iSub,iMod) = arrayfun(@(i) mean(modratioPD( bins == i)), 1:ncb);
    end
end

% Center the MSE from the alternative models relatively to the normative
% two-system model
MSE = MSE(:,1:end-1) - MSE(:,end);

% Perform paired t-tests
[h,pval,tci,stats] = ttest(MSE);
Emergence_PrintTstats(pval,tci,stats);

% In case of models with 2 dimensions
nCrv = 1;
if dispcont && contains(SimuType, 'Predictability', 'IgnoreCase', true)
   nCrv = numel(unique(options(4,1:end-1)));
   nMod = (nMod-1) / nCrv;
   MSE = reshape(MSE, [nSub,nCrv,nMod]);
end

% Display correlation between subjects' residual ratio and ideal observer's
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Create a new window
if dispfit
    figure('Position', [1 631 1920 200]);
    
    % For each model
    for iMod = 1:nMod
        
        % Average over subjects
        xm = mean(modavg(:,:,iMod), 2)';
        xs = sem( modavg(:,:,iMod), 2)';
        ym = mean(subavg(:,:,iMod), 2)';
        ys = sem( subavg(:,:,iMod), 2)';
        
        % Create colormap used to color bins according to IO belief
        prec = 100;
        gege = [flipud(Emergence_Colormap('Blues', prec)); Emergence_Colormap('Reds', prec)];
        grid = linspace(-1, 1, prec*2);
        [~,colidx] = min(abs(xm'-grid), [], 2);
        
        % Show correlation between subjects' P/D residual ratio and IO P/D ratio
        subplot(1,nMod,iMod); hold('on');
        plot(xm+[-xs;xs], repmat(ym,2,1), 'k-'); % horizontal error bars
        plot(repmat(xm,2,1), ym+[-ys;ys], 'k-'); % vertical   error bars
        plot(xm, ym, 'k-', 'LineWidth', 1);
        scatter(xm, ym, 100, gege(colidx,:), 'filled', 'MarkerEdgeColor', 'k');
        
        % Customize the axes
        set(gca, 'Box', 'Off', 'Layer', 'Bottom');
        axis('square');
        
        % Add some text labels
        xlabel('Beliefs from IO');
        ylabel('Residual beliefs from subjects');
    end
end

% Display the mean squared error
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Average MSE over subjects
m = squeeze(mean(MSE));
s = squeeze(sem( MSE));

% Prepare a new window
figure('Position', [1 905 300 200]); hold('on');

% If continuous models, use lines
if dispcont
    x = contparam';
    for i = 1:nCrv
        plotMSEM(x, m(i,:), s(i,:), 0.15, modc(i,:), modc(i,:), 2);
    end
    xlim(sort(x([1,end]))); xlabel('Parameters');
    
% If discrete models, use dots instead
elseif ~dispcont
    x = 1:size(MSE,2);
    plot(x, m, 'k-', 'LineWidth', 2);
    plot(repmat(x, [2,1]), m+s.*[-1;1], 'k-');
    scatter(x, m, 100, modc(x,:), 'Filled', 'MarkerEdgeColor', 'k');
    xlim([0,nMod]); xlabel('Models');
end

% Overlap p-values on top of the plot
if nCrv == 1
    lima = [0.001 0.01 0.05];
    sz = [6,4,2];
    for i = 1:numel(lima)
        idx = pval < lima(i);
        limy = ylim;
        plot(x(idx), repmat(max(limy), [1,sum(idx)]), 'ks', ...
           'MarkerEdgeColor', 'None', ...
           'MarkerFaceColor', 'k', 'MarkerSize', sz(i));
       ylim(limy);
    end
end

% Customize the axes
set(gca, 'Box', 'Off', 'Layer', 'Top');
lab = {'Continuous', 'Predictability'};
if dispcont && useXlog, set(gca, 'XScale', 'log'); end

% Add some text labels
ylabel('Mean squared error');
title(sprintf('%s - %s', SimuType, mes));
