% This script shows the error of regressions between belief difference in
% subjects and different models.
% 
% Copyright (c) 2018 Maxime Maheu

% Define options
% ~~~~~~~~~~~~~~

% Get model simulations
SimuType = 'IndependentDifferenceDiscrete';
Emergence_MC_ModelSimulations;

% Whether models are discrete or continuous
dispcont = contains(SimuType, 'Continuous', 'IgnoreCase', true);

% Define the number of bins to use when plotting the results
nBin = 10;

% Define the type of binning method
binmeth = 'equil';

% Select the sequences to look at
% 1: all fully-stochastic parts
% 2: only fully-stochastic sequences
% 3: only fully-stochastic sequences that were correctly labeled
% 4: only first-part of stochastic-to-regular sequences
restopt = 1;

% Perform subject-specific regressions
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Select observations
randidx = Emergence_SelectFullyStochSeq(G, filter, restopt);

% Get subjects' P/D ratio
pHpgYpHdgYsub = cellfun(@(x,i) ...
    (x.BarycCoord(i,2) - x.BarycCoord(i,1)) ./ (x.BarycCoord(i,2) + x.BarycCoord(i,1)), ...
    G, randidx, 'uni', 0);

% Prepare output variable
coef    = NaN(     nSub,nMod);
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
        iobelinS  = cell2mat(pHsgYmod(:,iSub));
        iointerac = modratioPD .* iobelinS;
        offset    = ones(numel(modratioPD), 1);    
        desmat    = [offset, modratioPD, iobelinS, iointerac];
        
        % Center predictors
        desmat = desmat - [0, mean(desmat(:,2:end), 1)];
        
        % Get probability bins
        regofint = 2;
        if strcmpi(binmeth, 'unif') % bins of the same amplitude
            edges = linspace(-1, 1, nBin+1);
        elseif strcmpi(binmeth, 'equil') % bins with the same number of observations
            edges = prctile(desmat(:,regofint), linspace(0, 100, nBin+1));
        else, error('Please check the binnig method that is provided');
        end
        
        % In the case the number of bins is larger than the number of
        % unique values predicted by the model (e.g. in the case of the max
        % version of the independent two-system model), use less bins
        unval = unique(desmat(:,regofint));
        ncb = numel(unval);
        if nBin > ncb
            edges = movmean(unval, 2);
            edges = [edges(1)-1; edges(2:end); edges(end)+1];
        else, ncb = nBin;
        end
        
        % Regress subject's P/D ratio against several IO predictors
        beta = regress(subratioPD, desmat);
        coef(iSub,iMod) = beta(2);
        
        % Compute error of the regression model
        nullmse = (subratioPD - mean(subratioPD)) .^2;
        modmse  = (subratioPD - (desmat * beta))  .^2;
        MSE(iSub,iMod) = mean(modmse - nullmse);
        
        % Remove variance from all the predictiors except the P/D ratio
        pidx = setdiff(1:numel(beta), regofint);
        pred = desmat(:,pidx) * beta(pidx); % predictions based on beta values
        subratioPD2 = subratioPD - pred; % pseudo-residuals
        
        % Average IO P/D ratio and subject's residual P/D ratio
        [~,~,bins] = histcounts(desmat(:,regofint), edges);
        subavg(1:ncb,iSub,iMod) = arrayfun(@(i) mean(subratioPD2(bins == i)),     1:ncb);
        modavg(1:ncb,iSub,iMod) = arrayfun(@(i) mean(desmat(bins == i,regofint)), 1:ncb);
    end
end

% Perform paired t-tests
if ~dispcont
    [h,pval,tci,stats] = ttest(MSE(:,1:end-1) - MSE(:,end));
    Emergence_PrintTstats(pval,tci,stats);
end

% Display correlation between subjects' residual ratio and ideal observer's
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Create a new window
if ~dispcont
    figure('Position', [302 905 800 200]);
    
    % For each model
    for iMod = 1:nMod
        
        % Average over subjects
        xm = mean(modavg(:,:,iMod), 2)';
        xs = sem( modavg(:,:,iMod), 2)';
        ym = mean(subavg(:,:,iMod), 2)';
        ys = sem( subavg(:,:,iMod), 2)';
        
        % Regress average P/D ratio from subjects and the ideal observer
        B = Emergence_Regress(ym, xm, 'TLS', {'beta0', 'beta1'});
        confint = Emergence_Regress(ym, xm, 'TLS', 'confint');
        xval = Emergence_Regress(ym, xm, 'TLS', 'confintx');
        
        % Create colormap used to color bins according to IO belief
        prec = 100;
        gege = [flipud(Emergence_Colormap('Blues', prec)); Emergence_Colormap('Reds', prec)];
        grid = linspace(-1/2, 1/2, prec*2);
        [~,colidx] = min(abs(xm'-grid), [], 2);
        
        % Display origin (because variables are centered)
        subplot(1,nMod,iMod); hold('on');
        
        % Display the regression line
        fill([xval, fliplr(xval)], [confint(1,:), fliplr(confint(2,:))], 'k', ...
               'EdgeColor', 'none', 'FaceColor', g, 'FaceAlpha', 1/2);
        plot(xval, xval.*B(2) + B(1), 'k-', 'LineWidth', 3);
        
        % Show correlation between subjects' P/D residual ratio and IO P/D ratio
        plot(xm+[-xs;xs], repmat(ym,2,1), 'k-'); % horizontal error bars
        plot(repmat(xm,2,1), ym+[-ys;ys], 'k-'); % vertical   error bars
        plot(xm, ym, 'k-', 'LineWidth', 1);
        scatter(xm, ym, 100, gege(colidx,:), 'filled', 'MarkerEdgeColor', 'k');
        
        % Display help lines
        ylim(max(abs(ylim)) .* [-1,1]);
        limy = ylim;
        plot(zeros(1,2), [-1,1], '-', 'Color', g);
        plot([-1,1], zeros(1,2), '-', 'Color', g);
        plot(limy, limy, '-', 'Color', g);
        ylim(limy);
        
        % Customize the axes
        set(gca, 'Box', 'Off'); axis('square');
        
        % Add some text labels
        xlabel('Centered beliefs from IO');
        ylabel('Residual centered beliefs from subjects');
    end
end

% Display mean squared error of the regression
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Get data to plot
modidx = 1:nMod-1;

% Get x-values
if dispcont, x = sigmparam';
else, x = modidx';
end

% Get models' color
c = modc(modidx,:);

% Average MSE over subjects
m = mean(MSE(:,modidx))';
s = sem( MSE(:,modidx))';

% Prepare a new window
figure('Position', [1 905 300 200]); hold('on');

% Display error made by the fully deterministic hypothesis
plotMSEM(x([1,end]), repmat(mean(MSE(:,end)),1,2), ...
    repmat(sem(MSE(:,end)),1,2), 1/2, modc(end,:));

% If continuous models, use lines
if dispcont
    y = (1:numel(x))';
    fill([x', flipud(x)'], [(m+s)', flipud(m-s)'], 'k', 'EdgeColor', 'None', ...
        'CData', [y;flipud(y)], 'FaceColor', 'Interp'); alpha(1/10);
    surface('XData', [x x], 'YData', [m m], 'ZData', zeros(numel(x),2), ...
            'CData', [y y], 'FaceColor', 'None', 'EdgeColor', 'Interp', ...
            'Marker', 'none', 'LineWidth', 2);
    colormap(c); caxis(y([1,end]));
    
% If discrete models, use dots instead
elseif ~dispcont
    plot(repmat(x', [2,1]), m'+s'.*[-1;1], 'k-');
    scatter(x, m, 100, c, 'Filled', 'MarkerEdgeColor', 'k');
end

% Customize the axes
set(gca, 'YScale', 'Log', 'Box', 'Off', 'Layer', 'Bottom');
if dispcont, set(gca, 'XScale', 'log');
else, xlim([0,nMod]); end
ylim([-0.02,-0.004]);

% Add some text labels
xlabel('Models');
ylabel('Mean squared error');
