% This script shows mean squared error between different alternaive models
% and the subjects measured either on the abruptness of the detection
% dynamics of deterministic rules (pattern-specific) or on the average
% probability of the deterministic rule hypothesis in the fully-random
% parts of sequences
% 
% Copyright (c) 2018 Maxime Maheu

% Initialization
% ~~~~~~~~~~~~~~

% On which measure to quantify the MSE
mes = 'abruptness'; % 'abruptness' or 'falsealarm'

% Get model simulations
SimuType = 'IndependentDifferenceDiscrete';
Emergence_MC_ModelSimulations;

% Whether models are discrete or continuous
dispcont = contains(SimuType, 'Continuous', 'IgnoreCase', true);

% Compute abruptness of deterministic detection dynamics
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if strcmpi(mes, 'abruptness')
    
    % Get positions of the change point
    cp = cellfun(@(x) x.Jump+1/2, G(cidx{2},:), 'uni', 0);
    
    % Measure of abruptness between change and end points of deterministic
    % regularities
    abruptness = @(p) mean(abs(diff(p,2,1)));
    modbeh = cellfun(@(p,c) abruptness(p(c:end,2)), pMgY(cidx{2},:,:), repmat(cp, [1 1 nMod]));
    subbeh = cellfun(@(p,c) abruptness(p.BarycCoord(c:end,2)), G(cidx{2},:), cp);
    
    % Remove deterministic regularities that were not identified by subjects
    detecmask = (filter{2} == 3);
    subbeh(~detecmask) = NaN;
    modbeh(~detecmask) = NaN;
    
    % Compute the error between models and subjects
    MSE = (subbeh - modbeh) .^ 2;
    
    % Average over sequences
    MSE = squeeze(mean(MSE, 1, 'OmitNaN'));

% Compute average false alarm position
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
elseif strcmpi(mes, 'falsealarm')
    
    % Select fully-stochastic parts of sequences
    idxtrimap = Emergence_SelectFullyStochSeq(G, filter, 1);
    
    % Get average positions in the (pseudo-)deterministic hypothesis from
    % the models
    modbeh = cellfun(@(x,i) x(i,2), pMgY, repmat(idxtrimap, [1,1,nMod]), 'uni', 0);
    modbeh = squeeze(mat2cell(modbeh, nSeq, ones(1,nSub,1), ones(1,1,nMod)));
    modbeh = cellfun(@(x) mean(cell2mat(x)), modbeh);
    
    % Get average positions in the (pseudo-)deterministic hypothesis from
    % the subjects
    subbeh = cellfun(@(x,i) x.BarycCoord(i,2), G, idxtrimap, 'uni', 0);
    subbeh = squeeze(mat2cell(subbeh, nSeq, ones(1,nSub)))';
    subbeh = cellfun(@(x) mean(cell2mat(x)), subbeh);
    
    % Compute the error between models and subjects
    MSE = (subbeh - modbeh) .^ 2;
end

% Perform paired t-tests
if ~dispcont
    [h,pval,tci,stats] = ttest(MSE(:,1:end-1) - MSE(:,end));
    Emergence_PrintTstats(pval,tci,stats);
end

% Display the results
% ~~~~~~~~~~~~~~~~~~~

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
if     strcmpi(mes, 'abruptness'), ylim([0.01,1].*1e-3);
elseif strcmpi(mes, 'falsealarm'), ylim([2,200] .*1e-3);
end
set(gca, 'YScale', 'Log', 'Box', 'Off', 'Layer', 'Bottom');
if dispcont, set(gca, 'XScale', 'log');
else, xlim([0,nMod]); end

% Add some text labels
xlabel('Models');
ylabel('Mean squared error');
