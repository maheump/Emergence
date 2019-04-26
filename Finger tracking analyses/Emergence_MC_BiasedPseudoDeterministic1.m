% This script shows the impact of pseudo-deterministic hypothesis learning
% higher-order transitions on the likelihood of false alarms, in
% particular, the likelihood of deterministic false alarms.
% 
% Copyright (c) 2018 Maxime Maheu

% Load simulations
% ~~~~~~~~~~~~~~~~

SimuType = 'PseudoDeterministic'; % 'PseudoDeterministic' or 'BiasedPseudoDeterministic'
Emergence_MC_ModelSimulations;

% Perform the analysis
% ~~~~~~~~~~~~~~~~~~~~

% Select fully-stochastic parts of sequences
idxtrimap = Emergence_SelectFullyStochSeq(G, filter, 1);

% Prepare the output variables
trajmap = cell(1,nMod);
CorCoef = [];

% For each model
for iMod = 1:nMod
    
    % Compute the 2D histogram
    points = cellfun(@(x,y) x(y,:), pMgY(:,:,iMod), idxtrimap, 'uni', 0);
    [trajmap{iMod}, ~, xgrid, ygrid, mask] = Emergence_TriHist(points);
end

% Display triangular histrograms
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Prepare the window
figure('Position', [201 453 1520 250]);

% For each density map
for iMod = 1:nMod
    sp = subplot(1,nMod,iMod);
    
    % Display the density maps of finger's position
    imagesc(xgrid, ygrid, trajmap{iMod}, 'AlphaData', mask); hold('on');
    contour(xgrid, ygrid, trajmap{iMod}, 5, 'k-', 'LineWidth', 1/2);
    image(xgrid, ygrid, repmat(~mask, [1,1,3]), 'AlphaData', ~mask);
    
    % Display an empty triangle
    fill(tricc(:,1), tricc(:,2), 'k', 'FaceColor', 'None', 'LineWidth', 2);
    
    % Overlap the grid corresponding to the resoution of the marginal
    % histrograms
    Emergence_PlotGridOnTri(10);
    
    % Customize the axes
    set(gca, 'LineWidth', 1, 'FontSize', 15);
    axis('xy'); axis('off'); axis('equal');
    axis([0,1,0,sqrt(3)/2]);
    
    % Customize the colormap
    colormap(sp, Emergence_Colormap('Greens'));
    caxis([0,max(trajmap{end}(:))]);
    
    % Add some text labels
    title(options{4,iMod});
end

% Compute correlation between deterministic "false alarms" and the order of
% the transitions that are monitored
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Prepare the output variable
avgdetfa = NaN(nSub,nMod);

% For each model
for iMod = 1:nMod
    
    % Compute the 2D histogram
    avgbel = cellfun(@(x,y) mean(x(y,2)), pMgY(:,:,iMod), idxtrimap);
    
    % Average over sequences
    avgdetfa(:,iMod) = mean(avgbel, 1);
end

% Compare "false alarm" likelihood across order of transitions
[~,pval,tci,stats] = ttest(avgdetfa(:,1:end-1) - avgdetfa(:,end));
Emergence_PrintTstats(pval,tci,stats);

% Correlate ?false alarm? likelihood with order of transitions
corcoef = cellfun(@(x) Emergence_Regress(x, 2:nMod, 'CC', 'r'), ...
    mat2cell(avgdetfa(:,1:nMod-1), ones(nSub,1), nMod-1)); 
[~,pval,tci,stats] = ttest(corcoef);
Emergence_PrintTstats(pval,tci,stats);
