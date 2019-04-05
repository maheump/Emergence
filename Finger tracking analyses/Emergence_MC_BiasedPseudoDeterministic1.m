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
    points = cellfun(@(x,y) x(y,:), pMgY(:,:,iMod), idxtrimap, 'UniformOutput', 0);
    [trajmap{iMod}, ~, xgrid, ygrid, mask] = Emergence_TriHist(points);
end

% Display triangular histrograms
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Prepare the window
figure('Position', [201 453 1520 250]);

% For each density map
for iMod = 1:nMod
    subplot(1,nMod,iMod);
    
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
    colormap(cbrewer2('Greens', 2000));
    caxis([0,max(trajmap{end}(:))]);
    
    % Add some text labels
    title(options{4,iMod});
end
