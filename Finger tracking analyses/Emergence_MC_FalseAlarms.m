% This script shows the impact of pseudo-deterministic hypothesis learning
% higher-order transitions on the likelihood of false alarms, in
% particular, the likelihood of deterministic false alarms.
% 
% Copyright (c) 2018 Maxime Maheu

% Define options
% ~~~~~~~~~~~~~~

% Get model simulations
SimuType = 'PseudoDeterministic';
Emergence_MC_ModelSimulations;

% Define the number of bins to use
nBin = 10;

% Select the sequences to look at
% 1: all fully-stochastic parts
% 2: only fully-stochastic sequences
% 3: only fully-stochastic sequences that were correctly labeled
% 4: only first-part of stochastic-to-regular sequences
restopt = 1;

% Select observations
idxtrimap = Emergence_SelectFullyStochSeq(G, filter, restopt);

% Compute triangular histograms
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Prepare output variables
avgtrihist = cell(1,nMod);
trihist    = cell(1,nSub);

% For each density map to be created
for iMod = 1:nMod
    
    % For each subject
    points = cell(1,nSub);
    for iSub = 1:nSub
        
        % Get the finger's positions in the current condition
        points{iSub} = cellfun(@(x,y) x(y,:), ...
            pMgY(:,iSub,iMod), idxtrimap(:,iSub), 'uni', 0);
        
        % Compute triangular histogram
        [trihist{iSub}, ~, xgrid, ygrid, mask] = ...
            Emergence_TriHist(points{iSub}, 1./(nBin*10));
    end
    
    % Average over subjects
    avgtrihist{iMod} = mean(cell2mat(reshape(trihist, [1,1,nSub])), 3);
end

% Display triangular histograms
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Prepare the window
figure('Position', [201 453 1520 250]);

% For each density map
for iMod = 1:nMod
    sp = subplot(1,nMod+1,iMod);
    
    % Display the density maps of finger's position
    imagesc(xgrid, ygrid, avgtrihist{iMod}, 'AlphaData', mask); hold('on');
    contour(xgrid, ygrid, avgtrihist{iMod}, 11, 'k-', 'LineWidth', 1/2);
    image(xgrid, ygrid, repmat(~mask, [1,1,3]), 'AlphaData', ~mask);
	
    % Overlap a grid
    Emergence_PlotGridOnTri(nBin);
    
    % Display an empty triangle
    Emergence_PlotTriInfo(tricc, tricol);
    
    % Customize the axes
    set(gca, 'LineWidth', 1, 'FontSize', 15);
    axis('xy'); axis('off'); axis('equal');
    
    % Add a colorbar
    ax = get(gca, 'Position');
    cbrpos = [ax(1)+ax(3)*3/4 ax(2) 0.01 ax(4)*1/3];
    cbr = colorbar('Position', cbrpos, 'LineWidth', 1);
    
    % Customize the colormap
    colormap(sp, Emergence_Colormap('Rainbow', 1e5));
    caxis([0,1/2*max(caxis)]);
    set(gca,'ColorScale','log');
    cbr.Label.String = {'Log-density'};
    cbr.Ruler.Scale = 'log';
    cbr.Ruler.MinorTick = 'on';
end
ScaleAxis('c');

% Display overlap
% ~~~~~~~~~~~~~~~

% Overlap the top 10% of the histogram
sp = subplot(1,nMod+1,nMod+1); hold('on');
toto = cellfun(@(x) x ./ max(x(:)), avgtrihist, 'uni', 0);
for iMod = 1:nMod
    contour(xgrid, ygrid, toto{iMod}, [0.9,1], ...
        'EdgeColor', modc(iMod,:), 'LineWidth', 1); hold('on');
end

% Overlap a grid
Emergence_PlotGridOnTri(2);

% Display an empty triangle
Emergence_PlotTriInfo(tricc, tricol);

% Customize the axes
set(gca, 'LineWidth', 1, 'FontSize', 15);
axis('xy'); axis('off'); axis('equal');
axis([0.2,0.8,0,0.5]);
