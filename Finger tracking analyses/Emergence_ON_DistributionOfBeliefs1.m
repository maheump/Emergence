% This scripts computes histograms of finger's positions inside the
% triangular arena for the different types of sequences. Moreover, the
% marginal histograms (along the relevant dimension) for the two types of
% regularity (probabilistic and deterministic ones) are compared.
% 
% Copyright (c) 2018 Maxime Maheu

% Define option
% ~~~~~~~~~~~~~

% Whether to restrict to sequences that were accurately labeled
restodet = true;

% Get beliefs in regular sequences
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Prepare the output variable
idxtrimap = repmat({cellfun(@(x) false(N,1), G, 'UniformOutput', 0)}, 1, 5);

% For sequences entailing regularities
for iHyp = 1:2
    
    % Get change point position
    cp = cellfun(@(x) x.Jump, G(cidx{iHyp},:));
    
    % Get logical indices of post-change-point observations
    idx = arrayfun(@(c) [false(c-1/2,1); true(N-c+1/2,1)], cp, 'UniformOutput', 0);
    
    % Keep only indices for sequences with the current type of regularity
    idxtrimap{iHyp}(cidx{iHyp},:) = idx;
    
    % Keep only sequences that have been correctly identified
    if restodet
        detecmask = false(size(G));
        detecmask(cidx{iHyp},:) = (filter{iHyp} == 1 | filter{iHyp} == 3);
        idxtrimap{iHyp}(~detecmask) = {false(N,1)};
    end
end

% Get beliefs in fully-stochastic sequences
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% First part of stochastic-to-regular sequences
idxtrimap{3} = Emergence_SelectFullyStochSeq(G, filter, 4);

% Fully-stochastic sequences that were correctly labelled
if    ~restodet, idxtrimap{4} = Emergence_SelectFullyStochSeq(G, filter, 2);
elseif restodet, idxtrimap{4} = Emergence_SelectFullyStochSeq(G, filter, 3);
end

% All fully-stochastic parts (first part of stochastic-to-regular sequences
% and fully-stochastic ones)
idxtrimap{5} = Emergence_SelectFullyStochSeq(G, filter, 4);

% Compute triangular histograms
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Get the number of triangular histogram maps we have to compute
nMap = numel(idxtrimap);

% Prepare output variables
trajmap = cell(1,nMap+1);
gppoints = cell(1,nMap);

% For each density map to be created
for iMap = 1:nMap

    % Get the finger's positions in the current condition
    points = cellfun(@(x,y) x.BarycCoord(y,:), ...
        D, idxtrimap{iMap}, 'UniformOutput', 0);
    
    % Compute triangular histogram
    [trajmap{iMap+1}, gppoints{iMap}, xgrid, ygrid, mask] = Emergence_TriHist(points);
end

% Also add a map of the difference between probabilistic and deterministic
% regularities
trajmap{1} = (trajmap{3} ./ max(trajmap{3}(:))) ...
           - (trajmap{2} ./ max(trajmap{2}(:)));
       
% Display triangular and marginal histograms
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Prepare the window
figure('Position', [461 253 1000 650]);

% Number of bins for the marginal histrograms
nBin = 30;

% Create colormaps
cmaps = {[flipud(Emergence_Colormap('Blues')); (Emergence_Colormap('Reds'))], ...
         Emergence_Colormap('Blues'), Emergence_Colormap('Reds'), ...
         Emergence_Colormap('Greens'), Emergence_Colormap('Greens'), Emergence_Colormap('Greens')};

% For each density map
for iMap = 1:nMap+1
    sp = subplot(2,3,iMap);
    
    % Display the density maps of finger's position
    imagesc(xgrid, ygrid, trajmap{iMap}, 'AlphaData', mask); hold('on');
    contour(xgrid, ygrid, trajmap{iMap}, 5, 'k-', 'LineWidth', 1/2);
    image(xgrid, ygrid, repmat(~mask, [1,1,3]), 'AlphaData', ~mask);
    
    % Display an empty triangle
    tr = Emergence_PlotTriInfo(tricc, tricol);
    
    % Draw marginal histograms on the limits of the triangle
	if iMap > 1, Emergence_PlotMargHist(gppoints{iMap-1}, nBin); end
    
    % Overlap the grid corresponding to the resoution of the marginal
    % histrograms
    Emergence_PlotGridOnTri(10);
    
    % Customize the axes
    set(gca, 'LineWidth', 1, 'FontSize', 15);
    axis('xy'); axis('off'); axis('equal');
    
    % Customize the colormap
    colormap(sp, cmaps{iMap});
    if iMap == 1, caxis([-abs(max(caxis)), abs(max(caxis))]);
    else, caxis([0,max(cellfun(@(x) max(x(:)), trajmap([2,3,end])))]);
    end
    ax = get(gca, 'Position');
    cbrpos = [ax(1)+ax(3)*3/4 ax(2) 0.01 ax(4)*1/3];
    cbr = colorbar('Position', cbrpos, 'LineWidth', 1);
    cbr.Label.String = {'Normalized density map', '(in log scale)'};
end

% Use the same zoom level for all the maps
ScaleAxis('x'); ScaleAxis('y');

% Save the figure
if isfield(D{1}, 'Seq'), save2pdf(fullfile(ftapath, 'figs', 'F_D_DensityMapsS.pdf'));
else, save2pdf(fullfile(ftapath, 'figs', 'F_D_DensityMapsIO.pdf'));
end
