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

% Get the different conditions
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Length of the sequences
N = numel(G{1}.Seq);

% Create indices...
idxtrimap = cell(1,4);

% First part of stochastic-to-regular sequences
idxtrimap{1} = cellfun(@(x) find(x.Gen == 1), G, 'UniformOutput', 0);
idxtrimap{1}(setdiff(1:size(G,1), cell2mat(cidx(1:2))),:) = {[]};

% Get subjects' responses
detecmask = cellfun(@(x) x.Questions(2), G);
detecmask(isnan(detecmask)) = 3;

% For sequences entailing regularities
for iHyp = 1:3
    
    % Get indices of post-change-point observations
    if iHyp < 3
        idxtrimap{iHyp+1} = cellfun(@(x) (x.Jump-1/2):N, G, 'UniformOutput', 0);
    elseif iHyp == 3
        idxtrimap{iHyp+1} = repmat({1:N}, size(G));
    end
    
    % Keep only indices for sequences with the current type of regularity
    idxtrimap{iHyp+1}(setdiff(1:size(G,1), cidx{iHyp}),:) = {[]};
    
    % Keep only sequences that have been correctly identified
    if restodet, idxtrimap{iHyp+1}(detecmask ~= iHyp) = {[]}; end
end

% Define properties of the triangular histogram
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Triangle's specifications
tcn = [0, sqrt(3)/2; 1, sqrt(3)/2; 1/2, 0];

% Options for creating the density map
j = 0.01;
xgrid = tcn(1,1):j:tcn(2,1);
ygrid = tcn(3,2):j:tcn(1,2);

% Create a mask that specifies which cells of the matrices computed below
% are outside the triangle
ng = 1/j+1;
T = tril(ones(ng));
T = [flip(T,2),T];
mask = round(imresize(T, [round(ng/2*tan(pi/3)), ng]));

% Options for the gaussian smooth
filtSigma = 2;
filtWidth = 2*ceil(2*filtSigma)+1;
imageFilter = fspecial('Gaussian', filtWidth, filtSigma);

% Define the degree of interpolation that we apply on the trajectory
dtint = 1/3; % must be <= 1, 1 is no interpolation

% Perform the analysis
% ~~~~~~~~~~~~~~~~~~~~

% Get the number of triangular histogram maps we have to compute
nMap = numel(idxtrimap);

% Prepare output variables
trajmap = cell(1,nMap+1);
subpoints = cell(nSub,nMap);
gppoints = cell(1,nMap);

% For each density map to be created
for iMap = 1:nMap
    
    % Get the finger's positions in the current condition
    points = cellfun(@(x,y) x.BarycCoord(y,:), ...
        D, idxtrimap{iMap}, 'UniformOutput', 0);
    
    % Get finger's positions of each subjects
    for iSub = 1:nSub
        subpoints{iSub,iMap} = cell2mat(points(:,iSub));
    end
    
    % Interpolate the trajectories
    interppoints = cellfun(@(x) interp1(1:size(x,1), x, 1:dtint:size(x,1)), ...
        points(~cellfun(@isempty, points)), 'UniformOutput', 0);
    gppoints{iMap} = cell2mat(interppoints);
    
    % Compute the density map
    cc = gppoints{iMap}*tcn;
    density = hist3(cc, 'Edges', {xgrid, ygrid})';
    
    % Smooth the density map
    density = convn(density, imageFilter, 'same');
    
    % Convert it to log scale
    density = log(density+1);
    
    % Normalize it
    trajmap{iMap+1} = density ./ sum(density(:));
end

% Also add a map of the difference
trajmap{1} = (trajmap{4} ./ max(trajmap{4}(:))) ...
           - (trajmap{3} ./ max(trajmap{3}(:)));
       
% Display triangular and marginal histrograms
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Prepare the window
figure('Position', [1 600 1920 500]);

% Number of bins for the marginal histrograms
nBin = 30;

% Create colormaps
prec = 2000;
cmaps = {[flipud(cbrewer2('Blues', prec)); (cbrewer2('Reds', prec))], ...
         cbrewer2('Greens', prec), ...
         cbrewer2('Blues',  prec), cbrewer2('Reds', prec), ...
         cbrewer2('Greens', prec)};

% For each density map
for iMap = 1:nMap+1
    sp = subplot(1, 5, iMap);
    
    % Display the density maps of finger's position
    imagesc(xgrid, ygrid, trajmap{iMap}, 'AlphaData', mask); hold('on');
    contour(xgrid, ygrid, trajmap{iMap}, 5, 'k-', 'LineWidth', 1/2);
    image(xgrid, ygrid, repmat(~mask, [1,1,3]), 'AlphaData', ~mask);
    
    % Display an empty triangle
    tr = Emergence_PlotTriInfo(tcn, tricol);
    
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
    else, caxis([0, max(caxis)]);
    end
    ax = get(gca, 'Position');
    cbrpos = [ax(1)+ax(3)*3/4 ax(2) 0.01 ax(4)*1/3];
    cbr = colorbar('Position', cbrpos, 'LineWidth', 1);
    cbr.Label.String = {'Normalized density map', '(in log scale)'};
end

% Use the same zoom level for all the maps
ScaleAxis('c', 1:3); ScaleAxis('x'); ScaleAxis('y');

% Save the figure
if isfield(D{1}, 'Seq'), save2pdf(fullfile(ftapath, 'figs', 'F_D_DensityMapsS.pdf'));
else, save2pdf(fullfile(ftapath, 'figs', 'F_D_DensityMapsIO.pdf'));
end
