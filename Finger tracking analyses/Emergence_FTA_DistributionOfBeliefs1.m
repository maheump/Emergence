% This scripts computes histograms of finger's positions inside the
% triangular arena for the different types of sequences. Moreover, the
% marginal histograms (along the relevant dimension) for the two types of
% regularity (probabilistic and deterministic ones) are compared.
%
% Copyright (c) 2018 Maxime Maheu

%% COMPUTE HISTOGRAMS OF FINGER'S POSITIONS
%  ========================================

% Get the different conditions
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Length of the sequences
N = numel(G{1}.Seq);

% Create indices...
idxtrimap = cell(1,3);

% For sequences entailing regularities
for iHyp = 1:2
    
    % Get indices of post-change-point observations
    idxtrimap{iHyp} = cellfun(@(x) (x.Jump-1/2):N, G, 'UniformOutput', 0);
    
    % Keep only indices for sequences with the current type of regularity
    idxtrimap{iHyp}(setdiff(1:size(G,1), cidx{iHyp}),:) = {[]};
    
    % Keep only sequences that have been correctly identified
    detecmask = zeros(size(G));
    detecmask(cidx{iHyp},:) = cellfun(@(xval) xval.Questions(2) == iHyp, G(cidx{iHyp},:));
    idxtrimap{iHyp}(~detecmask) = {[]};
end

% For stochastic parts (both the fully-stochastic sequences and the
% stochastic parts at the beginning of sequences with a regularity)
idxtrimap{end} = cellfun(@(x) find(x.Gen == 1), G, 'UniformOutput', 0);

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
trajmap{1} = (trajmap{3} ./ max(trajmap{3}(:))) ...
           - (trajmap{2} ./ max(trajmap{2}(:)));
       
% Display triangular and marginal histrograms
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Prepare the window
figure('Position', [1 600 1920 500]);

% Number of bins for the marginal histrograms
nBin = 30;

% Create colormaps
cmaps = {[flipud(cbrewer2('Blues', 2000)); (cbrewer2('Reds', 2000))], ...
         cbrewer2('Blues', 2000), cbrewer2('Reds', 2000), cbrewer2('Greens', 2000)};

% For each density map
for iMap = 1:nMap+1
    sp = subplot(1, 4, iMap);
    
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

%% COMPARE MARGINAL HISTOGRAMS FOR PROBABILISTIC AND DETERMINISTIC REGULARITIES
%  ============================================================================

% Create bins of beliefs
% ~~~~~~~~~~~~~~~~~~~~~~

% Sample uniformaly the probability continuum
nBin = 51;
pg = linspace(0, 1, nBin)';

% Prepare the output variable
bins = NaN(2,nBin,nSub);

% For each type of regularity, get marginal histrograms of finger's
% position along the relevant dimension
for iHyp = 1:2
    y = cellfun(@(x) histc(x(:,iHyp), pg)', subpoints(:,iHyp), 'UniformOutput', 0);
    y = cellfun(@(x) x./sum(x), y, 'UniformOutput', 0); % normalize it
    bins(iHyp,:,:) = cell2mat(y)';
end

% Average over subjects
avgbins = mean(bins, 3);
sembins = sem(bins, 3);

% Compare marginal histograms
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Prepare a new window
figure('Position', [1 225 200 300]); lgd = NaN(1,2);

% Useful variables
x = pg+1/(nBin-1)/2; % x-axis
d = [1,-1]; % to plot it symmetrically

% For each type of regularity, plot the histogram in a symmetrical manner
for iHyp = 1:2
    lgd(iHyp) = bar(x, d(iHyp).*avgbins(iHyp,:), 0.9, 'FaceColor', ...
        tricol(iHyp,:), 'EdgeColor', 'none'); hold('on');
    plot(repmat(x', 2, 1), d(iHyp).*avgbins(iHyp,:) + ...
        sembins(iHyp,:).*[-1;1], 'k-', 'LineWidth', 1/2);
    plot([0,1], (-2*(iHyp==2)+ones(1,2))/nBin, '-', 'Color', g, 'LineWidth', 1);
end

% Overlap a line along zero 
plot([0,1], zeros(1,2),  'k-');

% Customize the axes
ylim([-max(abs(ylim)), max(abs(ylim))]);
set(gca, 'Box', 'Off', 'XLim', [0,1]);
set(gca, 'YTick', get(gca, 'YTick'), 'YTickLabel', cellfun(@(x) ...
    sprintf('%1.2f', x), num2cell(abs(get(gca, 'YTick'))), 'UniformOutput', 0));
view([-90, 90]);

% Display bins for which the difference is significant
dif = squeeze(bins(1,:,:) - bins(2,:,:));
h = ttest(dif');
h(isnan(h)) = 0;
plot(pg(logical(h)), max(ylim)-0.05*diff(ylim), 'kp', ...
    'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');

% Add some text labels
xlabel('Strength of beliefs');
ylabel('Frequency');
legend(lgd, cellfun(@(x) x(1:5), proclab(1:2), 'UniformOutput', 0), ...
    'Orientation', 'Horizontal', 'Location', 'NorthOutside');

% Save the figure
if isfield(D{1}, 'Seq'), save2pdf(fullfile(ftapath, 'figs', 'F_D_MargHistS.pdf'));
else, save2pdf(fullfile(ftapath, 'figs', 'F_D_MargHistIO.pdf'));
end
