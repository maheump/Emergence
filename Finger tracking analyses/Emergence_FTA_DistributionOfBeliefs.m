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
    idxtrimap{iHyp} = cellfun(@(x) (x.Jump-1/2):N, G, 'UniformOutput', 0);
    idxtrimap{iHyp}(setdiff(1:size(D,1), cidx{iHyp}),:) = {[]};
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
ng = numel(xgrid);
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
    idx = find(sum(~cellfun(@isempty, points), 2) > 0);
    lists  = idxtrimap{iMap}(idx,:);
    points = points(idx,:);
    interppoints = cellfun(@(x,v) interp1(x', v, x(1):dtint:x(end)), ...
        lists, points, 'UniformOutput', 0);
    gppoints{iMap} = cell2mat(interppoints(:));
    
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
figure('Position', [1 805 800 300]);

% Number of bins for the marginal histrograms
nBin = 30;

% For each density map
for iMap = 1:nMap+1
    sp = subplot(1, 4, iMap);

    % Display the density maps of finger's position
    imagesc(xgrid, ygrid, trajmap{iMap}, 'AlphaData', mask); hold('on');
    contour(xgrid, ygrid, trajmap{iMap}, 5, 'k-', 'LineWidth', 1/2);
    image(xgrid, ygrid, repmat(~mask, [1,1,3]), 'AlphaData', ~mask);

    % Display an empty triangle
    tr = Emergence_PlotTriInfo(tcn, tricol);
    
    % Customize the colormap
    if iMap == 1
        colormap(sp, [flipud(cbrewer2('YlGnBu', 2000)); ...
                            (cbrewer2('YlOrRd', 2000))]);
        caxis([-abs(max(caxis)), abs(max(caxis))]);
    else, colormap(sp, LoadCWcbr); caxis([0, max(caxis)]);
    end
    cbr = colorbar('Location', 'SouthOutside', 'LineWidth', 1);
    cbr.Label.String = {'Normalized density map', '(in log scale)'};
    
    % Draw marginal histograms on the limits of the triangle
    %if iMap > 1, Emergence_PlotMargHist(gppoints{iMap-1}, nBin); end
    
    % Overlap the grid corresponding to the resoution of the marginal
    % histrograms
    Emergence_PlotGridOnTri(nBin/2);
    
    % Customize the axes
    set(gca, 'LineWidth', 1, 'FontSize', 15);
    axis('xy'); axis('off'); axis('equal');
end

% Use the same zoom level for all the maps
ScaleAxis('c', 1:3); ScaleAxis('x'); ScaleAxis('y');

% Save the figure
if isfield(D{1}, 'Seq'), save2pdf('figs/F_D_DensityMapsS.pdf');
else, save2pdf('figs/F_D_DensityMapsIO.pdf');
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
figure('Position', [802 805 200 300]); lgd = NaN(1,2);

% Useful variables
x = pg+1/(nBin-1)/2; % x-axis
d = [1,-1]; % to plot it symmetrically

% For each type of regularity, plot the histogram in a symmetrical manner
for iHyp = 1:2
    lgd(iHyp) = bar(x, d(iHyp).*avgbins(iHyp,:), 0.9, 'FaceColor', ...
        tricol(iHyp,:), 'EdgeColor', 'none'); hold('on');
    plot(repmat(x', 2, 1), d(iHyp).*avgbins(iHyp,:) + ...
        sembins(iHyp,:).*[-1;1], 'k-', 'LineWidth', 1/2);    
end

% Overlap a line along zero 
plot([0,1], zeros(1,2),  'k-');

% Customize the axes
SymAxis('y')
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
if isfield(D{1}, 'Seq'), save2pdf('figs/F_D_MargHistS.pdf');
else, save2pdf('figs/F_D_MargHistIO.pdf');
end
