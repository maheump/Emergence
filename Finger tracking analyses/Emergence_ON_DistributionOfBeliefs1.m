% This scripts computes histograms of finger's positions inside the
% triangular arena for the different types of sequences. Moreover, the
% marginal histograms (along the relevant dimension) for the two types of
% regularity (probabilistic and deterministic ones) are compared.
% 
% Copyright (c) 2018 Maxime Maheu

% Define option
% ~~~~~~~~~~~~~

% Define the number of bins to use
nBin = 10;

% Whether to restrict to sequences that were accurately labeled
restodet = true;

% Get beliefs in regular sequences
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Prepare the output variable
idxtrimap = repmat({cellfun(@(x) false(N,1), G, 'uni', 0)}, 1, 5);

% For sequences entailing regularities
for iHyp = 1:2
    
    % Get change point position
    cp = cellfun(@(x) x.Jump, G(cidx{iHyp},:));
    
    % Get logical indices of post-change-point observations
    idx = arrayfun(@(c) [false(c-1/2,1); true(N-c+1/2,1)], cp, 'uni', 0);
    
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
idxtrimap{5} = Emergence_SelectFullyStochSeq(G, filter, 1);

% Compute triangular histograms
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Get the number of triangular histogram maps we have to compute
nMap = numel(idxtrimap);

% Prepare output variables
avgtrihist  = cell(1,nMap+1);
avgmarghist = cell(1,nMap+1);

% For each density map to be created
for iMap = 1:nMap
    trihist  = cell(1,nSub);
    marghist = cell(1,nSub);
    
    % For each subject
    for iSub = 1:nSub
        
        % Get the finger's positions in the current condition
        points = cell2mat(cellfun(@(x,y) x.BarycCoord(y,:), ...
            D(:,iSub), idxtrimap{iMap}(:,iSub), 'uni', 0));
        
        % Compute triangular histogram
        gridprec = 0.01;
        [trihist{iSub}, ~, xgrid, ygrid, mask] = Emergence_TriHist(points, 1./(nBin*10));

        % Compute the marginal histograms along each dimension/hypothesis
        Bins = linspace(0, 1, nBin);
        mh = NaN(nBin, 3);
        for iDim = 1:3, mh(:,iDim) = hist(points(:,iDim), Bins); end
        marghist{iSub} = mh./sum(mh,1);
    end
    
    % Average over subjects
    avgtrihist{iMap+1}  = mean(cell2mat(reshape(trihist,  [1,1,nSub])), 3);
    avgmarghist{iMap+1} = mean(cell2mat(reshape(marghist, [1,1,nSub])), 3);
end

% Also add a map of the difference between probabilistic and deterministic
% regularities
avgtrihist{1} = avgtrihist{3} - avgtrihist{2};

% Display triangular and marginal histograms
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Prepare the window
figure('Position', [461 253 1000 650]);

% For each density map
for iMap = 1:nMap+1
    sp = subplot(2,3,iMap);
    
    % Display the density maps of finger's position
    imagesc(xgrid, ygrid, avgtrihist{iMap}, 'AlphaData', mask); hold('on');
    contour(xgrid, ygrid, avgtrihist{iMap}, 11, 'k-', 'LineWidth', 1/2);
    image(xgrid, ygrid, repmat(~mask, [1,1,3]), 'AlphaData', ~mask);
    
    % Draw marginal histograms on the limits of the triangle
	if iMap > 1, Emergence_PlotMargHist(avgmarghist{iMap}, nBin); end
    
    % Overlap the grid corresponding to the resoution of the marginal
    % histrograms
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
    if iMap == 1
        colormap(sp, [flipud(Emergence_Colormap('Blues')); ...
                             Emergence_Colormap('Reds')]);
        caxis([-abs(max(caxis)), abs(max(caxis))]);
    elseif iMap > 1
        colormap(sp, Emergence_Colormap('Rainbow', 1e5));
        caxis([0,1/2*max(caxis)]);
        set(gca,'ColorScale','log');
        cbr.Label.String = {'Log-density'};
        cbr.Ruler.Scale = 'log';
        cbr.Ruler.MinorTick = 'on';
    end
end

% Use the same zoom level for all the maps
ScaleAxis('x'); ScaleAxis('y'); ScaleAxis('c', 1:nMap);

% Save the figure
if isfield(D{1}, 'Seq'), save2pdf(fullfile(ftapath, 'figs', 'F_D_DensityMapsS.pdf'));
else, save2pdf(fullfile(ftapath, 'figs', 'F_D_DensityMapsIO.pdf'));
end
