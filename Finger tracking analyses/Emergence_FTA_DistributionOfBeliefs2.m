% This scripts computes histograms of finger's positions inside the
% triangular arena from the subjects/ideal observer conditioned on the
% finger's positions from the ideal observer/subjects. It allows to see if
% there is a correlation between subjects and the ideal observer. Or, in
% other terms, is there agreement between human subjects and the ideal
% observer. For instance, is it true that when subjects are at the bottom
% of the triangle, the ideal observer is also at the bottom of the
% triangle? And is the reverse also true?
% 
% Copyright (c) 2018 Maxime Maheu

% Bin the barycentric dimensions
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Define the number of bins (same number in each of the 3 dimensions)
nBin = 10; % number of bins in each barycentric dimension
nt = nBin^2; % total number of vertices

% Build the barycentric grids based on the number of bins
pgrid = linspace(0, 1, nBin+1);
pPgrid = cell2mat(arrayfun(@(x,y) repmat(x,y,1), pgrid, nBin+1:-1:1, 'UniformOutput', 0)');
pDgrid = cell2mat(arrayfun(@(x) (0:1/nBin:x)', 1-pgrid, 'UniformOutput', 0)');

% Triangle's specifications
tcn = [0, sqrt(3)/2; 1, sqrt(3)/2; 1/2, 0];

% Deduce vertices
% ~~~~~~~~~~~~~~~

% Create short vectors for the griding based on the delaunay triangulation
% of the set of points previously defined
x = (1/2) - pPgrid*cos(pi/3) + pDgrid/2;
y = (sqrt(3)/2) - ((sqrt(3)/2) - pPgrid*sin(pi/3) - pDgrid*cot(pi/6)/2);
tri = delaunay(x,y);

% Exclude small vertices 
triareas = polyarea(x(tri)', y(tri)');
triidx = find(triareas > mean(triareas));
tri = tri(triidx,:);

% Compute limits of the bins in barycentric coordinates
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Get the (2) cartesian coordinates of the (3) points (the vertices)
% defining each small triangle
rsptri = reshape(tri', [1,3,nt]);
intermvar = cat(1, x(rsptri), y(rsptri));
tricartescoord = squeeze(mat2cell(intermvar, 2, 3, ones(1,1,nt)));

% Convert cartesian coordinates of vertices to barycentric ones
tribaryccoord = cellfun(@(x) cartes2baryc(x, tcn), tricartescoord, 'UniformOutput', 0);

% Get limits in the 3 barycentric dimensions of each triangular bin
barycbinlim = cellfun(@(x) [min(x,[],2)'; max(x,[],2)'], tribaryccoord, 'UniformOutput', 0);

% Plot the results
% ~~~~~~~~~~~~~~~~

% Prepare a new figure
figure('Position', [1 805 600 300]);

% Display 2 triangular histograms: one reflecting subjects' beliefs
% conditioned on the ideal observer beliefs, and one showing the reverse
lab = {'G', 'IO'};
lab = {lab, fliplr(lab)};
for iMap = 1:2
    subplot(1,2,iMap);
    
    % Get beliefs from the subjects/IO (represented as color level based on
    % a mix of 3 colors) when the IO/subjects is in one part of the
    % triangular arena (based on bins defined earlier)
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    % Get reported subjective probabilities from the subjects/IO
    condpoints = cellfun(@(x) x.BarycCoord, eval(lab{iMap}{2}), 'UniformOutput', 0);
    condpoints = cell2mat(condpoints(:));
    
    % Bin those reported subjective probabilities according to the
    % barycentric grid defined earlier
    obsidx = cellfun(@(x) all(condpoints > x(1,:)  & ...
                              condpoints < x(2,:), 2), barycbinlim, 'UniformOutput', 0);
    
    % Count the number of observations in each bin
    nobs = cellfun(@sum, obsidx);
    facealpha = log(nobs+1); % log number of observations
    facealpha = facealpha./ max(facealpha); % normalize
    
    % Get reported subjective probabilities from the IO/subjects
    obspoints = cellfun(@(x) x.BarycCoord, eval(lab{iMap}{1}), 'UniformOutput', 0);
    
    % Select only sequences that were correctly labelled
    obspoints = cell2mat(obspoints(:));
    
    % Get and average data from the IO/subjects in those bins  
    obsavgbel = cellfun(@(x) mean(obspoints(x,:), 1), obsidx, 'UniformOutput', 0);
    
    % Convert averaged probabilities in each hypothesis into color level
    % (through weighted average)
    facecolor = cell2mat(cellfun(@(x) x*tricol, obsavgbel, 'UniformOutput', 0));
    
    % Display the triangular histrogram whose faces' color is indexed on
    % the other's averaged probabilistic beliefs
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    % Display histrogram color-coded faces
    patch('Faces', tri, 'Vertices', [x,y], ...
          'FaceColor', 'flat', 'FaceAlpha', 'flat', ...
          'FaceVertexCData', facecolor, 'FaceVertexAlphaData', facealpha, ...
          'AlphaDataMapping', 'None', ...
          'EdgeColor', g); hold('on');
    
    % Add info about the triangle
    tr = Emergence_PlotTriInfo(tcn, tricol);
    
    % Customize the axes
    set(gca, 'LineWidth', 1, 'FontSize', 15);
    axis('equal'); axis('off'); axis('square');
    
    % Add some text labels
    title(sprintf('Beliefs from %s conditioned on %s', lab{iMap}{:}));
end

% Save the figure
save2pdf(fullfile(ftapath, 'figs', 'F_D_CondMaps.pdf'));
