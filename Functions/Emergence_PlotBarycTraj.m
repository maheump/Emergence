function f = Emergence_PlotBarycTraj( pMgY, tricol, x )
% EMERGENCE_PLOTBARYCTRAJ displays beliefs in each hypothesis in a
% cumulative manner.
%   - "pMgY": a Nx3 matrix specifying the beliefs of each of the 3 possible
%       hypotheses as a function of observations.
%   - "tricol": a 3x3 matrix specifying the RGB (columns) colors to use for
%       each hypothesis.
%   - "x": a 1xN vector specifying the positions of observations.
% Copyright (c) 2018 Maxime Maheu

% By default, display cumulative probabilities starting at x = 1
nObs = size(pMgY, 1);
if nargin < 3, x = 1:nObs; end

% Colors to use 
if nargin < 2 || isempty(tricol)
    tricol = [cbrewer2('Blues', 1); cbrewer2('Reds', 1); cbrewer2('Greens', 1)];
end

% Compute cumulative probabilities
cum_pMs = cumsum(pMgY,2);

% Display cumulative probabilities as layers
f = NaN(1,3);
f(1) = fill([x,fliplr(x)], [cum_pMs(:,2); cum_pMs(:,3)]', ...
    'k', 'FaceColor', tricol(3,:), 'LineWidth', 1); hold('on');
f(2) = fill([x,fliplr(x)], [cum_pMs(:,1); flipud(cum_pMs(:,2))]', ...
    'k', 'FaceColor', tricol(2,:), 'LineWidth', 1);
f(3) = fill([x,fliplr(x)], [cum_pMs(:,1); zeros(nObs,1)]', ...
    'k', 'FaceColor', tricol(1,:), 'LineWidth', 1);
        
end