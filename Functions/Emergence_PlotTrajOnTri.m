function lgd = Emergence_PlotTrajOnTri( pMgY, J, tricol, ms, markers, lgd, fs )
% EMERGENCE_PLOTTRAJONTRI displays a given trajectory within an equilateral
% triangle with the following properties:
% 
%                       width = 1
%                       |-------|
%                       P ----- D   -----|
%                        \     /         | height =
%                         \   /          | sqrt(3)/2
%                           S       -----|
% 
%   - "pMgY": a Nx3 matrix specifying the beliefs of each of the 3 possible
%       hypotheses as a function of observations (default: NaN(1,3)).
%   - "J": a scalar specifying the position if the jump (default: NaN).
%   - "tricol": a 3x3 matrix specifiying colors to use (RGB in columns) for
%       each triangle's vertex (rows).
%   - "ms": a scalar specifying the size of dots to overlap on the sequence
%       that are used to see the amount of update between successive finger
%       positions (default = eps, i.e. the trajectory is simply drawn as a
%       line).
%   - "markers": a cell array specifying the markers to use for the first
%       observation of the sequence, the change point and the last
%       observation of the sequence (default = {'v', 'p', '^'}).
%   - "lgd": a logical boolean specifying wether to display the legend
%       or not (default = false).
%   - "fs": a scalar specifyong the size of the fonts (default = default
%       of "Emergence_PlotTriInfo").
% 
% Copyright (c) 2018 Maxime Maheu

% Define properties of the triangle
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Coordinates of the triangles limits (cartesian coordinates)
tricc = [0,   sqrt(3)/2; ... % top left (P)
         1,   sqrt(3)/2; ... % top right (D)
         1/2, 0        ];    % bottom (R)

% Default colors to use
if nargin < 3
    tricol = [066 146 198; ...     % top left (blue)
              239 059 033; ...     % top right (red)
              065 171 093] ./ 255; % bottom (green)
end

% Get barycentric coordinates of the trajectory
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Convert barycentric coordinates to cartesian coordinates
if nargin < 1, pMgY = NaN(1,3); end
cartcoor = pMgY*tricc;
cartcoor = cartcoor(any(~isnan(cartcoor), 2),:);

% Plot the triangle
% ~~~~~~~~~~~~~~~~~

% Clear the axes
cla;

% Plot the triangle's background
patch('Faces', 1:3, 'Vertices', tricc, 'FaceVertexCData', tricol, ...
      'FaceColor', 'interp', 'EdgeColor', 'None'); hold('on');

% Display information about the triangle
if nargin < 7, fs = []; end
Emergence_PlotTriInfo(tricc, tricol, fs);

% Get information about important moments of the sequence
lout = NaN(1,4);
lnames = {'Trajectory', 'Starting point', 'Ending point', 'Change point'};
if nargin < 4 || isempty(markers), markers = {'v', 'p', '^'}; end
if nargin < 4, ms = eps; end

% Plot the trajectory
if nargin >= 1 && ~isempty(cartcoor)
    lout(1) = plot(cartcoor(:,1), cartcoor(:,2), 'k.-', 'MarkerSize', ms, ...
        'MarkerEdgeColor', 'k', 'LineWidth', 1);
    lout(2) = plot(cartcoor(1,1), cartcoor(1,2), markers{1}, 'MarkerSize', 10, ...
        'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', 'LineWidth', 1);
    lout(3) = plot(cartcoor(end,1), cartcoor(end,2), markers{3}, 'MarkerSize', 10, ...
        'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', 'LineWidth', 1);
end

% Plot the position of the change point
if nargin < 2, J = NaN; end
if ~isnan(J)
    lout(4) = plot(cartcoor(J,1), cartcoor(J,2), markers{2}, 'MarkerSize', 15, ...
        'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', 'LineWidth', 1);
else, lout = lout(1:3); lnames = lnames(1:3);
end

% Customize axes
axis('equal'); axis('off');
axis([min(tricc(:,1)), max(tricc(:,1)), min(tricc(:,2)), max(tricc(:,2))]);

% Display the legend
if nargin < 6, lgd = false; end
if lgd, lgd = legend(lout, lnames, 'Location', 'SouthOutside'); end

end