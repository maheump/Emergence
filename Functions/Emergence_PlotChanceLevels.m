function Emergence_PlotChanceLevels( tricc )
% EMERGENCE_PLOTCHANCELEVELS displays some lines on the triangle that help
% to understand it better.
%   - "tricc": a 3x2 matrix specifying the x/y coordinates (columns) of
%       each vertex of the triangle. By default the width is 1, the height
%       sqrt(3)/2 and it is drawn at (0,0).
% Copyright (c) 2018 Maxime Maheu

% Coordinates of the triangles limits (cartesian coordinates)
if nargin < 1 || isempty(tricc)
    tricc = [0,   sqrt(3)/2; ... % top left (P)
             1,   sqrt(3)/2; ... % top right (D)
             1/2, 0        ];    % bottom (R)
end

% Get useful measures
th = tricc(1,2) - tricc(3,2); % height (pixels)
tw = tricc(2,1) - tricc(1,1); % width (pixels)
tx = tw/2 + tricc(1,1); % x coordinate of the triangle's center (pixels)
ty = tricc(1,2) - sqrt(3)/6 * tw; % y coordinate of the triangle's center (pixels)

% Draw help lines
g = repmat(0.6,1,3); % grey
plot([tx,tx],      [ty,tricc(1,2)],      '-', 'Color', g, 'LineWidth', 1); hold('on');
plot([tx,tx-tw/4], [ty,tricc(3,2)+th/2], '-', 'Color', g, 'LineWidth', 1);
plot([tx,tx+tw/4], [ty,tricc(3,2)+th/2], '-', 'Color', g, 'LineWidth', 1);

end
