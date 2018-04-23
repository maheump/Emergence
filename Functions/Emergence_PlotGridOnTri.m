function [ lines ] = Emergence_PlotGridOnTri( nGrid, dims, col, tricc )
% EMERGENCE_PLOTGRIDONTRI displays a grid on the triangle.
%   - "nGrid": a scalar number specifying the number of grid lines to
%       display.
%   - "dims": a scalar or array specifying along which of the 3 dimensions
%       to plot the grid lines.
%   - "col": and RGB array specifying which color to use for the lines.
%   - "tricc": a 3x2 matrix specifying the x/y coordinates (columns) of
%       each vertex of the triangle. By default the width is 1, the height
%       sqrt(3)/2 and it is drawn at (0,0).
%
% Copyright (c) 2018 Maxime Maheu

% Default number of sectors
if nargin < 1 || isempty(nGrid), nGrid = 10; end

% By default, display the grid over all the 3 dimensions of the triangle
if nargin < 2 || isempty(dims), dims = 1:3; end

% Use MATLAB standard for grids' color
if nargin < 3 || isempty(col), col = repmat(0.15,1,3); end
    
% Coordinates of the triangles limits (cartesian coordinates)
if nargin < 4 || isempty(tricc)
    tricc = [0,   sqrt(3)/2; ... % top left (P)
             1,   sqrt(3)/2; ... % top right (D)
             1/2, 0        ];    % bottom (R)
end

% The lines forming the grid are uniformly spaced
grids = linspace(1/nGrid, 1, nGrid);

% Prepare output variables for line objects
lines = NaN(nGrid-1, 3);

% Since the boundary is a straight line, we just need two points to draw it
% N.B. Use e.g. 100 if you want to convince you that it is a straight line!
nDots = 2;

% For each side of the triangle
for iDim = dims
	
    % For each boundary to display
    for iGrid = 1:nGrid
        
        % Specify that the first dimension should be equal to the currently
        % considered threshold
        x1 = grids(iGrid);
        
        % Deduce what degrees of freedom we have left for the other
        % dimensions
        rem = 1-x1;
        
        % Deduce the barycentric coordinates the other two dimensions can
        % take
        x1 = repmat(x1, [nDots,1]);
        x2 = linspace(0, rem, nDots)';
        x3 = flipud(x2);
        
        % Combine the three barycentric coordinates
        if     iDim == 1, bc = [x1, x2, x3];
        elseif iDim == 2, bc = [x2, x1, x3];
        elseif iDim == 3, bc = [x2, x3, x1];
        end
        
        % Convert barycentric to cartesian coordinates
        cc = bc*tricc;
        
        % Display the boundaries the same way MATLAB does for the grid
        lines(iGrid,iDim) = patchline(cc(:,1), cc(:,2), 'LineStyle', '-', ...
            'EdgeColor', col, 'LineWidth', 1/2, 'EdgeAlpha', 0.15); 
    end
end

end