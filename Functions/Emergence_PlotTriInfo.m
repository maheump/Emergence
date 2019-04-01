function tr = Emergence_PlotTriInfo( tricc, tricol, fs )
% EMERGENCE_PLOTTRIINFO displays information about the triangle.
%   - "tricc": a 3x2 matrix with x/y coordinates (columns) of each
%       triangle's vertex (rows): upper-left, upper-right and bottom one.
%   - "tricol": a 3x3 matrix specifiying colors to use (RGB in columns) for
%       each triangle's vertex (rows).
%   - "fs": a scalar specifying the size of the font.
% 
% Copyright (c) 2018 Maxime Maheu

% Coordinates of the triangles limits (cartesian coordinates)
if nargin < 1 || isempty(tricc)
    tricc = [0,   sqrt(3)/2; ... % top left (P)
             1,   sqrt(3)/2; ... % top right (D)
             1/2, 0        ];    % bottom (R)
end

% Default colors to use
if nargin < 2 || isempty(tricol)
    tricol = [049, 130, 189; 222, 045, 038; 049, 163, 084] ./ 255;
end

% Default font-size
if nargin < 3 || isempty(fs), fs = 15; end

% Draw triangle's limits
tr = fill(tricc(:,1), tricc(:,2), 'k', 'FaceColor', 'None', 'LineWidth', 2);
hold('on');

% Draw labels
proclab = {'P','D','S'};
hal = {'Right','Left','Center'};
val = {'Bottom','Bottom','Top'};
for k = 1:3
    text(tricc(k,1), tricc(k,2), ['$\mathcal{M}_{\mathrm{' proclab{k} '}}$'], ...
        'Interpreter', 'LaTeX', 'FontSize', fs, 'Color', tricol(k,:), ...
        'FontWeight', 'Bold', 'HorizontalAlignment', hal{k}, 'VerticalAlignment', val{k});
end

% Hide the axes
axis('equal');
axis('off');

end