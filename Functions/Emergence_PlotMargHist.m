function [ marghist, bars ] = Emergence_PlotMargHist( marghist, tricc, tricol )
% EMERGENCE_PLOTMARGHIST draws the 3 marginal histograms along the limits
% of the triangle. It assumes an equilateral triangle oriented with one tip
% to the bottom and the two remaing tips to the top, like that: ^
%   - "marghist": a Nx3 specifying the histogram to plot.
%   - "tricc": a 3x2 matrix specifying the x/y coordinates (columns) of
%       each vertex of the triangle. By default the width is 1, the height
%       sqrt(3)/2 and it is drawn at (0,0).
%   - "tricol": a 3x3 matrix specifying the RGB (columns) of each vertex of
%       the triangle.
% 
% Copyright (c) 2018 Maxime Maheu

% Fill in the inputs
% ~~~~~~~~~~~~~~~~~~

% Coordinates of the triangles limits (cartesian coordinates)
if nargin < 3 || isempty(tricc)
    tricc = [0,   sqrt(3)/2; ... % top left (P)
             1,   sqrt(3)/2; ... % top right (D)
             1/2, 0        ];    % bottom (R)
end

% Colors to use 
if nargin < 4 || isempty(tricol)
    tricol = [049, 130, 189; 222, 045, 038; 049, 163, 084] ./ 255;
end

% Prepare the variable that will entail bars' properties (colors, ...)
nBin = size(marghist,1);
bars = NaN(nBin,3);

% Prepare a variable
axlim = NaN(3,2);

% For each of the 3 side of the triangle
for iDim = 1:3
    
    % Get the histogram for the current side of the triangle
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    bins = marghist(:,iDim)';
    if iDim == 3, bins = fliplr(bins); end
    
    % Segment the limits of the triangle into equal parts
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    
    % Get x-coordinates of each segment
    if     iDim == 1, x = linspace(tricc(1,1), tricc(3,1), nBin+1); x = fliplr(x);
    elseif iDim == 2, x = linspace(tricc(1,1), tricc(2,1), nBin+1);
    elseif iDim == 3, x = linspace(tricc(3,1), tricc(2,1), nBin+1);
    end
    
    % Get y-coordinates of each segment
    if iDim == 2, y = repmat(sqrt(3)/2, 1, nBin+1);
    else,         y = linspace(0, sqrt(3)/2, nBin+1);
    end
    
    % Get coordinates of the bars
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %  2__3      __   => Example marginal histogram with nBin = 3
    %  |  |__   |  |  => with the location of points #1, #2, #3 and #4
    %  |  |  |__|  |  => whose coordinates are derived below
    % =1==4==|==|==|= => The "=" represent the limit of the triangle
    
    
    % Coordinates of point #1 (on the limit of the triangle)
    x1 = x(1:nBin);
    y1 = y(1:nBin);
    
    % For the horizontal dimension (top one)
    %   => regular (horizontal) histogram (as drawn above)
    if iDim == 2
        
        % Coordinates of point #2
        x2 = x1;
        y2 = y1 + bins;
        
        % Coordinates of point #3
        x3 = x(2:nBin+1);
        y3 = y1 + bins;
        
    % For the two tilted dimensions (to the left and to the right)
    %   => tilted histogram obtained thanks to trigonometry
    else
        
        % Coordinates of point #2
        o = bins*sind(30); % opposite side of the right triangle
        a = bins*cosd(30); % adjacent side of the right triangle
        if     iDim == 1, x2 = x1-a; % to the left
        elseif iDim == 3, x2 = x1+a; % to the right
        end
        y2 = y1-o; % to the top
        
        % Coordinates of point #3
        o = (1/nBin)*sind(60); % opposite side of the right triangle
        a = (1/nBin)*cosd(60); % adjacent side of the right triangle
        if     iDim == 1, x3 = x2-a; % to the left
        elseif iDim == 3, x3 = x2+a; % to the right
        end
        y3 = y2+o; % to the bottom
    end
    
    % Coordinates of point #4 (on the limit of the triangle)
    x4 = x(2:nBin+1);
    y4 = y(2:nBin+1);
    
    % Draw the histogram
    % ~~~~~~~~~~~~~~~~~~
    
    % Draw all the bars at the same time
    xcoord = [x1; x2; x3; x4]; % each row is one of the vertex of the bars
    ycoord = [y1; y2; y3; y4]; % and each column is a bin
    bars(:,iDim) = fill(xcoord, ycoord, 'k', 'FaceColor', ...
        tricol(iDim,:), 'LineWidth', 1); hold('on');
    
    % Identify x and y, min and max, limits of the histogram
    coord = [xcoord(:); ycoord(:)];
    axlim(iDim,:) = [min(coord), max(coord)];
end

% Tightly scale the axes such that (i) the highest bar of the histogram
% coincides with the limit of the axes and (ii) the 
axis(repmat([min(axlim(:,1)), max(axlim(:,2))], 1, 2));
axis('square');

end