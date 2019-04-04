function [ trihist, barc, xgrid, ygrid, mask ] = Emergence_TriHist( barc, gridprec, smooth, interp, tricc )
% EMERGENCE_TRIHIST returns a 2D triangular histogram of a series of
% barycentricpositions.
%   - "barc": a Nx3 matrix or a cell-array containing Nx3 matrices, where N
%       is the number of observations, from which to build the histogram.
%   - "gridprec": a scalar specifying the size in both dimensions of each
%       square parcels of the histrogram.
%   - "smooth": a scalar specifying the width of the 2D gaussian filter to
%       apply.
%   - "interp": a scalar specifying how much the trajectories (within each
%       cell of "barc) should be interpolated (1 is no interpolation,
%       smaller than 1 is interpolation).
%   - "tric": a 3x2 matrix specifying the cartesian coordintes of the
%       triangle vertices.
% 
% Copyright (c) 2018 Maxime Maheu

% Initialization
% ~~~~~~~~~~~~~~

% Defaut length 1 edge triangle
if nargin < 5, tricc = [0, sqrt(3)/2; 1, sqrt(3)/2; 1/2, 0]; end

% Make sure it is a column vector and it is free from empty cells
if ~iscell(barc), barc = {barc}; end
barc = barc(~cellfun(@isempty, barc));

% Interpolate the trajectories
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Define the degree of interpolation that we apply on the trajectory
if nargin < 4, interp = 1/3; end % must be <= 1, 1 is no interpolation

% Interpolate the trajectories
barc = cellfun(@(p) interp1(1:size(p,1), p, 1:interp:size(p,1)), barc, 'UniformOutput', 0);
barc = cell2mat(barc);

% Make a 2D histogram
% ~~~~~~~~~~~~~~~~~~~

% Convert to cartesian coordinates
carc = barc*tricc;

% Define the precision of the mesh
if nargin < 2, gridprec = 0.01; end
xgrid = tricc(1,1):gridprec:tricc(2,1);
ygrid = tricc(3,2):gridprec:tricc(1,2);

% Compute the density map
density = hist3(carc, 'Edges', {xgrid, ygrid})';

% Smooth the histogram
% ~~~~~~~~~~~~~~~~~~~~

% Options for the gaussian smooth
if nargin < 3, smooth = 2; end
filtWidth = 2*ceil(2*smooth)+1;
imageFilter = fspecial('Gaussian', filtWidth, smooth);

% Smooth the density map
density = convn(density, imageFilter, 'same');

% Scale the histogram
% ~~~~~~~~~~~~~~~~~~~

% Convert it to log scale
density = log(density+1);

% Normalize it
trihist = density ./ sum(density(:));
    
% Create a mask
% ~~~~~~~~~~~~~

% Create a mask that specifies which cells of the matrix are outside the
% triangle
if nargout > 2
    ng = 1/gridprec+1;
    T = tril(ones(ng));
    T = [flip(T,2),T];
    mask = round(imresize(T, [round(ng/2*tan(pi/3)), ng]));
end

end