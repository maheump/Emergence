function [ traj, corr ] = Emergence_ProjOnTri( traj, tripxl )
% EMERGENCE_PROJONTRI projects finger's positions onto the triangle borders
% (based on minimal, orthogonal, Euclidean distance) when those
% positions are (rarely) outside the triangle (mostly because of finger's
% imprecision).
%   - traj: a Nx2 matrix specifying the cartesian coordinates of the
%       finger's positions over time.
%   - "tripxl": a binary matrix specifying the screen's pixels in which the
%       triangle was displayed.
%
% Copyright (c) 2020 Maxime Maheu

% Check wether the position is not within the triangle
tocorrect = find(~ismember(traj, tripxl, 'rows'));
nc = numel(tocorrect); % number of corrections

% Prepare the output variable
corr = NaN(nc,2);
corr(:,1) = tocorrect;

% Loop over finger's positions to correct
for i = 1:nc
    I = tocorrect(i);
    
    % Compute the Euclidean distance between its position and all the
    % pixels that constitute the triangle on the 2D touchscreen
    EucDist = sqrt(sum((traj(I,:) - tripxl).^2,2));
    
    % Project it on the closest point in the triangle (i.e. on the
    % borders)
    [corr(i,2),idx] = min(EucDist);
    traj(I,:) = tripxl(idx,:);
end

end
