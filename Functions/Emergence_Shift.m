function [ xp, yp ] = Emergence_Shift( x, y, d )
% EMERGENCE_SHIFT shifts two arrays/matrices, one with respect to the
% other.
%   - "x": the first  matrix of Nobs*Ndim observations (e.g. subject)
%   - "y": the second matrix of Nobs*Ndim observations (e.g. model)
%   - "d": a scalar specifying the amount of shift (in terms of number of
%       observations) of one with respect to the other.
% 
% Copyright (c) 2020 Maxime Maheu

% Number of observations
n = size(x,1);

% If the subject is slower than the ideal observer
if d >= 0
    xp = x(d+1:n,:); % to be explained: subject's trajectory
    yp = y(1:n-d,:); % explaining variable: IO's belief
    
% If the subject is faster than the ideal observer
elseif d < 0
    xp = x(1:n+d,:);  % to be explained: subject's trajectory
    yp = y(-d+1:n,:); % explaining variable: IO's belief
end

end