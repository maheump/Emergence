function trajlocked = Emergence_LockOnPoint( traj, lockobs, win )
% EMERGENCE_LOCKONPOINT locks a given trajectory on a particular
% observation number and return the trajectory in a window around that
% locked observation.
%   - "traj": a NxK matrix where N is the number of observations.
%   - "lockobs": a scalar specifying which observation number on which to
%       lock the trajectory.
%   - "win": an array specifying either the limit of the window (e.g.
%       [-20,20]) or the index of the window (e.g. -20:20).
% 
% Copyright (c) 2018 Maxime Maheu

% By default, use all the data
if nargin < 3, win = [0, max(size(traj))-lockobs]; end

% If only the limits of the window are specified
if numel(win) == 2, win = win(1):win(2); end

% Get the length of the window
winlen = numel(win);

% Number of observations
[n,k] = size(traj);

% Get the observations indices to sample from
obsidx = lockobs + win;

% Make sure the observation indices stay within the limits
nanidx = obsidx >= 1 & obsidx <= n;
obsidx = obsidx(nanidx);

% Sample from the trajectory
trajlocked = NaN(winlen,k);
trajlocked(nanidx,:) = traj(obsidx,:);

end