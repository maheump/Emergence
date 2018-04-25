% Convenient tiny script to setup the MATLAB environment and figures'
% preferences required to run Emergence_FAT_* analysis scripts.
% 
% Copyright (c) 2018 Maxime Maheu

%% LOAD DATA
%  =========

% Locate data
homedir = mfilename('fullpath');
addpath(genpath(homedir));

% Load the data
load('Emergence_Behaviour_GroupData.mat');

% Define bad subjects
badsubs = [1, 7, 13, 14, 19];
goodsubs = setdiff(1:nSub, badsubs);

% Remove data from bad subjects
G = G(:,goodsubs);
S = S(goodsubs);
IO = IO(:,goodsubs);

% By default, analyse subjects' data
D = G;

% Get the number of subjects and the number of conditions
[nCond, nSub] = size(G);

%% CUSTOMIZE DEFAULT FIGURES
%  =========================

% Set default properties for the figure
set(groot, 'DefaultAxesLineWidth', 1); % axes
set(groot, 'DefaultBarLineWidth', 1); % bars
set(groot, 'DefaultColorbarLineWidth', 1); % colorbar
set(groot, 'DefaultLineLineWidth', 1); % line plots
set(groot, 'DefaultFigureColor', ones(1,3)); % Remove the uggly gray surrounding
set(groot, 'DefaultAxesColor', 'None'); % transparent background
set(groot, 'DefaultTextColor', zeros(1,3)); % black text
set(groot, 'DefaultAxesLineWidth', 1); % axes
set(groot, 'DefaultBarLineWidth', 1); % bars
set(groot, 'DefaultColorbarLineWidth', 1); % colorbar
set(groot, 'DefaultLineLineWidth', 1); % line plots
set(groot, 'DefaultAxesFontName', 'Helvetica'); % Helvetica font
set(groot, 'DefaultTextFontName', 'Helvetica'); % Helvetica font
set(groot, 'DefaultAxesFontUnits', 'Points'); % unit
set(groot, 'DefaultAxesFontSize', 10); % default font-size
set(groot, 'DefaultTextVerticalAlignment', 'Middle');
set(groot, 'DefaultTextHorizontalAlignment', 'Center');
set(groot, 'DefaultAxesBox', 'off'); % no upper and right axes
set(groot, 'DefaultAxesLayer', 'Top'); % axes on top of the plot
set(groot, 'DefaultAxesTickDir', 'In'); % ticks oriented outside the plot
set(groot, 'DefaultAxesTickLength', repmat(0.01, [1 2])); % length of ticks
set(groot, 'DefaultAxesTickDirMode', 'Auto'); % otherwise the previous commands do not work

% Some useful colour
g = repmat(0.7,1,3); % grey
