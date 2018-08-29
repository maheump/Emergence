% Setup default options and add required functions and scripts to the
% MATLAB path.
%
% Copyright (c) 2018 Maxime Maheu

% Clear the place
clear;
close('all');

% Locate
scriptpath = mfilename('fullpath');
folderpath = scriptpath(1:strfind(scriptpath,'Emergence')+8);

% Unzip folder containing various general purpose functions
unzip(fullfile(folderpath, 'Functions', 'dependencies.zip'), ...
      fullfile(folderpath, 'Functions'));

% Add scripts and functions to the MATLAB path
addpath(genpath(folderpath));

% Add the VBA toolbox to the path
try
    VBA_Setup;
catch
    warning('The VBA toolbox is not available. Some scripts will not run');
end

% Set default properties for the figures
Emergence_DefaultFigureProperties;
