% Setup default options and add required functions and scripts to the
% MATLAB path.
%
% Copyright (c) 2018 Maxime Maheu

% Clear the place
clear;
close('all');

% Locate
folderpath = fileparts(mfilename('fullpath'));

% Unzip folder containing various general purpose functions
unzip(fullfile(folderpath, 'Functions', 'dependencies.zip'), ...
      fullfile(folderpath, 'Functions'));

% Add scripts and functions to the MATLAB path
addpath(genpath(folderpath));

% Set default properties for the figures
Emergence_DefaultFigureProperties;

% Add the VBA toolbox to the path
try
    VBAsetupscript = 'VBA_setup.m';
    cd(fileparts(which(VBAsetupscript)));
    run(VBAsetupscript);
    cd(folderpath);
catch
    warning('The VBA toolbox is not available. Some scripts will not run');
end
