% Setup default options and add required functions and scripts to the
% MATLAB path.
%
% Copyright (c) 2018 Maxime Maheu

% Welcome message
fprintf('\n');
fprintf('#######################################################\n');
fprintf('# This is the setup routine for the Emergence toolbox #\n');
fprintf('#######################################################\n\n');

% Info message
fprintf(['Thank you for having downloaded our toolbox.\nWe will now go ', ...
    'through a series of setup steps that ensure that the scripts ', ...
    'featured in that toolbox can run.\n\n']);

% Check the version of MATLAB that is running
fprintf('Checking MATLAB version... ');
ismatlabold = verLessThan('matlab', '9.5');
fprintf('Done!\n');
if ismatlabold
    fprintf(['Your version of MATLAB is older than the one on which ', ...
        'the toolbox was developed. If it is too old, some functions ', ...
        'might not run. Consider to update it if it is the case.\n\n']);
elseif ~ismatlabold
    fprintf(['Your version of MATLAB is at least as recent than the one ', ...
        'on which the toolbox was developed. All scripts should work ', ...
        'properly.\n\n']);
end

% Locate the toolbox directory
fprintf('Locating the current directory ...');
folderpath = fileparts(mfilename('fullpath'));
if isempty(folderpath), folderpath = pwd; end
fprintf('Done!\n\n');

% Check if the dependencies have are available
fprintf('Checking that the dependencies are available... ');
depenfolder = fullfile(folderpath, 'Dependencies');
availdep = exist('cbrewer2.m', 'file') & ... % toolbox with colormaps
           exist('pat2str.m', 'file') & ...  % toolbox with home-made functions
           exist('VBA_setup.m', 'file');     % toolbox for variational inference
fprintf('Done! \n');
if availdep
    fprintf('All dependencies are available.\n\n');
elseif ~availdep
    fprintf('It appears that at least one dependency is missing.\n');
    instaldep = logical(input('Do you wish to install the dependencies? (0 for no, 1 for yes) '));
    if instaldep
        fprintf('Installing dependencies.\n');
        system('git submodule init', '-echo');
        system('git submodule update', '-echo');
        fprintf('Done!\n\n');
    elseif ~instaldep
        fprintf('Note that some, if not all, scripts will not run because of missing dependencies.\n');
    end
end

% Add scripts and functions to the MATLAB path
fprintf('Adding the toolbox functions to the MATLAB path... ');
addpath(genpath(folderpath));
fprintf('Done!\n\n');

% Set default properties for the figures
fprintf('Changing the default figures'' options... ');
Emergence_DefaultFigureProperties;
fprintf('Done!\n\n');

% Display some warning messages
fprintf('!!! Nothing that has been done here will survive the restart of MATLAB !!!\n');
fprintf('!!! Do not forget to rerun this script if you want to use the toolbox again !!!\n\n');

% Try to run a script in order to check that everything is working well
runtestscript = input(['Do you wish to try running an example script in ', ...
    'order to test that everything works fine? (0 for no, 1 for yes) ']);
fprintf('\n');
if runtestscript
    Emergence_IO_ToyExampleFullIO;
    fprintf('\n\nEverything works perfectly fine!\n');
end

% Say goodbye
fprintf('You can start using the toolbox. Goodbye!\n\n');
