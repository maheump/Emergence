% Setup default options and add required functions and scripts to the
% MATLAB path.
%
% Copyright (c) 2018 Maxime Maheu

% Welcome message
fprintf('\n\n');
fprintf('#######################################################\n');
fprintf('# This is the setup routine for the Emergence toolbox #\n');
fprintf('#######################################################\n\n');

% Info message
fprintf(['Thank you for having downloaded our toolbox.\nWe will now go ', ...
    'through a series of setup steps that ensure that the scripts ', ...
    'featured in that toolbox can run.\n\n']);
pause(4); % pause for 4 seconds

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
availdep = exist(fullfile(depenfolder, 'cbrewer', 'cbrewer2.m'), 'file') & ...
           exist(fullfile(depenfolder, 'matlab', 'sequences', 'pat2str.m'), 'file');
fprintf('Done! ');
if availdep
    fprintf('All dependencies are available.\n\n');
elseif ~availdep
    fprintf('There are missing dependencies. Getting them...\n');
    system('git submodule init', '-echo');
    system('git submodule update', '-echo');
    fprintf('Done!\n\n');
end

% Add the VBA toolbox to the path
fprintf('The VBA toolbox is required for some analysis scripts to run.\n');
fprintf('Checking that the VBA toolbox is available... ');
VBAsetupscript = 'VBA_setup';
addpath(genpath(fullfile(folderpath, 'Dependencies', 'VBA-toolbox')));
VBAfolder = fileparts(which(VBAsetupscript));
if ~isempty(VBAfolder)
    fprintf('Done! The VBA toolbox is located in %s.\n', VBAfolder);
    try
        infos = VBA_version;
        fprintf('We assume it has been correctly installed.\n');
    catch
        fprintf(['It has not been correctly installed. Please run the %s ', ...
            'command from the root of %s after.\n'], VBAsetupscript, VBAfolder);
    end
    fprintf('Also make sure you have the version 1.9.1 of that toolbox ');
    fprintf('(if it has been installed by the SETUP script, it is the case).\n\n')
elseif isempty(VBAfolder)
    fprintf('The VBA toolbox has not been found.\n');
    install = logical(input(sprintf(['Do you want us to download and ', ...
        'install the VBA toolbox for you? (0 for no, 1 for yes) '])));
    if install
        fprintf('Installing the VBA toolbox from Github...\n');
        cd('Dependencies'); % go to the dependencies directory
        system('git clone --branch v1.9.1 https://github.com/MBB-team/VBA-toolbox');
        cd('..'); % come back to the root directory
        fprintf('The VBA toolbox has been installed.\n\n');
    elseif ~install
        fprintf(['The scripts "Emergence_FTA_QuantitativeDeviations1" ', ...
            'and "Emergence_FTA_Dynamics1" will not run properly without ', ...
            'that toolbox. Everything else should work properly.\n\n']);
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
fprintf('!!! Do not forget to rerun this script if you want to use the toolbox again!!!\n\n');

% Try to run a script in order to check that everything is working well
runtestscript = input(['Do you wish to try running an example script in ', ...
    'order to test that everythink works fine? (0 for no, 1 for yes) ']);
fprintf('\n\n');
if runtestscript
    Emergence_IO_ToyExampleFullIO;
    fprintf('\n\nEverything works perfectly fine!\n\n');
end

% Say goodbye
fprintf('You can start using the toolbox. Goodbye!\n\n');
