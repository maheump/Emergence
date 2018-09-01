% This script setups the MATLAB environment and figures' preferences
% required to run Emergence_FTA_* analysis scripts.
%
% Copyright (c) 2018 Maxime Maheu

% Display info
fprintf('Loading data... ');

% Add functions to the MATLAB path
scriptpath = mfilename('fullpath');
folderpath = scriptpath(1:strfind(scriptpath,'Emergence')+8);
addpath(genpath(folderpath));

% Set default properties for the figures
Emergence_DefaultFigureProperties;

% Load the data
load('Emergence_Behaviour_GroupData.mat');

% Define bad subjects
badsubs = [1, 7, 13, 14, 19];
goodsubs = setdiff(1:nSub, badsubs);

% Remove data from bad subjects
if ~exist('rmbadsub', 'var') || rmbadsub
    G = G(:,goodsubs);
    S = S(goodsubs);
    IO = IO(:,goodsubs);
    rmbadsub = true;
else, rmbadsub = false;
end

% Get the number of subjects and the number of conditions
[nCond, nSub] = size(G);

% By default, look at subjects' data
D = G;

% Get the length of each sequence
N = numel(G{1}.Seq);

% Display info
fprintf('Done! ');
if      rmbadsub, fprintf('Bad subjects are excluded.\n');
elseif ~rmbadsub, fprintf('Bad subjects are NOT excluded.\n');
end
