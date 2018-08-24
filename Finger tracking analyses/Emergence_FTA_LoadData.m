% This script setups the MATLAB environment and figures' preferences
% required to run Emergence_FTA_* analysis scripts.
% 
% Copyright (c) 2018 Maxime Maheu

%% LOAD DATA
%  =========

% Add functions to the MATLAB path
scriptpath = mfilename('fullpath');
folderpath = scriptpath(1:strfind(scriptpath,'Emergence')+8);
addpath(genpath(folderpath));

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
end

% Get the number of subjects and the number of conditions
[nCond, nSub] = size(G);

% By default, look at subjects' data
D = G;

% Get the length of each sequence
N = numel(G{1}.Seq);
