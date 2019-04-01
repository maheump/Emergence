% This script setups the MATLAB environment and figures' preferences
% required to run Emergence_FTA_* analysis scripts.
% 
% Copyright (c) 2018 Maxime Maheu

% Display info
fprintf('Loading preprocessed data... ');

%% ADD FUNCTIONS TO THE MATLAB PATH
%  ================================

% Add functions to the MATLAB path
scriptpath = mfilename('fullpath');
ind = strfind(scriptpath,'Emergence');
folderpath = scriptpath(1:ind(end-1)+8);
addpath(genpath(folderpath));

% Make sure there is a folder to receive figure files
ftapath = fullfile(folderpath, 'Finger tracking analyses');
figdir = fullfile(ftapath, 'figs');
if ~exist(figdir, 'dir'), mkdir(figdir); end

%% DEFINE FIGURE PROPERTIES AND LOAD DATA
%  ======================================

% Set default properties for the figures
Emergence_DefaultFigureProperties;

% Load the data
load('Emergence_Behaviour_GroupData.mat');

%% EXCLUDE BAD SUBJECTS
%  ====================

% Define bad subjects
badsubs = [1, 7, 13, 14, 19];
goodsubs = setdiff(1:nSub, badsubs);

% Remove data from bad subjects
if ~exist('rmbadsub', 'var'), rmbadsub = true; end
if rmbadsub
    G = G(:,goodsubs);
    S = S(goodsubs);
    IO = IO(:,goodsubs);
    rmbadsub = true;
end

% Get the number of subjects and the number of conditions
[nSeq, nSub] = size(G);

% By default, look at subjects' data
D = G;

% Get the length of each sequence
N = numel(G{1}.Seq);

%% DEFINE SEQUENCE SELECTION FILTERS
%  =================================
%  N.B. We restrict some analyses to some of the sequences (e.g. sequences
%  accurately classified by subjects, with a regular detection point from
%  the subjects and the ideal observer). Define the matrix selection
%  filters here.

% Prepare the output variable
filter = cellfun(@(x) NaN(numel(x),nSub), cidx, 'UniformOutput', 0);

% For each type of sequence
for iHyp = 1:3
    
    % For sequences with a regularity
    if iHyp < 3
        
        % CRITERION 1: keep only sequences that were correctly labelled by
        % subjects at the post-sequence questions
        subrp = cellfun(@(x) x.Questions(2) == iHyp,  G(cidx{iHyp},:));
        offlinedetecmask = subrp;
        
        % CRITERION 2: keep only sequences for which we have a regular
        % detection point for the subjects AND the ideal observer
        cp = cellfun(@(x) x.Jump+1/2, G(cidx{iHyp},:), 'UniformOutput', 0);
        subdp = cellfun(@(x,c) Emergence_FindDetecPoint(x.BarycCoord(c:end,iHyp)),  G(cidx{iHyp},:), cp);
        iodp  = cellfun(@(x,c) Emergence_FindDetecPoint(x.BarycCoord(c:end,iHyp)), IO(cidx{iHyp},:), cp);
        onlinedetecmask = ~isnan(subdp) & ~isnan(iodp);
        
    % For sequences without regularities
    else
        
        % CRITERION: keep only sequences that were correctly labelled by
        % subjects at the post-sequence questions
        subrp = cellfun(@(x) isnan(x.Questions(2)), G(cidx{iHyp},:));
        offlinedetecmask = subrp;
        
        % No other criterion: we don't look for detection points in
        % fully-stochastic sequences
        onlinedetecmask = false(size(G(cidx{iHyp},:)));
    end
    
    % Apply those 2 criterions to select the sequences
    filter{iHyp}(~offlinedetecmask & ~onlinedetecmask) = 0; % all sequences
    filter{iHyp}( offlinedetecmask & ~onlinedetecmask) = 1; % accurately classified sequences
    filter{iHyp}(~offlinedetecmask &  onlinedetecmask) = 2; % regular detection point
    filter{iHyp}( offlinedetecmask &  onlinedetecmask) = 3; % both conditions
end

%% RETURN TO THE COMMAND WINDOW
%  ============================

fprintf('Done! ');
if      rmbadsub, fprintf('Bad subjects are excluded.\n');
elseif ~rmbadsub, fprintf('Bad subjects are NOT excluded.\n');
end
