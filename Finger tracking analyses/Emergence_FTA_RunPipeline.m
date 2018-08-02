% Run the entire analysis pipeline for the finger tracking experiment.
%
% Copyright (c) 2018 Maxime Maheu

%% PREPROCESS AND LOAD DATA
%  ========================

% Clear the place
clear;
close('all');

% Preprocess behavioural data and run the ideal observer on the sequences
% that were presented to the subjects
Emergence_FTA_PreprocData;

%% PERFORM SUBJECT-LEVEL ANALYSES
%  ==============================

% Load data from all the subjects
rmbadsub = false;
Emergence_FTA_LoadData;

% Get trajectories of all the subjects in all the conditions
Emergence_FTA_SubjectAnalysis;

% Display example trajectories (Figure 1)
Emergence_FTA_ExampleTrajectories;
Emergence_FTA_TrajMovie;

%% PERFORM GROUP-LEVEL ANALYSES
%  ============================

% Load data and exclude bad subjects
rmbadsub = true;
Emergence_FTA_LoadData;

% Analyse data only from the subjects
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Emergence_FTA_FalseAlarms;

% Analyse data only from the ideal observer
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Emergence_FTA_ChangePoint1;

% Analyse data separately from both the subjects and from the ideal observer
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for iLoop = 1:2
    if     iLoop == 1, D = G;
    elseif iLoop == 2, D = IO;
    end

    Emergence_FTA_Detection1;
    Emergence_FTA_Detection2;
    Emergence_FTA_DistributionOfBeliefs;
    Emergence_FTA_Dynamics1;
    Emergence_FTA_HypothesisWeighting1;
    Emergence_FTA_ChangePoint2;
end

% Correlation between subjects and ideal observer
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Emergence_FTA_Dynamics2;
Emergence_FTA_HypothesisWeighting2;
Emergence_FTA_QuantitativeDeviations1;
Emergence_FTA_QuantitativeDeviations2;