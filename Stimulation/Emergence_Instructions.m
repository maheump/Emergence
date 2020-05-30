% This script presents the different conditions to the subjects.
% 
% Copyright (c) 2018 Maxime Maheu

%% PREPARE THE AUDITORY STIMULI
%  ============================

% Add the stimulation folder to the MATLAB path
addpath(genpath('~/Emergence/Stimulation'));

% Initialize the PTB module
InitializePsychSound(1);
sndforced = 1;

% Close eventually open buffers
PsychPortAudio('Close');

% Open an audio buffer
dv = PsychPortAudio('GetDevices');
AudioFreq = dv(1).DefaultSampleRate; % Hertz
AudioChan = 2;
pahandle = PsychPortAudio('Open', [], [], 1, AudioFreq, AudioChan);

% Define properties of the auditory stimuli
StimFreqs = [350, 700, 1400; 500, 1000, 2000]; % Herzt
StimRise  = 0.007; % second
StimDur   = 0.050; % second

% Choose sounds
soundsname = {'bip', 'boup'};
nSnd = numel(soundsname);

% Create the sounds
SoundFreq = zeros(1, nSnd);
SoundMatrix = cell(1, nSnd);
for k = 1:nSnd
    SoundMatrix{k} = Emergence_Stimulation_CreateSound(StimFreqs(k,:), ...
        StimRise, StimDur, AudioFreq); % options for the sound
    SoundMatrix{k} = repmat(SoundMatrix{k}, [2,1]); % force stereo in case of mono
    
    % Play the sound and measure the delay
    PsychPortAudio('FillBuffer', pahandle, SoundMatrix{k});
    wantedt = GetSecs+1;
    obtainedt = PsychPortAudio('Start', pahandle, sndforced, wantedt, 1);
    fprintf('\nSound delay %s (%i/%i) = %1.4fs\n', ...
        soundsname{k}, k, nSnd, obtainedt - wantedt);
end
WaitSecs(1);

%% PLAY SEQUENCES THAT ARE USED AS EXAMPLES IN THE INSTRUCTIONS
%  ============================================================

% Select which sequence to present
t = 1;

% Define example sequences
seqs = {[2 1 1 1 2 1 1 2 2 1 1 1 1 2 1 2 2 1 1 1 1 1 1 2], ...
        [2 1 1 1 1 1 1 2 2 2 2 1 1 2 2 2 2 2 2 1 1 1 1 2], ...
        [2 1 1 2 1 2 1 2 1 2 1 2 1 2 2 1 2 1 1 2 1 2 1 2], ...
        [1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2], ...
        [2 2 2 1 2 2 2 1 2 2 2 1 2 2 2 1 2 2 2 1 2 2 2 1]};
    
% Define SOA
SOA = 0.3;

% Play sequentially each observation in the sequence
xg = GetSecs;
for tt = 1:numel(seqs{t})
    PsychPortAudio('FillBuffer', pahandle, SoundMatrix{seqs{t}(tt)});
    xg = PsychPortAudio('Start', pahandle, sndforced, xg + SOA, 1);
end

%% DISPLAY EXAMPLE SEQUENCES
%  =========================

% Load the demo
load('Emergence_Stimulation_Demo.mat');

% Define colors
col = [255 081 117; 049 132 182] ./ 255;

% Display the sequences
figure('Color', ones(1,3), 'Position', [1 456 1920 200]);
for i = 1:numel(D)
    subplot(numel(D), 1, i);
    for j = 1:2
        plot(find(D{i}.Seq == j), j, 'ko', 'MarkerFaceColor', ...
            col(j,:), 'LineWidth', 1); hold('on');
    end
    plot(repmat(D{i}.Jump,1,2), [-1,4], 'k-', 'LineWidth', 2);
    axis([1, numel(D{i}.Seq), -1, 4]); axis('off');
end
