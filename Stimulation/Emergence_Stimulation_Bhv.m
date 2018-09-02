% This is the stimulation script. It presents auditory sequences and record
% finger motions.
% 
% Copyright (c) 2018 Maxime Maheu

%% INITIALIZE THE MATLAB ENVIRONMENT
%  =================================

% Clear the MALAB environment
clear; close('all'); clc; tic;

% Welcome message
wm = 51;
fprintf('%s\n##%sEmergence experiment%s##\n%s\n\n', repmat('#',1,wm), ...
    repmat(' ',1,(wm-17-8)/2), repmat(' ',1,(wm-17-6)/2), repmat('#',1,wm));

% Get project's home directory
user = char(java.lang.System.getProperty('user.name'));
if strcmpi(user, 'Maxime')
    homepath = '/Users/Maxime/Documents/My projects/Emergence';
elseif strcmpi(user, 'marie')
    homepath = 'C:\Users\marie\Documents\MATLAB\Emergence';
else
    homepath = mfilename('fullpath');
end

% Add it to the path
try rmpath(homepath); addpath(genpath(homepath));
catch, end
homepath = fullfile(homepath, 'Stimulation', 'bhv', 'Data');
cd(homepath);

% Launch the PsychToolbox
addpath(genpath('/usr/share/psychtoolbox-3/'));
PsychtoolboxVersion; sca;

% Seeds the random number generator based on the current time
rng('shuffle');

%% PREPARE SUBJECT'S DATA FILES AND DATA FOLDER
%  ============================================

% Get the number of the last subject
subjectcat = 'Subject'; % 'Pilot' or 'Subject'
folders = dir(homepath);
folders = {folders.name};
folders = folders(cellfun(@(x) ~any(strfind(x, '.')), folders)); % keep only folders
folders = folders(cellfun(@(x) strncmpi(x, subjectcat, numel(subjectcat)), folders));
if isempty(folders)
    lastsubject = 0;
else
    lastsubject = str2double(folders{end}(end-1:end));
end

% Get subject's number
S = [];
S.SubNum = input(sprintf('%s(%i) : .......: ', subjectcat, lastsubject+1));

% Define subject's name
S.SubLab = sprintf('%s%02i', subjectcat, S.SubNum);

% Create a folder for that subject
mkdir(S.SubLab); cd(S.SubLab);

% Ask the type od session to run
S.Session = '';
while ~any(strcmpi(S.Session, {'D', 'T', 'E'}))
    S.Session = input('Session : ..........: ', 's');
end
if     strcmpi(S.Session, 'D'), S.Session = 'Demo';
elseif strcmpi(S.Session, 'T'), S.Session = 'Training';
elseif strcmpi(S.Session, 'E'), S.Session = 'Experiment';
end

% Ask demographic questions about the subject
S.Sex = input('Sexe : .............: ', 's');
S.Age = input('Âge : ..............: ');
S.SDL = input('NSD : ..............: ');

% Wait a bit
WaitSecs(5);

%% PREPARE SCREEN
%  ==============

% Get screen specifications
screens = Screen('Screens');
screennum = screens(1);

% Define colors
white = ones(1,3).*255;
black = zeros(1,3);
grey  = repmat(128,1,3);
red   = [255, 000, 000];
blue  = [000, 000, 255];

% Skip test if the experiment is run on my Mac
if strcmpi(user, 'Maxime'), Screen('Preference', 'SkipSyncTests', 1); end

% Open a window
fullscreen = false;
[scw, sch] = Screen('WindowSize', screennum);
if      fullscreen, [windowPtr, rect] = Screen('OpenWindow', screennum, grey); HideCursor; 
elseif ~fullscreen, [windowPtr, rect] = Screen('OpenWindow', screennum, grey, ...
                                        round([1, 1, 1+0.16*3*scw, 1+0.09*3*sch]));
end
w_px = rect(3);
h_px = rect(4);

% Center of the screen
stimxc = (1/2) * w_px; % pixels
stimyc = (1/2) * h_px; % pixels

% Enable transparency
Screen('BlendFunction', windowPtr, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

% Text options
Screen('TextSize', windowPtr, 50);
Screen('TextFont', windowPtr, 'Arial');
Screen('TextStyle', windowPtr, 1);

% Keyboard specifications
KbName('UnifyKeyNames');

% Display home message
Screen('Flip', windowPtr);
DisplayCenteredText(windowPtr, '~~~~~~~~~~ EMERGENCE ~~~~~~~~~~', stimxc, h_px*0.4, white);
DisplayCenteredText(windowPtr, 'Préparation de l''expérience...', stimxc, h_px*0.6, white);
Screen('Flip', windowPtr);
WaitSecs(0.8);

%% PREPARE THE TRIANGLE-LIKE RESPONSE MODE
%  =======================================

% Get dimension of the equilateral triangle
trglw = 800; % pixels
trglh = (sqrt(3) / 2) * trglw; % pixels

% Centered position of the triangle on the screen
trglxc = (1/2) * w_px; % pixels
trglyc = (3/5) * h_px; % pixels

% Coordinates of the triangles sommet (x and y coordinates)
trglc = [trglxc - (trglw/2), trglyc - (trglh/2); ... % top left
         trglxc + (trglw/2), trglyc - (trglh/2); ... % top right
         trglxc            , trglyc + (trglh/2)];    % bottom
     
% Gravitational center of the triangle
trglxgc = trglc(3,1); % pixels
trglygc = trglc(1,2) + (1/3) * trglh; % pixels

% Define the main colors of the triangle
leftcol  = [066 146 198]; % blue
rightcol = [065 171 093]; % green
downcol = [239 059 033]; % red
tricol = [leftcol; rightcol; downcol] ./ 255;

% Draw the triangle
figure('Visible', 'Off');
patch('Faces', 1:3, 'Vertices', trglc, 'FaceVertexCData', tricol, ...
      'FaceColor', 'interp', 'EdgeColor', 'None');
axis('equal'); axis('off'); set(gca, 'YDir', 'Reverse');

% Save it as a *.png file
print('Triangle.png', '-dpng');
close(gcf);

% Read the *.png file
trlg = imread('Triangle.png');

% Remove the *.png file
delete('Triangle.png');

% Make a transparency matrix
alpha = repmat(255, size(trlg(:,:,1)));
alpha(trlg(:,:,1) == 255) = 0;
trlg(:,:,4) = alpha;

% Remove white borders
colson = sum(alpha, 1) == 0;
rowson = sum(alpha, 2) == 0;
alpha = alpha(~rowson,~colson);
trlgTex = NaN([size(alpha), 4]);
for i = 1:3, trlgTex(:,:,i) = trlg(~rowson,~colson,i); end
trlgTex(:,:,4) = alpha;

% Convert it into a PTB texture
trlgTex = Screen('MakeTexture', windowPtr, trlgTex);
trlgPos = round(CenterRectOnPoint([0, 0, trglw, trglh], trglxc, trglyc));

%% POSITIONS OF ELEMENTS TO DISPLAY
%  ================================

% Default vertical distance to center
vertdecal = 300;

% Define the press button for starting each trial
go = trglc(3,:);

% Define positions for questions
VertPosQuest = [h_px*0.2, h_px*0.7];

% Define positions for answers
VertPosAnswr = [h_px*0.3, h_px*0.8];
HorzPosAnswr = [w_px*0.25, w_px*0.75];

% Define positions for axes
HorzPosAxis = [w_px*0.15, w_px*0.65];
TickHeight = 10;
CircleRad = 20; % pixels

% Define positions for the confidence scale
HorzPosCart = w_px*0.75;
CartHeight = 400; % pixels
CartWidth  = 100; % pixels

% Size of the feedback bar
FBbarW = 0.7 * w_px; % pixels
FBbarH = 0.2 * h_px; % pixels

%% PREPARE THE AUDITORY STIMULI
%  ============================

% Initialize the PTB module
InitializePsychSound(1);
sndforced = 1;

% Close eventually open buffers
PsychPortAudio('Close');

% Open an audio buffer
dv = PsychPortAudio('GetDevices');
AudioFreq = dv(1).DefaultSampleRate;% Hertz
AudioChan = 2;
pahandle = PsychPortAudio('Open', [], [], 1, AudioFreq, AudioChan);

% Define properties of the auditory stimuli
StimFreqs = [350, 700, 1400; 500, 1000, 2000]; % Hz
StimRise  = 0.007; % ms
StimDur   = 0.05; % ms

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
    fprintf('\nSound delay %s (%i/%i) = %1.4fs', soundsname{k}, k, nSnd, obtainedt - wantedt);
end
WaitSecs(1);

%% PREPARE DESIGN
%  ==============

% Random function
RandFun = @(len) (rand(1, len) < 0.5) + 1;

% Save triangle specifications
S.Tri.w = trglw;
S.Tri.h = trglw;
S.Tri.x = trglxgc;
S.Tri.y = trglygc;
S.Tri.Coord = trglc;
S.Tri.Col   = tricol;
S.Tri.CoordRev = [S.Tri.Coord(:,1), h_px - S.Tri.Coord(:,2)];

% Define timings
S.SensMod = 'Audio';
S.SOA     = 0.3; % seconds
S.StimDur = max(StimDur); % seconds
S.ISI     = S.SOA - S.StimDur; % seconds
fprintf('\nThere will be %1.2f stimuli per second.', 1 / (S.StimDur + S.ISI));

% Define sequences' parameters
S.LenNoise.avg  = 100; % stimuli
S.LenNoise.std  = 15;  % stimuli
S.LenNoise.clip = [S.LenNoise.avg - 3*S.LenNoise.std, ... % stimuli
                   S.LenNoise.avg + 3*S.LenNoise.std];    % stimuli

%% CREATE SEQUENCES
%  ================

% Load already sequences
if strcmpi(S.Session, 'Demo')
    
%     % To produce the file:
%     RtT = {[1/5 1/5], [1 2], NaN};
%     save('Emergence_Stimulation_Demo.mat', 'RtT', 'D');
    
    % Load the demo file
    load('Emergence_Stimulation_Demo.mat');
    
    % Number of trials and number of stimuli in each sequence
    S.Ntrials = numel(RtT);
    S.Nstims  = numel(D{1}.Seq);
    S.SeqDur  = S.Nstims * (S.StimDur + S.ISI);
    
% Define regularities
else
    % Regularities to test
    if strcmpi(S.Session, 'Training')
        RtT = {NaN, NaN, ... % fully random
               [1 1 2], [1 1 2 2 2 2], ... % probabilistic regularities
               [1/6 1/6], [1/6 5/6], [5/6 5/6]}; % deterministic regularities
    elseif strcmpi(S.Session, 'Experiment')
        RtT = {NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, ... % fully random
               [1 1 2 2], ...             % AABB       => [An,Bn] (n = 2)
               [1 1 1 2 2 2], ...         % AAABBB     => [An,Bn] (n = 3)
               [1 1 1 1 2 2 2 2], ...     % AAAABBBB   => [An,Bn] (n = 4)
               [1 1 1 1 1 2 2 2 2 2], ... % AAAAABBBBB => [An,Bn] (n = 5)
               [1 1 2 1 2 2], ...         % AABABB     => [A2,BAn,B2] (n = 1)
               [1 1 2 1 2 1 2 2], ...     % AABABABB   => [A2,BAn,B2] (n = 2)
               [1 1 2 1 2 1 2 1 2 2], ... % AABABABABB => [A2,BAn,B2] (n = 3)
               [1 1 1 2 1 2], ...         % AAABAB     => not compressible (n = 6)
               [1 1 1 2 2 1 2 2], ...     % AAABBABB   => not compressible (n = 8)
               [1 1 1 2 1 1 2 2 1 2], ... % AAABAABBAB => not compressible (n = 10)
               [1/3 1/3], ... % p(A|B) = 1/3 & p(B|A) = 1/3 => Rep. freq (low)
               [1/4 1/4], ... % p(A|B) = 1/4 & p(B|A) = 1/4 => Rep. freq (med)
               [1/5 1/5], ... % p(A|B) = 1/5 & p(B|A) = 1/5 => Rep. freq (high)
               [2/3 2/3], ... % p(A|B) = 2/3 & p(B|A) = 2/3 => Alt. freq (low)
               [3/4 3/4], ... % p(A|B) = 3/4 & p(B|A) = 3/4 => Alt. freq (med)
               [4/5 4/5], ... % p(A|B) = 4/5 & p(B|A) = 4/5 => Alt. freq (high)
               [3/4 1/2], ... % p(A|B) = 3/4 & p(B|A) = 1/2 => Square (right)
               [1/2 1/4], ... % p(A|B) = 1/2 & p(B|A) = 1/4 => Square (bottom)
               [1/2 3/4], ... % p(A|B) = 1/2 & p(B|A) = 3/4 => Square (top)
               [1/4 1/2], ... % p(A|B) = 1/4 & p(B|A) = 1/2 => Square (left)
               [1/3 2/3], ... % p(A|B) = 1/3 & p(B|A) = 2/3 => Item freq. (low)
               [1/4 3/4], ... % p(A|B) = 1/4 & p(B|A) = 3/4 => Item freq. (medium)
               [1/5 4/5]};    % p(A|B) = 1/5 & p(B|A) = 4/5 => Item freq. (high)
    end

    % Shuffle the order of the regularities to test
    RtT = RtT(randperm(numel(RtT)));

    % Number of trials and number of stimuli in each sequence
    S.Ntrials = numel(RtT);
    S.Nstims  = 200;
    S.SeqDur  = S.Nstims * (S.StimDur + S.ISI);
    
    % Prepare sequences to be presented
    D = cell(1, S.Ntrials);
    for t = 1:S.Ntrials

        % Trial's specifications
        D{t}.Rule = RtT{t};
        if     isnan(D{t}.Rule), D{t}.Cond = 'Stochastic';
        elseif any(RtT{t} < 1),  D{t}.Cond = 'Probabilistic';
        elseif all(RtT{t} >= 1), D{t}.Cond = 'Deterministic';
        else,                    D{t}.Cond = 'ERROR';
        end

        % Which stimulus was used as As and Bs
        StimRand = randi([1,2]);
        D{t}.Stims = {soundsname{StimRand}, ...
                      soundsname{abs(StimRand-3)}};

        % Define the position of the change point
        LenNoise = 0;
        while LenNoise < S.LenNoise.clip(1) || LenNoise > S.LenNoise.clip(2)
            LenNoise = round(normrnd(S.LenNoise.avg, S.LenNoise.std));
        end
        if strcmpi(D{t}.Cond, 'Stochastic'), D{t}.Jump = NaN;
        else,                                D{t}.Jump = LenNoise + 0.5; end

        % Define the hidden generative process(es)
        if strcmpi(D{t}.Cond, 'Stochastic')
            D{t}.Gen = ones(1, S.Nstims);
        else
            D{t}.Gen = NaN(1, S.Nstims);
            D{t}.Gen(1:LenNoise)     = 1;
            D{t}.Gen(LenNoise+1:end) = 2;
        end

        % Create the sequence S.SubLab will be presented with
        if     strcmpi(D{t}.Cond, 'Deterministic')
            D{t}.Seq = [RandFun(LenNoise), repmat(D{t}.Rule, 1, S.Nstims)];
        elseif strcmpi(D{t}.Cond, 'Probabilistic')
            D{t}.Seq = [RandFun(LenNoise), GenRandSeq(S.Nstims, D{t}.Rule)];
        elseif strcmpi(D{t}.Cond, 'Stochastic')
            D{t}.Seq = RandFun(S.Nstims);
        end

        % Make sure the sequence is exactly the expected length
        D{t}.Seq = D{t}.Seq(1:S.Nstims);
    end

    % Save the design and the settings
    save(sprintf('%s_%s_Settings.mat', S.SubLab, S.Session), '-v7.3');
end

%% MAIN PRESENTATION LOOP
%  ======================

% Display some info
fprintf('\nThe sequences will last %i seconds (%i stimuli).\n', S.SeqDur, S.Nstims);
fprintf('There will be %i sequences.\n', S.Ntrials);
fprintf('The experiment will thus last ~%i minutes.\n', (S.SeqDur * S.Ntrials) / 60 * 2);
fprintf('\nExperiment started...');

% For each trial
for t = 1%:S.Ntrials
    
    % Prepare output variables
    D{t}.StimOnTime  = NaN(1, S.Nstims);
    D{t}.StimOffTime = NaN(1, S.Nstims);
    D{t}.DiscreteMouseCoord = NaN(2, S.Nstims);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%             DISPLAY THE TRIANGLE ON THE BACKGROUND              %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Clean the screen
    Screen('Flip', windowPtr);
    [xbeg, ybeg, xend, yend] = DisplayCenteredText(windowPtr, 'Go!', go(1), go(2), black);
    
    % Wait until the subject press on the go button
    begin = false;
    while ~begin
        
        % If the click is within the "Go!" area, break the loop
        [x,y,b] = GetMouse(windowPtr);
        if any(b == 1) && x >= xbeg && x <= xend && y >= ybeg && y <= yend
            begin = true;
        end
        
        % Display some information
        if ~begin
            
            % Display the position of the serie within the experiment
            DisplayCenteredText(windowPtr, sprintf('%s - Série %i/%i', datestr(now, ...
                'HH:MM'), t, S.Ntrials), stimxc, stimyc - vertdecal*1.5, black);

            % Display instructions
            DisplayCenteredText(windowPtr, 'Prêt ?', stimxc, stimyc - vertdecal, black);
            [xbeg, ybeg, xend, yend] = DisplayCenteredText(windowPtr, 'Go!', go(1), go(2), black);
        end
        
        % Draw the colored background of the triangle
        Screen('DrawTexture', windowPtr, trlgTex, [], trlgPos);

        % Draw inner lines
        for zzz = 1:3
            Screen('DrawLine', windowPtr, grey, trglxgc,      trglygc, ...
                                                trglc(zzz,1), trglc(zzz,2));
        end

        % Draw triangle's limits
        Screen('FramePoly', windowPtr, black, round(trglc), 3);

        % Display the vertices' labels
        DisplayCenteredText(windowPtr, 'Entièrement', ...
            trglc(3,1), trglc(3,2)+50, downcol);
        DisplayCenteredText(windowPtr, 'aléatoire', ...
            trglc(3,1), trglc(3,2)+100, downcol);
        if mod(S.SubNum, 2) == 0 % pair
            DisplayCenteredText(windowPtr, 'Événements', ... % left
                trglc(1,1)-200, trglc(1,2)-50, leftcol);
            DisplayCenteredText(windowPtr, 'plus fréquents', ... % left
                trglc(1,1)-200, trglc(1,2)-00, leftcol);
            DisplayCenteredText(windowPtr, 'Motif', ... % right
                trglc(2,1)+100, trglc(2,2)-50, rightcol);
            DisplayCenteredText(windowPtr, 'répété', ... % right
                trglc(2,1)+100, trglc(2,2)-00, rightcol);
        elseif mod(S.SubNum, 2) == 1 % impair
            DisplayCenteredText(windowPtr, 'Événements', ... % right
                trglc(2,1)+200, trglc(1,2)-50, rightcol);
            DisplayCenteredText(windowPtr, 'plus fréquents', ... % right
                trglc(2,1)+200, trglc(1,2)-00, rightcol);
            DisplayCenteredText(windowPtr, 'Motif', ... % left
                trglc(1,1)-100, trglc(2,2)-50, leftcol);
            DisplayCenteredText(windowPtr, 'répété', ... % left
                trglc(1,1)-100, trglc(2,2)-00, leftcol);
        end
        
        % Flip the screen
        Screen('Flip', windowPtr);
    end
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%                     PREPARE A MOUSE CALLBACK                    %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
    % Prepare a text file for the recording of the mouse coordinates
    mousecoordtxtfile = sprintf('%s%s%s_MouseCoord.txt', homepath, filesep, S.SubLab);
    fid = fopen(mousecoordtxtfile, 'w');
    fprintf(fid, 'Trial,Time,Xcoord,Ycoord\n');
    fclose(fid);

    % Prepare the mouse recording as a MATLAB timer
    MouseCoordObj = timer('TimerFcn', {@Emergence_Stimulation_MouseCoordCallback, ...
                          t, mousecoordtxtfile}, ... % which text file
                          'StartDelay', 0, ... % delay in execution
                          'TasksToExecute', +Inf, ... % execute until it is manually stopped
                          'Period', 0.01, ... % the rate (hertz) at which record mouse position
                          'ExecutionMode', 'fixedRate');
    
    % Set the mouse to the bottom of the triangle
    SetMouse(go(1), go(2), windowPtr);
    
    % Start the timer for the current sequence
    start(MouseCoordObj);
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%                      STIMULUS PRESENTATION                      %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % For each stimulus
    fprintf('\n\t- Sequence %03i/%03i: ', t, S.Ntrials);
    xg = GetSecs;
    for tt = 1:S.Nstims
        
        % Get the sound to play
        stimulus = D{t}.Seq(tt);
        soundname = D{t}.Stims{stimulus};
        fprintf('%s, ', soundname);
        sndix = find(cellfun(@(x) any(strfind(x, soundname)), soundsname));
        
        % Play the sound
        PsychPortAudio('FillBuffer', pahandle, SoundMatrix{sndix});
        xg = PsychPortAudio('Start', pahandle, sndforced, xg + S.SOA, 1);
        
        % Get the presentation timing
        D{t}.StimOnTime(tt)  = xg;
        D{t}.StimOffTime(tt) = xg+StimDur;
        
        % Get the position of the finger
        [D{t}.DiscreteMouseCoord(1,tt), y] = GetMouse(windowPtr);
        D{t}.DiscreteMouseCoord(2,tt) = abs(y - h_px);
    end
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%                       IMPORT TRAJECTORY                         %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Import mouse coordinates data
    stop(MouseCoordObj);
    mouse = importdata(mousecoordtxtfile, ',', 1);
    delete(mousecoordtxtfile);
    
    % Save that into a MATLAB variable
    D{t}.ContinuousMouseTime  = mouse.data(:,2)'; % trial numnber and time
    D{t}.ContinuousMouseCoord = mouse.data(:,end-1:end)'; % (0,0) => bottom-left
    
    % Keep only measurements between the first and the last stimulus
    [~,begidx] = min(abs(D{t}.ContinuousMouseTime - D{t}.StimOnTime(1)));
    [~,endidx] = min(abs(D{t}.ContinuousMouseTime - D{t}.StimOffTime(end)));
    D{t}.ContinuousMouseTime  = D{t}.ContinuousMouseTime(   begidx:endidx);
    D{t}.ContinuousMouseCoord = D{t}.ContinuousMouseCoord(:,begidx:endidx);
    
    % Mirror the trajectory such that deterministic regularity is always on
    % the right and probabilistic regularity on the left
    if mod(S.SubNum, 2) == 1 % impair
        D{t}.DiscreteMouseCoord(1,:)   = w_px - D{t}.DiscreteMouseCoord(1,:);
        D{t}.ContinuousMouseCoord(1,:) = w_px - D{t}.ContinuousMouseCoord(1,:);
    end
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%             QUESTION ON THE PRESENCE OF REGULARITY              %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Clear the sceen
    Screen('Flip', windowPtr);
    
    % Prepare the variable that will receive answers to the questions
    D{t}.Questions = NaN(1,6);
    
    % Display question
    DisplayCenteredText(windowPtr, 'Avez-vous détecté une régularité ?', ...
        stimxc, VertPosQuest(1), black);
    
    % Display possible answers
    xbeg = NaN(1,2); ybeg = NaN(1,2); xend = NaN(1,2); yend = NaN(1,2);
    [xbeg(1), ybeg(1), xend(1), yend(1)] = DisplayCenteredText(windowPtr, ...
        'Oui', HorzPosAnswr(1), VertPosAnswr(1), black);
    [xbeg(2), ybeg(2), xend(2), yend(2)] = DisplayCenteredText(windowPtr, ...
        'Non', HorzPosAnswr(2), VertPosAnswr(1), black);
    Screen('Flip', windowPtr, GetSecs, 1);
    
    % Look for answer
    while isnan(D{t}.Questions(1))
        
        % Get mouse position
        [x,y,b] = GetMouse(windowPtr);
        
        % Check if the finger is in the response 1 box
        if any(b == 1) && x >= xbeg(1) && x <= xend(1) && y >= ybeg(1) && y <= yend(1)
            D{t}.Questions(1) = 1;
            DisplayCenteredText(windowPtr, 'Oui', HorzPosAnswr(1), VertPosAnswr(1), red);
        end
        
        % Check if the finger is in the response 2 box
        if any(b == 1) && x >= xbeg(2) && x <= xend(2) && y >= ybeg(2) && y <= yend(2)
            D{t}.Questions(1) = 2;
            DisplayCenteredText(windowPtr, 'Non', HorzPosAnswr(2), VertPosAnswr(1), red);
        end
        Screen('Flip', windowPtr, GetSecs, 1); % do not erase the question
    end
    
    % Remaining questions are ask only if the answer to the first one is yes
    if D{t}.Questions(1) == 1
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%         QUESTION ON THE IDENTITY OF THE REGULARITY          %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Display question
        DisplayCenteredText(windowPtr, 'De quel type était cette régularité ?', ...
            stimxc, VertPosQuest(2), black);
        
        % Display possible answers
        xbeg = NaN(1,2); ybeg = NaN(1,2); xend = NaN(1,2); yend = NaN(1,2);
        if mod(S.SubNum, 2) == 0 % pair
            [xbeg(1), ybeg(1), xend(1), yend(1)] = DisplayCenteredText(windowPtr, ... % left
                'Événements plus fréquents', HorzPosAnswr(1), VertPosAnswr(2), black);
            [xbeg(2), ybeg(2), xend(2), yend(2)] = DisplayCenteredText(windowPtr, ... % right
                'Motif répété', HorzPosAnswr(2), VertPosAnswr(2), black);
        elseif mod(S.SubNum, 2) == 1 % impair
            [xbeg(1), ybeg(1), xend(1), yend(1)] = DisplayCenteredText(windowPtr, ... % right
                'Événements plus fréquents', HorzPosAnswr(2), VertPosAnswr(2), black);
            [xbeg(2), ybeg(2), xend(2), yend(2)] = DisplayCenteredText(windowPtr, ... % left
                'Motif répété', HorzPosAnswr(1), VertPosAnswr(2), black);
        end
        Screen('Flip', windowPtr, GetSecs, 1); % do not erase the question
        
        % Look for answer
        while isnan(D{t}.Questions(2))
            
            % Get mouse position
            [x,y,b] = GetMouse(windowPtr);
            
            % Check if the finger is in the response 1 box
            if any(b == 1) && x >= xbeg(1) && x <= xend(1) && y >= ybeg(1) && y <= yend(1)
                D{t}.Questions(2) = 1;
                if mod(S.SubNum, 2) == 0 % pair
                    DisplayCenteredText(windowPtr, 'Événements plus fréquents', ...
                        HorzPosAnswr(1), VertPosAnswr(2), red); % left
                elseif mod(S.SubNum, 2) == 1 % impair
                    DisplayCenteredText(windowPtr, 'Événements plus fréquents', ...
                        HorzPosAnswr(2), VertPosAnswr(2), red); % right
                end
            end
            
            % Check if the finger is in the response 2 box
            if any(b == 1) && x >= xbeg(2) && x <= xend(2) && y >= ybeg(2) && y <= yend(2)
                D{t}.Questions(2) = 2;
                if mod(S.SubNum, 2) == 0 % pair
                    DisplayCenteredText(windowPtr, 'Motif répété', ...
                        HorzPosAnswr(2), VertPosAnswr(2), red); % right
                elseif mod(S.SubNum, 2) == 1 % impair
                    DisplayCenteredText(windowPtr, 'Motif répété', ...
                        HorzPosAnswr(1), VertPosAnswr(2), red); % left
                end
            end
            Screen('Flip', windowPtr, GetSecs, 1); % do not erase the question
        end
        
        % Clear the screen
        WaitSecs(1);
        Screen('Flip', windowPtr);
        
        % Present 2 different questions with a scale and a confidence
        % judgment
        quest = {'est-elle apparue', 'l''avez vous détecté'};
        idx = [2,4];
        for sci = 1:2
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%            QUESTION WITH A SEQUENCE-LIKE SCALE          %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Display question
            DisplayCenteredText(windowPtr, sprintf('A quel moment %s ?', ...
                quest{sci}), stimxc, VertPosQuest(sci), black);

            % For the first question, diaplsy the entire axis but, for the
            % second one, only the right part of it
            if     sci == 1, begax = HorzPosAxis(1);
            elseif sci == 2, begax = lastaxpos;
            end
            
            % Display the axis
            Screen('DrawLine', windowPtr, black, begax, VertPosAnswr(sci), ...
                HorzPosAxis(2), VertPosAnswr(sci), 3);
            
            % Add ticks on the scale
            majorticks = linspace(HorzPosAxis(1), HorzPosAxis(2), 11);
            if sci == 2, majorticks = majorticks(majorticks >= begax); end
            for tick = majorticks
                Screen('DrawLine', windowPtr, black, tick, VertPosAnswr(sci)...
                    -TickHeight, tick, VertPosAnswr(sci)+TickHeight, 3);
                ticklab = (tick - HorzPosAxis(1)) / (HorzPosAxis(2) - HorzPosAxis(1)); % pct
                Screen('TextSize', windowPtr, Screen('TextSize', windowPtr) / 2);
                DisplayCenteredText(windowPtr, sprintf('%i', round(ticklab * S.Nstims)), ...
                    tick, VertPosAnswr(sci) + TickHeight*4, black);
                Screen('TextSize', windowPtr, Screen('TextSize', windowPtr) * 2);
            end
            
            % On the second axis, display the answer to the previous question 
            if sci == 2
                Screen('FillOval', windowPtr, red, CenterRectOnPoint([0, ...
                    0, CircleRad, CircleRad], begax, VertPosAnswr(sci)));
            end
            
            % Flip the screen
            Screen('Flip', windowPtr, GetSecs, 1);
            
            % Look for answer
            while isnan(D{t}.Questions(idx(sci)+1))
                
                % Get mouse position
                [x,y,b] = w(windowPtr);
                
                % Check if the finger position is on the axis
                if any(b == 1) && y >= VertPosAnswr(sci) - TickHeight ...
                               && y <= VertPosAnswr(sci) + TickHeight ...
                               && x >= begax && x <= HorzPosAxis(2)
                    
                    % Record answer
                    lastaxpos = x;
                    pct = (x - HorzPosAxis(1)) / (HorzPosAxis(2) - HorzPosAxis(1)); % percent
                    D{t}.Questions(idx(sci)+1) = round(pct*S.Nstims); % in percent
                    
                    % Display chosen answer
                    Screen('FillOval', windowPtr, red, CenterRectOnPoint(...
                        [0, 0, CircleRad, CircleRad], x, VertPosAnswr(sci)));
                    Screen('Flip', windowPtr, GetSecs, 1);
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%           QUESTION WITH A CONFIDENCE CARTOUCHE          %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Display the cartouche
            CartPos = CenterRectOnPoint([0, 0, CartWidth, CartHeight], ...
                HorzPosCart, VertPosAnswr(sci));
            Screen('FrameRect', windowPtr, black, CartPos, 3);
            Screen('Flip', windowPtr, GetSecs, 1);
            
            % Look for answer
            while isnan(D{t}.Questions(idx(sci)+2))
                
                % Get mouse position
                [x,y,b] = GetMouse(windowPtr);
                
                % Check if the finger position is within the rectangle
                if any(b == 1) && x >= CartPos(1) && x <= CartPos(3) ...
                               && y >= CartPos(2) && y <= CartPos(4)
                    
                    % Record answer
                    pct = (y - CartPos(2)) / (CartPos(4) - CartPos(2)); % percent
                    D{t}.Questions(idx(sci)+2) = 100 - round(pct*100); % in percent
                    
                    % Display chosen answer
                    Screen('FillRect', windowPtr, red, [CartPos(1), y, CartPos(3), CartPos(4)]);
                    Screen('FrameRect', windowPtr, black, CartPos, 3);
                    Screen('Flip', windowPtr, GetSecs, 1);
                end
            end
        end
    end
    WaitSecs(1);
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%                     DISPLAY THE FEEDBACK                        %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
	% Clear the sceen
    Screen('Flip', windowPtr);
    
    % Do not display feedback during the real experiment
    if ~strcmpi(S.Session, 'Experiment')

        % Convert the change point position in percent
        j = D{t}.Jump/S.Nstims;
        if isnan(j), j = 1; end

        % Display a title
        DisplayCenteredText(windowPtr, 'Composition de la précédente série :', ...
            stimxc, stimyc - FBbarH, white);

        % Display the generative process of the first part of the sequence
        Screen('FillRect', windowPtr, downcol, [stimxc - FBbarW/2, ...
                                                stimyc - FBbarH/2, ...
                                                stimxc - FBbarW/2 + FBbarW*j, ...
                                                stimyc + FBbarH/2]);
        DisplayCenteredText(windowPtr, 'Entièrement', ...
            stimxc - FBbarW/2 + FBbarW*j/2, stimyc - FBbarH/6, white);
        DisplayCenteredText(windowPtr, 'aléatoire', ...
            stimxc - FBbarW/2 + FBbarW*j/2, stimyc + FBbarH/6, white);

        % Use the same colors as the ones in the triangle
        if mod(S.SubNum, 2) == 0 % pair
            fbcols = [leftcol; rightcol];
        elseif mod(S.SubNum, 2) == 1 % impair
            fbcols = [rightcol; leftcol];
        end

        % If there was a second non-stochastic part in the sequence
        if ~strcmpi(D{t}.Cond, 'Stochastic')

            % Display the generative process of the second part of the sequence
            if strcmpi(D{t}.Cond, 'Probabilistic')
                Screen('FillRect', windowPtr, fbcols(1,:), [stimxc - FBbarW/2 + FBbarW*j, ...
                                                            stimyc - FBbarH/2, ...
                                                            stimxc + FBbarW/2, ...
                                                            stimyc + FBbarH/2]);
                DisplayCenteredText(windowPtr, 'Événements', ...
                    stimxc - FBbarW/2 + FBbarW*j + FBbarW*(1-j)/2, stimyc - FBbarH/6, white);
                DisplayCenteredText(windowPtr, 'plus fréquents', ...
                    stimxc - FBbarW/2 + FBbarW*j + FBbarW*(1-j)/2, stimyc + FBbarH/6, white);

            % Display the generative process of the second part of the sequence
            elseif strcmpi(D{t}.Cond, 'Deterministic')
                Screen('FillRect', windowPtr, fbcols(2,:), [stimxc - FBbarW/2 + FBbarW*j, ...
                                                            stimyc - FBbarH/2, ...
                                                            stimxc + FBbarW/2, ...
                                                            stimyc + FBbarH/2]);
                DisplayCenteredText(windowPtr, 'Motif', ...
                    stimxc - FBbarW/2 + FBbarW*j + FBbarW*(1-j)/2, stimyc - FBbarH/6, white);
                DisplayCenteredText(windowPtr, 'répété', ...
                    stimxc - FBbarW/2 + FBbarW*j + FBbarW*(1-j)/2, stimyc + FBbarH/6, white);
            end
        end

        % Add a line around the bar
        Screen('FrameRect', windowPtr, black, [stimxc - FBbarW/2, stimyc - FBbarH/2, ...
                                               stimxc + FBbarW/2, stimyc + FBbarH/2], 3);

        % Add major ticks        
        majorticks = linspace(stimxc - FBbarW/2, stimxc + FBbarW/2, 11);
        for tick = majorticks
            Screen('DrawLine', windowPtr, black, tick, stimyc + FBbarH/2, ...
                                                 tick, stimyc + FBbarH/2 - TickHeight, 3);
        end

        % Flip the screen
        Screen('Flip', windowPtr);

        % Wait for a press somewhere in the screen
        begin = false;
        while ~begin
            [x,y,b] = GetMouse(windowPtr);
            if any(b == 1), begin = true; end
        end
    end
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%                   SAVE DATA IN A SEPARATE FILE                  %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Save data
    s = S;
    d = D{t};
    save(sprintf('%s_%s_Seq%i.mat', S.SubLab, S.Session, t), '-v7.3', 's', 'd');
end

%% END OF THE EXPERIMENT
%  =====================

% Get the time required for the experiment
totaltime = toc;

% Thanks for participating message
DisplayCenteredText(windowPtr, 'Merci d''avoir participé !', ...
    stimxc, stimyc, white);
Screen('Flip', windowPtr);

% Save data and remove temporary files
try save(sprintf('%s_%s_All.mat', S.SubLab, S.Session), '-v7.3');
catch, warning('!!! SAVING HAS FAILED !!!'); end

% Close the PsychToolbox
WaitSecs(3);
sca; fprintf('\n');
PsychPortAudio('Close', pahandle);
if fullscreen, ShowCursor; end

%% DISPLAY THE TRAJECTORIES
%  ========================

figure('Color', ones(1,3), 'Units', 'Normalized', 'Position', [0 0 1 1]);
for t = 1:S.Ntrials
    subplot(4,9,t);
    fill(S.Tri.CoordRev(:,1), S.Tri.CoordRev(:,2), 'k-', 'FaceAlpha', 0); hold('on');
    try
        plot(D{t}.DiscreteMouseCoord(1,:), D{t}.DiscreteMouseCoord(2,:), '.-');
        plot(D{t}.ContinuousMouseCoord(1,:), D{t}.ContinuousMouseCoord(2,:), '.--');
    catch
    end
    text(S.Tri.CoordRev(1,1), S.Tri.CoordRev(1,2), 'P', ...
        'HorizontalAlignment', 'Center', 'FontWeight', 'Bold');
    text(S.Tri.CoordRev(2,1), S.Tri.CoordRev(2,2), 'D', ...
        'HorizontalAlignment', 'Center', 'FontWeight', 'Bold');
    text(S.Tri.CoordRev(3,1), S.Tri.CoordRev(3,2), 'S', ...
        'HorizontalAlignment', 'Center', 'FontWeight', 'Bold');
    axis('equal');
    axis([1, w_px, 1, h_px]);
    set(gca, 'XTick', [], 'YTick', []);
    title(D{t}.Cond);
end
