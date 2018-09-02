% Generate an animated gif or an audiovisual movie of a single example
% sequence with progressively updating beliefs from an example subject and
% the ideal observer.
%
% Copyright (c) 2018 Maxime Maheu

% N.B. To be run without subject exclusion. To do so, turn the boolean
% variable "rmbadsub" in "Emergence_FTA_LoadData" to false.

%% INITIALIZATION
%  ==============

% Make sure data from bad subjects is available
if rmbadsub
    rmbadsub = false;
    Emergence_FTA_LoadData;
end

% Save the file
tosave = 'gif'; % can be '', 'gif' or 'avi'

% Which trajectory to use
seq = 28; % subject number
sub = 22; % sequence number

% Get corresponding data
X = cell(1,2);
X{1} = G{sub,seq};
X{2} = IO{sub,seq};
X{2}.Seq = G{sub,seq}.Seq;
N = numel(X{1}.Seq);

%% PREPARE THE AUDITORY STIMULI
%  ============================

if strcmpi(tosave, 'avi')

% Define properties of the auditory stimuli
AudioFreq = 44100; % hertz
StimFreqs = [350, 700, 1400; 500, 1000, 2000]; % hertz
StimRise  = 0.007; % second
StimDur   = 0.050; % second
SOA       = 0.300; % second
StimOff   = SOA - StimDur; % seconds

% Choose sounds
soundsname = {'beep', 'boop'};
nSnd = numel(soundsname);

% Create the sounds
SoundFreq = zeros(1, nSnd);
SoundMatrix = cell(1, nSnd);
for k = 1:nSnd
    SoundMatrix{k} = Emergence_Stimulation_CreateSound(StimFreqs(k,:), ...
        StimRise, StimDur, AudioFreq); % options for the sound
    SoundMatrix{k} = repmat(SoundMatrix{k}, [2,1]); % force stereo in case of mono
    silence = zeros(2,round(StimOff*AudioFreq));
    SoundMatrix{k} = [SoundMatrix{k}, silence];
end

% Create the entire auditory sequence
audseq = cell2mat(SoundMatrix(X{1}.Seq));
seqdur = size(audseq,2) / AudioFreq;
audfile = fullfile(ftapath, 'figs', 'aud.wav');
audiowrite(audfile, audseq(1,:), AudioFreq);

end

%% GET FINGER'S TRAJECTORY
%  =======================

% Prepare the window
ff = figure('Color', ones(1,3), 'Units', 'Pixels', 'Position', [1 806 400 300]);
fs = 8;
ms = 8;

% Type of observer
obslab = {{'Example', 'subject'}, {'Ideal','observer'}};

% Prepare the video file
if strcmpi(tosave, 'gif')
    filename = fullfile(ftapath, 'figs', 'F_M.gif');
    delete(filename);
elseif strcmpi(tosave, 'avi')
    vidfile = fullfile(ftapath, 'figs', 'vid.avi');
    vidObj = VideoWriter(vidfile);
    vidObj.FrameRate = N/seqdur;
    open(vidObj);
end

% After each observation
for iObs = 1:N
    
    % For subject's and IO's inference
    for ind = 1:2
        shift = (2*3)*(ind-1);
        
        % Triangle
        % ~~~~~~~~
        subplot(4,3,[1,4]+shift); cla;
        
        % Check whether the change point already happened or not
        J = X{ind}.Jump+1/2;
        if iObs < J, J = NaN; end
        
        % Plot the trajectory (up to sample "i") on the triangle
        Emergence_PlotTrajOnTri(X{ind}.BarycCoord(1:iObs,:), J, ...
            tricol, eps, {'', '', ''}, false, fs);
        Emergence_PlotGridOnTri(3);
        cc = X{ind}.BarycCoord(iObs,:)*tricc;
        plot(cc(1), cc(2), 'k.', 'MarkerSize', 10)
        
        % Display the type of observer
        text(min(xlim), min(ylim), obslab{ind}, 'HorizontalAlignment', ...
            'Left', 'FontSize', fs, 'Interpreter', 'LaTeX');
        
        % Customize the axes
        set(gca, 'FontSize', fs);
        
        % Add some labels
        txt = condlab{sub};
        spc = strfind(txt, ' ');
        
        % Auditory sequence
        % ~~~~~~~~~~~~~~~~~
        subplot(4,3,(2:3)+shift); cla;
        
        % Display the sequence (up to sample "i")
        for j = 1:2
            idx = find(X{ind}.Seq(1:iObs) == j);
            plot(1.2+(X{ind}.Seq(1:iObs)-1).*0.6, '-', 'Color', g, ...
                'LineWidth', 1/2); hold('on');
            plot(idx, repmat(j, [1,numel(idx)]), 'k.', 'MarkerSize', 5); 
            text(0, j, {letters{j},''}, 'Color', 'k', 'FontSize', fs, ...
                'HorizontalAlignment', 'Right', 'VerticalAlignment', 'Top', ...
                'Interpreter', 'LaTeX');
        end
        
        % Display the change point's position if it already occured
        if iObs >= J
            plot(repmat(X{ind}.Jump,1,2), [0.8,2.2], 'k-', 'LineWidth', 1/2);
        end
        
        % Customize the axes
        axis([1,N,-1,3]); axis('off');
        set(gca, 'XTick', [1,20:20:N], 'XTickLabel', {}, 'YDir', 'reverse');
        set(gca, 'Box', 'Off', 'Layer', 'Bottom', 'FontSize', fs, ...
            'TickLabelInterpreter', 'LaTeX', 'LineWidth', 1/2);
        if ind == 1
            title(sprintf('%s sequence (%s)', txt(1:spc(1)-1), ...
                txt(spc(1)+1:end)), 'Interpreter', 'LaTeX');
        end
        
        % Barycentric coordinates
        % ~~~~~~~~~~~~~~~~~~~~~~~
        if iObs > 1
            subplot(4,3,(5:6)+shift); cla;
            
            % Display the cumulative Barycentric coordinates
            Emergence_PlotBarycTraj(X{ind}.BarycCoord(1:iObs,:), tricol);
            
            % Draw help lines (change levels)
            plot([1,iObs], repmat(1/3,1,2), '-', 'Color', g);
            plot([1,iObs], repmat(2/3,1,2), '-', 'Color', g);
            
            % Display the change point's position if it already occured
            if iObs >= J
                plot(repmat(X{ind}.Jump,1,2), [0,1], 'k-', 'LineWidth', 1/2);
            end
            
            % Customize the axes
            axis([1,N,0,1]);
            set(gca, 'XTick', [1,20:20:N]);
            set(gca, 'YTick', 0:1/3:1, 'YTickLabel', {'0','1/3','2/3','1'});
            set(gca, 'Box', 'Off', 'Layer', 'Bottom', 'FontSize', fs, ...
                'TickLabelInterpreter', 'LaTeX', 'LineWidth', 1/2);
            
            % Add some labels
            if ind == 2, xlabel('Observation ($K$)', 'Interpreter', 'LaTeX'); end
            ylabel('$p(\mathcal{M}|y_{1:K})$', 'Interpreter', 'LaTeX');
        end
    end
    
    % Draw iteratively and get the frame
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    % Update the figure
    drawnow;
    
    % Get the image
    frame = getframe(ff);
    
    % Append to the file
    if strcmpi(tosave, 'gif')
        im = frame2im(frame);
        [imind, cm] = rgb2ind(im, 256);
        if     iObs == 1, imwrite(imind, cm, filename, 'gif', 'Loopcount', inf);
        elseif iObs  > 1, imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append');
        end
    elseif strcmpi(tosave, 'avi')
    	writeVideo(vidObj, frame);
    end
end

% Close the file
if strcmpi(tosave, 'avi'), close(vidObj); end

%% GENERATE AN AUDIO-VISUAL AVI FILE
%  =================================

% Combine auditory and visual files into a movie
if strcmpi(tosave, 'avi')
    % The avconv plugin works on Mac OS and allow to combine visual and auditory
    % files together. Instructions for installation can be found at:
    % https://libav.org/avconv.html
    allfile = fullfile(ftapath, 'figs', 'F_M.avi');
    cmd = sprintf('avconv -i ''%s'' -i ''%s'' -c copy ''%s''', vidfile, audfile, allfile);
    try system(cmd);
    catch
    end
    fprintf(['If the *.avi file has not been generated. First, make ', ...
        'sure you have installed ''avconv'' from ''https://libav.org/avconv.html''. ', ...
        'Then directly run the following command in the terminal:\n%s\n'], cmd);
end
