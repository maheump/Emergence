% Generate an animated gif or an audiovisual movie of a single example
% sequence with progressively updating beliefs from an example subject and
% the ideal observer.
%
% Copyright (c) 2018 Maxime Maheu

%% INITIALIZATION
%  ==============

% Save the file
tosave = 'gif'; % can be '', 'gif' or 'avi'

% Which trajectory to use
sub = 23; % subject number
seq = cidx{2}(1); % sequence number

% Get corresponding data
X = cell(1,2);
X{1} = G{seq,sub};
X{2} = IO{seq,sub};
X{2}.Seq = G{seq,sub}.Seq;
N = numel(X{1}.Seq);

% Colors
g = repmat(0.6,1,3); % grey

%% PREPARE THE AUDITORY STIMULI
%  ============================

if strcmpi(tosave, 'avi')

% Define properties of the auditory stimuli
AudioFreq = 44100; % hertz
StimFreqs = [350, 700, 1400; 500, 1000, 2000]; % hertz
StimRise  = 0.007; % seconds
StimDur   = 0.05; % seconds
SOA       = 0.3; % seconds
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
audiowrite('aud.wav', audseq(1,:), AudioFreq);

end

%% GET FINGER'S TRAJECTORY
%  =======================

% Coordinates of the triangles limits (cartesian coordinates)
normtrglc = [0, sqrt(3)/2; ... % top left (P)
             1, sqrt(3)/2; ... % top right (D)
             1/2, 0];          % bottom (R)

% Prepare the window
figure('Color', ones(1,3), 'Units', 'Pixels', 'Position', [1 806 400 300]);
markers = {'.', 'p', 'o'};
stimcol = {'b', 'r'};
fs = 8;
ms = 8;

% Type of observer
obslab = {'Subject', {'Ideal','observer'}};

% Prepare the video file
if strcmpi(tosave, 'gif')
    filename = 'Emergence_FTA_ExampleSequence.gif';
    delete(filename);
elseif strcmpi(tosave, 'avi')
    vidObj = VideoWriter('vid.avi');
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
            tricol, eps, markers, false, fs);
        
        % Display the type of observer
        text(min(xlim), min(ylim), obslab{ind}, 'FontWeight', 'Bold', ...
            'HorizontalAlignment', 'Left', 'FontSize', fs);
        
        % Customize the axes
        set(gca, 'FontSize', fs);

        % Add some labels
        txt = condlab{seq};
        spc = strfind(txt, ' ');

        % Auditory sequence
        % ~~~~~~~~~~~~~~~~~
        subplot(4,3,(2:3)+shift); cla;

        % Display the sequence (up to sample "i")
        for j = 1:2
            idx = find(X{ind}.Seq(1:iObs) == j);
            plot(idx, repmat(j, [1,numel(idx)]), '.', 'Color', stimcol{j}); hold('on');
            text(0, j, {letters{j},''}, 'Color', stimcol{j}, 'FontSize', fs, ...
                'HorizontalAlignment', 'Left', 'VerticalAlignment', 'Top');
        end

        % Display the change point's position if it already occured
        if iObs >= J
            plot(repmat(X{ind}.Jump,1,2), [0,3], 'k-');
            plot(X{ind}.Jump, 3, markers{2}, 'MarkerSize', ms, ...
                'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w');
        end    

        % Customize the axes
        axis([1,N,0,3]);
        set(gca, 'XTick', [1,20:20:N], 'YColor', 'none');
        set(gca, 'Box', 'Off', 'Layer', 'Bottom', 'FontSize', fs);
        if ind == 1
            title(sprintf('%s sequence (%s)', txt(1:spc(1)-1), txt(spc(1)+1:end)));
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
                plot(repmat(X{ind}.Jump,1,2), [0,1], 'k-');
                plot(X{ind}.Jump, 1, markers{2}, 'MarkerSize', ms, ...
                    'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w');
            end

            % Customize the axes
            axis([1,N,0,1]);
            set(gca, 'XTick', [1,20:20:N]);
            set(gca, 'YTick', 0:1/3:1, 'YTickLabel', {'0','1/3','2/3','1'});
            set(gca, 'Box', 'Off', 'Layer', 'Bottom', 'FontSize', fs);

            % Add some labels
            if ind == 2, xlabel('Observation (#)'); end
            ylabel({'Barycentric','coordinates'});
        end
    end
    
    % Draw iteratively and get the frame
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    % Update the figure
    drawnow;
    
    % Get the image
    frame = getframe(gcf);
    
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
    system('avconv -i vid.avi -i aud.wav -c copy figs/F_M.avi');
end
