% Generate animated gif or audiovisual movie from a single example
% sequence.
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
figure('Color', ones(1,3), 'Units', 'Pixels', 'Position', [1 655 700 400]);
markers = {'.', 'p', 'o'};

% Prepare the video file
if strcmpi(tosave, 'gif')
    filename = 'images/Emergence_FTA_ExampleSequence.gif';
elseif strcmpi(tosave, 'avi')
    vidObj = VideoWriter('vid.avi');
    vidObj.FrameRate = N/seqdur;
    open(vidObj);
end

% After each observation
for i = 1:N
    
    % For subject's and IO's inference
    for ind = 1:2
        shift = (2*3)*(ind-1);

        % Triangle
        % ~~~~~~~~
        subplot(4,3,[1,4]+shift); cla;

        % Check whether the change point already happened or not
        J = X{ind}.Jump+1/2;
        if i < J, J = NaN; end

        % Plot the trajectory (up to sample "i") on the triangle
        Emergence_PlotTrajOnTri(X{ind}.BarycCoord(1:i,:), J, tricol, eps, markers);

        % Customize the axes
        set(gca, 'LineWidth', 1, 'FontSize', 15);

        % Add some labels
        txt = condlab{seq};
        spaces = strfind(txt, ' ');

        % Auditory sequence
        % ~~~~~~~~~~~~~~~~~
        subplot(4,3,(2:3)+shift); cla;

        % Display the sequence (up to sample "i")
        col = {'b', 'r'};
        for j = 1:2
            idx = find(X{ind}.Seq(1:i) == j);
            plot(idx, repmat(j, [1,numel(idx)]), '.', 'Color', col{j}); hold('on');
            text(0, j, {soundsname{j},''}, 'Color', col{j}, ...
                'HorizontalAlignment', 'Left', 'VerticalAlignment', 'Top');
        end

        % Display the change point's position if it already occured
        if i >= J
            plot(repmat(X{ind}.Jump,1,2), [0,3], 'k-', 'LineWidth', 1);
            plot(X{ind}.Jump, 3, markers{2}, 'MarkerSize', 15, ...
                'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', 'LineWidth', 1);
        end    

        % Customize the axes
        axis([1,N,0,3]);
        set(gca, 'XTick', [1, get(gca, 'XTick')], 'YColor', 'none');
        set(gca, 'Box', 'Off', 'Layer', 'Bottom', 'LineWidth', 1, 'FontSize', 12);
        if ind == 1, title({sprintf('%s sequence (%s)', txt(1:spaces(1)-1), ...
                txt(spaces(1)+1:end)), ''}); end

        % Barycentric coordinates
        % ~~~~~~~~~~~~~~~~~~~~~~~
        if i > 1
            subplot(4,3,(5:6)+shift); cla;

            % Display the cumulative Barycentric coordinates
            Emergence_PlotBarycTraj(X{ind}.BarycCoord(1:i,:), tricol);

            % Draw help lines (change levels)
            plot([1,i], repmat(1/3,1,2), '-', 'Color', g, 'LineWidth', 1);
            plot([1,i], repmat(2/3,1,2), '-', 'Color', g, 'LineWidth', 1);

            % Display the change point's position if it already occured
            if i >= J
                plot(repmat(X{ind}.Jump,1,2), [0,1], 'k-', 'LineWidth', 1);
                plot(X{ind}.Jump, 1, markers{2}, 'MarkerSize', 15, ...
                    'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', 'LineWidth', 1);
            end

            % Customize the axes
            axis([1,N,0,1]);
            set(gca, 'XTick', [1, get(gca, 'XTick')]);
            set(gca, 'YTick', 0:1/3:1, 'YTickLabel', {'0','1/3','2/3','1'});
            set(gca, 'Box', 'Off', 'Layer', 'Bottom', 'LineWidth', 1, 'FontSize', 12);

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
        if     i == 1, imwrite(imind, cm, filename, 'gif', 'Loopcount', inf); 
        elseif i  > 1, imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append'); 
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
    system('avconv -i vid.avi -i aud.wav -c copy Emergence_FTA_ExampleSequence.avi');
end
