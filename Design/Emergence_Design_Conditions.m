% This script shows the different types of regularities that have been used
% (both the probabilistic and the deterministic ones) as well as the
% different types of stimuli.
% 
% Maxime Maheu, March 2017

%% PROBABILISTIC REGULARITIES
%  ==========================

% Compute an entropy map
t = 0:0.001:1;
H = @(p) -(p .* log2(p) + (1-p) .* log2(1-p));
entmap = H(t)' + H(t);

% Prepare the window
figure('Position', [1 765 255 340]);

% Display entropy map
imagesc(t, t, entmap); hold('on');
colormap(plasma(2000)); caxis([0,2]);
cbr = colorbar('Location', 'SouthOutside', 'LineWidth', 1);
cbr.Label.String = 'Entropy';

% Display the different entropy levels of the rules
entlev = unique(sum(H(cell2mat(prob')), 2)); 
contour(t, t, entmap, entlev, 'k-');

% Display diagonals
plot([0,1], [0,1], 'k-', 'LineWidth', 1);
plot([0,1], [1,0], 'k-', 'LineWidth', 1);

% Get theoretical generative probabilities
p = 1 - cell2mat(prob');

% Plot used probabilities
x = unique(p(:,1));
y = unique(p(:,2));

% For each condition
for i = 1:size(p,1)
    plot(p(i,1), p(i,2), 'k.', 'MarkerSize', 25);
end

% Customize the axes
axis('square'); axis('xy');
set(gca, 'XTick', [0;x;1], 'YTick', [0;y;1]);
set(gca, 'FontSize', 15, 'LineWidth', 1);

% Add some labels
xlabel('p(A|B)'); ylabel('p(B|A)');

%% DETERMINISTIC REGULARITIES
%  ==========================

% Get the rows and columns in the table
nrows = 4;
ncols = 3;

% Prepare a window
figure('Position', [257 939 630 166]); j = 1;

% For each deterministic regularity
for i = 1:numel(dr)
    
    % Create a subplot
    subplot(nrows, ncols, j);
    if i == 1, j = j + 3;
    else, j = j + 1;
    end
    
    % Display the deterministic regularity
    d = abs(10 - numel(det{i})) / 2;
    x = find(det{i} == 1) + d;
    plot(x, 1, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 15, 'LineWidth', 1); hold('on');
    text(x, ones(1,numel(x)), 'A', 'Color', 'w', 'FontWeight', 'Bold', ...
        'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle');
    x = find(det{i} == 2) + d;
    plot(x, 1, 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 15, 'LineWidth', 1);
    text(x, ones(1,numel(x)), 'B', 'Color', 'k', 'FontWeight', 'Bold', ...
        'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle');
    
    % Customize the axes
    axis([1,max(cellfun(@numel, dr)),0,2]);
    axis('off');
    set(gca, 'FontSize', 15, 'LineWidth', 1);
end

%% STIMULI
%  =======

% Define properties of the auditory stimuli
AudioFreq = 44100; % Hertz
StimFreqs = [350, 700, 1400; 500, 1000, 2000]; % Hz
StimRise  = 0.007; % ms
StimDur   = 0.05; % ms

% Choose sounds
soundsname = {'beep', 'boop'};
nSnd = size(StimFreqs, 1);

% Create the sounds
SoundMatrix = NaN(nSnd, AudioFreq*StimDur);
for iSnd = 1:nSnd
    SoundMatrix(iSnd,:) = Emergence_Stimulation_CreateSound(StimFreqs(iSnd,:), ...
        StimRise, StimDur, AudioFreq); % options for the sound
end

% Display the waveforms
figure('Position', [257 765 630 100]);
for iSnd = 1:nSnd
    subplot(2,1,iSnd);
    plot(linspace(0, StimDur, AudioFreq*StimDur) .* 1e3, SoundMatrix(iSnd,:));
    set(gca, 'YLim', [-1.25,1.25], 'YColor', 'None', 'Box', 'Off');
    if     iSnd == 1, set(gca, 'XColor', 'None');
    elseif iSnd == 2, set(gca, 'XTick', [0,StimDur*1e3], 'XTickLabel', {});
    end
end
