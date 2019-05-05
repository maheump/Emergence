% This script shows how different metrics can be used to quantify the
% similarity between posterior beliefs from the ideal observer model and
% finger trajectory from the subjects. We also investigate how taking into
% account the temporal shift between the two may help to maximize the
% similarity between the two.
% 
% Copyright (c) 2018 Maxime Maheu

% Initialization
% ~~~~~~~~~~~~~~

% Define similarity measures to test
mes = {'P', 'E', 'D', 'R', 'L', 'M'};
coord = {'C', 'B'};
nMes = numel(mes);
nCrd = numel(coord);

% Define the shifting grid
shiftgrid = 0:100; % in number of observations

% Select an example sequence
iSeq = cidx{2}(8);
iSub = 10;

% Get trajectories from an example subject and the ideal observer
subtraj = G{iSeq,iSub}.BarycCoord; % subject
modtraj = IO{iSeq,iSub}.BarycCoord; % ideal observer

% Display metric as a function of shift
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Prepare output variable for the optimal shift
optimshift = NaN(nMes,nCrd);

% Prepare a new window
figure('Position', [399.5000 300 560 800]);

% For each measure applied on each coordinate system
for i = 1:nMes
    for j = 1:nCrd
        subplot(nMes,nCrd,(nCrd*i-1)+1*(j-1));
        
        % Compute similarity measures between subjects' and IO's beliefs by
        % shifting one with respect to the other
        lab = strcat(mes{i}, coord{j});
        [out, allmet, optimshift(i,j)] = Emergence_Similarity(...
            subtraj, modtraj, lab, shiftgrid);
        
        % Plot metrics as a function of the shift
        plot(shiftgrid, allmet, 'k-', 'LineWidth', 2); hold('on');
        
        % Display shift inducing best value of the metric
        plot(shiftgrid(optimshift(i,j)), out, 'k.', 'MarkerSize', 15);
        
        % Add some text labels
        xlabel('Shift (# obs.)');
        ylabel('Metrics');
        title(lab);
    end
end

% Display trajectories of subjects and IO after optimaly shifting one with
% respect to the other
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Prepare a new window
figure('Position', [959.5000 300 560 800]);

% For each measure applied on each coordinate system
for i = 1:nMes
    for j = 1:nCrd
        subplot(nMes,nCrd,(nCrd*i-1)+1*(j-1));
        
        % Shift the two trajectories according to the best shift identified
        % previously based on that particular metric
        optimlag = shiftgrid(optimshift(i,j));
        [a, b] = Emergence_Shift(subtraj, modtraj, optimlag);
        pos = optimlag+1:N;
        
        % Display trajectories from the model and the subjects
        l = cell(1,2);
        l{1} = plot(pos, a, 'LineWidth', 2); hold('on');
        l{2} = plot(pos, b, 'LineWidth', 1/2);
        for t = 1:2, for h = 1:3, set(l{t}(h), 'Color', tricol(h,:)); end; end
        
        % Add some text labels
        xlabel('Observation number');
        ylabel('p(H|y)');
        lab = strcat(mes{i}, coord{j});
        title(lab);
    end
end
