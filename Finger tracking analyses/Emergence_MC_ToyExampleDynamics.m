% This script shows example metrics use to quantify different aspects of
% the trajectories
% 
% Copyright (c) 2018 Maxime Maheu

% Load simulations
SimuType = 'PseudoDeterministic';
Emergence_MC_ModelSimulations;

% Define functions
abruptness = @(p) mean(diff(p,1,1));
detecpoint = @Emergence_FindDetecPoint;
staircaseness = @(p) mean(abs(diff(p,2,1)));

% Get change point positions
cp = cell(1,2);
for iHyp = 1:2
    cp{iHyp} = cellfun(@(x) x.Jump+1/2, G(cidx{iHyp},:), 'UniformOutput', 0);
end

% Define colormaps to use
redcmap = flipud(cbrewer2('Reds', 4));
bluecmap = flipud(cbrewer2('Blues', 4));

% Define simulations
simu = {abruptness,         abruptness,         detecpoint,         staircaseness; ...
        [1 2],              [1 1 1],            [2 2 2],            [2 2]; ...
        [3 1],              [1 2 3],            [3 6 9],            [8 8]; ...
        [nMod nMod],        [nMod nMod nMod]    [nMod nMod nMod]    [nMod 5]; ...
        tricol(1:2,:)       bluecmap(3:-1:1,:)  redcmap(1:3,:)      [tricol(2,:); 0.4588 0.4196 0.6941]; ...
        'Belief update',    'Belief update',    'Delay',            'Staircaseness'};
nSimu = size(simu,2);

% Define window for the trajectories
win = [0,50];
nWin = diff(win)+1;

% Prepare a new window
figure('Position', [811 272 300 612]);

% For each simulation
for iSimu = 1:nSimu
    
    % Define regularities to explore
    iHyp = simu{2,iSimu};
    iReg = simu{3,iSimu};
    iMod = simu{4,iSimu};
    nExp = numel(iReg);
    
    % Display example trajectories
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    % Lock trajectories onto the change point
    traj = NaN(nWin,nSub,nExp);
    for i = 1:nExp
        traj(:,:,i) = cell2mat(cellfun(@(p,c) ...
            Emergence_LockOnPoint(p(:,iHyp(i)), c, win), ...
            pMgY(cidx{iHyp(i)}(iReg(i)),:,iMod(i)), cp{iHyp(i)}(iReg(i),:), ...
            'UniformOutput', 0));
    end
    
    % Display trajectories
    col = simu{end-1,iSimu};
    subplot(nSimu,3,(1:2)+3*(iSimu-1));
    for i = 1:nExp
        plotMSEM(win(1):win(end), mean(traj(:,:,i), 2), sem(traj(:,:,i), 2), ...
            1/5, col(i,:), col(i,:), 3);
    end
    
    % Customize the axes
    axis([win,0,1]);
    set(gca, 'YTick', 0:1/4:1, 'Box', 'Off');
    
    % Add some text labels
    xlabel('Position w.r.t. change point');
    ylabel('p(H|y)');
    
    % Display corresponding metrics
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    % Compute metrics on the trajectories
    mattraj = mat2cell(traj, nWin, ones(1,nSub,1), ones(1,1,nExp));
    fun = simu{1,iSimu};
    index = squeeze(cellfun(fun, mattraj));
    
    % Display metrics
    subplot(nSimu,3,3+3*(iSimu-1));
    Emergence_PlotSubGp(index, col);
    
    % Customize the axes
    xlim([0,nExp+1]);
    set(gca, 'XTick', [], 'XColor', 'none');
    
    % Add some text labels
    ylabel(simu{end,iSimu});
end
