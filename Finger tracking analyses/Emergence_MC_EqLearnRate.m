% This script shows variations of equivalent learning rates from the model.
% In particular, an increase of equivalent learning rates following
% change-points.
% 
% Copyright (c) 2020 Maxime Maheu

%% DISPLAY MEMORY DECAYS AS A FUNCTION OF LEARNING RATE
%  ====================================================

% Compute the relative contribution of sequence observation as a function
% of the (equivalent) learning rate using the delta rule
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Define parameters
no = 51; % number of osbervations
na = 20; % number of alpha values to test
alim = [10^-2.5, 1/2]; % limits of alpha values to test

% Create the grid of (equivalent) learning rate (i.e. alpha)
agrid = exp(linspace(log(alim(1)), log(alim(2)), na))';

% Prepare output variable
rc = NaN(na,no);

% Compute the relative contribution of past observations
rc(:,1) = 1;
for i = 2:no
    rc(:,i) = rc(:,i-1) - (agrid .* abs(0 - rc(:,i-1))); % delta rule
end

% Display the relative contribution of past observation depending on the
% (equivalent) learning rate considered
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Create a new window
figure('Position', [1 905 260 200]);

% Display relative contribution as a function of sequence observation for
% each value of learning rate
l = plot(1-(1:no), rc, 'LineWidth', 1.5);
col = flipud(spring(na));
for i = 1:na, set(l(i), 'Color', col(i,:)); end

% Display a colorbar
colormap(col);
cbr = colorbar;
cbr.Label.String = 'Equivalent learning rate';
caxis(alim);
set(gca,'ColorScale', 'log');

% Customize the axes
set(gca, 'Box', 'Off');
axis([-no,0,0,1]);

% Add text labels
xlabel('Observation # in the past');
ylabel('Relative contribution');

%% SHOW THE DISTRIBUTION OF EQUIVALENT LEARNING RATE USED BY THE MODEL
%  ===================================================================

% Compute distribution of (equivalent) learning rate
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Get distribution of equivalent learning rate
EqAlpha = cellfun(@(x) x.EqLearnRate, IO, 'uni', 0);
EqAlpha = mat2cell(EqAlpha, nSeq, ones(1,nSub));
EqAlpha = cellfun(@cell2mat, EqAlpha, 'uni', 0);

% Get corresponding kernel densities
nkgrid = 1000;
kgrid = linspace(alim(1), alim(2), nkgrid);
dist = cellfun(@(x) ksdensity(x, kgrid, 'BandWidth', 0.025, 'Support', ...
    [-eps,5], 'BoundaryCorrection', 'Reflection'), EqAlpha, 'uni', 0);

% Normalize such that the distributions sum to 1
dist = cell2mat(dist');
dist = dist ./ sum(dist, 2);

% Average over pseudo-subjects
avgdist = mean(dist, 1);
semdist = sem(dist, 1);

% Display the distribution of (equivalent) learning rate
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Create a new window
figure('Position', [262 905 260 200]);

% Display the kernel distribution
inflim = min(avgdist-semdist);
fill([kgrid, fliplr(kgrid)], [avgdist, repmat(inflim, 1, nkgrid)], ...
    'k'); alpha(1/10); hold('on');
plotMSEM(kgrid, avgdist, semdist, 1/2, 'k', 'k', 1.5, []);

% Display the colorbar
scatter(agrid, repmat(inflim, na,1), 100.*ones(na,1), col, 's', 'filled');

% Customize the axes
set(gca, 'Box', 'Off', 'XScale', 'Log', 'XLim', alim([1,2]), 'YScale', 'log');
ScaleAxis('y');

% Add text labels
xlabel('Equivalent learning rate (log-scale)');
ylabel('Frequency (log-scale)');

%% (EQUIVALENT) LEARNING RATES VARY, IN PARTICULAR FOLLOWING CHANGE-POINTS
%  =======================================================================

% Get changes in (equivalent) learning rate locked on the change-point
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Define grid  of osbervations around the change-point
obsgrid = -80:80; % grid of observations around the change-point to test

% Define the strength of temporal smoothing
smoothparam = 20; % in number of observations

% Prepare output variables
LockEqAlpha = NaN(numel(obsgrid),3,nSub);

% For each type of generative process
for iReg = 1:2
    
    % Get 
    genproc = cellfun(@(x) x.Gen, G(cidx{iReg},:), 'uni', 0);
    detecmask = (filter{iReg} == 3);
    
    % For each subject
    for iSub = 1:nSub
        
        % Get equivalent learning rate in the entire sequence
        eqAlpha1 = cell2mat(cellfun(@(x) x.EqLearnRate', IO(cidx{iReg},iSub), 'uni', 0));
        
        % Restrict to generative processes accurately labeled
        curgenproc = genproc(:,iSub);
        curgenproc(~detecmask(:,iSub)) = {NaN(1,N)};
        curgenproc = cell2mat(curgenproc);
        eqAlpha2 = eqAlpha1(curgenproc == 2);
            
        % Get the position of the change-point
        cp = cellfun(@(x) x.Jump, G(cidx{iReg},iSub), 'uni', 0);
        
        % Get equivalent learning rate aroung change-points (with a
        % temporal smoothing average)
        trajlocked = cellfun(@(x,y) Emergence_LockOnPoint(movmean(x',20), y, obsgrid), ...
            mat2cell(eqAlpha1, ones(numel(cidx{iReg}),1), 200), cp, 'uni', 0);
        
        % Averaged over pseudo-subjects
        LockEqAlpha(:,iReg,iSub) = mean(cell2mat(trajlocked'), 2, 'OmitNaN')';
    end
end


% Display the change in (equivalent) learning rate induced by change-points
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Create a new window
figure('Position', [523 905 260 200]);

% For each type of regularity
for iReg = 1:2
    
    % Display the dynamic of (equivalent) learning rate
    plotMSEM(obsgrid, mean(LockEqAlpha(:,iReg,:), 3), ...
                      sem( LockEqAlpha(:,iReg,:), 3), ...
        1/2, tricol(iReg,:), tricol(iReg,:), 1.5, []);
end

% Display position of change point
plot(zeros(1,2), alim([1,2]), 'k--');

% Display the colorbar
scatter(repmat(obsgrid(1),na,1), agrid, 100.*ones(na,1), col, 's', 'filled');

% Customize the axes
set(gca, 'Box', 'Off', 'YScale', 'Log', 'XLim', obsgrid([1,end]), ...
    'XTick', obsgrid(1):20:obsgrid(end), 'YLim', alim([1,2]));

% Add text labels
xlabel('Observation w.r.t. change-point');
ylabel('Equivalent learning rate (log-scale)');
