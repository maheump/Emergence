% This scripts displays error in the (pseudo-)deterministic hypothesis
% between ideal observer models and the subjects regarding two main aspects
% in which they might fail: the likelihood of that hypothesis in
% fully-stochastic parts (i.e. deterministic false alarms) and the
% staircaseness of the detection dynamics in deterministic parts.
% 
% Copyright (c) 2018 Maxime Maheu

% Load simulations
% ~~~~~~~~~~~~~~~~

SimuType = 'PseudoDeterministic'; % 'PseudoDeterministic' or 'BiasedPseudoDeterministic'
Emergence_MC_ModelSimulations;

% Whether to center on the error 
centeronD = true;

% Axis #1: belief in p(Hd|y) in fully-stochastic parts of sequences
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Look at fully-stochastic parts of sequences
idxtrimap = Emergence_SelectFullyStochSeq(G, filter, 1);

% Get average positions in the (pseudo-)deterministic hypothesis from the
% ideal observer model
iopos = cellfun(@(x,i) x(i,2), pMgY, repmat(idxtrimap, [1,1,nMod]), 'UniformOutput', 0);
iopos = squeeze(cellfun(@mean, iopos));

% Get average positions in the (pseudo-)deterministic hypothesis from the
% subjects
subpos = cellfun(@(x,i) x.BarycCoord(i,2), G, idxtrimap, 'UniformOutput', 0);
subpos = squeeze(cellfun(@mean, subpos));

% Compute the error between the two
MSEx = (subpos - iopos) .^ 2;

% Center on the error of the truly deterministic hypothesis
if centeronD, MSEx = MSEx - MSEx(:,:,end); end

% Average over sequences
MSEx = squeeze(mean(MSEx, 1));

% Perform paired t-tests
[h,pval,tci,stats] = ttest(MSEx);
Emergence_PrintTstats(pval,tci,stats);

% Axis #2: abruptness of detection of deterministic regularities
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Get positions of the change point
cp = cellfun(@(x) x.Jump+1/2, G(cidx{2},:), 'UniformOutput', 0);

% Define the measure of abruptness
staircaseness = @(p) mean(abs(diff(p,2,1)));

% Measure abruptness of detection in the ideal observer model
ioabr = cellfun(@(p,i) staircaseness(p(i:end,2)), pMgY(cidx{2},:,:), repmat(cp, [1 1 nMod]));

% Measure abruptness of detection in the subjects
subabr = cellfun(@(p,i) staircaseness(p.BarycCoord(i:end,2)), G(cidx{2},:,:), cp);

% Remove deterministic regularities that were not identified by subjects
detecmask = (filter{2} == 3);
subabr(~detecmask) = NaN;

% Compute the error between the two
MSEy = (subabr - ioabr) .^ 2;

% Center on the error of the truly deterministic hypothesis
if centeronD, MSEy = MSEy - MSEy(:,:,end); end

% Average over sequences
MSEy = squeeze(mean(MSEy, 1, 'OmitNaN'));

% Perform paired t-tests
[h,pval,tci,stats] = ttest(MSEy);
Emergence_PrintTstats(pval,tci,stats);

% Display the results
% ~~~~~~~~~~~~~~~~~~~

% Prepare a new window
figure('Position', [1 885 300 220]); hold('on');

% Add text labels
xlabel('Order of transitions');
title(['Error in $p(\mathcal{H}_{\mathrm{D''}}|y)$ relative to ', ...
    '$\mathcal{H}_{\mathrm{D}}$'], 'Interpreter', 'LaTeX');

% Axis #1
% ~~~~~~~

% Choose which models to display results from
modidx = 2:nMod-1;

% Get data to plot aand average over subjects
x = modidx-1/6;
m = mean(MSEx(:,modidx));
s = sem(MSEx(:,modidx));

% Create the axis
yyaxis('left');
set(gca, 'YColor', tricol(3,:));
ylabel('Average in fully-stochastic parts');

% Display error made by the fully deterministic hypothesis
plotMSEM([0,nMod+1], repmat(mean(MSEx(:,end)),1,2), ...
    repmat(sem(MSEx(:,end)),1,2), 1/2, tricol(3,:));

% Display error bars
plot(repmat(x, [2,1]), m+s.*[-1;1], 'k-');

% Display average
plot(x, m, 'k-', 'LineWidth', 3);
plot(x, m, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', tricol(3,:));

% Customize the axis
ylim([0,0.02]);
ylim([-diff(ylim)/10, max(ylim)]);

% Axis #2
% ~~~~~~~

% Get data to plot aand average over subjects
x = modidx+1/6;
m = mean(MSEy(:,modidx));
s = sem(MSEy(:,modidx));

% Create the axis
yyaxis('right');
set(gca, 'YColor', tricol(2,:));
ylabel('Abruptness in deterministic parts');

% Display error made by the fully deterministic hypothesis
plotMSEM([0,nMod+1], repmat(mean(MSEy(:,end)),1,2), ...
    repmat(sem(MSEy(:,end)),1,2), 1/2, tricol(2,:));

% Display error bars
plot(repmat(x, [2,1]), m+s.*[-1;1], 'k-');

% Display average
plot(x, m, 'k-', 'LineWidth', 3);
plot(x, m, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', tricol(2,:));

% Customize the aces
axis([min(modidx)-1,max(modidx)+1,-diff(ylim)/10, max(ylim)]);
set(gca, 'XTick', modidx);
