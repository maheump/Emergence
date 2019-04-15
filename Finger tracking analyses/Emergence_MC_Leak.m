% This script shows how a leaky integration (implemented as the probability
% of doing a substitution error when remembering past observations) impacts
% the probabilistic hypothesis in particular regarding the detection
% sensitivity of probabilistic regularities (versus concluding that the
% sequence is fully-stochastic) and the detection dynamics of probabilistic
% regularities.
% 
% Copyright (c) 2018 Maxime Maheu

%% INITIALIZATION
%  ==============

% Load simulations
SimuType = 'Leak';
Emergence_MC_ModelSimulations;

% Look at probabilistic regularities
iHyp = 1;

%% DETECION SENSITIVITY
%  ====================

% Compute detection sensitivity index
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% MAP generative hypothesis
[~,HypMAPreg] = cellfun(@(p) max(p(end,:)), pMgY(cidx{iHyp},:,:));
[~,HypMAPrnd] = cellfun(@(p) max(p(end,:)), pMgY(cidx{end},:,:));
signal = (HypMAPreg ~= 3);
noise  = (HypMAPrnd ~= 3);

% Get frequencies for each possible answer
H  = sum(signal == 1); % hit
M  = sum(signal == 0); % miss
CR = sum(noise  == 0); % correct rejection
FA = sum(noise  == 1); % false alarm

% Pad extreme probabilities
pad = N/(2*N);
H(H == 0) = pad;
M(M == 0) = pad;
CR(CR == 0) = pad;
FA(FA == 0) = pad;

% Compute the proportions of hits and false alarms
pH  = H ./ (H + M);
pFA = FA ./ (FA + CR);

% Compute sensitiviy index
dprime = squeeze(norminv(pH) - norminv(pFA));

% Average over subjects
m = mean(dprime);
s = sem(dprime);

% Display d' as a function of models
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

figure('Position', [1 905 250 200]);
x = cell2mat((options(2,:)));
plotMSEM(x, m, s, 1/20, 'k', 'k', 2, 1, '-');
scatter(x, m, 50, modc, 'filled');

% Customize the axes
set(gca, 'Box', 'Off', 'Layer', 'Bottom');

% Add some text labels
xlabel('Substitution error');
ylabel('Detection sensitivity');

% Display d' as a function of models
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Get most likely hypothesis at the end of each sequence under each model
[~,HypMAPreg] = cellfun(@(p) max(p(end,:)), pMgY(cidx{iHyp},:,:));

% Average proportion of each type of respose over subjects and sequences
avg = [squeeze(mean(HypMAPreg == 1, 1:2)), ...
       squeeze(mean(HypMAPreg == 2, 1:2)), ...
       squeeze(mean(HypMAPreg == 3, 1:2))];
err = [squeeze(sem(HypMAPreg == 1, 1:2)), ...
       squeeze(sem(HypMAPreg == 2, 1:2)), ...
       squeeze(sem(HypMAPreg == 3, 1:2))];

% Prepare a new window
figure('Position', [1 531 722 300]);

% For each model
for iMod = 1:nMod
    m = avg(iMod,:);
    s = err(iMod,:);

    % Display correct labeling
    lgd = NaN(1,3);
    lgd(iHyp) = bar(iMod, m(iHyp), 1, 'FaceColor', tricol(iHyp,:)); hold('on');

    % Display error bar for correct reports
    plot(repmat(iMod,1,2), m(iHyp)+s(iHyp).*[-1,1], 'k-');

    % Display incorrect labeling
    idx = setdiff(1:3, iHyp); % other responses than the real one
    [~,nd] = sort(m(idx), 2, 'descend');
    idx = idx(nd);
    h = bar([repmat(iMod,1,2); repmat(iMod-1,1,2)], ...
        [-m(idx); NaN(1,2)], 1, 'Stacked');
    for i = 1:2, set(h(i), 'FaceColor', tricol(idx(i),:)); end

    % Display error bars for wrong reports
    cm = cumsum(-m(idx));
    for j = 1:2
        plot(repmat(iMod,1,2), cm(j)+s(j).*[-1,1], 'k-');
    end
end

% Display chance level
plot([0,nMod+1],  ones(1,2)./3, '--', 'Color', g);
plot([0,nMod+1], -ones(1,2)./3, '--', 'Color', g);

% Customize the axes
axis([1/2,nMod+1/2,-1,1]);
set(gca, 'XTick', 1:nMod, 'XTickLabel', cellfun(@num2str, options(2,:), ...
    'uni', 0), 'XTickLabelRotation', 90, 'Box', 'Off');

% Add some text labels
xlabel('Substitution error');
ylabel('Proportion of responses');

% Compute belief differentce between regularities that were detected by
% subjects versus those that were missed
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Get regularities detected by subjects
detecmask = (filter{iHyp} == 1 | filter{iHyp} == 3);

% Get the belief in the correct hypothesis at the end of the sequence from
% the different models
ultimbel = cellfun(@(p) p(end,iHyp), pMgY(cidx{iHyp},:,:));

% Compute belief difference for each model and each subject
difdet = NaN(nSub,nMod);
for iMod = 1:nMod
    for iSub = 1:nSub
        modbel = ultimbel(:,iSub,iMod);
        subdet = detecmask(:,iSub);
        difdet(iSub,iMod) = mean(modbel(subdet)) - mean(modbel(~subdet));
    end
end

% Average over subjects
m = mean(difdet);
s = sem(difdet);

% Display belief difference as a function of substitution error
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Prepare a new figure
figure('Position', [252 905 250 200]);

% Display difference in beliefs as a function of substitution error
plotMSEM(x, m, s, 1/20, 'k', 'k', 2, 1, '-');
scatter(x, m, 50, modc, 'filled');

% Customize the axes
set(gca, 'Box', 'Off', 'Layer', 'Bottom');

% Add some text labels
xlabel('Substitution error');
ylabel('Belief difference btw un/detec reg');

%% DETECION DYNAMICS
%  =================

% Lock trajectories
% ~~~~~~~~~~~~~~~~~

% Get change point positions
cp = cellfun(@(x) x.Jump+1/2, G(cidx{iHyp},:), 'uni', 0);

% Get detection point positions
dp = cellfun(@(p,c) c + Emergence_FindDetecPoint(p(c:end,iHyp)), ...
    pMgY(cidx{iHyp},:,:), repmat(cp, [1,1,nMod]), 'uni', 0);

% Lock trajectories on detection point
win = 0:60;
traj = cellfun(@(p,d) Emergence_LockOnPoint(p,d,win), ...
    pMgY(cidx{iHyp},:,:), repmat(cp, [1,1,nMod]), 'uni', 0);

% Remove sequences entailing regularities that were missed by subjects and
% for which we could not find a detection point in subjects or the full
% ideal observer model. That allows to directly compare the resulting
% trajectories with thoses that are generated by Emergence_FTA_Dynamics1
% (they are therefore based *exactly* on the same sequences)
detecmask = (filter{iHyp} == 3);
traj(~repmat(detecmask, [1,1,nMod])) = {NaN(numel(win),3)};

% Transform the cell matrix into a double matrix
traj = cell2mat(reshape(traj, [1,1,size(traj)]));

% Average over sequences
seqavg = mean(traj(:,iHyp,:,:,:), 3, 'OmitNaN');

% Average over subjects
m = squeeze(mean(seqavg, 4));
s = squeeze(sem( seqavg, 4));

% Display the trajectories
% ~~~~~~~~~~~~~~~~~~~~~~~~

% Prepare a new window
figure('Position', [503 905 220 200]);

% Display detection threhsold
plot(win([1,end]), ones(1,2)./2, '--', 'Color', g); hold('on');

% Display trajectory under each substitution error parameter
for iMod = fliplr(1:nMod)
    plotMSEM(win, m(:,iMod), s(:,iMod), 1/20, ...
        modc(iMod,:), modc(iMod,:));
end

% Customize the axes
axis([win([1,end]),0,1]);
set(gca, 'Box', 'Off');

% Add some text labels
xlabel('Position w.r.t. change point');
ylabel('Posterior probability');
    