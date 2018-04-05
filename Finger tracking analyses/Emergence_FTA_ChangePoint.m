% This script analyses inferred positions of the change point. First, we
% derive predictions from the ideal observer and show that beliefs
% regarding the position of the change point are more precise in the case
% of the emergence of a deterministic compared to a probabilistic
% regularity. Second, we show that in both cases, estimates of change
% point's position are correlated with the real positions of the change
% point. Third, we show that the confidence in that estimate is stronger
% in the case of deterministic than probabilistic regularities.
%
% Copyright (c) 2018 Maxime Maheu

%% PREDICTIONS FROM THE IDEAL OBSERVER
%  ===================================

% Choose the size of the window to look into
winwidth = 50;

% Get the beliefs from the IO regarding the position of the change point
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Get a vector of position centered around the change point
x = -winwidth/2:winwidth/2;

% Prepare output variables
cppost = NaN(nSub, winwidth+1, 2);
varcppost = NaN(nSub, 2);

% For each type of regularity
for iHyp = 1:2
    
    % Get the position of the change point
    cp = cellfun(@(x) x.Jump-1/2, G(cidx{iHyp},:), 'UniformOutput', 0);
    
    % Get beliefs about the change point's position at the end of each
    % sequence
    lab = sprintf('pJkgYMs%s', lower(proclab{iHyp}(1)));
    bel = cellfun(@(x) x.(lab)(end,:), IO(cidx{iHyp},:), 'UniformOutput', 0);
    
    % Look around the change point's position in the window defined earlier
    belwrtcp = cellfun(@(b,c) b(c-winwidth/2:c+winwidth/2), bel, cp, 'UniformOutput', 0);
    
    % Average over sequences for each subject
    belwrtcp = mat2cell(belwrtcp, numel(cidx{iHyp}), ones(1,nSub));
    belwrtcp = cellfun(@(x) mean(cell2mat(x), 1), belwrtcp, 'UniformOutput', 0);
    cppost(:,:,iHyp) = cell2mat(belwrtcp');
    
    % Measure the precision of these posterior distributions
    varcppost(:,iHyp) = cellfun(@(y) 1/std(x,y), belwrtcp);
end

% Display precision of posterior distributions
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Get the length of the sequence
N = S{1}.Nstims;

% Prepare a new window
figure('Position', [1 906 340 200]);

% For each type of regularity
lgd = NaN(1,2);
for iHyp = fliplr(1:2)
    
    % Average posterior distributions
    d = cppost(:,:,iHyp);
    m = mean(d,1);
    
    % Display the posterior distribution of change point's positions
    % centered around the real position of the change point
    lgd(iHyp) = fill([x(1), x, x(end)], [0, m, 0], 'k', 'FaceColor', ...
        tricol(iHyp,:), 'EdgeColor', 'None', 'FaceAlpha', 1/3); hold('on');
    plot(x, m, '-', 'Color', tricol(iHyp,:), 'LineWidth', 2)
    
    % Perform a t-test against the uniform (prior) distribution at each
    % position around the true change point's position
    h = ttest(cppost(:,:,iHyp), 1/N, 'tail', 'right');
    plot(x(logical(h)), repmat(0.34+0.005*(iHyp-1),1,sum(h)), ...
        '-', 'Color', tricol(iHyp,:), 'LineWidth', 2);
end

% Display the position of the change point (i.e. 0)
plot(zeros(1,2), ylim, 'k--');

% Display the uniform scenario over change point's position
plot(x([1,end]), repmat(1/N,1,2), 'k:');

% Customize the axes
xlim(x([1,end]));
set(gca, 'Box', 'Off');

% Add some text labels
legend(lgd, proclab(1:2), 'Location', 'NorthWest');
xlabel('Position w.r.t. change point');
ylabel('p(j_k|H_i,y)');

% Save the figure
save2pdf('figs/F_CP_PredCorr.pdf');

% Display precision of posterior distributions
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Prepare a new window
figure('Position', [342 906 120 200]);

% Plot individual
plot((1:2)'+xshift', varcppost', '-', 'Color', g, 'LineWidth', 1/2); hold('on');
plot(1:2, mean(varcppost), 'k-', 'LineWidth', 3);

% For each type of regularity
for iHyp = 1:2
    
    % Get the precision levels of posterior distribution over change
    % point's positions
    d = varcppost(:,iHyp);
    m = mean(d);
    e = m+sem(d)*[-1,1];
    
    % Display individual data points
    scatter(iHyp+xshift, d, 'o', 'MarkerEdgeColor', 'None', ...
        'MarkerFaceColor', tricol(iHyp,:), 'MarkerFaceAlpha', 1/3);
    
    % Display group-average and related error
    plot(repmat(iHyp,1,2), e, 'k-');
    plot(iHyp, m, 'ko', 'MarkerFaceColor', tricol(iHyp,:), 'MarkerSize', 10);
end

% Customize the axes
xlim([0,3]);
set(gca, 'XTick', [], 'XColor', 'None', 'Box', 'Off');

% Add some text labels
ylabel('1/std');

% Save the figure
save2pdf('figs/F_CP_PredConf.pdf');

%% CORRELATION BETWEEN ESTIMATED AND REAL CHANGE POINT'S POSITION
%  ==============================================================

% Look at the position of the CHANGE point or of the DETECTION point?
xp = 1; % 1 for change point and 2 for detection point
ptname = {'change', 'detection'};

% Define the number of bins to use for the correlation at the group-level
nBin = 8;

% Correlate real and inferred positions of the change point
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Create the bins in which the positions will be averaged
pgrid = linspace(0, 100, nBin+1);
allcp = cellfun(@(x) x.Jump-1/2, D(cat(2, cidx{1:2}),:));
pct = percentile(allcp, pgrid);

% Prepare output variable
corcoef = NaN(nSub,2);
dotsm = NaN(nBin,2);
dotss = NaN(nBin,2);

% For each type of regularity
for iHyp = 1:2
    
    % Get the number of sequences with this type of regularity
    nR = numel(cidx{iHyp});
    
    % Get objective and subjective change/detection point's positions
    cp = cellfun(@(x) x.Jump-1/2, D(cidx{iHyp},:));
    if xp == 1 % for question about change point's position
        op = cp;
        sp = cellfun(@(x) x.Questions(3), D(cidx{iHyp},:));
    elseif xp == 2 % for question about detection point's position
        op = cellfun(@(x,c) c + find(x.BarycCoord(c:end,iHyp) > 1/2, 1, 'first'), ...
            D(cidx{iHyp},:), num2cell(cp), 'UniformOutput', 0);
        op(cellfun(@isempty, op)) = {NaN};        
        op = cell2mat(op);
        sp = cellfun(@(x) x.Questions(5), D(cidx{iHyp},:));
    end
    
    % Remove sequences in which regularities were not detected
    detecmask = logical(cellfun(@(x) x.Questions(2) == iHyp, G(cidx{iHyp},:)));
    op(~detecmask) = NaN;
    sp(~detecmask) = NaN;
    
    % Subject-specific regression between
    corcoef(:,iHyp) = cellfun(@(y,x) nancorr(y, x), ...
        mat2cell(sp, nR, ones(nSub,1)), ... % inferred change point's position
        mat2cell(op, nR, ones(nSub,1)), ... % true change point's position
        'UniformOutput', 1);
    
    % For each bin
    for iBin = 1:nBin
        
        % Find the data points that are in the range of the current bin
        idx = op > pct(iBin) & op <= pct(iBin+1);
        
        % Average over sequences
        dotsm(iBin,1,iHyp) = mean(op(idx), 'omitnan');
        dotsm(iBin,2,iHyp) = mean(sp(idx), 'omitnan');
        
        % Compute the SEM over sequences
        dotss(iBin,1,iHyp) = sem(op(idx));
        dotss(iBin,2,iHyp) = sem(sp(idx));
    end
end

% Display group-level correlation between real and estimated change point's position
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Prepare the figure
figure('Position', [1 632 340 200]);
subplot(1,3,1:2);

% Display the generative process underlying the position of the change points
bot = 60;
if xp == 1
    x = linspace(1, N, N*10);
    prior = normpdf(x, 100, 15);
    fill([0,x,N], [0,(prior-min(prior))*500+bot,bot], 'k-', ...
        'FaceColor', g, 'EdgeColor', 'None'); hold('on');
end

% Display the identity relation (i.e. a diagonal)
plot([bot, 200-bot], [bot, 200-bot], 'k--');

% For each type of regularity
for iHyp = 1:2
    
    % Shortcuts to avoid errors
    ocp_m = dotsm(:,1,iHyp);
    scp_m = dotsm(:,2,iHyp);
    ocp_s = dotss(:,1,iHyp);
    scp_s = dotss(:,2,iHyp);
    
    % Regression on binned data between objective and inferred position
    B = regress(scp_m, [ones(nBin,1), ocp_m]);
    [top_int, bot_int, xval] = regression_line_ci(0.05, B, ocp_m, scp_m, 100);
    fill([xval, fliplr(xval)], [top_int, fliplr(bot_int)], 'k', ...
        'EdgeColor', 'None', 'FaceColor', tricol(iHyp,:), 'FaceAlpha', 1/3);
    plot(xval, B(1)+B(2).*xval, '-', 'Color', tricol(iHyp,:), 'LineWidth', 3);
    
    % For each bin
    for iBin = 1:nBin
        
        % Display horizontal and vertical error bars
        plot(repmat(ocp_m(iBin),1,2), scp_m(iBin)+scp_s(iBin)*[-1,1], 'k-');
        plot(ocp_m(iBin)+ocp_s(iBin)*[-1,1], repmat(scp_m(iBin),1,2), 'k-');
        
        % Display the average value in that bin
        plot(ocp_m(iBin), scp_m(iBin), 'ko', 'MarkerFaceColor', ...
            tricol(iHyp,:), 'MarkerSize', 7);
    end
    
    % Compute the correlation coefficient on binned data
    rho = nancorr(scp_m, ocp_m);
    text(xval(end), B(1)+B(2)*xval(end), sprintf(' r = %0.2f', rho), ...
        'HorizontalAlignment', 'Left', 'VerticalAlignment', 'Middle', ...
        'Color', tricol(iHyp,:));
end

% Customize the axes
axis('square');
axis(repmat([bot, 200-bot], 1, 2));
set(gca, 'XTick', get(gca, 'YTick'), 'Box', 'Off');

% Add some labels
xlabel(sprintf('Objective %s point''s position', ptname{xp}));
ylabel(sprintf('Inferred %s point''s position', ptname{xp}));

% Display individual correlation coefficient from the correlation between
% real and estimated change point's position
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Display the scenario in which no correlation is observed
subplot(1,3,3);
plot([0,3], zeros(1,2), 'k--'); hold('on');

% Display difference in correlation coefficients between the two types of regularity
plot((1:2)'+xshift', corcoef', 'Color', g, 'LineWidth', 1/2);
plot(1:2, mean(corcoef), 'k-', 'LineWidth', 3);

% For each type of regularity
for iHyp = 1:2
    
    % Get distribution of correlation coefficients over subjects
    d = corcoef(:,iHyp);
    m = mean(d);
    e = m+sem(d)*[-1,1];
    
    % Display individual correlation coefficients
    scatter(iHyp+xshift, d, 'MarkerFaceColor', tricol(iHyp,:), ...
        'MarkerFaceAlpha', 1/3, 'MarkerEdgeColor', 'None');
    
    % Display group-average and related error
    plot(repmat(iHyp,1,2), e, 'k-');
    plot(iHyp, m, 'ko', 'MarkerFaceColor', tricol(iHyp,:), 'MarkerSize', 10);
end

% Perform a paired t-test between correlation coefficients
[~,pval,tci,stats] = ttest(diff(corcoef, 1, 2)); % between regularities
disptstats(pval,tci,stats);
[~,pval,tci,stats] = ttest(corcoef - 1); % against ideal scenario
disptstats(pval,tci,stats);
[~,pval,tci,stats] = ttest(corcoef); % against an absence of correlation
disptstats(pval,tci,stats);

% Customize the axes
axis([0,3,-1,1]);
ylim([-max(abs(ylim)),max(abs(ylim))]);
set(gca, 'XColor', 'None', 'Box', 'Off');

% Save the figure
if isfield(D{1}, 'Seq'), save2pdf('figs/F_CP_CorrS.pdf');
else, save2pdf('figs/F_CP_CorrIO.pdf');
end

%% CONFIDENCE IN THE ESTIMATE
%  ==========================

% Get confidence in the change point's estimation question
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Get confidence levels
if     xp == 1, conf = cellfun(@(x) x.Questions(4), D);
elseif xp == 2, conf = cellfun(@(x) x.Questions(6), D);
end

% Average over regularities for each condition (proba. & deter.)
avgConf = NaN(nSub,2);
for iHyp = 1:2, avgConf(:,iHyp) = nanmean(conf(cidx{iHyp},:),1)'; end
if any(avgConf(:) > 1), avgConf = avgConf ./ 100; end

% Display reported confidence in the estimation for the two types of regularities
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%  Prepare the window
figure('Position', [342 632 120 200]);

% Display difference in confidence ratings between the two types of regularity
plot((1:2)'+xshift', avgConf', 'Color', g, 'LineWidth', 1/2); hold('on');
plot(1:2, mean(avgConf), 'k-', 'LineWidth', 3);

% For each type of regularity
for iHyp = 1:2
    
    % Get distribution of confidence rations over subjects
    d = avgConf(:,iHyp);
    m = mean(d);
    e = m+sem(d)*[-1,1];
    
    % Display individual confidence ratings
    scatter(iHyp+xshift, d, 'MarkerFaceColor', tricol(iHyp,:), ...
        'MarkerFaceAlpha', 1/3, 'MarkerEdgeColor', 'None');
    
    % Display group-average and related error
    plot(repmat(iHyp,1,2), e, 'k-');
    plot(iHyp, m, 'ko', 'MarkerFaceColor', tricol(iHyp,:), 'MarkerSize', 10);
end

% Perform a paired t-test between confidence ratings
[h,pval,ci,stats] = ttest(diff(avgConf,1,2));
disptstats(pval,tci,stats);

% Customize the axes
axis([0,3,0,1]);
set(gca, 'XColor', 'None', 'Box', 'Off');

% Add some labels
ylabel('Confidence in the estimate');

% Save the figure
if isfield(D{1}, 'Seq'), save2pdf('figs/F_CP_ConfS.pdf');
else, save2pdf('figs/F_CP_ConfIO.pdf');
end
