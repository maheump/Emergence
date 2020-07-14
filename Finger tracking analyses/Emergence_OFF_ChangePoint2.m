% This script analyses inferred positions of the change point. We show
% that for both probabilistic and deterministic regularities, estimates of
% change point's position are correlated with the real positions of the
% change point. Moreover, we show that the confidence in that estimate is
% higher in the case of deterministic than probabilistic regularities.
%
% Copyright (c) 2020 Maxime Maheu

%% DISTRIBUTION OF CHANGE-POINTS
%  =============================

% Compute the theoretical distribution over a change point (i.e. a Gaussian
% distribution)
obs = 1:1:N;
pCP = normpdf(obs, S{1}.LenNoise.avg, S{1}.LenNoise.std);
pCP(obs < S{1}.LenNoise.clip(1)) = 0;
pCP(obs > S{1}.LenNoise.clip(2)) = 0;

% Prepare a new figure
figure('Position', [463 906 220 200]);

% Display the real distribution of change points
for iReg = 1:2
    
    % Get empirical distribution of change points
    cp = cellfun(@(x) x.Jump, D(cidx{iReg},:));
    nCP = histcounts(cp, [obs(1)-1,obs(2:end),obs(end)+1]);
    nCP(nCP == 0) = NaN;
    
    % Display properties of the generative statistics
    subplot(2,1,iReg); hold('on');
    m = S{1}.LenNoise.avg;
    plot(repmat(m, [1,2]), [0,max(pCP)], 'k-');
    c = S{1}.LenNoise.clip;
    plot(repmat(c, [2,1]), repmat([0;max(pCP)], [1,2]), 'k-');
    s = S{1}.LenNoise.std;
    plot(m+s.*[-1,1], repmat(max(pCP), [1,2]), 'k-');
    
    % Display theoretical distribution of change points
    plot(obs, pCP, 'Color', tricol(iReg,:), 'LineWidth', 2);
    
    % Display empirical distribution of change points
    scatter(obs, zeros(1,N), nCP.*10, ...
        tricol(iReg,:), 'filled', 'MarkerEdgeColor', 'k');
    
    % Customize the axes
    xlim(obs([1,end]));
    set(gca, 'Box', 'Off', 'Layer', 'Bottom', 'YTick', [], ...
        'XTick', [1, get(gca, 'XTick')]);
    
    % Add some text labels
    xlabel('Observation number');
    ylabel('Density');
end

%% CORRELATION BETWEEN ESTIMATED AND REAL CHANGE POINT'S POSITION
%  ==============================================================

% Define the number of bins to use for the correlation at the group-level
nBin = 7;

% Correlate real and inferred positions of the change point
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Create the bins in which the positions will be averaged
pgrid = linspace(0, 100, nBin+1);
allcp = cellfun(@(x) x.Jump, D(cat(2, cidx{1:2}),:));
pct = prctile(allcp(:), pgrid);

% Prepare output variable
corrcoef = NaN(nSub,2);
cpbinsm  = NaN(nBin,2);
cpbinss  = NaN(nBin,2);

% For each type of regularity
for iReg = 1:2
    
    % Get the number of sequences with this type of regularity
    nR = numel(cidx{iReg});
    
    % Get objective and subjective change/detection point's positions
    cp = cellfun(@(x) x.Jump, D(cidx{iReg},:));
    sp = cellfun(@(x) x.Questions(3), D(cidx{iReg},:));
    
    % Remove sequences in which regularities were not detected
    detecmask = (filter{iReg} == 1 | filter{iReg} == 3);
    cp(~detecmask) = NaN;
    sp(~detecmask) = NaN;
    
    % Subject-specific regression
    corrcoef(:,iReg) = cellfun(@(toexplain,explainingvar) ...
        Emergence_Regress(toexplain, explainingvar, 'CC', 'r'), ...
        mat2cell(sp, nR, ones(nSub,1)), ... % inferred change point's position
        mat2cell(cp, nR, ones(nSub,1)));    % true change point's position
    
    % Find the data points that correspond to each bin
    [~,~,idx] = histcounts(cp, pct);
    
    % For each bin
    for iBin = 1:nBin
        
        % Average over sequences
        cpbinsm(iBin,1,iReg) = mean(cp(idx == iBin), 'OmitNaN');
        cpbinsm(iBin,2,iReg) = mean(sp(idx == iBin), 'OmitNaN');
        
        % Compute the SEM over sequences
        cpbinss(iBin,1,iReg) = sem(cp(idx == iBin));
        cpbinss(iBin,2,iReg) = sem(sp(idx == iBin));
    end
end

% Display group-level correlation between real and estimated change point's position
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Prepare a new figure
figure('Position', [684 905 220 200]);

% Display the generative process underlying the position of the change points
bot = 60;
x = linspace(1, N, N*10);
prior = normpdf(x, 100, 15);
fill([0,x,N], [0,(prior-min(prior))*500+bot,bot], 'k-', ...
    'FaceColor', g, 'EdgeColor', 'None'); hold('on');

% Display the identity relation (i.e. a diagonal)
plot([bot, 200-bot], [bot, 200-bot], '-', 'Color', g);
text(bot, bot, '     Identity', 'Color', g, 'HorizontalAlignment', 'Left', ...
    'VerticalAlignment', 'Bottom', 'Rotation', 45);

% For each type of regularity
for iReg = 1:2
    
    % Shortcuts to avoid errors
    ocp_m = cpbinsm(:,1,iReg)';
    scp_m = cpbinsm(:,2,iReg)';
    ocp_s = cpbinss(:,1,iReg)';
    scp_s = cpbinss(:,2,iReg)';
    
    % Group-level regression between (binned) objective and subjective
    % positions of the change point
    beta    = Emergence_Regress(scp_m, ocp_m, 'OLS', {'beta0', 'beta1'});
    confint = Emergence_Regress(scp_m, ocp_m, 'OLS', 'confint');
    xval    = Emergence_Regress(scp_m, ocp_m, 'OLS', 'confintx');
    
    % Regression on binned data between objective and inferred position
    fill([xval, fliplr(xval)], [confint(1,:), fliplr(confint(2,:))], 'k', ...
        'EdgeColor', 'None', 'FaceColor', tricol(iReg,:), 'FaceAlpha', 1/3);
    plot(xval, beta(1)+beta(2).*xval, '-', 'Color', tricol(iReg,:), 'LineWidth', 3);
    
    % Display horizontal and vertical error bars
    plot(repmat(ocp_m, [2,1]), scp_m+scp_s.*[-1;1], 'k-');
    plot(ocp_m+ocp_s.*[-1;1], repmat(scp_m, [2,1]), 'k-');
    
    % Display the average value in that bin
    plot(ocp_m, scp_m, 'ko', 'MarkerFaceColor', tricol(iReg,:), 'MarkerSize', 7);
end

% Customize the axes
axis('square');
axis(repmat([bot, N-bot], 1, 2));
set(gca, 'XTick', get(gca, 'YTick'), 'Box', 'Off');

% Add some labels
xlabel('True change point''s position');
ylabel('Inferred changes point''s position');

% Display the real distribution of change points
yyaxis('right');
dist = linspace(bot, N-bot, 13);
[n,x] = histcounts(allcp(:), dist, 'Normalization', 'Probability');
bar((x(1:end-1) + x(2:end)) / 2, n, 'EdgeColor', 'None', 'FaceColor', 'k', 'FaceAlpha', 1/4);
set(gca, 'YColor', 'None');
axis([x([1,end]), 0, 1]);

% Display individual correlation coefficient from the correlation between
% real and estimated change point's position
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Prepare a new figure
figure('Position', [905 905 120 200]);

% Display difference in correlation coefficients between the two types of regularity
Emergence_PlotSubGp(corrcoef, tricol(1:2,:));

% Perform a paired t-test between correlation coefficients
[~,pval,tci,stats] = ttest(diff(corrcoef, 1, 2)); % between regularities
Emergence_PrintTstats(pval,tci,stats);
[~,pval,tci,stats] = ttest(corrcoef - 1); % against ideal scenario
Emergence_PrintTstats(pval,tci,stats);
[~,pval,tci,stats] = ttest(corrcoef); % against an absence of correlation
Emergence_PrintTstats(pval,tci,stats);

% Customize the axes
axis([0,3,-1,1]);
set(gca, 'XColor', 'None', 'Box', 'Off');

% Display whether the difference is significant or not
Emergence_DispStatTest(corrcoef);

% Add text labels
ylabel('Correlation coefficient');

%% CONFIDENCE IN THE ESTIMATE
%  ==========================

% Get confidence in the change point's estimation question
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Get confidence levels
conf = cellfun(@(x) x.Questions(4), D) ./ 100;

% Remove sequences in which regularities were not detected
detecmask = (filter{iReg} == 1 | filter{iReg} == 3);
conf(~detecmask) = NaN;

% Average over regularities for each condition (proba. & deter.)
avgConf = NaN(nSub,2);
for iReg = 1:2, avgConf(:,iReg) = mean(conf(cidx{iReg},:), 1, 'OmitNaN')'; end

% Prepare the output variable
prec = 101;
xgrid = linspace(0, 1, prec);
confkern = NaN(2,prec,nSub);

% Measure a kernel density from subjects' confidence ratings
for iReg = 1:2
    for iSub = 1:nSub
        dist = conf(cidx{iReg},iSub);
        [confkern(iReg,:,iSub), u] = ksdensity(dist, ... % which distribution to plot
            xgrid, ...                                   % grid of confidence levels
            'BandWidth',            0.15, ...            % bandwidth of the kernel smoothing window 
            'Support',              [0-eps,1+eps], ...   % restrict the kernel to a certain range of values
            'BoundaryCorrection',   'Reflection');       % type of correction for the boundaries
    end
end

% Normalize the distribution such that it can be directly compared
% between both conditions (because there are not necessarily the same
% number of accurately detected regularities in both conditions).
confkern = confkern ./ sum(confkern, 2);

% Mirror the second distribution
confkern(1,:) = -confkern(1,:);

% Get average confidence for each type of regularity
avgConf = NaN(nSub,2);
for iReg = 1:2
    avgConf(:,iReg) = mean(conf(cidx{iReg},:), 1, 'OmitNaN');
end

% Display the entire distribution (over subjects and sequences) of
% confidence ratings in the two types of regularities
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Prepare the window
figure('Position', [1026 905 220 200]);

% Display origin
plot([0,1], [0,0], 'k-'); hold('on');

% For each type of regularity
for iReg = 1:2
    
    % Average over subjects
    md = mean(confkern(iReg,:,:), 3);
    sd = sem(confkern(iReg,:,:), 3);
    
    % Display the kernel distribution
    plotMSEM(u, md, sd, 1/5, tricol(iReg,:), tricol(iReg,:), 2);
    
    % Display the average confidence rating
    ma = mean(avgConf(:,iReg));
    sa = sem(avgConf(:,iReg));
    [~,mi] = min(abs(u-ma));
    plot(repmat(u(mi), [1,2]), [0,md(mi)], '-', 'Color', tricol(iReg,:));
    [~,mi] = min(abs(u-ma+sa.*[-1;1]), [], 2);
    plot(repmat(u(mi), [2,1]), [0,0;md(mi)], '-', 'Color', tricol(iReg,:), 'LineWidth', 1/2);
end

% Display the difference between the two distributions
d = squeeze(sum(confkern, 1));
md = mean(d,2);
sd = sem(d,2);
plotMSEM(u, md, sd, 1/5, g, g, 2);

% Customize the axes
ylim([-1,1].*max(abs(ylim)));
set(gca, 'YAxisLocation', 'Origin', 'Box', 'Off');

% Add some text labels
xlabel('Confidence (a.u.)');
ylabel('Density');

% Display reported confidence in the estimation for the two types of regularities
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Prepare the window
figure('Position', [1247 905 120 200]);

% Display difference in confidence ratings between the two types of regularity
Emergence_PlotSubGp(avgConf, tricol(1:2,:));

% Perform a paired t-test between confidence ratings
[h,pval,tci,stats] = ttest(diff(avgConf,1,2));
Emergence_PrintTstats(pval,tci,stats);

% Customize the axes
axis([0,3,0,1]);
set(gca, 'XColor', 'None', 'Box', 'Off');

% Display whether the difference is significant or not
Emergence_DispStatTest(avgConf);

% Add some labels
ylabel('Average confidence in the estimate');
