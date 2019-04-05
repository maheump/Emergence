% This script analyses inferred positions of the change point. We show
% that for both probabilistic and deterministic regularities, estimates of
% change point's position are correlated with the real positions of the
% change point. Moreover, we show that the confidence in that estimate is
% higher in the case of deterministic than probabilistic regularities.
%
% Copyright (c) 2018 Maxime Maheu

%% CORRELATION BETWEEN ESTIMATED AND REAL CHANGE POINT'S POSITION
%  ==============================================================

% Define the number of bins to use for the correlation at the group-level
nBin = 7;

% Correlate real and inferred positions of the change point
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Create the bins in which the positions will be averaged
pgrid = linspace(0, 100, nBin+1);
allcp = cellfun(@(x) x.Jump-1/2, D(cat(2, cidx{1:2}),:));
pct = prctile(allcp(:), pgrid);

% Prepare output variable
corrcoef = NaN(nSub,2);
cpbinsm  = NaN(nBin,2);
cpbinss  = NaN(nBin,2);

% For each type of regularity
for iHyp = 1:2
    
    % Get the number of sequences with this type of regularity
    nR = numel(cidx{iHyp});
    
    % Get objective and subjective change/detection point's positions
    cp = cellfun(@(x) x.Jump-1/2, D(cidx{iHyp},:));
    sp = cellfun(@(x) x.Questions(3), D(cidx{iHyp},:));
    
    % Remove sequences in which regularities were not detected
    detecmask = (filter{iHyp} == 1 | filter{iHyp} == 3);
    cp(~detecmask) = NaN;
    sp(~detecmask) = NaN;
    
    % Subject-specific regression
    corrcoef(:,iHyp) = cellfun(@(toexplain,explainingvar) ...
        Emergence_Regress(toexplain, explainingvar, 'CC', 'r'), ...
        mat2cell(sp, nR, ones(nSub,1)), ... % inferred change point's position
        mat2cell(cp, nR, ones(nSub,1)));    % true change point's position
    
    % Find the data points that correspond to each bin
    [~,~,idx] = histcounts(cp, pct);
    
    % For each bin
    for iBin = 1:nBin
        
        % Average over sequences
        cpbinsm(iBin,1,iHyp) = mean(cp(idx == iBin), 'OmitNaN');
        cpbinsm(iBin,2,iHyp) = mean(sp(idx == iBin), 'OmitNaN');
        
        % Compute the SEM over sequences
        cpbinss(iBin,1,iHyp) = sem(cp(idx == iBin));
        cpbinss(iBin,2,iHyp) = sem(sp(idx == iBin));
    end
end

% Display group-level correlation between real and estimated change point's position
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Prepare a new figure
figure('Position', [463 906 220 200]);

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
for iHyp = 1:2
    
    % Shortcuts to avoid errors
    ocp_m = cpbinsm(:,1,iHyp)';
    scp_m = cpbinsm(:,2,iHyp)';
    ocp_s = cpbinss(:,1,iHyp)';
    scp_s = cpbinss(:,2,iHyp)';
    
    % Group-level regression between (binned) objective and subjective
    % positions of the change point
    beta    = Emergence_Regress(scp_m, ocp_m, 'OLS', {'beta0', 'beta1'});
    confint = Emergence_Regress(scp_m, ocp_m, 'OLS', 'confint');
    xval    = Emergence_Regress(scp_m, ocp_m, 'OLS', 'confintx');
    
    % Regression on binned data between objective and inferred position
    fill([xval, fliplr(xval)], [confint(1,:), fliplr(confint(2,:))], 'k', ...
        'EdgeColor', 'None', 'FaceColor', tricol(iHyp,:), 'FaceAlpha', 1/3);
    plot(xval, beta(1)+beta(2).*xval, '-', 'Color', tricol(iHyp,:), 'LineWidth', 3);
    
    % Display horizontal and vertical error bars
    plot(repmat(ocp_m, [2,1]), scp_m+scp_s.*[-1;1], 'k-');
    plot(ocp_m+ocp_s.*[-1;1], repmat(scp_m, [2,1]), 'k-');
    
    % Display the average value in that bin
    plot(ocp_m, scp_m, 'ko', 'MarkerFaceColor', tricol(iHyp,:), 'MarkerSize', 7);
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

% Save the figure
if isfield(D{1}, 'Seq'), save2pdf(fullfile(ftapath, 'figs', 'F_CP_Corr1S.pdf'));
else, save2pdf(fullfile(ftapath, 'figs', 'F_CP_Corr1IO.pdf'));
end

% Display individual correlation coefficient from the correlation between
% real and estimated change point's position
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Prepare a new figure
figure('Position', [684 905 120 200]);

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

% Save the figure
if isfield(D{1}, 'Seq'), save2pdf(fullfile(ftapath, 'figs', 'F_CP_Corr2S.pdf'));
else, save2pdf(fullfile(ftapath, 'figs', 'F_CP_Corr2IO.pdf'));
end

%% CONFIDENCE IN THE ESTIMATE
%  ==========================

% Get confidence in the change point's estimation question
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Get confidence levels
conf = cellfun(@(x) x.Questions(4), D) ./ 100;

% Remove sequences in which regularities were not detected
detecmask = (filter{iHyp} == 1 | filter{iHyp} == 3);
conf(~detecmask) = NaN;

% Average over regularities for each condition (proba. & deter.)
avgConf = NaN(nSub,2);
for iHyp = 1:2, avgConf(:,iHyp) = mean(conf(cidx{iHyp},:), 1, 'OmitNaN')'; end

% Display the entire distribution (over subjects and sequences) of
% confidence ratings in the two types of regularities
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Prepare the window
figure('Position', [804 906 120 200]);

% Prepare the output variable
prc = 1e4;
xgrid = linspace(0, 1, prc);
f = NaN(2,prc);

% For each type of regularity
for iHyp = 1:2
    
    % Measure a kernel density from subjects' confidence ratings
    dist = conf(cidx{iHyp},:);
    [f(iHyp,:), u] = ksdensity( dist(:), ...        % which distribution to plot
        xgrid, ...                                  % grid of confidence levels
        'BandWidth',            0.05, ...           % bandwidth of the kernel smoothing window 
        'Support',              [0-eps,1+eps], ...  % restrict the kernel to a certain range of values
        'BoundaryCorrection',   'Reflection');      % type of correction for the boundaries
    
    % Normalize the distribution such that it can be directly compared
    % between both conditions (because there are not necessarily the same
    % number of accurately detected regularities in both conditions).
    f(iHyp,:) = f(iHyp,:) ./ sum(f(iHyp,:));
    
    % Measure average confidence rating
    avg = mean(dist(:), 'OmitNaN');
    [~,i] = min(abs(avg - u));
    
    % Mirror the second distribution
    if iHyp == 1, f(iHyp,:) = -f(iHyp,:); end
    
    % Make sure the distribution does not exceeds the limits of the scale
    f(iHyp, u >= 1) = 0;
    f(iHyp, u <= 0) = 0;
    
    % Display the kernel distribution
    plot(f(iHyp,:), u, 'Color', tricol(iHyp,:)); hold('on');
    
    % Display the average confidence rating
    plot([f(iHyp,i),0], repmat(u(i),1,2), '-', 'Color', tricol(iHyp,:), 'LineWidth', 2)
end

% Compute the difference between the two kernel distributions
kd = f(1,:) + f(2,:);

% Detect parts of the scale that are more used in one condition or another
kdidx = {kd < 0, kd > 0};

% For each type of regularity
for iHyp = 1:2
    
    % Get windows containing confidence ratings that are more used in that
    % condition
    on = find(diff(kdidx{iHyp}) == 1) + 1;
    off = find(diff(kdidx{iHyp}) == -1);
    
    % Fill the area under the difference curve with the appropriate color
    for i = 1:numel(on)
        idx = on(i):off(i);
        fill([0,kd(idx),0], [u(idx(1)-1),u(idx),u(idx(end))], ...
            'k', 'FaceColor', tricol(iHyp,:), 'EdgeColor', 'None');
    end
end

% Plot the difference curve 
plot(kd, u, 'k-', 'LineWidth', 1);

% Customize the axes
axis([-max(abs(xlim)), max(abs(xlim)), 0, 1]);
set(gca, 'YAxisLocation', 'Origin', 'Box', 'Off');

% Save the figure
if isfield(D{1}, 'Seq'), save2pdf(fullfile(ftapath, 'figs', 'F_CP_ConfDistS.pdf'));
else, save2pdf(fullfile(ftapath, 'figs', 'F_CP_ConfDistIO.pdf'));
end

% Display reported confidence in the estimation for the two types of regularities
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Prepare the window
figure('Position', [925 905 120 200]);

% Display difference in confidence ratings between the two types of regularity
Emergence_PlotSubGp(avgConf, tricol(1:2,:));

% Perform a paired t-test between confidence ratings
[h,pval,ci,stats] = ttest(diff(avgConf,1,2));
Emergence_PrintTstats(pval,tci,stats);

% Customize the axes
axis([0,3,0,1]);
set(gca, 'XColor', 'None', 'Box', 'Off');

% Display whether the difference is significant or not
Emergence_DispStatTest(avgConf);

% Add some labels
ylabel('Confidence in the estimate');

% Save the figure
if isfield(D{1}, 'Seq'), save2pdf(fullfile(ftapath, 'figs', 'F_CP_ConfAvgS.pdf'));
else, save2pdf(fullfile(ftapath, 'figs', 'F_CP_ConfAvgIO.pdf'));
end
