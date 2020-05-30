% This script implements analyses looking at the inferred generative
% process of the sequence once the sequence has been entirely presented
% (post-sequence questions). We show that probabilistic regularity that are
% missed by subjects are associated with lower beliefs from the ideal
% observer. Finally, we investigate the determinants of missed and
% accurately detected regularities. In particuler, we show that, for
% probabilistic regularities, it depends on the entropy level and the type
% of bias (alternation are more often biased) while, for deterministic
% regularities, it depends on the length of the repeating rule and its
% inner organisation.
% 
% Copyright (c) 2020 Maxime Maheu

% The following script must be run beforehand
Emergence_OFF_Detection1;

%% GET THE IO'S POTERIOR BELIEFS IN DETECTED AN UNDETECTED PROBABILISTIC REGULARITIES
%  ==================================================================================

% Focus on probabilistic regularities because these are those that are most
% often missed by the ideal observer
iHyp = 1;

% Get the average beliefs of the ideal observer for probabilistic
% regularities that have been detected versus missed by the subjects
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Get the hidden generative process
truegenproc = cellfun(@(x) find(strncmpi(x.Cond(1), proclab, 1)), G);

% Deduce wether the estimation was correct
corest = (truegenproc == estgenproc);
corest = mat2cell(corest(cidx{iHyp},:), numel(cidx{iHyp}), ones(nSub,1));

% Get the ultimate value of belief in the probabilistic hypothesis from the
% ideal observer
ultibel = cellfun(@(x,y) x.BarycCoord(end,y), IO, num2cell(truegenproc));
ultibel = mat2cell(ultibel(cidx{iHyp},:), numel(cidx{iHyp}), ones(nSub,1));

% Check whether missed regularities can be explained by lower posterior
% beliefs in the ideal observer model
avgbel = cat(2, cellfun(@(x,y) mean(x(y == 0)), ultibel, corest)', ...
                cellfun(@(x,y) mean(x(y == 1)), ultibel, corest)');
[~,pval,tci,stats] = ttest(diff(avgbel, 1, 2));
Emergence_PrintTstats(pval, tci, stats);

% Check whether missed regularities can be explained by a latter change
% point position
subcp = mat2cell(cp{iHyp} + 1/2, numel(cidx{iHyp}), ones(1,nSub));
avgcp = cat(2, cellfun(@(x,y) mean(x(y == 0)), subcp, corest)', ...
               cellfun(@(x,y) mean(x(y == 1)), subcp, corest)');
[~,pval,tci,stats] = ttest(diff(avgcp, 1, 2));
Emergence_PrintTstats(pval, tci, stats);
[~,bf01] = BF_ttest(diff(avgcp, 1, 2));
fprintf('BF in favour of the null: %1.2f\n', bf01);

% Display average IO's beliefs following sequences entailing a
% probabilistic regularity that was either detected or missed by subjects
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Prepare a new window
figure('Position', [835 806 240 200]);

% Display middle of the sequence
subplot(1,2,1);
plot([0,3], 100.*ones(1,2), '-', 'Color', g, 'LineWidth', 1); hold('on');

% Display the change point position in missed versus detected regularities
Emergence_PlotSubGp(avgcp, tricol(1,:));

% Customize the axes
axis([0,3,60,140]);
set(gca, 'XTick', 1:2, 'XTickLabel', {'Undetected', 'Detected'}, ...
    'XTickLabelRotation', 30);
set(gca, 'Box', 'Off');

% Display whether the difference is significant or not
Emergence_DispStatTest(avgcp);

% Add some text labels
xlabel({'Probabilistic', 'regularities'});
ylabel('Change point position');

% Display chance level
subplot(1,2,2);
plot([0,3], ones(1,2)./3, '-', 'Color', g, 'LineWidth', 1); hold('on');

% Display difference in IO's beliefs between missed versus detected
% regularities
Emergence_PlotSubGp(avgbel, tricol(1,:));

% Customize the axes
axis([0,3,0,1]);
set(gca, 'XTick', 1:2, 'XTickLabel', {'Undetected', 'Detected'}, ...
    'XTickLabelRotation', 30);
set(gca, 'Box', 'Off');

% Display whether the difference is significant or not
Emergence_DispStatTest(avgbel);

% Add some text labels
xlabel({'Probabilistic', 'regularities'});
ylabel('Belief p(M_P|y) from the IO');

% Save the figure
if isfield(D{1}, 'Seq'), save2pdf(fullfile(ftapath, 'figs', 'F_D_BelUndetS.pdf'));
else, save2pdf(fullfile(ftapath, 'figs', 'F_D_BelUndetIO.pdf'));
end

%% DETERMINANTS OF DETECTION FOR PROBABILISTIC REGULARITIES
%  ========================================================

% Depending on the entropy
% ~~~~~~~~~~~~~~~~~~~~~~~~

% Get theoretical transition probabilities
pXgY = cell2mat(prob');
pAgB = pXgY(:,1);
pBgA = pXgY(:,1);

% Compute entropy levels
TPent = arrayfun(@(x,y) Emergence_MarkovEntropy(x, y), pAgB, pBgA);

% Define bins
entlim = 1.85;
entgpidx = {find(TPent <= entlim), find(TPent > entlim)};
entgpttl = {sprintf('Low (< %1.2f)', entlim), sprintf('High (> %1.2f)', entlim)};

% Get the average detection rate in each entropy bins
entdetecrate = cell2mat(cellfun(@(x) mean(estgenproc(cidx{1}(x),:) == 1, 1)', ...
    entgpidx, 'uni', 0));

% Prepare a window
figure('Position', [1076 806 240 200]);

% Display chance level
subplot(1,2,1);
plot([0,3], ones(1,2)./3, '-', 'Color', g, 'LineWidth', 1); hold('on');

% Display difference in detection rate according to entropy levels
Emergence_PlotSubGp(entdetecrate, tricol(1,:));

% Display whether the difference is significant or not
Emergence_DispStatTest(entdetecrate);

% Customize the axes
axis([0,3,0,1]);
set(gca, 'XTick', 1:2, 'XTickLabel', entgpttl, 'XTickLabelRotation', 30);
set(gca, 'Box', 'Off');

% Add some labels
xlabel('Entropy levels');
ylabel('Average detection rate');

% Depending on the type of bias
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Create groups
pAlt = cellfun(@(x) (2.*x(1).*x(2)) ./ (x(1) + x(2)), prob);
biasgpidx = {find(pAlt < 1/2), find(pAlt > 1/2), 7:9};
biasgpttl = {'Repetition', 'Alternation', 'Frequency',};
nGp = numel(biasgpidx);

% Get the average detection rate in each group of regularities
biasdetecrate = cell2mat(cellfun(@(x) mean(estgenproc(cidx{1}(x),:) == 1, 1)', ...
    biasgpidx, 'uni', 0));

% Display chance level
subplot(1,2,2);
plot([0,4], ones(1,2)./3, '-', 'Color', g, 'LineWidth', 1); hold('on');

% Display difference in detection rate according to biased dimensions
Emergence_PlotSubGp(biasdetecrate, tricol(1,:));

% Display whether the difference is significant or not
Emergence_DispStatTest(biasdetecrate);
[~,pval,tci,stats] = ttest(biasdetecrate(:,2)-biasdetecrate(:,1));
Emergence_PrintTstats(pval,tci,stats);
[~,pval,tci,stats] = ttest(biasdetecrate(:,2)-biasdetecrate(:,3));
Emergence_PrintTstats(pval,tci,stats);

% Customize the axes
axis([0,4,0,1]);
set(gca, 'XTick', 1:nGp, 'XTickLabel', biasgpttl, 'XTickLabelRotation', 30);
set(gca, 'Box', 'Off');

% Add some labels
xlabel('Biased dimension');
ylabel('Average detection rate');

% Save the figure
if isfield(D{1}, 'Seq'), save2pdf(fullfile(ftapath, 'figs', 'F_D_ProbRegS.pdf'));
else, save2pdf(fullfile(ftapath, 'figs', 'F_D_ProbRegIO.pdf'));
end

%% DETERMINANTS OF DETECTION FOR DETERMINISTIC REGULARITIES
%  ========================================================

% Depending on the rules' length
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Create groups
len = cellfun(@numel, dr(1:end));
lengpidx = cellfun(@(x) find(x == len), num2cell(unique(len)), 'uni', 0)';
nGp = numel(lengpidx);

% Get the average detection rate in each group of regularities
lendetecrate = cell2mat(cellfun(@(x) mean(estgenproc(cidx{2}(x),:) == 2, 1)', ...
    lengpidx, 'uni', 0));

% Prepare a new window
figure('Position', [1317 806 240 200]);

% Display chance level
subplot(1,2,1);
plot([0,5], ones(1,2)./3, '-', 'Color', g, 'LineWidth', 1); hold('on');

% Display difference in detection rate according to rules' length
Emergence_PlotSubGp(lendetecrate, tricol(2,:));

% Customize the axes
axis([0,nGp+1,0,1]);
set(gca, 'XTick', 1:nGp, 'XTickLabel', num2cell(unique(len)));
set(gca, 'Box', 'Off');

% Display whether the difference is significant or not
Emergence_DispStatTest(lendetecrate);

% Add some labels
xlabel('Rules'' length');
ylabel('Average detection rate');

% Depending on the rules' difficulty
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Create groups
cplxidx = {[1,2,5,8], [3,6,9], [4,7,10]};
cplxttl = {'Easy', 'Moderate', 'Difficult'};
nGp = numel(cplxidx);

% Get the average detection rate in each group of regularities
cplxetecrate = cell2mat(cellfun(@(x) mean(estgenproc(cidx{2}(x),:) == 2, 1)', ...
    cplxidx, 'uni', 0));

% Display chance level
subplot(1,2,2);
plot([0,4], ones(1,2)./3, '-', 'Color', g, 'LineWidth', 1); hold('on');

% Display difference in detection rate according to rules' difficulty
Emergence_PlotSubGp(cplxetecrate, tricol(2,:));

% Customize the axes
axis([0,4,0,1]);
set(gca, 'XTick', 1:numel(cplxidx), 'XTickLabel', cplxttl, 'XTickLabelRotation', 30);
set(gca, 'Box', 'Off');

% Display whether the difference is significant or not
Emergence_DispStatTest(lendetecrate);

% Add some labels
xlabel('Rule''s difficulty');
ylabel('Average detection rate');

% Save the figure
if isfield(D{1}, 'Seq'), save2pdf(fullfile(ftapath, 'figs', 'F_D_DetRegS.pdf'));
else, save2pdf(fullfile(ftapath, 'figs', 'F_D_DetRegIO.pdf'));
end
