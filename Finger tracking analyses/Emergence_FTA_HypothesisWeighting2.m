% With this script, we seek to push further the analysis of trajectories at
% the time of the detection of deterministic regularities. In particular,
% we show that the extent to which subjects transiently increase their
% beliefs in the probabilistic hypothesis correlates with the extent to
% which the ideal observer does. As previously shown, we confirmed that
% this effect is subtended by how biased in terms of transition
% probabilities the rules (i.e. their entropies).
% 
% Copyright (c) 2018 Maxime Maheu

%% RUN THE PREVIOUS SCRIPT SEPARATELY FOR SUBJECTS AND THE IO
%  ==========================================================

% Run the analysis focusing on pre-detection periods using subjects' data
D = G;
Emergence_FTA_HypothesisWeighting1;
subtranspbel = transpbel;

% Run the analysis focusing on pre-detection periods using IO's data
D = IO;
Emergence_FTA_HypothesisWeighting1;
iotranspbel = transpbel;

% Combine the beliefs in the probabilistic hypothesis, before the detection
% of deterministic regularities, from the subjects and the IO
transpbel = cat(3, iotranspbel, subtranspbel);
transpbel = permute(transpbel, [1,3,2]);

%% CORRELATE SUBJECTS' AND IO'S BELIEFS IN THE PROBABILISTIC HYPOTHESIS
%  BEFORE THE DETECTION OF DETERMINISTIC REGULARITIES
%  ====================================================================

% Measure the correlation coefficient between IO's and subjects' beliefs
% separately for each subject
sub = mat2cell(subtranspbel, nR, ones(nSub, 1));
io  = mat2cell(iotranspbel,  nR, ones(nSub, 1));
coef = cellfun(@(toexplain,explainingvar) ...
    Emergence_Regress(toexplain, explainingvar, 'CC', 'r'), sub, io)';

% Prepare the window
figure('Position', [1 632 170 200]);

% Display average regression line across rules
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Average beliefs over subjects
subdim = ndims(transpbel);
m = mean(transpbel, subdim, 'OmitNaN');
s = sem(transpbel, subdim);

% Regress, across the different rules, averaged beliefs in the
% probabilistic hypotheses from the subjects and the ideal observer
B = Emergence_Regress(m(:,2), m(:,1), 'TLS', {'beta0', 'beta1'});
confint = Emergence_Regress(m(:,2), m(:,1), 'TLS', 'confint');
xval = Emergence_Regress(m(:,2), m(:,1), 'TLS', 'confintx');

% Display identity line
plot([0,1]./2, [0,1]./2, '-', 'Color', g, 'LineWidth', 1); hold('on');
text(0.15, 0.15, 'Identity', 'Color', g, 'VerticalAlignment', 'Top', 'Rotation', 45);

% Display the regression line
fill([xval, fliplr(xval)], [confint(1,:), fliplr(confint(2,:))], 'k', ...
       'EdgeColor', 'none', 'FaceColor', g, 'FaceAlpha', 1/2);
plot(xval, xval.*B(2) + B(1), 'k-', 'LineWidth', 3);

% Display each individual rule, with its corresponding SEM, and colored
% according to its entropy level (computed on transition probabilities)
for iR = 1:nR
    plot(m(iR,1)+[-s(iR,1),s(iR,1)], repmat(m(iR,2), 1, 2), '-', 'Color', 'k');
    plot(repmat(m(iR,1), 1, 2), m(iR,2)+[-s(iR,2),s(iR,2)], '-', 'Color', 'k');
    plot(m(iR,1), m(iR,2), 'o', 'MarkerEdgeColor', 'k', ...
        'MarkerFaceColor', EntCol(iR,:), 'MarkerSize', 7, 'LineWidth', 1);
end

% Customize the axes
set(gca, 'Box', 'Off');
axis('equal');
axis([0,0.3,0,0.5]);

% Add some text labels
xlabel('Ideal observer');
ylabel('Subjects');

% Save the figure
save2pdf(fullfile(ftapath, 'figs', 'F_HW_GpCorr.pdf'));

% Display individual regression coefficients
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Prepare a new window
figure('Position', [172 632 120 200]);
plot([0,3], zeros(1,2), '-', 'Color', g, 'LineWidth', 1); hold('on');

% Display the dispersion of correlation coefficients across subjects
Emergence_PlotSubGp(coef, zeros(1,3));

% Customize the axes
axis([0,2,-1,1]);
set(gca, 'Box', 'Off', 'XTick', [], 'XColor', 'None');

% Display whether the difference is significant or not
Emergence_DispStatTest(coef);

% Add some text labels
ylabel('Correlation coefficients');

% Compare the distribution of correlation coefficients against zero
[~,pval,tci,stats] = ttest(coef);
disptstats(pval,tci,stats);

% Save the figure
save2pdf(fullfile(ftapath, 'figs', 'F_HW_SubCorr.pdf'));
