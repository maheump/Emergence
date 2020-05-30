% This script implements analyses of finger tracking that aim at showing
% quantitative deviations that subjects exhibit compared to the ideal
% inference scenario. In particular, we show that subjects are delayed
% regarding the update of their beliefs compared to the ideal observer (at
% least partly because of the rapid auditory presentation).
% 
% Copyright (c) 2020 Maxime Maheu

% Correlate subjects' with ideal observer's beliefs
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Create the parameter grid, i.e. the value of shift between the subjects
% and the ideal observer to be considered in the following parameter
% estimation loop
shiftlist = -20:100;
nParam = numel(shiftlist);

% Prepare output variable
[~,coef] = cellfun(@(p,q) Emergence_Similarity(p, q, 'PC', shiftlist), ...
    G, IO, 'uni', 0);
coef = cell2mat(reshape(coef, [1,nSeq,nSub]));

% Average over sequences (i.e. this is the same as fitting all the
% sequences together at once)
avgcoef = squeeze(mean(coef, 2, 'OmitNaN'));

% Display averaged R2 and corresponding best shift parameter for each subject
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Create a new window
figure('Position', [1 906 200 200]);

% Display group-level R2 curves according to shift parameter
gpavgcoef = mean(avgcoef, 2);
gpsemcoef = sem(avgcoef, 2);
f = plotMSEM(shiftlist, gpavgcoef, gpsemcoef, 0.15, 'k', 'k', 3, 1, '-', 'none'); hold('on');

% Display group-level best shift parameter
[gpmaxavgR2, gpmaxidxavgR2] = max(gpavgcoef, [], 1);
plot(shiftlist(gpmaxidxavgR2), gpmaxavgR2, 'o', 'LineWidth', 1, 'MarkerSize', 8, ...
    'MarkerEdgeColor', get(f, 'Color'), 'MarkerFaceColor', g);
text(shiftlist(gpmaxidxavgR2), gpmaxavgR2, sprintf('   delay = %1.0f', ...
    shiftlist(gpmaxidxavgR2)), 'VerticalAlignment', 'Bottom', 'HorizontalAlignment', 'Left');

% Customize the axes
xlim([shiftlist(1), 50]);
set(gca, 'Box', 'Off');

% Add some text labels
xlabel('Delay btw subjects and the IO (# obs.)');
ylabel('Correlation coefficient');

% Save the figure
save2pdf(fullfile(ftapath, 'figs', 'F_QD_CorrSubIO.pdf'));

% Display dispersion of lag over subjects
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Create a new window
figure('Position', [202 906 120 200]);

% Display the ideal scenario: no delay
plot([0,2], zeros(1,2), '-', 'Color', g); hold('on');
text(0, 0, ' Identity', 'Color', g, 'VerticalAlignment', 'Top', ...
    'HorizontalAlignment', 'Left');

% Display subject-level best shift parameters
[submaxavgR2, submaxidxavgR2] = max(avgcoef, [], 1);
Emergence_PlotSubGp(shiftlist(submaxidxavgR2), 'k');

% Customize the axes
axis([0, 2, -5, 20]);
set(gca, 'Box', 'Off', 'XColor', 'None');

% Display whether the difference is significant or not
inddelay = shiftlist(submaxidxavgR2);
Emergence_DispStatTest(inddelay);

% Perform a t-test on integration delay against zero
[~,pval,tci,stats] = ttest(inddelay');
Emergence_PrintTstats(pval,tci,stats);

% Add some text labels
ylabel('Estimated delay (# obs.)');

% Save the figure
save2pdf(fullfile(ftapath, 'figs', 'F_QD_Delay.pdf'));

% Correlate average lag and quality of fit with the IO over subjects
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Create a new window
figure('Position',  [323 906 220 200]);

% Display confidence intervals
x = shiftlist(submaxidxavgR2);
y = submaxavgR2;
confint  = Emergence_Regress(y, x, 'OLS', 'confint');
confintx = Emergence_Regress(y, x, 'OLS', 'confintx');
fill([confintx, fliplr(confintx)], [confint(1,:), fliplr(confint(2,:))], ...
    'k', 'EdgeColor', 'none', 'FaceColor', 'k', 'FaceAlpha', 0.15); hold('on');

% Display regression line
b = Emergence_Regress(y, x, 'OLS', {'beta0', 'beta1'});
plot(confintx([1,end]), confintx([1,end]).*b(2)+b(1), 'k-', 'LineWidth', 3);

% For each subject, display the dispersion of correlation coefficients over
% sequences
s = NaN(1,nSub);
for iSub = 1:nSub, s(iSub) = sem(coef(submaxidxavgR2(iSub),:,iSub), 2); end
plot([x; x], y + s.*[-1;1], 'k');

% Display individual subjects
plot(x, y, 'ko', 'MarkerFaceColor', g, 'MarkerSize', 8);

% Customize the axes
%axis([-2, 20, 0.5, 1]);
set(gca, 'Box', 'Off', 'Layer', 'Bottom');

% Add some text labels
xlabel('Estimated delay (# obs.)'); ylabel('Correlation with the IO');

% Test the correlation
[rho,pval] = corr(x', y');
disp(pval);

% Save the figure
save2pdf(fullfile(ftapath, 'figs', 'F_QD_Corr.pdf'));
