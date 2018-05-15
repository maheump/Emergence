% This script implements analyses of finger tracking that aim at showing
% quantitative deviations that subjects exhibit compared to the ideal
% inference scenario. In particular, we show that subjects are delayed
% regarding the update of their beliefs compared to the ideal observer (at
% least partly because of the rapid auditory presentation).
% 
% Copyright (c) 2018 Maxime Maheu

%% COMPUTE SUBJECT-SPECIFIC AVERAGE DETECTION LAG
%  ==============================================

% Correlate subjects' with ideal observer's beliefs
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Create the parameter grid, i.e. the value of shift between the subjects
% and the ideal observer to be considered in the following parameter
% estimation loop
param = -140:140;
nParam = numel(param);

% Length of the sequence
N = numel(G{1}.Seq);

% Prepare output variable
coef = NaN(nParam, nCond, nSub); % R2 of interest
bchm = NaN(nParam, nCond, nSub); % benchmark R2

% For each subject
for iSub = 1:nSub
    disp(iSub);

    % For each considered value of the shift (in # of observations)
    % between subject and ideal observer beliefs
    for iParam = 1:nParam
        i = param(iParam);
        
        % For each sequence
        for iCond = 1:nCond

            % Get relevant trajectories
            subbel = G{iCond,iSub}.BarycCoord; % from the subject
            iobel  = IO{iCond,iSub}.BarycCoord; % from the ideal observer
            
            % Make sure probabilities do evolve between 0 and 1
            subbel(subbel < 0) = 0; subbel(subbel > 1) = 1;
            iobel(iobel   < 0) = 0; iobel(iobel   > 1) = 1;
            
            % Shift the subject's trajectory with respect to the ideal
            % observer's trajectory
            if i >= 0 % if the subject is slower than the ideal observer
                subtraj = subbel(i+1:N,:); % to be explained: subject's trajectory
                iotraj  = iobel(1:N-i,:);  % explaining variable: IO's belief
            elseif i < 0 % if the subject is faster than the ideal observer
                subtraj = subbel(1:N+i,:); % to be explained: subject's trajectory
                iotraj  = iobel(-i+1:N,:); % explaining variable: IO's belief
            end
            
            % Measure the correlation between subject's and ideal
            % observer's trajectories
            coef(iParam,iCond,iSub) = Emergence_Regress(...
                subtraj(:), iotraj(:), 'CC', 'r');
            
            % Run a benchmark regression in which both the explained
            % variable and the one to be explained come from the ideal
            % observer. However, one can be shifted compared to the other.
            % This provides us a symmetric R2 behaviour that help quantify
            % the effect of the shift on two perfectly correlated
            % variables.
            if i >= 0 % if the subject is slower than the ideal observer
                iotraj1 = iobel(i+1:N,:);
                iotraj2 = iobel(1:N-i,:);
            elseif i < 0 % if the subject is faster than the ideal observer
                iotraj1 = iobel(1:N+i,:);
                iotraj2 = iobel(-i+1:N,:);
            end
            bchm(iParam,iCond,iSub) = Emergence_Regress(iotraj1(:), iotraj2(:), 'CC', 'r');
        end
    end
end

% Average over sequences (i.e. this is the same as fitting all the
% sequences together at once)
avgR2 = squeeze(mean(coef, 2));

% Display averaged R2 and corresponding best shift parameter for each subject
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Create a new window
figure('Position', [1 906 200 200]);

% Display a no-shift line
plot(zeros(1,2), [0,1], 'k--'); hold('on');

% Display the benchmark scenario (fitting the IO against itself)
avgBM = squeeze(mean(bchm, 3));
plot(param, mean(avgBM, 2), 'k-', 'LineWidth', 1/2);

% Display group-level R2 curves according to shift parameter
gpavgR2 = mean(avgR2, 2);
gpsemR2 = sem(avgR2, 2);
f = plotMSEM(param, gpavgR2, gpsemR2, 0.15, 'k', 'k', 3, 1, '-', 'none');

% Display group-level best shift parameter
[gpmaxavgR2, gpmaxidxavgR2] = max(gpavgR2, [], 1);
plot(param(gpmaxidxavgR2), gpmaxavgR2, 'o', 'LineWidth', 1, 'MarkerSize', 8, ...
    'MarkerEdgeColor', get(f, 'Color'), 'MarkerFaceColor', g);
text(param(gpmaxidxavgR2), gpmaxavgR2, sprintf('   delay = %1.0f', ...
    param(gpmaxidxavgR2)), 'VerticalAlignment', 'Bottom', 'HorizontalAlignment', 'Left');

% Customize the axes
axis([param(1), param(end), 0, 1]);
set(gca, 'Box', 'Off');

% Add some text labels
xlabel('Delay btw subjects and the IO (# obs.)');
ylabel('Correlation coefficient');

% Save the figure
save2pdf('figs/F_QD_CorrSubIO.pdf');

% Display dispersion of lag over subjects
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Create a new window
figure('Position', [202 906 120 200]);

% Display the ideal scenario: no delay
plot([0,2], zeros(1,2), 'k--'); hold('on');

% Display subject-level best shift parameters
[submaxavgR2, submaxidxavgR2] = max(avgR2, [], 1);
Emergence_PlotSubGp(param(submaxidxavgR2), 'k');

% Customize the axes
axis([0, 2, -10, 50]);
set(gca, 'Box', 'Off', 'XColor', 'None');

% Display whether the difference is significant or not
Emergence_DispStatTest(param(submaxidxavgR2));

% Add some text labels
ylabel('Estimated delay (# obs.)');

% Save the figure
save2pdf('figs/F_QD_Delay.pdf');

% Correlate average lag and quality of fit with the IO over subjects
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Create a new window
figure('Position',  [323 906 220 200]);

% Display confidence intervals
x = param(submaxidxavgR2);
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
axis([-2, 20, 0.5, 1]);
set(gca, 'Box', 'Off');

% Add some text labels
xlabel('Estimated delay (# obs.)'); ylabel('Correlation with the IO');

% Save the figure
save2pdf('figs/F_QD_Corr.pdf');
