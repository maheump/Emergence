% This script implements analyses of finger tracking that aim at showing
% quantitative deviations that subjects exhibit compared to the ideal
% inference scenario. In particular, we show that subjects are biased in
% the way they report their probabilistic beliefs in a similar manner of 
% what has been shown by economic theory.
% 
% Copyright (c) 2018 Maxime Maheu

%% EXPLAIN SUBJECTS' TRAJECTORIES USING IO'S BELIEFS WITH A MIXED EFFECT APPROACH
%  ===============================================================================

% Define the number of bins to use
nBin = 11;

% Define the type of binning method
binmeth = 'unif';

% Create bins of probability values
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Get probability bins
if strcmpi(binmeth, 'unif') % bins uniformed on the probability space
    pgrid = linspace(0, 1, nBin+1);
elseif strcmpi(binmeth, 'equil') % bins uniformed on the probability space
    iobel = cellfun(@(x) x.BarycCoord(:), IO, 'UniformOutput', 0);
    iobel = cell2mat(iobel(:));
    pgrid = prctile(iobel, linspace(0, 100, nBin));
else, error('Please check the binnig method that is provided');
end

% Prepare the output variable
binbel = NaN(nBin,2,nSub);

% For each subject
for iSub = 1:nSub
    
    % Get beliefs of the subject and those of the ideal observer in all the
    % sequences
    subbel = cellfun(@(x) x.BarycCoord(:), G(:,iSub), 'UniformOutput', 0);
    subbel = cell2mat(subbel(:));
    iobel = cellfun(@(x) x.BarycCoord(:), IO(:,iSub), 'UniformOutput', 0);
    iobel = cell2mat(iobel(:));
    
    % Average subject's and IO's beliefs in each bin
    for iBin = 1:nBin
        idx = iobel >= pgrid(iBin) & iobel < pgrid(iBin+1);
        binbel(iBin,1,iSub) = mean(iobel(idx));
        binbel(iBin,2,iSub) = mean(subbel(idx));
    end
end

% Fit subjects' trajectories against IO's
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% The observation functions to use (see at the end of this script)
g_fname = {@g_PW, ... % Model 1: probabilistic weighting function
           @g_LM};    % Model 2: linear regression

% Specify the dimensions of the problem
dim         = [];   % empty variable
dim.n_phi   = 2;    % number of observation parameters
dim.n       = 0;    % number of hidden states
dim.n_theta = 0;    % number of evolution parameters
dim.p       = nBin; % number of obervations per time sample
dim.n_t     = 1;    % number of time samples

% Define minimal options
options             = [];
options.DisplayWin  = 0;
options.verbose     = 0;
options             = repmat({options}, [nSub,1]);
optiongp            = [];
optiongp.DisplayWin = 0;
optiongp.verbose    = 0;

% Explaining variable: IO's binned beliefs
for iSub = 1:nSub, options{iSub}.inG.p = binbel(:,1,iSub); end

% Variable to explain
y = mat2cell(squeeze(binbel(:,2,:)), nBin, ones(nSub,1))';

% Try to load results from the previous analysis
try
filename = fullfile(homedir, 'Finger tracking analyses', 'ppdata', 'QD1_MFX.mat');
load(filename);
catch

% Prepare output variables
p_sub = cell(2,nSub); p_gp = cell(2,1); % posterior over parameters
o_sub = cell(2,nSub); o_gp = cell(2,1); % quality of fit

% Run the mixed-effect fitting scheme (it deduces priors through empirical
% Bayes) using a variational procedure (it uses Laplace approximation) in
% order to estimate the best parameters for each model and each subject
for m = 1:2
    [p_sub(m,:), o_sub(m,:), p_gp{m}, o_gp{m}] = ...
        VBA_MFX(y, [], [], g_fname{m}, dim, options, [], optiongp);
end

% Save the result of this analysis in a MATLAB file
save(filename, 'p_sub', 'o_sub', 'p_gp', 'o_gp');
    
end

% Perform model comparison
% ~~~~~~~~~~~~~~~~~~~~~~~~

% Perform Bayesian model selection
L = cellfun(@(x) x.F, o_sub);
options = [];
options.DisplayWin = 0;
[posterior,out] = VBA_groupBMC(L, options);

% Get individual model frequencies
pmf = posterior.r';

% Get labels of the models that have been fitted
funlab = cellfun(@func2str, g_fname, 'UniformOutput', 0);
funlab = cellfun(@(x) x(3:end), funlab, 'UniformOutput', 0);

% Prepare a new window
figure('Position', [1 906 120 200]);

% Display chance level
plot([0,3], ones(1,2)/2, '-', 'Color', g); hold('on');

% Display estimated model frequencies
avgpmf = mean(pmf);
bar(1:2, avgpmf, 'EdgeColor', 'k', 'FaceColor', repmat(1/2,1,3));

% Display estimated model frequencies
sempmf = sem(pmf);
plot(repmat(1:2, 2, 1), avgpmf + sempmf' .* [-1;1], 'k-');

% Customize the axes
axis([0,3,0,1]);
set(gca, 'XTick', 1:2, 'XTickLabel', funlab);
set(gca, 'Box', 'Off');

% Add some text labels
ylabel('Model frequencies');

% Save the figure
save2pdf('figs/F_QD_BMS.pdf');

% Display the fit at the group level
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Prepare window
figure('Position', [122 906 220 200]);

% Display identity line and help lines
plot([0,1], [0,1], '-', 'Color', g); hold('on');
text(0.15, 0.15, 'Identity', 'Color', g, 'VerticalAlignment', 'Top', 'Rotation', 45);

% Average best parameters over subjects
PWparams = cell2mat(cellfun(@(x) x.muPhi, p_sub(1,:), 'UniformOutput', 0));
mP = mean(PWparams,2);
sP = sem(PWparams,2);

% Compute upper and lower positions of the probability weighting function
P = [mP(1) + sP(1) * [-1;-1;1;1], mP(2) + sP(2) * [-1;1;-1;1]];
P = mat2cell(P, ones(1,4), 2);
p = linspace(0, 1, 1001);
y = cell2mat(cellfun(@(x) g_PW([], x, [], p), P, 'UniformOutput', 0));
fill([p, fliplr(p)], [max(y, [], 1), fliplr(min(y, [], 1))], ...
    'k', 'FaceAlpha', 0.15, 'EdgeColor', 'none');

% Display probability weighting function obtained with group-average
% parameters
fp = g_PW([], mP, [], p);
plot(p, fp, 'k-', 'LineWidth', 3);

% Display averaged bins and related error
avgbinbel = mean(binbel, 3);
sembinbel = sem(binbel, 3);
plot(repmat(avgbinbel(:,1), 1, 2)', avgbinbel(:,2)'+sembinbel(:,2)'.*[-1,1]', 'k-')
plot(avgbinbel(:,1), avgbinbel(:,2), 'ko', 'MarkerFaceColor', g, 'MarkerSize', 8);

% Customize the axes
axis(repmat([0,1], 1, 2)); axis('square');
set(gca, 'XTick', 0:0.2:1, 'YTick', 0:0.2:1);
set(gca, 'Box', 'Off');

% Add some text labels
xlabel('Beliefs from the ideal observer'); ylabel('Beliefs from subjects');

% Save the figure
save2pdf('figs/F_QD_AvgPWFit.pdf');

% Display best parameters
% ~~~~~~~~~~~~~~~~~~~~~~~

% Prepare a new window
figure('Position', [343 906 200 200]);

% Define parameters' name and null values
pwparamname = {'\gamma', 'p_{0}'};
nullparamval = [1, 1/2];

% For each parameter (gamma and p0) of the probability weighting function
for iParam = 1:2
    subplot(1,2,iParam); hold('on');
    
    % Display null value of the parameter aginst which it will be compared
    plot([0,2], repmat(nullparamval(iParam),1,2), 'k--');
    
    % Display distribution of parameter value
    Emergence_PlotSubGp(PWparams(iParam,:), 'k');
    
    % Customize the axes
    xlim([0,2]);
    if     iParam == 1, ylim([0,2]);
    elseif iParam == 2, ylim([0,1]);
    end
    set(gca, 'Box', 'Off', 'XTick', [], 'XColor', 'none');
    
    % Display whether the difference is significant or not
    Emergence_DispStatTest(PWparams(iParam,:));
    
    % Add some text labels
    ylabel(['$', pwparamname{iParam}, '$'], 'Interpreter', 'LaTeX', ...
        'Rotation', 0, 'VerticalAlignment', 'Bottom');
    
    % Display the output of the statistical comparison against null value
    [~,pval,tci,stats] = ttest(PWparams(iParam,:)' - nullparamval(iParam));
    disptstats(pval,tci,stats);
end

% Save the figure
save2pdf('figs/F_QD_GpPWFit.pdf');

% Probability weighting function
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% See: https://en.wikipedia.org/wiki/Prospect_theory
% See: Gonzalez, R., & Wu, G. (1999). On the shape of the probability
%   weighting function. Cognitive psychology, 38(1), 129-166.
function [g, dgdx, dgdP] = g_PW(x, P, u, in)

% Make sure that the values at which to evaluate the function are within a structure
if ~isstruct(in)
    p = in;
    in = [];
    in.p = p;
end

% Get parameters
gamma = P(1);
p0 = P(2);

% Compute predictions of the probabilistic weighting function
g = logitinv(gamma .* logit(in.p) + (1 - gamma) * logit(p0));

% Rely on numerical derivatives
dgdx = [];
dgdP = [];
end

% General linear model
% ~~~~~~~~~~~~~~~~~~~~
% See: https://en.wikipedia.org/wiki/General_linear_model
function [g,dgdx,dgdP] = g_LM(x, P, u, in)

% Add an offset to the design matrix 
X = [in.p, ones(numel(in.p),1)];

% Compute predictions of the linear regression
g = X*P;

% Use analytical derivatives
dgdx = [];
dgdP = X';
end
