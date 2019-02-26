function Emergence_PrintTstats( pval, ci, stats )
% EMERGENCE_PRINTTSTATS displays in the MATLAB command window the results
% of a two-tailed Student t-test(s) with alpha level = 0.05 (as computed
% with ttest.m).
%	- "pval": a 1xN array specifying the p-values.
%   - "ci": a 2xN matrix specifying the confidence intervals.
%   - "stats": a 1xN cell array specifying several metrics from the t-test.
% 
% Copyright (c) 2018 Maxime Maheu

% Get the number of tests' results to display
I = numel(pval);

% Define significiance thresholds
lima  = [0.001 0.01 0.05 1];
stars = {'***' '**' '*' 'ns'};

% Get group averages from the t-based confidence intervals
avgx = mean(ci);

% Get standard deviations from the t-based info
alpha = 0.05; % significiance level
n = stats.df + 1; % number of samples
stdx = ((avgx - ci(1,:)) ./ tinv((1 - alpha / 2), stats.df)) .* sqrt(n);

% Compute an estimate of effect size (Cohen's d)
cohend = avgx ./ stdx;

% Display statistical results in the command window
for i = 1:I
    idx = find(pval(i) < lima, 1, 'first');
    fprintf('m = %1.2f ± [%1.2f,%1.2f], d = %1.2f, t(%2.0f) = %1.2f, ', ...
        avgx(i), ci(:,i), cohend(i), stats.df(i), stats.tstat(i));
    if     pval(i) >= 0.001, fprintf('p = %1.3f %s\n', pval(i), stars{idx});
    elseif pval(i) <  0.001, fprintf('p = %1.2d %s\n', pval(i), stars{idx});
    end
end

end
