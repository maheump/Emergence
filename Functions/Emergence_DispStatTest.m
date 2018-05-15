function pval = Emergence_DispStatTest( y, ref )
% EMERGENCE_DISPSTATTEST displays on top of the current plot 
%   - "y": a NxK matrix where N is the number of subjects and K the number
%       of conditions. If K = 1, a t-test against zero (or what is
%       specified in "ref) is performed; if K = 2, a paired t-test is
%       performed; if K > 2 a 1-way repeated measures ANOVA is performed.
%   - "ref": the reference value against which "y" must be compared when
%       there is only 1 column (default is 0).
% 
% Copyright (c) 2018 Maxime Maheu

% Get the size of the input data matrix
[nSub,nCond] = size(y);

% If there is only 1 condition and the array is misoriented, take care of
% making it a column vector
if nSub == 1
    y = y(:);
    nSub = nCond;
    nCond = 1;
end

% By default, compare to zero
if nargin < 2, ref = 0; end

% Perform the appropriate statistical test
if nCond == 1 % against a mean value
    [~,pval] = ttest(y - ref);
elseif nCond == 2 % paired t-test
    [~,pval] = ttest(diff(y, 1, 2));
elseif nCond > 2 % 1-way repeated measures ANOVA
    RMtbl = rmANOVA(y);
    pval = RMtbl.pValue(3); % interaction intercept * factor 1
end

% Define significiance thresholds
lima  = [0.001 0.01 0.05 0.1 1];
stars = {'***' '**' '*' 'o' 'ns'};

% Find the lowest threshold that is reached
pidx = find(pval < lima, 1, 'first');

% Display the sifnificant stars on the current plot
limy = ylim;
hold('on');
text(mean(xlim), max(limy), stars{pidx}, ...
    'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle');
ylim(limy);

end