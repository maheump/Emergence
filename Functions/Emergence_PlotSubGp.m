function Emergence_PlotSubGp( y, col, xon, ms, aval )
% EMERGENCE_PLOTSUBGP plots subjects- and group-average values of a given
% variable in several conditions.
%   - "y": a data matrix with rows specifying the subjects and columns
%       specifying the conditions.
%   - "col": a Nx3 matrix specifying the RGB colors to use. N can be equal
%       to 1. In that case, all conditions are plotted using the same
%       color. N can be equal to the number of conditions (i.e. of columns
%       of "y"). In that case, each condition is plotted using its own
%       dedicated color.
%   - "xon": a scalar value specifying where to plot the data on the
%       x-axis.
%   - "ms": the marker size for group-average dots.
%   - "aval": the level of transparency for subject-level data.
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

% Fill in the inputs
if nargin < 2 || isempty(col), col = 'k'; end % black
if nargin < 3 || isempty(xon), xon = 0; end % no global shift
if nargin < 4 || isempty(ms), ms = 10; end % 10 points markers
if nargin < 5 || isempty(aval), aval = 0.15; end % transparency

% Wether all conditions should be displayed using the same colors or not
ncol = size(col,1);

% Create a vector of horizontal shifts to display group distribution more
% conveniently
rng(1, 'twister'); % for reproducibility
xshift = (rand(nSub,1)-1/2)./4;

% Vector of x-axis positions
o = (1:nCond) + xon;
x = o + xshift;

% Average over subjects
m = mean(y, 1, 'OmitNaN');
s = sem(y, 1);

% Do not erase what is already displayed on the plot
hold('on');

% Display individual differences
for iSub = 1:nSub
    p = plot(x(iSub,:), y(iSub,:), 'ko', ...
        'MarkerEdgeColor', 'none', 'LineStyle', '-', 'LineWidth', 1/2);
    set(p, 'Color', [get(p, 'Color'), aval]);
end

% Display individual data points
if ncol == 1
    scatter(x(:), y(:), 'MarkerFaceColor', col, ...
        'MarkerFaceAlpha', aval, 'MarkerEdgeColor','none');
elseif ncol > 1
    for i = 1:ncol
        scatter(x(:,i), y(:,i), 'MarkerFaceColor', col(i,:), ...
            'MarkerFaceAlpha', aval, 'MarkerEdgeColor','none');
    end 
end

% Display group difference
plot(o, m, 'k-', 'LineWidth', 3);
 
% Display dispersion over subjects
plot(repmat(o,2,1), m+s.*[-1;1], 'k-');

% Display group average
if ncol == 1
    plot(o, m, 'ko', 'MarkerFaceColor', col, 'MarkerSize', ms);
elseif ncol > 1
    for i = 1:ncol
        plot(o(i), m(i), 'ko', 'MarkerFaceColor', col(i,:), 'MarkerSize', ms);
    end
end

end