function out = Emergence_Regress( y, x, method, outvar )
% EMERGENCE_REGRESS performs regression/correlation between x and y
% variables.
%   - "y": the variable to explain.
%   - "x": the explaining variable.
%   - "method": can be "OLS" (ordinary least squares), "TLS" (total least
%       squares) or "CC" (correlation coefficient).
%   - "out" (not required): can be "beta0", "beta1", "R2", "confint", ...
%       depending on the type of analysis that has been chosen and
%       specified in the "method" field.
%
% Copyright (c) 2018 Maxime Maheu

% Check that the explaining variable and the variable to explain have the
% same size
if numel(x) ~= numel(y), error('Inputs must have the same size'); end

% If only one output variable is provided, make sure it is in a cell array
if nargin > 3 && ~iscell(outvar), outvar = {outvar}; end

% Remove NaNs
ok = ~isnan(x) & ~isnan(y);
x = x(ok);
y = y(ok);

% Make sure these are column vectors
x = x(:);
y = y(:);

% Number of observations
n = numel(y);

% Ordinary least squares: predictor variables are measured exactly, and
% only the response variable has an error component
if strcmpi(method, 'OLS')

    % Create the design matrix
    X = [x, ones(n,1)];

    % General linear model
    [b,~,~,~,stats] = regress(y, X);
    beta0 = b(2);
    beta1 = b(1);
    R2 = stats(1);

    % 95% confidence interval
    [confint, confintx] = regerr(beta0, beta1, x, y);

    % Export relevant metrics in the output structure
    out = struct('beta0', beta0, 'beta1', beta1, 'R2', R2, ...
        'confint', confint, 'confintx', confintx);

% Orthogonal regression / total least squares mathod / Deming regression:
% both predictor and response variables are measured with error
elseif strcmpi(method, 'TLS')

    % 2D matrix
    X = [x,y];

    % For conveniency, we use PCA as an implementation of TLS regression
    [coeff,score,~,~,explained,~] = pca(X);

    % Compute projection onto the boundary
    meanX = mean(X, 1);
    dirVect = coeff(:,1);
    t = [min(score(:,1)), max(score(:,1))];
    endpts = [meanX + t(1)*dirVect'; meanX + t(2)*dirVect'];

    % Deduce beta coefficients
    beta1 = (endpts(2,2) - endpts(1,2)) / (endpts(2,1) - endpts(1,1));
    beta0 = endpts(1,2) - beta1*endpts(1,1);

    % Compute the determinant coefficient
    R2 = explained(1)/100;

    % 95% confidence interval
    [confint, confintx] = regerr(beta0, beta1, x, y);

    % Export relevant metrics in the output structure
    out = struct('beta0', beta0, 'beta1', beta1, 'R2', R2, ...
        'confint', confint, 'confintx', confintx);

% Simple correlation
elseif strcmpi(method, 'CC')

    % Correlation coefficient
    r = corr(x,y);

    % Output the asked variable
    out = struct('r', r, 'r2', r^2);
end

% Output only one or several variable (not the entire structure
out = cellfun(@(x) out.(x), outvar, 'UniformOutput', 0);
try out = cell2mat(out); % if possible return a scalar or an array
catch
end

end

% Nested function to compute error bars of regression lines
function [ yval, xval ] = regerr(beta0, beta1, x, y, ng, alpha)

if nargin < 5, ng = 100; end % precision of the grid
if nargin < 6, alpha = 0.05; end % significance threshold

n = length(x);
X = linspace(min(x), max(x), ng);
Y = ones(size(X))*beta0 + beta1*X;

SE_y_cond_x = sum((y - beta0*ones(size(y))-beta1*x).^2)/(n-2);
SSX = (n-1)*var(x);
SE_Y = SE_y_cond_x*(ones(size(X))*(1/n + (mean(x)^2)/SSX) + (X.^2 - 2*mean(x)*X)/SSX);

Yoff = (2*finv(1-alpha,2,n-2)*SE_Y).^0.5;
yval = Y + Yoff.*[-1;1];
xval = X;

end
