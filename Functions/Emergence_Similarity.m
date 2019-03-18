function [ out, all, idx ] = Emergence_Similarity( pb, qb, metric, optim )
% EMERGENCE_SIMILARITY returns similarity index between 2 matrices of
% posterior probability.
%   - "pb": first N*3 matrix of posterior probability.
%   - "qb": second N*3 matrix of posterior probability.
%   - "metric": a string of characters specifying what to return:
%       - "P*": Pearson correlation coefficient...
%       - "E*": Mean squared error...
%       - "D*": Euclidean distance...
%       - "U*": Univariate log-likelihood...
%       - "M*": Multivariate log-likelihood...
%       - "*B": ...on Barycentric coordinates
%       - "*C": ...on Cartesian coordinates
%   - "optim": either a boolean specifying whether to optimize for temporal
%       delay or an array of integers specifying the temporal delays over
%       which to optimize (in terms of observation #, see Emergence_Shift).
% 
% Copyright (c) 2018 Maxime Maheu

% Make sure the variables are N*3 matrices
if isstruct(pb), pb = pb.BarycCoord; end
if isstruct(qb), qb = qb.BarycCoord; end
if any(size(pb) ~= size(qb)), error('Inputs must have the same size'); end
N = size(pb,1);

% Which metric to use
if nargin < 3 || isempty(metric), metric = 'MC'; end
checkfun = @(x) contains(metric, x, 'IgnoreCase', true);

% Convert to cartesian coordinates
if checkfun('C')
    tcn = [0, sqrt(3)/2; 1, sqrt(3)/2; 1/2, 0];
    pc = pb * tcn;
    qc = qb * tcn;
end

% Whether to optimize for temporal delay
if nargin < 4, optim = false; end
if islogical(optim)
    if      optim, dlist = 0:N;
    elseif ~optim, dlist = 0; 
    end
elseif isnumeric(optim)
    dlist = optim;
    if any(dlist < -N | dlist > N), error('Check the shift array.'); end
end
nd = numel(dlist);

% Get the number of dimensions
k = size(pc, 2);

% Prepare output variable
all = NaN(nd,1);

% Shift one with respect to the other
for d = 1:nd
    curd = dlist(d);
    if     checkfun('C'), [p,q] = Emergence_Shift(pc, qc, curd);
    elseif checkfun('B'), [p,q] = Emergence_Shift(pb, qb, curd);
    end
    n = N - abs(curd); % number of remaining observations
    
    % Pearson correlation coefficient
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if checkfun('P')
        p = p(:); % column vector
        q = q(:); % column vector
        mp = sum(p) / (k*n); % average (over observations and dimensions)
        mq = sum(q) / (k*n); % average (over observations and dimensions)
        all(d) = sum((p - mp) .* (q - mq)) / ... % correlation coefficient
            (sqrt(sum((p - mp) .^ 2)) * sqrt(sum((q - mq) .^ 2)));
        
    % Mean squared error
    % ~~~~~~~~~~~~~~~~~~
    elseif checkfun('E')
        e = p - q; % signed error
        se = e.^2; % squared error
        sse = sum(se, 2); % sum of squared error (over dimensions)
        all(d) = sum(sse) / n; % average squared error (over observations)
        
    % Euclidean distance
    % ~~~~~~~~~~~~~~~~~~
    elseif checkfun('D')
        e = p - q; % signed error
        se = e.^2; % squared error
        sse = sum(se, 2); % sum of squared error (over dimensions)
        ed = sqrt(sse); % euclidean distance
        all(d) = sum(ed) / n; % average euclidean distance (over observations)
        
    % Log-likelihood under univariate normally distributed errors
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    elseif checkfun('U')
        p = p(:); % column vector
        q = q(:); % column vector
        e = p - q; % signed error
        s2 = (1/n) * sum(e.^2); % variance with mu = 0
        all(d) = - n/2 * log(2*pi) - n/2 * log(s2) - 1/(2*s2) * sum(e.^2);
        % log-likelihood of the univariate distance
        
    % Log-likelihood under multivariate normally distributed errors
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    elseif checkfun('M')
        e = p - q; % signed error
        S = (1/n) * (e' * e); % covariance matrix
        all(d) = -1/2 * sum(log(det(S)) + sum(e * inv(S) .* e, 2) + k * log(2*pi));
        % log-likelihood of the bivariate distance
        % N.B. It is equivalent to the univariate case when e is a vector, 
        % and thus k = 1.
    end    
end

% Whether to find the maximum or the minimum over shifts
if     checkfun('P'), funtooptim = @max; % maximum correlation coefficient
elseif checkfun('E'), funtooptim = @min; % minimum error
elseif checkfun('D'), funtooptim = @min; % minimum distance
elseif checkfun('U'), funtooptim = @max; % maximum log-likelihood
elseif checkfun('M'), funtooptim = @max; % maximum log-likelihood
end
    
% Find the best temporal shift
[out,idx] = funtooptim(all);

end
