function [ out, gro, idx ] = Emergence_Similarity( pb, qb, metric, optim )
% EMERGENCE_SIMILARITY returns similarity index between 2 matrices of
% posterior probability.
%   - "pb": first N*3 matrix of posterior probability.
%   - "qb": second N*3 matrix of posterior probability.
%   - "metric": a string of characters specifying what to return:
%       - "P*": Pearson correlation coefficient...
%       - "E*": Mean squared error...
%       - "D*": Euclidean distance...
%       - "R*": Mean squared error from univariate regressions...
%       - "M*": Mean squared error from multivariate regression...
%       - "L*": log-likelihood from multivariate regression...
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
[N,k] = size(pb);

% Which metric to use
if nargin < 3 || isempty(metric), metric = 'EC'; end
metric = upper(metric);
checkfun = @(x) any(x == metric);

% Convert to cartesian coordinates
if checkfun('C')
    tricc = [0, sqrt(3)/2; 1, sqrt(3)/2; 1/2, 0];
    pc = pb * tricc;
    qc = qb * tricc;
    k = 2; % 2 dimensions in that case
end

% Whether to optimize for temporal delay
if nargin < 4, optim = false; end
if islogical(optim)
    if      optim, dlist = 0:N;
    elseif ~optim, dlist = 0; 
    end
elseif isnumeric(optim)
    dlist = optim;
    if any(dlist < -N | dlist > N)
        dlist = dlist(dlist >= -N & dlist <= N);
    end
end
nd = numel(dlist);

% Prepare output variable
gro = NaN(nd,1);

% Shift one with respect to the other
for d = 1:nd
    m = dlist(d);
    if     checkfun('C'), [p,q] = Emergence_Shift(pc, qc, m);
    elseif checkfun('B'), [p,q] = Emergence_Shift(pb, qb, m);
    end
    n = N - abs(m); % number of remaining observations
    
    % Pearson correlation coefficient
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if checkfun('P')
        p = p(:); % column vector
        q = q(:); % column vector
        mp = sum(p) / (k*n); % average (over observations and dimensions)
        mq = sum(q) / (k*n); % average (over observations and dimensions)
        gro(d) = sum((p - mp) .* (q - mq)) / ... % correlation coefficient
            (sqrt(sum((p - mp) .^ 2)) * sqrt(sum((q - mq) .^ 2)));
        gro(d) = gro(d) / abs(m);
        
    % Mean squared error
    % ~~~~~~~~~~~~~~~~~~
    elseif checkfun('E')
        e = p - q; % signed error
        se = e.^2; % squared error
        sse = sum(se, 2); % sum of squared error (over dimensions)
        gro(d) = sum(sse) / n; % average squared error (over observations)
        
    % Euclidean distance
    % ~~~~~~~~~~~~~~~~~~
    elseif checkfun('D')
        e = p - q; % signed error
        se = e.^2; % squared error
        sse = sum(se, 2); % sum of squared error (over dimensions)
        ed = sqrt(sse); % euclidean distance
        gro(d) = sum(ed) / n; % average euclidean distance (over observations)
        
	% Univariate linear regression
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    elseif checkfun('R')
        e = NaN(size(p));
        for i = 1:k % for each dimension
            Y = p(:,i); % variable to explain
            X = [ones(n,1), q(:,i)]; % design matrix
            B = pinv(X) * Y; % regression coefficients (slower than X \ Y but handles singular values)
            Yhat = X * B; % predictions
            e(:,i) = Y - Yhat; % signed error
        end
        se = e.^2; % squared error
        sse = sum(se, 2); % sum of squared error (over dimensions)
        gro(d) = sum(sse) / n; % average squared error (over observations)
        
	% Multivariate linear regression
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % See http://users.stat.umn.edu/~helwig/notes/mvlr-Notes.pdf
    elseif checkfun('M') || checkfun('L')
        Y = p; % variable to explain
        X = [ones(n,1), q]; % design matrix
        B = inv(X' * X) * X' * Y; % regression coefficients
        Yhat = X * B; % predictions
        E = Y - Yhat; % signed error
        if checkfun('M')
            gro(d) = mean(E(:).^2); % average squared error (over observations)
        elseif checkfun('L')
            S = ((Y'-Yhat') * (Y'-Yhat')') ./ (n-k-1); % covariance matrix
            LLH = NaN(1,n); % log likelihood of each observation
            for i = 1:n
                LLH(i) = (Y(i,:)' - B'*X(i,:)')' * S * (Y(i,:)' - B'*X(i,:)');
            end
            gro(d) = -1/2 * sum(LLH);
        end
    end
end

% Whether to find the maximum or the minimum over shifts
if     checkfun('P'), funtooptim = @max; % maximum correlation coefficient
elseif checkfun('E'), funtooptim = @min; % minimum error
elseif checkfun('D'), funtooptim = @min; % minimum distance
elseif checkfun('R'), funtooptim = @min; % minimum error
elseif checkfun('M'), funtooptim = @min; % minimum error
elseif checkfun('L'), funtooptim = @max; % maximum likelihood
end
    
% Find the best temporal shift
[out,idx] = funtooptim(gro);

end
