function [ D_JS, D_KL ] = Emergence_IO_DistDist( varargin )
%EMERGENCE_IO_DISTDIST measures the distance between distributions. It
%returns both the Jensen-Shannon (symmetrical) and Kullback-Leibler
%(non-symmetrical) divergences.
%   - Emergence_IO_DistDist(dist1, dist2): measures the distance between
%       the two distributions "dist1" and "dist2". They must have the same
%       size. It can be joint distributions (i.e. n-dimensional matrices).
%   - Emergence_IO_DistDist(dist): measures iteratively the distance
%       between the same distribution at time t and t+1. "dist" must be a
%       N-dimensional matrix whose dimensions correspond to:
%       - dim1 = grid (histogram-like distribution) over which the
%       distribution is defines
%       - dim2 = dimension over which to iterate (e.g. time-points),
%       - dimN (N>2) = other non-relevant dimensions that are not used by the
%       function but that must be conserved.
% 
% Copyright (c) 2018 Maxime Maheu

%% Initialization
%  ==============

% If two distribution of the same size are provided, compute the distance between them
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if numel(varargin) == 2 && ...
        numel(size(varargin{1})) == numel(size(varargin{2})) && ...
        all(size(varargin{1}) == size(varargin{2}))
    
    % Get (joint) distributions to be column vectors
    dist1 = varargin{1}(:);
    dist2 = varargin{2}(:);

% If a single distribution is provided, compute iterative update
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
else
    
    % Get the single input distribution
    dist = varargin{1};
    
    % Make 2 distributions out of it
    dist1 = dist(1:end-1,:); % previous beliefs
    dist2 = dist(2:end,:);   %  current beliefs
    flag = true;
end

% Remove zeros
if any(dist1 == 0), dist1(dist1 == 0) = NaN; end
if any(dist2 == 0), dist2(dist2 == 0) = NaN; end

%% Compute update
%  ==============

% Define Kullback-Leibler divergence
% !!! It is NOT symmetrical: D_KL(p,q) =/ D_KL(q,p))
% see https://en.wikipedia.org/wiki/Kullback%E2%80%93Leibler_divergence
KL = @(p,q) sum(p .* log2(p./q), 2, 'OmitNaN')';

% Define Jensen-Shannon divergence
% !!! It is symmetrical: D_JS(p,q) = D_JS(q,p))
% see https://en.wikipedia.org/wiki/Jensen%E2%80%93Shannon_divergence
JS = @(p,q) (1/2) .* KL(p, (p+q)/2) + (1/2) .* KL(q, (p+q)/2);

% Compute the Jensen-Shannon divergence
D_JS = JS(dist1, dist2);

% Compute the Kullback-Leibler divergence
if nargout > 1, D_KL = KL(dist1, dist2); end

% Add a first NaN such that the output variables have the same dimension as
% the input
if flag
    D_JS = cat(2, NaN, D_JS);
    if nargout > 1, D_KL = cat(2, NaN, D_KL); end
end

end
