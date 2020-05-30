function H = Emergence_IO_Entropy( dist, isindep, corr )
% EMERGENCE_IO_ENTROPY estimates Shannon's entropy from different
% distributions.
%   - "dist": a scalar, array or N-dimensional matrix specifying the
%       distribution for which the entropy must be estimated
%   - "isindep": (default: true) a boolean specifying whether, in case of
%       multidimensional distributions (i.e. N-dimensional matrix),
%       marginal distributions are independent from one another. Note that
%       if the marginal distributions are indeed independent from one
%       another, then Emergence_IO_Entropy(dist, true) should be equal to
%       Emergence_IO_Entropy(dist, false).
%   - "corr": (default: true) a boolean specifying whether to correct
%       pathological values such as 0s or Inf.
% 
% Copyright (c) 2020 Maxime Maheu

% By default, take care of pathological values.
if nargin < 3, corr = true; end

% Remove infinite values
infloc = isinf(dist);
if corr && any(infloc), dist(infloc) = NaN; end

% Remove zeros and ones
if corr
    zeroloc = (dist <= 0);
    oneloc  = (dist >= 1);
    if any(zeroloc(:)), dist(zeroloc) =   eps; end
    if any(oneloc(:)),  dist(oneloc)  = 1-eps; end
end

% Get the number of dimensions of the input
nDim = ndims(dist);
Sz = size(dist);

% Unidimensional distributions
if nDim <= 2 && any(Sz == 1)
    
    % 1) Bernoulli distribution
    % ~~~~~~~~~~~~~~~~~~~~~~~~~
    if all(Sz == 1) % scalar
        p = dist;
        H = -(p .* log2(p) + (1-p) .* log2(1-p));
    
    % 2) Unidimensional distribution
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    else % array
        dist = dist ./ sum(dist, 'OmitNaN'); % normalization
        H = -sum(dist .* log2(dist), 'OmitNaN');
    end
    
% Multidimensional distributions
else % N-D matrix
    if nargin < 2, isindep = false; end % by default, assume non-independency
    dist = dist ./ sum(dist(:), 'OmitNaN'); % normalization
    
    % 3) Independent multidimensional distribution
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if isindep
        H = NaN(1,nDim);
        for iDim = 1:nDim % for each (independent) dimension
            margdist = marg(dist, iDim); % marginalisation
            H(iDim) = -sum(margdist .* log2(margdist), 'OmitNaN');
        end
        H = sum(H);
        
    % 4) Non-independent multidimensional distribution
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    elseif ~isindep
        H = -sum(dist(:) .* log2(dist(:)), 'OmitNaN');
    end
end

% If there are only NaNs, do not return H = 0
if all(isnan(dist(:))), H = NaN; end

end

% Marginalisation function
function margdist = marg(dist, dim)
    nd = ndims(dist);
    v = 1:nd;
    v([1,dim]) = v([dim,1]);
    ditprime = reshape(permute(dist,v), size(dist,dim), []);
    margdist = squeeze(sum(ditprime, 2, 'OmitNaN'));
end
