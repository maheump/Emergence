function [ seqlh, post, predA, surp, predent, update, pe, postent ] = ...
    Emergence_IO_RunMi( IOfun, seq, inputs )
% EMERGENCE_IO_RUNMI runs iteratively a given ideal observer algorithm on a
% sequential input. This allows to get updated beliefs of the ideal
% observer after each observation of the sequence, and not considering the
% entire sequence.
%   - "IOfun": the ideal observer function to be run iteratively
%   - "seq": a 1xN array specifying the sequence of binary observations (1s
%       and 2s).
%   - "inputs": a cell array specifying the options to pass to the ideal
%       observer function.
%
% Copyright (c) 2018 Maxime Maheu

% Avoid error when no inputs are provided
if nargin < 3, inputs = {}; end

% Number of observation in the sequence
N = numel(seq);

% Prepare the output variables
seqlh   = NaN(1,N);   % p(y|M)
post    = cell(1,N);  % p(theta|y,M)
predA   = NaN(1,N+1); % p(y_k+1=A|y,M)
surp    = NaN(1,N);   % -log p(y_k+1|y,M)
predent	= NaN(1,N);   % H(p(y_k+1=A|y,M))
update  = NaN(1,N);   % D(p(theta|y_1:K-1,M)||p(theta|y_1:K,M))
pe      = NaN(1,N);   % |y_k+1 - p(y_k+1|y,M)|
postent = NaN(1,N);   % H(p(theta|y,M))

% The initial prediction is simply the probability of observing one
% observation or the other by chance, which is 1/2 because sequences are
% binary
predA(1) = 1/2;

% Loop over observations of the sequence
for k = 1:N
    
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
    % Compute expectations-related quantities %
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

    % Identity of the current observation (1 for A, 2 for B)
    obs = seq(k);

    % Estimated probability of receiving this observation
    if     obs == 1, pObs =     predA(k); % p(y_k = A)
    elseif obs == 2, pObs = 1 - predA(k); % p(y_k = B) = 1 - p(y_k = A)
    end
    
    % Prediction error
    pe(k) = 1 - pObs;
    
    % Theoretical Shannon surprise (simply log pe)
    surp(k) = -log(pObs);

    % Entropy of the prediction
    predent(k)  = -(pObs .* log2(pObs) + (1-pObs) .* log2(1-pObs));

    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
    % Update the beliefs of the ideal observer %
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
    [seqlh(k), ...      % likelihood of the sequence up to the k-th observation
     predA(k+1), ...    % probability that the k+1-th observation will be a A
     post{k}, ...       % posterior probability distribution over model's parameters
     postent(k)] = ...  % entropy of the posterior distribution over model's parameters
     IOfun(seq(1:k), inputs{:});

    % ~~~~~~~~~~~~~~~~~~~~ %
    % Compute model update %
    % ~~~~~~~~~~~~~~~~~~~~ %
    if k > 1 && ~isempty(post{k})

    % Get probability distributions to compare
    cB = post{k}(:)';   % current beliefs
    pB = post{k-1}(:)'; % previous beliefs

    % Make sure there is no zeros or NaNs in the posterior distributions
    pB(isnan(pB)) = 1/numel(pB); % uniform prior
    cB(isnan(cB)) = 1/numel(cB); % uniform prior
    pB(pB == 0) = realmin; % smallest non-zero positive double
    cB(cB == 0) = realmin; % smallest non-zero positive double

    % Normalize the distributions
    pB = pB ./ sum(pB);
    cB = cB ./ sum(cB);

    % Compute the update between the current (updated) posterior and
    % the previous one (before receiving the current k-th observation)
    % using the Jensen-Shannon divergence
    update(k) = (1/2) .* sum(pB .* log2(pB ./ ((pB+cB)/2))) ...
              + (1/2) .* sum(cB .* log2(cB ./ ((pB+cB)/2)));
    end
end

% Discard the last prediction because we do not have any observation to
% compare with (the sequence is over)
predA = predA(1:end-1);

% Transform the posterior distribution into a matrix
post = cell2mat(reshape(post, [ones(1,ndims(post{end})), N]));
post = squeeze(post);

end
