function [ seqlh, post, predA, surp, predent, update, pe, postent ] = ...
    Emergence_IO_RunIO( IOfun, seq, inputs )
% EMERGENCE_IO_RUNIO runs iteratively a given ideal observer algorithm on a
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
pe      = NaN(1,N);   % |y_k+1 - p(y_k+1|y,M)|
postent = NaN(1,N);   % H(p(theta|y,M))

% The initial prediction is simply the probability of observing one
% observation or the other by chance, which is 1/2 because sequences are
% binary
predA(1) = 1/2;

% Loop over observations of the sequence
for k = 1:N
    
    % Get the estimated probability of receiving the current observation
    if     seq(k) == 1, pObs =     predA(k); % p(y_k = A)
    elseif seq(k) == 2, pObs = 1 - predA(k); % p(y_k = B) = 1 - p(y_k = A)
    end
    
    % Compute the entropy of the prediction
    predent(k)  = Emergence_IO_Entropy(pObs);
    
    % Compute theoretical Shannon surprise
    pe(k)   = 1 - pObs;   % prediction error
    surp(k) = -log(pObs); % Shannon surprise (simply log pe)
    
    % Update the beliefs of the ideal observer
    [seqlh(k), ...      % likelihood of the sequence up to the k-th observation
     predA(k+1), ...    % probability that the k+1-th observation will be a A
     post{k}, ...       % posterior probability distribution over model's parameters
     postent(k)] = ...  % entropy of the posterior distribution over model's parameters
     IOfun(seq(1:k), inputs{:});
end

% Discard the last prediction because we do not have any observation to
% compare with (the sequence is over)
predA = predA(1:end-1);

% Transform the posterior distribution into a matrix
post = cell2mat(reshape(post, [ones(1,ndims(post{end})), N]));
post = squeeze(post);

% Compute the iterative update induced by each observation. To do so, we
% measure the Jensen-Shannon divergence between posterior distributions
% before and after having received each observation:
% D(p(theta|y_1:K-1,M)||p(theta|y_1:K,M))
update = Emergence_IO_DistDist(post);

end
