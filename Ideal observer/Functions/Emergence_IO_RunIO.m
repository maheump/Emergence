function [ seqlh, predA, predent, pe, surp, margpost, jointpost, postent, update ] = ...
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
seqlh    = NaN(1,N);   % p(y|M)
margpost = cell(1,N);  % p(theta|y,M)
predA    = NaN(1,N+1); % p(y_k=A|y_1:k-1,M)
surp     = NaN(1,N);   % -log p(y_k|y_1:k-1,M)
predent	 = NaN(1,N);   % H(p(y_k=A|y,M))
pe       = NaN(1,N);   % |y_k - p(y_k|y,M)|
postent  = NaN(1,N);   % H(p(theta|y,M))

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
    predent(k) = Emergence_IO_Entropy(pObs);
    
    % Compute theoretical Shannon surprise
    pe(k)   = 1 - pObs;   % prediction error
    surp(k) = -log(pObs); % Shannon surprise (simply log pe)
    
    % Update the beliefs of the ideal observer
    [seqlh(k), ...      % likelihood of the sequence up to the k-th observation
     predA(k+1), ...    % probability that the k+1-th observation will be a A
     margpost{k}, ...   % posterior probability distribution over model's parameters
     postent(k)] = ...  % entropy of the posterior distribution over model's parameters
     IOfun(seq(1:k), inputs{:});
end

% Discard the last prediction because there is no observation to compare it
% with (the sequence is over)
if nargout > 1, predA = predA(1:end-1); end

% If posterior distributions are marginal distributions, we assume
% they are independend ones and derive the joint posterior distribution
% by making a product out of the different dimensions (i.e. the
% different marginal distributions)
nTdim = size(margpost{1},2);
if nTdim > 1
    if nargout > 5
        if nargout > 6
            jointpost = cellfun(@indepmarg2joint, margpost, 'UniformOutput', 0);
            jointpost = cell2mat(reshape(jointpost, [ones(1,nTdim), N]));
        end
        margpost  = permute(squeeze(cell2mat(reshape(margpost,  [ones(1,nTdim), N]))), [1,3,2]);
    end
elseif nTdim == 1
    if nargout > 5
        margpost  = cell2mat(margpost);
        if nargout > 6, jointpost = margpost; end
    end        
end

% Compute the iterative update induced by each observation. To do so,
% we measure the Jensen-Shannon divergence between posterior
% distributions before and after having received each observation:
% D(p(theta|y_1:K-1,M)||p(theta|y_1:K,M))
if nargout > 8
    update = Emergence_IO_DistDist(permute(jointpost, [nTdim+1, 1:nTdim]));
end

end

% Useful function to turn independent marginal distributions into a joint
% distribution
function joint = indepmarg2joint( marg )

[ngrid,ndim] = size(marg);
idim = 1:ndim;
step1 = arrayfun(@(x) repmat(marg(:,x), [1,repmat(ngrid,1,ndim-1)]), idim, 'UniformOutput', 0);
step2 = reshape(step1, [ones(1,ndim),ndim]);
joint = prod(cell2mat(step2), ndim+1);

end
