function [ seqLH, post, predA, surp, entr, update ] = ...
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

% Number of observation in the sequence
N = numel(seq);

% Prepare the output variables
seqLH  = NaN(1,N);
post   = cell(1,N);
predA  = NaN(1,N);
surp   = NaN(1,N);
entr   = NaN(1,N);
update = NaN(1,N);

% Loop over observations of the sequence
for k = 1:N
    
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
    % Compute expectations-related quantities %
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

    % Identity of the current observation (1 for A, 2 for B)
    obs = seq(k);

    % Estimated probability of receiving this observation
    if k == 1, pObs = 1/2; % chance
    elseif k > 1 % prediction derived from the most recent inference
        if     obs == 1, pObs =     predA(k-1); % p(y_k = A)
        elseif obs == 2, pObs = 1 - predA(k-1); % p(y_k = B) = 1 - p(y_k = A)
        end
    end

    % Theoretical Shannon surprise
    surp(k) = -log2(pObs);

    % Entropy of the prediction
    entr(k)  = -(pObs .* log2(pObs) + (1-pObs) .* log2(1-pObs));

    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
    % Update the beliefs of the ideal observer %
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
    [seqLH(k), ...  % likelihood of the sequence up to the k-th observation
     predA(k), ...  % probability that the k+1-th observation will be a A
     post{k}] = ... % posterior probability distribution over model's parameters
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
    pB(pB == 0) = eps; % smallest non-zero positive double
    cB(cB == 0) = eps; % smallest non-zero positive double

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

% Transform the posterior distribution into a matrix
post = cell2mat(reshape(post, [ones(1,ndims(post{end})), N]));
post = squeeze(post);

end
