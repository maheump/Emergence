function [ pY, pAgY, pRigY, HpRigY, R ] = Emergence_IO_Tree( y, nu, scaleme, usegrid, prior, decw, corout )
% EMERGENCE_IO_TREE implements an observer learning repeating patterns from
% a sequence of binary observations.
%   - "y": a 1xN array specifying the sequence of binary observations (1s
%       and 2s).
%   - "nu": a scalar specifying the depth of the tree to explore.
%   - "scaleme": a string ('lin' or 'log') specifying whether the model
%       evidence sould be computed on a linear or logarithmic scale.
%   - "usegrid": a boolean specifying whether to use grid-based or analytical
%       solutions.
%   - "prior": the name of the prior distribution over patterns to use,
%       either 'Bayes-Laplace' or based on the 'Size-principle'.
%   - "decw": a Nx1 or 1xN array of (decaying) weights that will weight
%       past observations of the sequence.
%   - "corout": a boolean specifying whether to correct the output
%       variables when the result of the inference is singular.
%
% Copyright (c) 2018 Maxime Maheu

% In order to understand the nomenclature used troughout the function,
% let's consider the following example
% EXAMPLE: if y = [ A A B ] and nu = 5, then the tree is the following one:
% ===================================================================|==> Here X = A (thus Y = B)
% |                 _______________[X]_______________                |==> K = 1
% |                 |                               |                |
% |         _______[X]_______               ________Y_______         | Entirely observed
% |         |               |               |              |         | patterns {Ro} = 2^3
% |     ____X____       ___[Y]___       ____X____      ____Y____     |
% |     |       |       |       |       |       |      |       |     |
% |=====|=======|=======|=======|=======|=======|======|=======|=====|==> K = 3
% |                     |       |                                    |
% |                   _(X)_   _(Y)_                                  | Partialy observed
% |                   |   |   |   |                                  | patterns: {Ru} = 2^3
% |                  (X) (Y) (X) (Y)                                 |
% ===================================================================|==> nu = 5

%% INITIALIZATION
%  ==============

% Make sure the sequence of observations is organised as a column vector.
% This allows to rely on linear algebra below, which speed up the
% computations.
y = y(:);

% Number of observation in the sequence
K = numel(y);
% N.B. K = K-k if the sequence is a part after a change point.

% By default, the depth of the tree (i.e. the maximum pattern length allowed)
% is the length of the sequence
if nargin < 2 || isempty(nu), nu = K; end
% N.B. Deep trees induce longer computation time

% By default, export the model evidence in log-scale
if nargin < 3, scaleme = 'log'; end
if     strcmpi(scaleme, 'lin'), islin = true;  islog = false;
elseif strcmpi(scaleme, 'log'), islin = false; islog = true;
end

% By default, use analytical solutions to speed up the computations
if nargin < 4, usegrid = false; end

% By default, use a prior distribution based on the size principle
if nargin < 5 || isempty(prior), prior = 'Size-principle'; end

% By default, use a perfect integration
if nargin < 6 || isempty(decw), decw = ones(1,K); % no decay
else, decw = decw(:); % make sure it is a column vector (for linear algebra below)
end
if numel(decw) > K, decw = decw(end-K+1:end); % take only the most recent weights
elseif numel(decw) == K % nothing to do in that case
elseif numel(decw)  < K, error('Not enough weights for the length of the sequence');
end
if     any(decw ~= 1), idealinteg = false;
elseif all(decw == 1), idealinteg = true;
end

% By default, make sure that the output is corrected when the result of
% the inference is singular. This happens only in some rare specific cases
% (e.g. when the tree has been fully explored and no patterns have been found)
if nargin < 7 || isempty(corout), corout = true; end

%% GET CHARACTERISTICS OF THE TREE
%  ===============================

% The number of possible patterns depends on the sequence's length and the
% depth of the tree to explore
nlRo = min([nu, K]); % number or levels in the   observed part of the tree
nlRu = nu - nlRo;    % number or levels in the unobserved part of the tree

% Get levels in both parts of the tree
ilRo = 1:nlRo; % levels' indices in the   observed part of the tree
ilRu = 1:nlRu; % levels' indiced in the unobserved part of the tree
% (!!! they are different from the true levels' indices which can be
% computed using (nlRo+1):nu).

% Get the total number of patterns
if usegrid
    nRo = 2.^(ilRo - 1);      % number of observed patterns
    nRu = 2.^ilRu;            % number of unobserved patterns
    nR = sum(nRo) + sum(nRu); % total number of patterns
elseif ~usegrid
    nRo = 2^nlRo - 1;         % number of observed patterns
    nRu = 2 * (2^nlRu - 1);   % number of unobserved patterns
    nR = nRo + nRu;           % total number of patterns
end

%% PATTERNS' LIKELIHOOD
%  ====================

% Prepare output variables for the patterns
if nargout > 3, R = cell(1,nlRo); end % if patterns have to be returned
pYgRio  = NaN(1,nlRo); % probability of the pattern given the sequence
pXgYRio = NaN(1,nlRo); % expected identity of the next item given the patterns

% For each pattern's length (i.e. each level of the tree)
for i = ilRo
    
    % There is only one correct pattern per level in the tree, get that
    % pattern for the current level
    Ri = y(1:i); % pattern Ri
    if nargout > 3, R{i} = Ri; end % return the pattern
    
    % Check whether the observed sequence could be a fixed (integer) number
    % of repetition(s) of that pattern. This allows to use the following
    % conditionals which speed up the computations.
    ncol = ceil(K/i);
    ismult = (mod(K/i, 1) == 0);
    
    % Get the observed sequence...
    % - If the length of the sequence is not a fixed number of the pattern's
    % length, add zeros such that it becomes the case
    if ~ismult
        j = (ncol*i)-K;
        ymat = zeros(K+j,1);
        ymat(1:K) = y;
    % - If the length of the sequence is a fixed number of the repetition
    % of the pattern, simply get the sequence
    elseif ismult
        ymat = y;
    end
    
    % Transform the array specifying the observed sequence into a matrix
    % whose number of rows correspond to the length of the current pattern
    % and the number of columns to the number of possible repetitions of
    % that pattern
    ymat = reshape(ymat, [i, ncol]);
    
    % Compare each observation of the sequence to each observation of a
    % sequence obtained from the repetition of the current pattern
    seqcheck = ymat == Ri; % a 2D matrix
    % N.B. This is more efficient than "repmat"ting the pattern but it
    % works only on the most recent versions of MATLAB
    
    % In case of a perfect integration
    if idealinteg
        
        % Check if the repetition of that current pattern entirely matches
        % the observed sequence
        isruletrue = sum(seqcheck(:)) == K;
        pYgRio(i) = (1/2) .* isruletrue;
        % N.B. the probability of observing the sequence given a particular
        % pattern. Because patterns are defined based on Xs/Ys, the first
        % item can be an A or a B with probability 1/2
        
    % In case of non-uniform integration, e.g. a memory leak implemented as
    % observations' weights indexed on an exponentialy decaying function
    elseif ~idealinteg
        
        % Get the probability (based on the memory leak) with which each
        % observation of the sequence matches the one predicted from the
        % repetition of that pattern
        pYkgRi = NaN(1,K);
        x = seqcheck(1:K);
        ok = (x == 1); %    match between pattern and sequence
        er = (x == 0); % mismatch between pattern and sequence
        pYkgRi(ok) =     decw(ok); %    match
        pYkgRi(er) = 1 - decw(er); % mismatch
        
        % Compute how much the repetition of that current pattern matches
        % the observed sequence
        % p(Ri|y) = (1/2) * prod_(y_k == yhat_k) w_k * prod_(y_k ~= yhat_k) (1-w_k)
        if     islin, pYgRio(i) = (1/2) * prod(pYkgRi);
        elseif islog, pYgRio(i) = -log(2) + sum(log(pYkgRi));
        end
    end
    
    % Derive the expected identity of the next observation based on that
    % pattern (be it correct or not)
    if    ~ismult, pXgYRio(i) = Ri(i-j+1);
    elseif ismult, pXgYRio(i) = Ri(1);
    end
end

%% PRIOR PROBABILITIES
%  ===================

% In the case of a uniform (Bayes-Laplace) prior distribution
if strcmpi(prior, 'Bayes-Laplace')
    
    % The prior probability is simply 1 over the total number of patterns
    % (i.e. both fully and partially observed ones)
    % p(Ri) = 1/{R}
    if     idealinteg || ~idealinteg && islin, p_pRi = ones(1, nu) ./ nR;
    elseif               ~idealinteg && islog, p_pRi = repmat(-log(nR), 1, nu);
    end
    
% In the case of a prior distribution based on the size-principle
elseif strcmpi(prior, 'Size-principle')
    
    % We want to satisfy 2 constraints:
    %   1) We want a pattern whose length equals i+1 to be 1/3 more likely
    %      (a priori) than a pattern whose length equals i.
    %   2) We want the sum of prior probabilities over all (entirely
    %      observed and partialy observed) patterns equal to 1.
    
    % To satisfy these 2 constraints, we define the prior probability of
    % patterns from the first level of the tree as follows
    if     idealinteg || ~idealinteg && islin, p_pR1 = 1 / ((9/2) * (1 - (2/3)^(nu+1)) - (3/2));
    elseif               ~idealinteg && islog, p_pR1 = -log((9/2) * (1 - (2/3)^(nu+1)) - (3/2));
    end
    
    % We then deduce the prior probability of all patterns
    if     idealinteg || ~idealinteg && islin, p_pRi =  (1/3).^(0:nu-1)   .* p_pR1;
    elseif               ~idealinteg && islog, p_pRi = -(0:nu-1) .* log(3) + p_pR1;
    end
end

% Split the prior distribution into two: a first distribution for entirely
% observed patterns and a second distribution for partially observed ones
p_pRio = p_pRi(1:nlRo);
p_pRiu = p_pRi((nlRo+1):nu);

%% SEQUENCE'S MARGINAL LIKELIHOOD
%  ==============================

% Entirely observed patterns
% ~~~~~~~~~~~~~~~~~~~~~~~~~~

% Posterior probability of each observed pattern
% p(Rio|y) propto p(y|Rio) * p(Ri)
if      idealinteg || ~idealinteg && islin, pRiogY =     pYgRio .* p_pRio;
elseif                ~idealinteg && islog, pRiogY = exp(pYgRio  + p_pRio);
end

% Sum of posterior probabilities
% p(y|Ro) = sum(i=1:min([K,nu])) p(y|Rio) * p(Ri)
pYgRo = sum(pRiogY);

% Partially observed patterns
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Posterior probability of each partially observed pattern
% p(y|Riu) = (1/2) * p(Ri)
if     idealinteg || ~idealinteg && islin, pRiugY =   (1/2) .* p_pRiu;
elseif               ~idealinteg && islog, pRiugY = -log(2)  + p_pRiu;
end

% Sum of posterior probabilities
% p(y|Ru) = sum(i=1:{Ru}) (1/2) * p(Ri) = sum(i=1:(nu-K)}) 2^i * (1/2) * p(Ri)
if     idealinteg || ~idealinteg && islin, pYgRu = sum(2.^(1:nlRu) .* pRiugY);
elseif               ~idealinteg && islog, pYgRu = sum(exp((1:nlRu).*log(2) + pRiugY));
end

% All patterns together
% ~~~~~~~~~~~~~~~~~~~~~

% Compute the probability of observing the sequence under a deterministic
% model
% p(y|Md) = p(y|Ro) * p(Ro) + p(y|Ru) * p(Ru)
pY = pYgRo + pYgRu;
if islog, pY = log(pY); end

% If none of the observed patterns can explain the sequence and the
% sequence is longer than the depth of the tree (i.e. the longest possible
% pattern that is considered), the likelihood of the sequence is null. When
% exported in log scale, we return the log of the smallest postitive
% floating point number.
if corout && islog && isinf(pY), pY = log(realmin); end

%% PREDICTION
%  ==========

% If the expectation has to be returned
if nargout > 1
    
    % Conditionaly on entirely observed patterns
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    % Get predictions given each observed patterns
    % p(A|y,Rio) E {0,1}
    pAgYRio = double(pXgYRio == 1);
    
    % Marginalize predictions over all observed patterns
    % p(A|y,Ro) = sum_(i=1:{Ro}) p(A|y,Rio) * p(Rio|y)
    pAgYRo = sum(pAgYRio .* pRiogY) / pYgRo;
    
    % Conditionaly on partialy observed patterns
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    % The likelihood that the next observation would be an A conditioned on the
    % yet unonbserved patterns is simply chance level
    % p(A|Riu) = 1/2
    pAgYRu = 1/2;
    
    % Combine both expectancies
    % ~~~~~~~~~~~~~~~~~~~~~~~~~
    
    % Compute the likelihood that the next observation will be a A by means of
    % Bayesian Model Averaging
    % p(A|y) = (p(A|y,Ro)*p(Ro|y) + p(A|y,Ru)*p(Ru|y)) / (p(Ro|y) + p(Ru|y))
    pAgY = (pAgYRo*pYgRo + pAgYRu*pYgRu) / (pYgRo + pYgRu);
    
    % If none of the observed patterns can explain the sequence and the
    % sequence is longer than the depth of the tree (i.e. the longest possible
    % pattern that is considered), we return chance level for the expectancy of
    % the next observation.
    if corout && (pYgRo + pYgRu) == 0, pAgY = 1/2; end
end

%% ENTROPY OF THE POSTERIOR DISTRIBUTION
%  =====================================

% If the posterior distribution over patterns has to be returned
if nargout > 2
    
    % Combine both the posterior distributions over entirely and partially
    % observed patterns
    if     idealinteg || ~idealinteg && islin, ppRigY = [pRiogY,     pRiugY ];
    elseif               ~idealinteg && islog, ppRigY = [pRiogY, exp(pRiugY)];
    end
    
    % Normalize the posterior distribution
    pRigY = ppRigY(:) ./ sum(ppRigY);
    
    % If the sequence is longer than the longest allowed pattern (i.e. the
    % depth of the tree) and that no patterns can reproduce the sequence,
    % (i.e. we cannot normalize the posterior because it would mean divide
    % it by zero), return an array of zero
    if all(isnan(pRigY)), pRigY = zeros(nlRo,1); end
end

% If the entropy of the posterior distribution has to be returned
if nargout > 3
    HpRigY = Emergence_IO_Entropy(pRigY);
end

% Fill the output variable (for the partially observed patterns) such that
% its size always equals the depth of the tree (i.e. nu)
if nargout > 4, R = [R, repmat({''}, 1, nlRu)]; end

end
