function [ pYgMd, pAgYMd, pRigY, HpRigY, R ] = Emergence_IO_Tree( y, nu, scaleme, usegrid, prior, decw, corout )
% EMERGENCE_IO_TREE implements an observer learning repeating patterns from
% a sequence of binary observations.
%   - "y": a 1xN array specifying the sequence of binary observations (1s
%       and 2s).
%   - "nu": a scalar specifying the depth of the tree to explore.
%   - "scaleme": a string ('lin' or 'log') specifying whether the model
%       evidence sould be computed on a linear or logarithmic scale.
%   - "usegrid": a boolean specifing whether use grid-based or analytical
%       solutions.
%   - "prior": the name of the prior distribution over patterns to use,
%       either 'Bayes-Laplace' or based on the 'Size-principle'.
%   - "decw": a Nx1 or 1xN array of (decaying) weights that will weight
%       past observations of the sequence.
%   - "corout": a boolean specofying whether to correct the output
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
    
%% Initialization
%  ==============

% Make sure the sequence of observations is organised as a column vector.
% This allows to rely on linear algebra below, which speed up the
% computations.
y = y(:);

% Number of observation in the sequence
K = numel(y);
% N.B. K = K-k if the sequence is a part after a change point.

% By default, make sure that the output is corrected when the result of
% the inference is singular. This happens only in some rare specific cases
% (e.g. when the tree has been fully explored and no patterns have been found)
if nargin < 7 || isempty(corout), corout = true; end

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

% By default, use a prior distribution based on the size principle
if nargin < 5 || isempty(prior), prior = 'Size-principle'; end

% By default, use analytical solutions to speed up the computations
if nargin < 4, usegrid = false; end

% By default, export the model evidence in log-scale
if nargin < 3, scaleme = 'log'; end

% By default, the depth of the tree (i.e. the maximum pattern length allowed)
% is the length of the sequence
if nargin < 2, nu = K; end
% N.B. Deep trees induce longer computation time

%% Get the dimensionality of the problem
%  =====================================

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

%% Patterns' likelihood
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
    
    % Check whether the observed sequence could be a fixe (integer) number
    % of repetition(s) of that pattern. This allows to use the following
    % conditionals which speed up the computations.
    ncol = ceil(K/i);
    ismult = (mod(K/i, 1) == 0);
    
    % Get the observed sequence...
    % - If the length of the sequence is not a fixe number of the pattern's
    % length, add zeros such that it becomes the case
    if ~ismult
        j = (ncol*i)-K;
        ymat = zeros(K+j,1);
        ymat(1:K) = y;
    % - If the length of the sequence is a fixe number of the repetition
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
        if     strcmpi(scaleme, 'lin'), pYgRio(i) =  (1/2) * prod(pYkgRi);
        elseif strcmpi(scaleme, 'log'), pYgRio(i) = -log(2) + sum(pYkgRi);
        end
   end
    
    % Derive the expected identity of the next observation based on that
    % pattern (be it correct or not)
    if    ~ismult, pXgYRio(i) = Ri(i-j+1);
    elseif ismult, pXgYRio(i) = Ri(1);
    end
end

%% Prior probabilities
%  ===================

% In the case of a uniform (Bayes-Laplace) prior distribution
if strcmpi(prior, 'Bayes-Laplace')
    
    % The prior probability is simply 1 over the total number of patterns
    % (i.e. both fully and partially observed ones)
    % p(Ri) = 1/{R}
    if     idealinteg || ~idealinteg && strcmpi(scaleme, 'lin'), p_pRi = ones(1, nu) ./ nR;
    elseif               ~idealinteg && strcmpi(scaleme, 'log'), p_pRi = repmat(-log(nR), 1, nu);
    end
    
% In the case of a prior distribution based on the size-principle
elseif strcmpi(prior, 'Size-principle')
    
    % The size principle states that, since these are binary patterns,
    % their prior simply depends on their lengths:
    % p(pRi) = (1/3) .^ |Ri|;
    if     idealinteg || ~idealinteg && strcmpi(scaleme, 'lin'), p_pRi = (1/3) .^ (1:nu);
    elseif               ~idealinteg && strcmpi(scaleme, 'log'), p_pRi = -(1:nu) .* log(3);
    end
end

% Split the prior distribution into two: for observed versus unobserved
% patterns
p_pRio = p_pRi(1:nlRo);
p_pRiu = p_pRi((nlRo+1):nu);

%% Sequence's marginal likelihood
%  ==============================

% Entirely observed patterns
% ~~~~~~~~~~~~~~~~~~~~~~~~~~

% Posterior probability of each observed pattern
% p(Rio|y) propto p(y|Rio) * p(Ri)
if      idealinteg || ~idealinteg && strcmpi(scaleme, 'lin'), pRiogY =     pYgRio .* p_pRio;
elseif                ~idealinteg && strcmpi(scaleme, 'log'), pRiogY = exp(pYgRio  + p_pRio);
end

% Sum of posterior probabilities
% p(Ro|y) = sum(i=1:min([K,nu])) p(y|Rio) * p(Ri)
pRogY = sum(pRiogY); 

% Partially observed patterns
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Posterior probability of each partially observed pattern
% p(Rio|y) propto (1/2) * p(Ri)
if     idealinteg || ~idealinteg && strcmpi(scaleme, 'lin'), pRiugY =    (1/2) .* p_pRiu;
elseif               ~idealinteg && strcmpi(scaleme, 'log'), pRiugY = - log(2)  + p_pRiu;
end

% Sum of posterior probabilities
% p(Ru|y) = sum(i=1:{Ru}) (1/2) * p(Ri) = sum(i=1:(nu-K)}) 2^i * (1/2) * p(Ri)
if     idealinteg || ~idealinteg && strcmpi(scaleme, 'lin'), pRugY = sum(2.^(1:nlRu) .* pRiugY);
elseif               ~idealinteg && strcmpi(scaleme, 'log'), pRugY = sum(exp((1:nlRu).* log(2) + pRiugY));
end

% All patterns together
% ~~~~~~~~~~~~~~~~~~~~~

% Compute the probability of observing the sequence under a deterministic
% model
% p(y|Md) = p(y|Ro) * p(Ro) + p(y|Ru) * p(Ru)
pYgMd = pRogY + pRugY;
if idealinteg && strcmpi(scaleme, 'log'), pYgMd = log(pYgMd); end

% If none of the observed patterns can explain the sequence and the
% sequence is longer than the depth of the tree (i.e. the longest possible
% pattern that is considered), the likelihood of the sequence is null. When
% exported in log scale, we return the log of the smallest postitive
% floating point number.
if corout && strcmpi(scaleme, 'log') && isinf(pYgMd), pYgMd = log(realmin); end

%% Predictions
%  ===========

% If the expectation has to be returned 
if nargout > 1
    
    % Conditionaly on entirely observed patterns
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    % Get predictions given each observed patterns
    % p(A|y,Rio) E {0,1}
    pAgYRio = double(pXgYRio == 1);
    
    % Marginalize predictions over all observed patterns
    % p(A|y,Ro) = sum_(i=1:{Ro}) p(A|y,Rio) * p(Rio|y)
    pAgYRo = sum(pAgYRio .* pRiogY) / pRogY;
    
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
    pAgYMd = (pAgYRo*pRogY + pAgYRu*pRugY) / (pRogY + pRugY);
    
    % If none of the observed patterns can explain the sequence and the
    % sequence is longer than the depth of the tree (i.e. the longest possible
    % pattern that is considered), we return chance level for the expectancy of
    % the next observation.
    if corout && (pRogY + pRugY) == 0, pAgYMd = 1/2; end
end

%% Wrap things up
%  ==============

% If the posterior distribution over patterns has to be returned 
if nargout > 2
    
    % Combine both the posterior distributions over entirely and partially
    % observed patterns
    if     idealinteg || ~idealinteg && strcmpi(scaleme, 'lin'), pRigY = [pRiogY,     pRiugY ];
    elseif               ~idealinteg && strcmpi(scaleme, 'log'), pRigY = [pRiogY, exp(pRiugY)];
    end
    
    % Normalize the posterior distribution
    pRigY = pRigY ./ sum(pRigY);
    
    % If the sequence is longer than the longest allowed pattern (i.e. the
    % depth of the tree) and that no patterns can reproduce the sequence,
    % (i.e. we cannot normalize the posterior because it would mean divide
    % it by zero), return an array of zero 
    if all(isnan(pRigY)), pRigY = zeros(1, nlRo); end
end

% If the entropy of the posterior distribution has to be returned
if nargout > 3
    HpRigY = Emergence_IO_Entropy(pRigY);
end

% Fill the output variable (for the partially observed patterns) such that
% its size always equals the depth of the tree (i.e. nu)
if nargout > 4, R = [R, repmat({''}, 1, nlRu)]; end

end
