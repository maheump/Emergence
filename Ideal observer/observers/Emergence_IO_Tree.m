function [ pYgMd, pAgYMd, pRigY, H_pRigY, R ] = Emergence_IO_Tree( y, nu, scaleme, usegrid, prior, decw, corout )
% EMERGENCE_IO_TREE implements an observer learning repeating patterns from
% a sequence of binary observations.
%   - "y": a 1xN array specifying the sequence of binary observations (1s
%       and 2s).
%   - "nu": a scalar specifying the depth of the tree to explore.
%   - "scaleme": a string ('lin' or 'log') specifying whether the model
%       evidence sould be computed on a linear or logarithmic scale.
%   - "usegrid": a boolean specifing whether use grid-based or analytical
%       solutions.
%   - "prior": the name of the prior distribution over rules to use, either
%       'Uniform' or based on the 'Size-principle'.
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
% |         |               |               |              |         | rules {Ro} = 2^3
% |     ____X____       ___[Y]___       ____X____      ____Y____     |
% |     |       |       |       |       |       |      |       |     |
% |=====|=======|=======|=======|=======|=======|======|=======|=====|==> K = 3
% |                     |       |                                    |
% |                   _(X)_   _(Y)_                                  | Partialy observed
% |                   |   |   |   |                                  | rules: {Ru} = 2^3
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
% (e.g. when the tree has been fully explored and no rules have been found)
if nargin < 7 || isempty(corout), corout = true; end

% By default, use a perfect integration
if nargin < 6 || isempty(decw), decw = ones(1,K); % no decay
else, decw = decw(:); % make sure it is a column vector (for linear algebra below)
end
if numel(decw) > K, decw = decw(end-K+1:end); % take only the most recent weights
elseif numel(decw) == K % nothing to do in that case
else, error('Not enough weights for the length of the sequence');
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

% By default, the depth of the tree (i.e. the maximum rule length
% allowed) is the length of the sequence
if nargin < 2, nu = K; end
% N.B. Deep trees induce longer computation time

%% Get the dimensionality of the problem
%  =====================================

% The number of possible rules depends on the sequence's length and the
% depth of the tree to explore
nlRo = min([nu, K]); % number or levels in the observed   part of the tree
nlRu = nu - nlRo;    % number or levels in the unobserved part of the tree

% Get levels in both parts of the tree
ilRo = 1:nlRo; % levels' indices in the   observed part of the tree
ilRu = 1:nlRu; % levels' indiced in the unobserved part of the tree
% (!!! they are different from the true levels' indices which can be
% computed using (nlRo+1):nu).

% Get the total number of rules
if usegrid
    nRo = 2.^(ilRo - 1);      % number of observed rules
    nRu = 2.^ilRu;            % number of unobserved rules
    nR = sum(nRo) + sum(nRu); % total number of rules
elseif ~usegrid
    nRo = 2^nlRo - 1;         % number of observed rules
    nRu = 2 * (2^nlRu - 1);   % number of unobserved rules
    nR = nRo + nRu;           % total number of rules
end

%% Rules' likelihood
%  =================

% Prepare output variables for the rules
if nargout > 3, R = cell(1,nlRo); end % if rules have to be returned
pYgRio  = NaN(1,nlRo); % probability of the rule given the sequence
pXgYRio = NaN(1,nlRo); % expected identity of the next item given the rules

% For each rule's length (i.e. each level of the tree)
for i = ilRo
    
    % There is only one correct rule per level in the tree, get that rule
    % for the current level
    Ri = y(1:i); % rule Ri
    if nargout > 3, R{i} = Ri; end % return the rule
    
    % Check whether the observed sequence could be a fixe (integer) number
    % of repetition(s) of that rule. This allows to use the following
    % conditionals which speed up the computations.
    ncol = ceil(K/i);
    ismult = (mod(K/i, 1) == 0);
    
    % Get the observed sequence...
    % - If the length of the sequence is not a fixe number of the rule's
    % length, add zeros such that it becomes the case
    if ~ismult
        j = (ncol*i)-K;
        ymat = zeros(K+j,1);
        ymat(1:K) = y;
    % - If the length of the sequence is a fixe number of the repetition
    % of the rule, simply get the sequence
    elseif ismult
        ymat = y;
    end
    
    % Transform the array specifying the observed sequence into a matrix
    % whose number of rows correspond to the length of the current rule and
    % the number of columns to the number of possible repetitions of that
    % rule
    ymat = reshape(ymat, [i, ncol]);
    
    % Compare each observation of the sequence to each observation of a
    % sequence obtained from the repetition of the current rule
    seqcheck = ymat == Ri; % a 2D matrix
    % N.B. This is more efficient than "repmat"ting the rule but it works
    % only on the most recent versions of MATLAB
    
    % In case of a perfect integration
    if idealinteg
    
        % Check if the repetition of that current rule entirely matches the
        % observed sequence
        isruletrue = sum(seqcheck(:)) == K;
        pYgRio(i) = (1/2) .* isruletrue;
        % N.B. the probability of observing the sequence given a particular
        % rule. Because rules are defined based on Xs/Ys, the first item
        % can be an A or a B with probability 1/2
        
    % In case of non-uniform integration, e.g. a memory leak implemented as
    % observations' weights indexed on an exponentialy decaying function.
    elseif ~idealinteg
        
        % Get the probability (based on the memory leak) with which each
        % observation of the sequence matches the one predicted from the
        % repetition of that rule
        pYkgRi = NaN(1,K);
        x = seqcheck(1:K);
        ok = (x == 1); %    match between rule and sequence
        er = (x == 0); % mismatch between rule and sequence
        pYkgRi(ok) =     decw(ok); % match
        pYkgRi(er) = 1 - decw(er); % mismatch
        
        % Compute how much the repetition of that current rule matches the
        % observed sequence
        % p(Ri|y) = (1/2) * prod_(y_k == yhat_k) w_k * prod_(y_k ~= yhat_k) (1-w_k)
        if     strcmpi(scaleme, 'lin'), pYgRio(i) =  (1/2) * prod(pYkgRi);
        elseif strcmpi(scaleme, 'log'), pYgRio(i) = -log(2) + sum(pYkgRi);
        end
   end
    
    % Derive the expected identity of the next observation based on that
    % rule (be it correct or not)
    if    ~ismult, pXgYRio(i) = Ri(i-j+1);
    elseif ismult, pXgYRio(i) = Ri(1);
    end
end

%% Prior probabilities
%  ===================

% For uniform (Bayes-Laplace) prior over rules
if strcmpi(prior, 'Bayes-Laplace')
    
    % The prior probability is simply 1 over the total number of rules
    % (i.e. both fully and partially observed ones)
    % p(Ri) = 1/{R}
    if     idealinteg || ~idealinteg && strcmpi(scaleme, 'lin'), p_pRi = ones(1, nlRo) ./ nR;
    elseif               ~idealinteg && strcmpi(scaleme, 'log'), p_pRi = repmat(-log(nR), 1, nlRo);
    end
    
% For prior based on the size-principle
elseif strcmpi(prior, 'Size-principle')
    
    % The size principle states that, since these are binary rules,
    % their prior simply depends on their lengths:
    % p(pRi) = (1/3) .^ |Ri|;
    if     idealinteg || ~idealinteg && strcmpi(scaleme, 'lin'), p_pRi = (1/3) .^ ilRo;
    elseif               ~idealinteg && strcmpi(scaleme, 'log'), p_pRi = -ilRo .* log(3);
    end
end

%% Sequence's marginal likelihood
%  ==============================

% Entirely observed rules
% ~~~~~~~~~~~~~~~~~~~~~~~

% Sequence's likelihood for each rule
% p(y|Rio) propto p(Rio|y) * p(Ri)
if      idealinteg || ~idealinteg && strcmpi(scaleme, 'lin'), pRiogY =     pYgRio .* p_pRi;
elseif                ~idealinteg && strcmpi(scaleme, 'log'), pRiogY = exp(pYgRio  + p_pRi);
end

% Sequence's likelihood for entirely observed rules
% p(Ro|y) = sum(i=1:{Ro}) p(Ri|y)
pRogY = sum(pRiogY); 

% Partially observed rules
% ~~~~~~~~~~~~~~~~~~~~~~~~

% Sequence's likelihood for partially observed rules
% p(Ru|y) = sum(i=1:{Ru}) {Rui} * p(y|Rui) * p(Rui)
if      usegrid, pRugY = sum(2.^(1:nlRu) .* (1/2) .* ((1/3).^(K+(1:nlRu))));
elseif ~usegrid, pRugY = max([0, 3^-K - 3^-nu * 2^(nu-K)]);
end

% All rules together
% ~~~~~~~~~~~~~~~~~~

% Compute the probability of observing the sequence under a 
% deterministic model
% p(y|Md) = p(y|Ro) * p(Ro) + p(y|Ru) * p(Ru)
pYgMd = pRogY + pRugY;
if idealinteg && strcmpi(scaleme, 'log'), pYgMd = log(pYgMd); end

% If none of the observed rules can explain the sequence and the sequence
% is longer than the depth of the tree (i.e. the longest possible rule that
% is considered), the likelihood of the sequence is null. When exported in
% log scale, we return the log of the smallest postitive floating point
% number.
if corout && strcmpi(scaleme, 'log') && isinf(pYgMd), pYgMd = log(realmin); end

%% Predictions
%  ===========

% Conditionaly on entirely observed rules
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Get predictions given
% p(A|y,Rio) E {0,1}
% p(A|y,Ro) = sum_(i=1:{Ro}) p(A|y,Rio) * p(Rio|y)
pAgYRio = double(pXgYRio == 1);
pAgYRo = sum(pRiogY .* pAgYRio) / pRogY;

% Conditionaly on partialy observed rules
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% The likelihood that the next observation would be an A conditioned on the
% yet unonbserved rules is simply chance level
% p(A|Riu) = 1/2
pAgYRu = 1/2;

% Combine both expectancies
% ~~~~~~~~~~~~~~~~~~~~~~~~~

% Compute the likelihood that the next observation will be a A by means of
% Bayesian Model Averaging
% p(A|y) = (p(A|y,Ro)*p(Ro|y) + p(A|y,Ru)*p(Ru|y)) / (p(Ro|y) + p(Ru|y))
pAgYMd = (pAgYRo*pRogY + pAgYRu*pRugY) / (pRogY + pRugY);

% If none of the observed rules can explain the sequence and the sequence
% is longer than the depth of the tree (i.e. the longest possible rule that
% is considered), we return chance level for the expectancy of the next
% observation.
if corout && (pRogY + pRugY) == 0, pAgYMd = 1/2; end

%% Wrap things up
%  ==============

% Fill output variables (for the unobserved rules) such that their
% size always equals the depth of the tree (i.e. nu)
if nargout > 2, pRigY = [pRiogY / pRogY, zeros(1,nlRu)]; end % NaN
if nargout > 4, R = [R, repmat({''}, 1, nlRu)]; end

% Measure the entropy of the posterior distribution
if nargout > 3, H_pRigY = sum(-pRigY .* log2(pRigY), 'OmitNaN'); end

end
