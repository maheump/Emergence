function [ pY, pAgY, mpTgY, H_pTgY ] = Emergence_IO_Markov( y, scaleme, usegrid, prior, decw, dt )
% EMERGENCE_IO_MARKOV implements an observer learning parameters (i.e.
% frequency of transitions) of a first-order binary Markov chain from a
% sequence of binary observations.
%   - "y": a 1xN array specifying the sequence of binary observations (1s
%       and 2s).
%   - "scaleme": a string ('lin' or 'log') specifying whether the model
%       evidence sould be computed on a linear or logarithmic scale.
%   - "usegrid": a boolean specifing whether use grid-based or analytical
%       solutions.
%   - "prior": can be a string ('Bayes-Laplace' or 'Jeffreys') or a 2x2
%       matrix specifying the prior knowledge regarding the frequency of
%       transitions (expressed in pseudo-counts).
%   - "decw": a Nx1 or 1xN array of (decaying) weights that will weight
%       past observations of the sequence.
%   - "dt": a scalar specifying the grid precision for the posterior
%       distributions over theta.
% 
% Copyright (c) 2018 Maxime Maheu

%% INITIALIZATION
%  ==============

% Get the number of observations in the sequence
K = numel(y);

% By default, export the model evidence in log-scale
if nargin < 2 || isempty(scaleme), scaleme = 'log'; end
if     strcmpi(scaleme, 'lin'), islin = true;  islog = false;
elseif strcmpi(scaleme, 'log'), islin = false; islog = true;
end

% By default, use a uniform, non-informative, prior distribution
if nargin < 4 || isempty(prior), prior = 'Bayes-Laplace'; end

% By default, use a perfect integration
if nargin < 5 || isempty(decw), decw = ones(K,1); % no decay
else, decw = decw(:); % make sure it is a column vector (for linear algebra below)
end
if numel(decw) > K, decw = decw(end-K+1:end); % take only the most recent weights
elseif numel(decw) == K % nothing to do in that case
elseif numel(decw)  < K, error('Not enough weights for the length of the sequence');
end
decw = cat(1,decw,NaN);
% N.B adding a NaN at the end allow the use of linear algebra when the
% sequence is composed of only one single observation

% By default, use analytical solutions to speed up the computations
if nargin < 3, usegrid = false; end

% Analytical solutions are not available for analytical solutions, use
% grid-based approximations instead
if any(decw(1:end-1) ~= 1), usegrid = true; end

% By default, use a small grid precision to speed-up the computation
if nargin < 6, dt = 0.1; end

% This requires to create a grid for theta
if ~isempty(dt) && nargout > 2, returnpost = true;
else, returnpost = false; mpTgY = []; % return empty vector
end

% Create a vector for theta values, if:
%   - we use non-analytical, grid-based, solutions, or...
%   - we use analytical solutions but are asked to return the posterior,
%   because in that case, the posterior is never explicitly computed. We
%   thus have to create it afterwards, based on analytical solutions, using
%   probability distribution functions.
if usegrid || ~usegrid && returnpost
    theta = (0+dt/2):dt:(1-dt/2); % grid for theta based on grid precision
    nt = 1/dt; % number of values for one dimension of theta
end

%% PRIOR PROBABILITIES
%  ===================

% In case of a perfect integration, we can resort on analytical solutions
% that require to define prior probabilities in terms of transitions'
% pseudo-counts (and not probabilities per se) such that a conjugate prior
% can be used.
if ~usegrid
    
    % Uniform and non-informative Bayes-Laplace prior
    if isa(prior, 'char') && strcmpi(prior, 'Bayes-Laplace')
        p_nAgB = 1; % pN(A|B)
        p_nBgA = 1; % pN(B|A)
        p_nAgA = 1; % pN(A|A)
        p_nBgB = 1; % pN(B|B)
        
    % Non-informative Jeffreys prior
    elseif isa(prior, 'char') && strcmpi(prior, 'Jeffreys')
        p_nAgB = 1/2; % pN(A|B)
        p_nBgA = 1/2; % pN(B|A)
        p_nAgA = 1/2; % pN(A|A)
        p_nBgB = 1/2; % pN(B|B)
        
    % Custom prior
    elseif isa(prior, 'double') && numel(prior) == 4
        p_nAgB = prior(2,1); % pN(A|B)
        p_nBgA = prior(1,2); % pN(B|A)
        p_nAgA = prior(1,1); % pN(A|A)
        p_nBgB = prior(2,2); % pN(B|B)
        
    % Catch possible mistakes
    else, error('The prior input has the wrong form.');
    end
    
% Prior beliefs over theta values    
elseif usegrid
    
    % Uniform and non-informative Bayes-Laplace prior
    if isa(prior, 'char') && strcmpi(prior, 'Bayes-Laplace')
        pT = ones(nt, nt); % same value over all the grid
        
    % Non-informative Jeffreys prior
    elseif isa(prior, 'char') && strcmpi(prior, 'Jeffreys')
        pTgAgY = Emergence_IO_BetaPDF(theta, 1/2, 1/2, nt); % p(X|A), x-axis
        pTgBgY = Emergence_IO_BetaPDF(theta, 1/2, 1/2, nt); % p(X|B), y-axis
        pT = pTgAgY'*pTgBgY;
        
    % Custom prior
    elseif isa(prior, 'double') && numel(prior) == 4
        pTgAgY = Emergence_IO_BetaPDF(theta, prior(1,2), prior(1,1), nt); % p(X|A), x-axis
        pTgBgY = Emergence_IO_BetaPDF(theta, prior(2,1), prior(2,2), nt); % p(X|B), y-axis
        pT = pTgAgY'*pTgBgY;
        
    % Catch possible mistakes
    else, error('The prior input has the wrong form.');
    end
    
    % Return log prior distribution when the model evidence must be
    % returned in log-scale
    if islog, pT = log(pT); end
end

%% SEQUENCE'S MARGINAL LIKELIHOOD
%  ==============================

% Initialization
% ~~~~~~~~~~~~~~

% The likelihood of the first event is simply 1 over the number of
% different possible stimuli in the sequence
% p(y_1) = 1/2
if     islin, pY1 = 1/2;
elseif islog, pY1 = -log(2);
end

% Get the identity of the previous observations of each observation in the
% sequence. This eases the frequency counts of transitions
cond = cat(2, NaN, y(1:end-1));

% Get events' positions
A = (y == 1); % observations that are A
B = (y == 2); % observations that are B
gA = cond(A); % observations that follow an A
gB = cond(B); % observations that follow a B

% Get positions of each transition in the sequence
% N.B we add a false at the end to allow linear algebra below
% (it has no effect on the counts)
AgA = false(1, K+1); AgB = AgA; BgA = AgA; BgB = AgA;
AgA(A) = (gA == 1); % A => A cases
BgA(B) = (gB == 1); % A => B cases
AgB(A) = (gA == 2); % A => B cases
BgB(B) = (gB == 2); % B => B cases

% Get observational counts
nAgA = sum(AgA);
nAgB = sum(AgB);
nBgA = sum(BgA);
nBgB = sum(BgB);

% Compute the model evidence in case of a perfect integration
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% N.B. We can resort on analytical solutions
if ~usegrid
    
    % Combined observational and prior (pseudo-) counts
    nXgA = [nAgA + p_nAgA, nBgA + p_nBgA]; % N(A|A) + pN(A|A) & N(B|A) + pN(B|A)
    nXgB = [nAgB + p_nAgB, nBgB + p_nBgB]; % N(A|B) + pN(A|B) & N(B|B) + pN(B|B)
    
    % The integral of a beta distribution whose parameters are transitions'
    % counts. This can be solved analytically by the means of gamma
    % distributions:
    % p(y|Msp) = p(y_1) int[p(y_2:K|t(A|B),t(B|A))] dt
    if islin
        pYgTgA = prod(gamma(nXgA)) ./ gamma(sum(nXgA));
        pYgTgB = prod(gamma(nXgB)) ./ gamma(sum(nXgB));
    elseif islog
        pYgTgA = sum(gammaln(nXgA)) - gammaln(sum(nXgA));
        pYgTgB = sum(gammaln(nXgB)) - gammaln(sum(nXgB));
    end
    
    % The likelihood is the product of the two integrals and the
    % likelihood of the first event:
    % p(y|t(A|B),t(B|A)) = p(y_1) * p(y_2:K|t(X|A)) * p(y_2:K|t(X|B))
    % <=> log(p(y|t(A|B),t(B|A))) = log(p(y_1)) + log(p(y_2:K|t(X|A))) + log(p(y_2:K|t(X|B)))
    if     islin, pY = pY1 * pYgTgA * pYgTgB;
    elseif islog, pY = pY1 + pYgTgA + pYgTgB;
    end
    
    % If asked, return the posterior distribution
    if returnpost
        pTgAgY = Emergence_IO_BetaPDF(theta, nXgA(2), nXgA(1), nt); % p(X|A), x-axis
        pTgBgY = Emergence_IO_BetaPDF(theta, nXgB(1), nXgB(2), nt); % p(X|B), y-axis
        pTgAgY = pTgAgY ./ sum(pTgAgY); % normalize the marginal distribution
        pTgBgY = pTgBgY ./ sum(pTgBgY); % normalize the marginal distribution
        mpTgY = [pTgAgY', pTgBgY']; % concatenate marginal distributions
    end
    
    % Compute the expected value of each transition
    pXgA = nXgA ./ sum(nXgA); % p(A|A) & p(B|A)
    pXgB = nXgB ./ sum(nXgB); % p(A|B) & p(B|B)

% Compute the model evidence in case of a leaky integration
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% N.B. We need to rely on grid-based approximation solutions
elseif usegrid
    
    % For each transition, compute the likelihood of each observation,
    % given the weight of its corresponding position in the sequence, and
    % for each value of theta that is considered (depending on the grid)
    % e.g. p(observing A|A at position #10 with the corresponding weight
    % of position #10 if estimated p(A|A) = 1/10)
    dAgA = decw(AgA) * (1-theta) + (1 - decw(AgA)) *    theta;
    dBgA = decw(BgA) *    theta  + (1 - decw(BgA)) * (1-theta);
    dAgB = decw(AgB) *    theta  + (1 - decw(AgB)) * (1-theta);
    dBgB = decw(BgB) * (1-theta) + (1 - decw(BgB)) *    theta ;
	
    % Compute a distribution that is proportional to the joint posterior
    % distribution over theta 
    if islin
        pYgtBgA = prod(dAgA, 1) .* prod(dBgA, 1); % marginal distribution
        pYgtAgB = prod(dAgB, 1) .* prod(dBgB, 1); % marginal distribution
        pYgT = pY1 .* pYgtBgA' * pYgtAgB;         % joint distribution
        ppTgY = pYgT .* pT;                       % numerator in Bayes' rule
    elseif islog
        pYgtBgA = sum(log(dAgA), 1) + sum(log(dBgA), 1); % marginal distribution
        pYgtAgB = sum(log(dAgB), 1) + sum(log(dBgB), 1); % marginal distribution
        pYgT = pY1 + pYgtBgA' + pYgtAgB;                 % joint distribution
        ppTgY = pYgT + pT;                               % numerator in Bayes' rule
    end
    
    % Compute the sequence's marginal likelihood (denominator in Bayes'
    % rule) by integrating over theta (the joint distribution)
    if     islin, pY =     sum(    ppTgY(:))   * 1/nt^2;
    elseif islog, pY = log(sum(exp(ppTgY(:)))) - 2*log(nt);
    end
    % N.B. We use sum(X) / N(X) instead of the MATLAB "mean" function
    % because it is much faster (it avoids checks that are useless in the
    % context of this function).
    
    % Compute the posterior distribution over theta using Bayes' rule
    if     islin, pTgY =    (ppTgY ./ pY .* 1/nt^2);
    elseif islog, pTgY = exp(ppTgY  - pY  - 2*log(nt));
    end
    
    % Return (independent) marginal distributions
    pTgAgY = sum(pTgY, 2)'; % marginal distribution (vector over theta values)
    pTgBgY = sum(pTgY, 1);  % marginal distribution (vector over theta values)
    if returnpost, mpTgY = [pTgAgY; pTgBgY]'; end % concatenate marginal distributions
    
	% Compute the expected value of each transition
    pAgB = (pTgBgY / sum(pTgBgY)) * theta'; % expected value of theta(X|B)
    pBgA = (pTgAgY / sum(pTgAgY)) * theta'; % expected value of theta(X|A)
    pXgA = [1 - pBgA, pBgA]; % p(A|A) & p(B|A)
    pXgB = [pAgB, 1 - pAgB]; % p(A|B) & p(B|B)
end

%% PREDICTION
%  ==========

% Compute the likelihood that the next observation will be a A (an
% analytical formula can be used, which is simply a ratio).
% N.B. This depends on the identity of the previously received observation.
if nargout > 1
    pAgXY = [pXgA(1), pXgB(1)]; % p(A|A) & p(A|B)
    X     = y(end); % get the identity of the last observation
    pAgY  = pAgXY(X); % returns p(A) conditionaly on the last observation
end

%% ENTROPY OF THE POSTERIOR DISTRIBUTION
%  =====================================

% If the entropy of the posterior distribution has to be returned
if nargout > 3

    % If analytical solutions have been used
    if ~usegrid
        
        % Use analytical formulation of the entropy of a beta distribution
        % => differential entropy
        anent = @(a,b) betaln(a,b) - (a-1) * psi(a) - (b-1) * psi(b) + (a+b-2) * psi(a+b);
        H_pTBgA = anent(nXgA(2), nXgA(1));
        H_pTAgB = anent(nXgB(1), nXgB(2));
        
    % If grid-based numerical solutions have been used
    elseif usegrid
        
        % Compute the entropy of the (histogram-like) distribution
        % => discrete entropy
        H_pTBgA = Emergence_IO_Entropy(pTgAgY);
        H_pTAgB = Emergence_IO_Entropy(pTgBgY);
    end
    
    % When marginal distributions are independant, we have:
    % H(X,Y) = H(X)+H(Y)
    H_pTgY = H_pTBgA + H_pTAgB;
end
    
end
