function [ pY, pAgY, pTgY, H_pTgY ] = Emergence_IO_Bernoulli( y, scaleme, usegrid, prior, decw, dt )
% EMERGENCE_IO_BERNOULLI implements an observer learning the frequency of
% items from a binary sequence.
%   - "y": a 1xN array specifying the sequence of binary observations (1s
%       and 2s).
%   - "scaleme": a string ('lin' or 'log') specifying whether the model
%       evidence sould be computed on a linear or logarithmic scale.
%   - "usegrid": boolean value specifing whether to use grid-based or 
%       analytical solutions.
%   - "prior": can be a string ('Bayes-Laplace' or 'Jeffreys') or a 1x2
%       array specifying the prior knowledge regarding the frequency of
%       transitions (expressed in pseudo-counts).
%   - "decw": a  Nx1 or 1xN array of (decaying) weights that will weight
%       past observations of the sequence).
%   - "dt": a scalar specifying the grid precision for the posterior
%       distributions over theta.
% 
% Copyright (c) 2020 Maxime Maheu

%% INITIALIZATION
%  ==============

% Number of observations in the sequence
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

% By default, use a coarse grid precision to speed-up the computation
if nargin < 6, dt = 0.1; end

% This requires to create a grid for theta
if ~isempty(dt) && nargout > 2, returnpost = true;
else, returnpost = false; pTgY = []; % return empty vector
end

% Create a vector for theta values, if:
%   - we use non-analytical, grid-based, solutions, or...
%   - we use analytical solutions but are asked to return the posterior,
%   because in that case, the posterior is never explicitly computed. We
%   thus have to create it afterwards, based on analytical solutions, using
%   probability distribution functions.
if usegrid || ~usegrid && returnpost
    theta = (0+dt/2):dt:(1-dt/2); % grid for theta based on grid precision
    nt = 1/dt; % number of values for theta
end

%% PRIOR PROBABILITIES
%  ===================

% In case of a perfect integration, we can resort to analytical solutions
% that require to define prior probabilities in terms of transitions'
% pseudo-counts (and not probabilities per se) such that a conjugate prior
% can be used.
if ~usegrid
    
    % Uniform and non-informative Bayes-Laplace prior
    if isa(prior, 'char') && strcmpi(prior, 'Bayes-Laplace')
        p_nA = 1; % pN(A)
        p_nB = 1; % pN(B)
    
	% Non-informative Jeffreys prior
    elseif isa(prior, 'char') && strcmpi(prior, 'Jeffreys')
        p_nA = 1/2; % pN(A)
        p_nB = 1/2; % pN(B)
    
	% Custom prior
    elseif isa(prior, 'double') && numel(prior) == 2
        p_nA = prior(1); % pN(A)
        p_nB = prior(2); % pN(B)
    
    % Catch possible mistakes
    else, error('The prior input has the wrong form.');
    end
    
% Prior beliefs over theta values
elseif usegrid
    
    % Uniform and non-informative Bayes-Laplace prior
    if isa(prior, 'char') && strcmpi(prior, 'Bayes-Laplace')
        pT = ones(1, nt); % same value over all the grid

    % Non-informative Jeffreys prior
    elseif isa(prior, 'char') && strcmpi(prior, 'Jeffreys')
        pT = Emergence_IO_BetaPDF(theta, 1/2, 1/2, nt);

    % Custom prior
    elseif isa(prior, 'double') && numel(prior) == 2
        pT = Emergence_IO_BetaPDF(theta, prior(1), prior(2), nt);

    % Catch possible mistakes
    else, error('The prior input has the wrong form.');
    end
    
    % Return log prior distribution when the model evidence must be
    % returned in log-scale
    if islog, pT = log(pT); end
end

%% SEQUENCE MARGINAL LIKELIHOOD
%  ============================

% Initialization
% ~~~~~~~~~~~~~~

% Get events' positions
A = [(y == 1), false]; % we add a false at the end to allow linear algebra when there
B = [(y == 2), false]; % is only one observation (it has no effect on the counts)

% Compute the model evidence in case of a perfect integration
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% N.B. We can resort on analytical solutions
if ~usegrid
    
    % Get observational counts
    nA = sum(A); % oN(A)
    nB = sum(B); % oN(B)
    
    % Combined observational and prior (pseudo-) counts
    nX = [nA + p_nA, nB + p_nB]; % N(A) & N(B)
    
    % Since we use a conjugate prior, the model evidence is a beta
    % distribution, and its integral can be analytically computed.
    if     islin, pY =   beta(nX(1), nX(2));
    elseif islog, pY = betaln(nX(1), nX(2));
    end
    
    % If asked, return the posterior distribution
    if returnpost
        pTgY = Emergence_IO_BetaPDF(theta, nX(1), nX(2), nt); % beta distribution
        pTgY = pTgY ./ sum(pTgY(:)); % normalize the posterior
    end
    
% Compute the model evidence in case of a leaky integration
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% N.B. We need to rely on grid-based approximation solutions
elseif usegrid
    
    % For each item, compute the likelihood of each observation, given the
    % weight of its corresponding position in the sequence, and for each
    % value of theta that is considered (depending on the grid)
    % e.g. p(observing A at position #10 with the corresponding weight
    % of position #10 if estimated p(A) = 1/10)
    dA = decw(A) *    theta  + (1 - decw(A)) * (1-theta);
    dB = decw(B) * (1-theta) + (1 - decw(B)) *    theta ;
    % Note that this is a different, but mathematicaly equivalent, way of
    % implementing the leaky integration we propose in Meyniel, Maheu &
    % Dehaene, PCB (2016). In that case, we formalize it using substitution
    % probability instead of an exponential leak. The free parameter is
    % thus the probability of making a memory substitution (at each
    % observation) and not the slope of the exponential leak.
    
    % Compute the sequence's likelihood
    if     islin, pYgT = prod(   dA,  1) .* prod(    dB,  1);
    elseif islog, pYgT = sum(log(dA), 1)  + sum( log(dB), 1);
    end
    
    % Compute a distribution that is proportional to the posterior
    % distribution over theta
    if     islin, ppTgY = pYgT .* pT; % numerator in Bayes' rule
    elseif islog, ppTgY = pYgT  + pT; % numerator in Bayes' rule
    end
    
    % Compute the sequence's marginal likelihood (denominator in Bayes'
    % rule) by integrating over theta
    if     islin, pY =     sum(    ppTgY)   * 1/nt;
    elseif islog, pY = log(sum(exp(ppTgY))) - log(nt);
    end
    % N.B. We use sum(X) / N(X) instead of the MATLAB "mean" function
    % because it is much faster (it avoids checks that are useless in the
    % context of this function).
    
    % Compute the posterior distribution over theta using Bayes' rule
    if     islin, pTgY =    (ppTgY ./ pY .* 1/nt);
    elseif islog, pTgY = exp(ppTgY  - pY - log(nt));
    end
end

%% PREDICTION
%  ==========

% Compute the likelihood that the next observation will be a A
if nargout > 1
    if    ~usegrid, pAgY = nX(1) / sum(nX); % analytical formula (a ratio)
    elseif usegrid, pAgY = pTgY * theta';   % based on the grid
    end
end

%% ENTROPY OF THE POSTERIOR DISTRIBUTION
%  =====================================

% Return posterior distribution as a column vector
if nargout > 2, pTgY = pTgY(:); end

% If the entropy of the posterior distribution has to be returned
if nargout > 3
    
    % If analytical solutions have been used
    if ~usegrid
        
        % Use analytical formulation of the entropy of a beta distribution
        % => differential entropy
        anent = @(a,b) betaln(a,b) - (a-1) * psi(a) - (b-1) * psi(b) + (a+b-2) * psi(a+b);
        H_pTgY = anent(nX(1), nX(2));
        
    % If grid-based numerical solutions have been used
    elseif usegrid
        
        % Compute the entropy of the (histogram-like) distribution
        % => discrete entropy
        H_pTgY = Emergence_IO_Entropy(pTgY);
    end
end

end
