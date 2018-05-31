function [ pYgTs, pAgTs, pTgY ] = Emergence_IO_Bernoulli( y, scaleme, usegrid, prior, decw, dt )
% EMERGENCE_IO_BERNOULLI implements an observer learning the frequency of
% items from a binary sequence.
%   - "y": a 1xN array specifying the sequence of binary observations (1s
%       and 2s).
%   - "scaleme": a string ('lin' or 'log') specifying whether the model
%       evidence sould be computed on a linear or logarithmic scale.
%   - "usegrid": a boolean specifing whether use grid-based or analytical
%       solutions.
%   - "prior": can ba a string ('Bayes-Laplace' or 'Jeffreys') or a 1x2
%       array specifying the prior knowledge regarding the frequency of
%       transitions (expressed in pseudo-counts).
%   - "decw": a  Nx1 or 1xN array of (decaying) weights that will weight
%       past observations of the sequence).
%   - "dt": a scalar specifying the grid precision for the posterior
%       distributions over theta.
% 
% Copyright (c) 2018 Maxime Maheu

%% Initialization
%  ==============

% Number of observations in the sequence
K = numel(y);

% By default, export the model evidence in log-scale
if nargin < 2 || isempty(scaleme), scaleme = 'log'; end

% By default, use a uniform, non-informative, prior distribution
if nargin < 4 || isempty(prior), prior = 'Bayes-Laplace'; end

% By default, use a perfect integration
if nargin < 5 || isempty(decw), decw = ones(K,1); % no decay
else, decw = decw(:); % make sure it is a column vector (for linear algebra below)
end
if numel(decw) > K, decw = decw(end-K+1:end); % take only the most recent weights
elseif numel(decw) == K % nothing to do in that case
else, error('Not enough weights for the length of the sequence');
end

% By default, use analytical solutions to speed up the computations
if nargin < 3, usegrid = false; end

% If the weights are uniform over the entire sequence, it means that there
% are no memory decay, therefore analytical solutions can be used.
% Otherwise, it means that there is a memory decay which requires to used
% non-analytical grid-based solutions.
if     any(decw ~= 1), usegrid = true;
elseif all(decw == 1), usegrid = false;
end

% By default, use a small grid precision to speed-up the computation
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
    theta = 0:dt:1; % grid for theta based on grid precision
    nt = (1/dt)+1;  % number of values for theta
end

%% Prior probabilities
%  ===================

% In case of a perfect integration, we can resort on analytical solutions
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
    elseif isa(prior, 'double') && numel(prior) == 4
        pT = Emergence_IO_BetaPDF(theta, prior(1), prior(2), nt);

    % Catch possible mistakes
    else, error('The prior input has the wrong form.');
    end

    % Make sure the prior distribution is normalized
    pT = pT ./ sum(pT(:));
end

%% Sequence's marginal likelihood
%  ==============================

% Initialization
% ~~~~~~~~~~~~~~

% Get events' positions
A = (y == 1);
B = (y == 2);

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
    if     strcmpi(scaleme, 'lin'), pYgTs =   beta(nX(1), nX(2));
    elseif strcmpi(scaleme, 'log'), pYgTs = betaln(nX(1), nX(2));
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
    
    % Compute the sequence's likelihood
    pYgT = prod(dA, 1) .* prod(dB, 1);

    % Compute the posterior distribution over theta
    BayesNum = pYgT .* pT; % likelihood times prior (numerator in Bayes' rule)
    pY = sum(BayesNum);    % marginal likelihood (denominator in Bayes' rule)
    pTgY = BayesNum ./ pY; % posterior distribution over theta using Bayes' rule

    % Derive (log-) model evidence
    if     strcmpi(scaleme, 'lin'), pYgTs = pY / nt;
    elseif strcmpi(scaleme, 'log'), pYgTs = log(pY) - log(nt);
    end
    % N.B. We use sum(X) / N(X) instead of the MATLAB "mean" function
    % because it is much faster (it avoids checks that are useless in the
    % context of this function).

    % Compute the likelihood that the next observation will be a A
    pAgTs = pTgY * theta';
end

%% Predictions
%  ===========

% Compute the likelihood that the next observation will be a A
if    ~usegrid, pAgTs = nX(1) / sum(nX); % analytical formula (a ratio)
elseif usegrid, pAgTs = pTgY * theta';   % based on the grid
end

end
