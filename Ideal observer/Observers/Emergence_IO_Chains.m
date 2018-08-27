function [ pYgMp, pAgYMp, mpTgY, H_pTgY, pXgX ] = Emergence_IO_Chains( y, scaleme, prior, i, dt )
% EMERGENCE_IO_CHAINS implements a N-order Markov chain that is estimated
% from a binary sequence.
%   - "y": a 1xN array specifying the sequence of binary observations (1s
%       and 2s).
%   - "scaleme": a string ('lin' or 'log') specifying whether the model
%       evidence sould be computed on a linear or logarithmic scale.
%   - "prior": can ba a string ('Bayes-Laplace' or 'Jeffreys') or a 2x(2^N)
%       matrix specifying the prior knowledge regarding the frequency of
%       transitions (expressed in pseudo-counts).
%   - "i": a scalar specifying the order of the Markov chain to estimate.
% 
% Copyright (c) 2018 Maxime Maheu

%% Initialization
%  ==============

% Get the number of observations in the sequence
K = numel(y);

% By default, export the model evidence in log-scale
if nargin < 2 || isempty(scaleme), scaleme = 'log'; end

% By default, use a uniform, non-informative, prior distribution
if nargin < 4 || isempty(prior), prior = 'Bayes-Laplace'; end

% By default, use a small grid precision to speed-up the computation
if nargin < 5, dt = 0.1; end
if ~isempty(dt) && ~isnan(dt), returnpost = true;
else, returnpost = false;
end

% Create a vector for theta values
if returnpost
    theta = (0+dt/2):dt:(1-dt/2); % grid for theta based on grid precision
    nt = 1/dt;  % number of values for theta
elseif ~returnpost, mpTgY = []; % return empty vector
end

%% Construct the chain
%  ===================

% Derive all the transitions in the chain
cond = ff2n(i) + 1;
ncond = 2^i;

% Count the number of As and Bs after each transition
o_nXgX = NaN(1,ncond);
for c = 1:ncond
    
    % Locate positions of the current transition
    condpos = strfind(y, cond(c,:));
    
    % Get the position of the observations coming just after
    following = condpos + i;
    following = following(following <= K);
    
    % Count the number of such observations separately for As and Bs
    o_nXgX(1,c) = sum(y(following) == 1); % for As
    o_nXgX(2,c) = sum(y(following) == 2); % for Bs
end

%% Prior probabilities
%  ===================

% Uniform and non-informative Bayes-Laplace prior
if isa(prior, 'char') && strcmpi(prior, 'Bayes-Laplace')
    p_nXgX = ones(2,ncond);
    
% Non-informative Jeffreys prior
elseif isa(prior, 'char') && strcmpi(prior, 'Jeffreys')
    p_nXgX = ones(2,ncond)./i;
    
% Custom prior
elseif isa(prior, 'double') && numel(prior) == ncond*2
    p_nXgX = prior;
    
% Catch possible mistakes
else, error('The prior input has the wrong form.');
end

%% Sequence's marginal likelihood
%  ==============================

% Get the likelihood of the first observations (chance) for which the chain
% cannot account for
% p(y_1:i) = (1/2)^i
% <=> log(p(y_1:i)) = -i * log(2)
x = min([K,i]); % number of beginning observations
if     strcmpi(scaleme, 'lin'), L1 = (1/2)^x;
elseif strcmpi(scaleme, 'log'), L1 = -x*log(2);
end

% Combined observational and prior (pseudo-) counts
nXgX = o_nXgX + p_nXgX;

% Get the likelihood of the chain with its current parameters (i.e. the
% counts). This equals to the integral of a beta distribution whose
% parameters are transitions' counts. This can be solved analytically by
% the means of gamma distributions.
if     strcmpi(scaleme, 'lin'), L2 = prod(gamma(nXgX), 1) ./ gamma(sum(nXgX, 1));
elseif strcmpi(scaleme, 'log'), L2 = sum(gammaln(nXgX), 1) - gammaln(sum(nXgX, 1));
end

% Combine both likelihoods together
% p(y|chain) = p(y_1:i) * prod_k_2^i (p(y_i+1:K|theta_k))
% <=> log(p(y|chain)) = log(p(y_1:i)) + sum_k_2^i (p(y_i+1:K|theta_k))
if     strcmpi(scaleme, 'lin'), pYgMp = L1 * prod(L2);
elseif strcmpi(scaleme, 'log'), pYgMp = L1 + sum(L2);
end

% Return marginal posterior distributions
if returnpost
    
    % Estimate marginal distributions
    % p(t_j|y) ~ Beta(N(A|xj), N(B|xj));
    mpTgY = NaN(nt,ncond);
    for j = 1:ncond
        mpTgY(:,j) = Emergence_IO_BetaPDF(theta, nXgX(1,j), nXgX(2,j), nt);
    end
    
    % Normalize the distribution
    mpTgY = mpTgY ./ sum(mpTgY);

% To speed-up the computations, avoid returning the distributions
else, mpTgY = [];
end

%% Predictions
%  ===========

% If the expectation has to be returned 
if nargout > 1
    
    % Compute the observed probabilities
    pXgX = nXgX ./ sum(nXgX,1);
    
    % Conditioned on what has been observed in the past
    if K > i
        
        % Get which history of observations (among those possible from the
        % chain) we are in
        prevobs = y(end-i+1:end);
        currtrans = (sum(cond == prevobs, 2) == i);
        
        % Get the corresponding probability of observing a A after this
        % history of observation, the trasition X => A
        pAgYMp = pXgX(1, currtrans);
    else
        
        % For the first observations, predictions are simply at chance level
        pAgYMp = 1/2;
    end
end

%% Entropy of the posterior distribution
%  =====================================
    
% If the entropy of the posterior distribution has to be returned
if nargout > 3

    % Use analytical formulation of the entropy of a beta distribution
    % => differential entropy
    anent = @(a,b) betaln(a,b) - (a-1) * psi(a) - (b-1) * psi(b) + (a+b-2) * psi(a+b);
    H_pTXgY = cellfun(@(x) anent(x(1),x(2)), mat2cell(nXgX, 2, ones(1,ncond)));
    
    % When marginal distributions are independant, we have:
    % H(X,Y) = H(X)+H(Y)
    H_pTgY = sum(H_pTXgY);
end

end

% % This is is a snipset of code that can be used to derive the joint
% % distributions from the marginal ones. It is commented because it
% % quickly becomes highly computationaly demanding
% joint = mat2cell(pTgY, ones(ncond, 1), nt);
% reshmat = mat2cell(eye(ncond)*(nt-1)+1, ones(ncond,1), ncond);
% joint = cellfun(@(x,y) reshape(x', y), joint, reshmat, 'UniformOutput', 0);
% joint = cellfun(@(x,y) repmat(x, -y+nt+1), joint, reshmat, 'UniformOutput', 0);
% joint = cell2mat(reshape(joint, [ones(1,ncond), ncond]));
% joint = prod(joint, N+2);
