function out = Emergence_IO_FullIO( ...
    s, ...          % (1)  the binary sequence with 1s and 2s (mandatory)
    pEd, ...        % (2)  the probability of an error in the deterministic observer
    pEp, ...        % (3)  the probability of an error in the probabilistic observer
    patlen, ...     % (4)  the depth of the tree used by the deterministic observer
    stat, ...       % (5)  the type of statistics learnt by the probabilistic obserber
    p_pRi, ...      % (6)  the prior for the patterns in the deterministic observer
    p_pTi, ...      % (7)  the prior for the statistics in the probabilistic observer
    p_pJk, ...      % (8)  the prior for the position of the change point
    computfor, ...  % (9)  whether to do iterative or final computations
    scaleme, ...    % (10) the scale for the model evidence
    postprec, ...   % (11) the precision of the posterior in the probabilistic observer
    verbose )       % (12) the type of command window output
% 
% Emergence_IO_FullIO is a function implementing the (Bayesian) ideal
% observer for the "Emergence" project.
% 
% In the related behavioral task, subjects are asked to report the
% generative (hidden) process of a binary sequence they are presented with.
% They are 3 possible generative processes: at a given point in time, the
% sequence can either be fully stochastic (i.e. p(A|B) = 1/2 & p(B|A) =
% 1/2, thus p(A) = 1/2), probabilistically biased (e.g. p(A|B) = 1/3, &
% p(B|A) = 1/3), or deterministically determined by the repetition of a
% pattern (e.g. [AAB]^n). The presented sequence is actually generated such
% that there is either (i) a change point in it separating a first part
% that is fully stochastically generated from a second one that can either
% be probabilistically or deterministically generated or (ii) no change
% point such that the sequence is simply a fully stochastic one. The
% possible sequences are thus:
%   - y = [Fully stochastic (1:N)]
%   - y = [Fully stochastic (1:Jk), Probabilistically generated (Jk+1:N)]
%   - y = [Fully stochastic (1:Jk), Deterministically generated (Jk+1:N)]
% where Jk is the position of the change point.
% 
% The present observer tries to infer the posterior probability of each of
% these three possible generative htpotheses:
%   - Mss: a fully stochastic hypothesis (where observations are generated
%       according to p(A) = 1).
%   - Mps: a probabilistic hypothesis (composed of a random part and a part
%       in which transition probabilities are biased).
%   - Msd: a deterministic hypothesis (composed of a random part and a part
%       in which a pattern is repeated).
% based solely on the received input (i.e. the sequence). The last two
% observers consider that there is a change point whose position is unknown
% separating the 2 parts. In these cases, the sequence's marginal likelhood
% (required to compure the hypotheses' posterior probability) is computed
% by combining the likelihood of the 2 parts of the sequence marginalized
% over all possible positions the change point can take.
% 
% The inputs of the function are the following ones:
%   - "s" (mandatory): a binary sequence.
%   - "pEd" (default: 0): the probability that the deterministic observer
%       will make a memory mistake at each received observation (this
%   - "pEp" (default: 0): the probability that the probabilistic observer
%       will make a memory mistake at each received observation (this
%       results in a leaky integration).
%       results in a leaky integration).
%   - "patlen" (default: length(s)): the depth of the patterns' tree to 
%       explore (i.e. the longest pattern length considered).
%   - "stat" (default: 'Transitions'): the type of statistics to be learnt
%       by the probabilistic observer. It can be either 'Items',
%       'Alternation', 'Transitions', 'Chain3', 'Chain10', ...
%   - "p_pRi" (default: 'Size-principle'): prior distribution over pattern
%       in the deterministic observer. It can be either 'Uniform' or based
%       a 'Size-principle' (i.e. smallest pattern are favoured).
%   - "p_pTi" (default: 'Bayes-Laplace'): prior distribution over cases in
%       the probabilistic observer. It can be either 'Bayes-Laplace',
%       'Jeffreys' or a custom prior in terms of pseudo-counts. In that
%       latter case, it must match the "stat" field: e.g. if stat = "Item",
%       then p_pTi should be a 1x2 array (e.g. [1 1] is a 'Bayes-Laplace'
%       prior).
%   - "p_pJk" (default: 'Uniform'): prior distribution over change point's
%       position. It can be either 'Uniform', [mu, sigma] for a Gaussian
%       prior or any custom prior distribution (a 1xN array where N is the
%       length of the sequence "s").
%   - "computfor" (default: 'all'): whether to return the iteratively
%       updated inference ('all') or simply its final status ('last').
%   - "scaleme" (default: 'log'): scale for the model evidences, either a
%       linear ('lin') or logarithmic ('log') scale.
%   - "postprec" (default: []): precision of the posterior distribution for
%       the probabilistic observer. Because it is computationaly demanding,
%       when it is empty, posterior distributions for the probabilistic
%       observer are not returned. Note also that posterior distributions
%       cannot be returned for Markov chain of order higher than 1.
%   - "verbose" (default: true): whether to print status of the inference
%       in the command window. It can be either false (only prints the
%       running time of the function), true (it prints chosen options) or
%       2 (it prints status of the inference after each observation).
% 
% The function returns many different quantities related to the inference
% of the input sequence, including (1) the posterior probability of each
% hypothesis, (2) the posterior over change point's position, (3)
% predictions, surprise and entropy levels and (4) update-related
% quantities.
% 
% Copyright (c) 2018 Maxime Maheu

%% Check the inputs
%  ================

% Start the timer
tic;

% By default, print informative messages
if nargin < 12, verbose = true; end

% Welcome message
eqnum = 60;
if verbose
    fprintf('%s\nThis is the Bayesian ideal observer for "Emergence"\n', ...
        repmat('=', 1, eqnum));
end

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Get information about the sequence %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

s = s(:)'; % make sure the sequence vector is a row vector
Ns = numel(unique(s)); % get the number of different observations in the sequence
if Ns > 2, error('The sequence must be binary'); end
N = numel(s); % get the total number of observations in the sequence
if verbose, fprintf('  * The sequence is composed of %i observations\n', N); end

% ~~~~~~~~~~~~~~~~~~~ %
% Recode the sequence %
% ~~~~~~~~~~~~~~~~~~~ %

if Ns == 2
    A = min(s);
    B = max(s);
    S = NaN(1,N);
    S(s == A) = 1; % As are denoted by 1s
    S(s == B) = 2; % Bs are denoted by 2s
elseif Ns == 1
    S = ones(1,N);
end

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Define the integration function %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% We define integration function as exponential/power leaks. Values of this
% function weight the received observations. One parameter controls the
% decay of these integration functions. It is simply the probability of
% making a "memory mistake" at each newly received observation.

% By default, use a perfect integration (i.e. no "memory mistake")
if nargin < 2 || isempty(pEd) || isnan(pEd), pEd = 0; end
if nargin < 3 || isempty(pEp) || isnan(pEp), pEp = 0; end

% Check that the error probabilities are between 0 and 1
if pEp < 0 || pEp > 1/2, error('The P-error probability is not defined between 0 and 1/2'); end
if pEd < 0 || pEd > 1/2, error('The D-error probability is not defined between 0 and 1/2'); end

% Get the vectors of decaying weights indexed on those probabilities
if pEp == 0, decayP = ones(1, N); % equivalent to >> Emergence_IO_Leak(pEp, N); but faster
else, decayP = Emergence_IO_Leak(pEp, N);
end
if pEd == 0, decayD = ones(1, N); % equivalent to >> Emergence_IO_Leak(pEd, N); but faster
elseif pEd > 0 && pEd == pEp, decayD = decayP; % faster than calling "Emergence_IO_Leak" again
elseif pEd > 0 && pEd ~= pEp, decayD = Emergence_IO_Leak(pEd, N);
end

% Display the chosen probability for each observer of making a "memory
% mistake"
if verbose
    fprintf(['  * Probability of making a memory P-mistake: ', ...
        '%1.0f%%\n'], ceil(pEp*100));
    fprintf(['  * Probability of making a memory D-mistake: ', ...
        '%1.0f%%\n'], ceil(pEd*100));
end

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Define the depth of the patterns' tree to explore %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

% By default, it is the length of the sequence
if nargin < 4, patlen = N; end

% If it is larger than the sequence, notify it
if patlen > N 
    warning(['The depth of patterns'' tree to explore is larger than the ', ...
        'number of observations in the sequence.']);
end
if verbose, fprintf('  * Longest possible pattern: %i observations\n', patlen); end

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Define the type of statistics to be learnt by the probabilistic model %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

% By default, it learns transition probabilities
if nargin < 5 || isempty(stat), stat = 'Transitions'; end

% Define a useful string comparison function
probamod = find(cellfun(@(x) contains(stat, x, 'IgnoreCase', true), ...
    {'Item', 'Alternation', 'Transition', 'Chain'}));

% Check that the specified statistic to learn is among the allowed ones
if isempty(probamod)
    error(['The type of statistic to be learnt by the probabilistic model is ', ...
        'not available']);
else
    if verbose, fprintf('  * Estimated statistic(s): %s\n', lower(stat)); end
end

% If a Markov chain has to be estimated, get its order
nOrd = str2double(stat(regexp(stat, '\d')));

% Define the number of (independant) theta distributions to estimate
nTdim = [1, 1, 2, 2^nOrd];
nTdim = nTdim(probamod);

% Return an error if a leaky integration is required
if probamod == 4 && pEp ~= 0
    error(['The higher-order Markov chains are not (yet)', ...
        'compatible with leaky integrations.']);
end

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Define the type of prior for the deterministic observer %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

% By default, use the size-principle for the prior over patterns
if nargin < 6 || isempty(p_pRi), p_pRi = 'Size-principle'; end

% The two available analytical prior distributions over patterns
if verbose, fprintf('  * Type of prior over patterns: %s\n', lower(p_pRi)); end
if ~any(strcmpi(p_pRi, {'Uniform', 'Size-principle'}))
    error(['Please provide a prior over patterns that is either "Uniform" ', ...
        'or based on the "Size-principle"']);
end

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Define the type of prior for the probabilistic observer %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

% By default, use a uniform non-informative prior
if nargin < 7 || isempty(p_pRi), p_pTi = 'Bayes-Laplace'; end

% If it is a custom prior distribution, check that its size matches the
% type of statistic that will be estimated
if ischar(p_pTi)
    if ~any(strcmpi(p_pRi, {'Uniform', 'Size-principle'}))
        error(['Please provide a prior over statistics that is either ', ...
            '"Bayes-Laplace" or based on the "Size-principle"']);
    end
    if verbose, fprintf('  * Type of prior over statistics: %s\n', p_pTi); end
elseif ~ischar(p_pTi)
    if any(size(p_pTi) ~= [2, nTdim])
        error('The custom prior has the wrong size');
    end
    if verbose, fprintf('  * Type of prior over statistics: custom\n'); end
end

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Define the type of prior belief regarding the change point's position %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

% By default, the prior probability regarding the position of the change
% point is uniform over the sequence
if     nargin <  8 || isempty(p_pJk) || ischar(p_pJk),     Jprior = 'Uniform';
elseif nargin >= 8 && ~ischar(p_pJk) && numel(p_pJk) == 2, Jprior = 'Gaussian';
elseif nargin >= 8 && ~ischar(p_pJk) && numel(p_pJk) == N, Jprior = 'Custom';
else, error(['Please provide a prior distribution over change point''s', ...
        'positions that is either ''Uniform'', [mu,sigma] or 1xN array']);
end
if verbose
    fprintf('  * Type of prior over change point''s position: %s\n', lower(Jprior));
end

% Create p_pJ, a row vector of length N, in which the value in each cell
% denotes the prior probability that a change point would appear after the
% corresponding observation.

% The prior over change point's positions can be a uniform uninformative
% Bayes-Laplace prior:
% for all k E [1,N]: p(Jk) = 1/N
if strcmpi(Jprior, 'Uniform')
    p_pJk = ones(1,N) ./ N;
    
% The prior over change point's positions can also be a non-uniform 
% informative prior that follows a Gaussian distribution:
% p(J) ~ Normal(1:N, mu, sigma)
elseif strcmpi(Jprior, 'Gaussian')
    mu    = p_pJk(1);                % mean
    sigma = p_pJk(2);                % variance
    gauss = normpdf(1:N, mu, sigma); % the gaussian distribution
    p_pJk  = gauss ./ sum(gauss);    % normalized

% If the prior distrubution over change point's position is a user-provided
% custom one, make sure it is a row vector to ensure that linear algebra
% below works
elseif strcmpi(Jprior, 'Custom')
    p_pJk = p_pJk(:)';
end

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Define whether to return iteratively updated inference or final one %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

% By default, compute the model's quantities updated following each
% observation received
if nargin < 9, computfor = 'all'; end

% If the computations have to be done after each observation
if strcmpi(computfor, 'all')
    obs = 2:N;
    if verbose, fprintf('  * Update type: iterative\n'); end

% If the computations have to be done only after the last observation
elseif strcmpi(computfor, 'last')
    obs = [N-1, N];
    if verbose, fprintf('  * Type of inference: final\n'); end
else, error('Please provide a computation method that is either "all" or "last"');
end

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Define the scale of the model evidence %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% The use of the log-scale is intended to improve the computation accuracy
% and to avoid numeric overflow.

% By default, use the logarithmic scale for the marginal likelihoods
if nargin < 10, scaleme = 'log'; end

% Check that the specified scale is either linear or logarithmic
if ~any(strcmpi(scaleme, {'lin', 'log'}))
    error('Please provide a scale that is either "lin" or "log"');
end
if verbose, fprintf('  * Scale of model evidences: %s-scale\n', scaleme); end

% Convert to boolean to speed up the computations
if     strcmpi(scaleme, 'lin'), islin = true;  islog = false;
elseif strcmpi(scaleme, 'log'), islin = false; islog = true;
end

% Turn the prior on change point positions to log scale if required
if islog, p_pJk = log(p_pJk); end

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Define whether to return posterior distribution over models' unknown parameters %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

% By default, do not return posterior distributions over theta (unknown
% parameter(s) of a given model) for the models because it is
% computationaly demanding
if nargin < 11, postprec = []; end
if numel(postprec) > 1 || any(postprec < 0 | postprec > 1)
    error('Check the precision of the probability grid ');
end
if ~isempty(postprec), fprintf(['  * Precision of the posterior over ', ...
        'models'' parameters = %1.0f%%\n'], ceil(postprec*100)); end
if ~isempty(postprec), nTprec = 1/postprec; end

%% Prepare output variables
%  ========================

% Output variables for the model evidence
% p(y|Msi)
pYgMsp = NaN(1,N);
pYgMsd = NaN(1,N);
pYgMss = NaN(1,N);

% Output variables for models' posterior distribution over theta (unknown
% model parameters)
% p(theta|y,Msi)
if ~isempty(postprec), pTgYMsp = NaN(nTprec, N, nTdim); end
pTgYMsd = NaN(patlen, N);

% Output variables for the posterior belief over change point's position
% p(Jk|y,Msi)
pJkgYMsp = NaN(N,N);
pJkgYMsd = NaN(N,N);

% Output variables for expectancies
% p(A|y,Msi)
pAgYMsp = NaN(1,N);
pAgYMsd = NaN(1,N);

%% Compute models evidence, beliefs regarding change point's position and predictions
%  ==================================================================================

% For detailed display of the computations' status
if strcmpi(computfor, 'all') && verbose == 2
    cols = ['|  K  | y(K) | p(y|Mss) | ', 'p(y|Msp) | p(y|Msd) |'];
    letters = {'A','B'};
    fprintf('%s\n', repmat('=', 1, eqnum));
end

% Following each observation
for K = obs

    %% Get ready
    %  =========
    
    % Display the status of the computations
    if strcmpi(computfor, 'all') && verbose == 1
        if K == obs(1), fprintf('Computing...     '); end
        fprintf('\b\b\b\b%3.0f%%', (K/N)*100);
    end
    
    % Display the status of the computations in a more detailed manner
    if strcmpi(computfor, 'all') && verbose == 2
        
        % Print values
        if K > obs(1)
            fprintf('\n| %3.0f |   %s  ', K, letters{S(K)});
            fprintf('|   %4.0f   |   %4.0f   |   %4.0f   |', ...
                pYgMss(K-1), pYgMsp(K-1), pYgMsd(K-1));
        end
        
        % Print labels
        if any(K == [2, 25:25:N-1])
            fprintf('\n%s\n%s\n%s', repmat('-', 1, numel(cols)), ...
                cols, repmat('-', 1, numel(cols)));
        end
    end
    
    % K is the number of observations received for now. It allows to get
    % the sequence from the first observation to the current one.
    Y = S(1:K);
    
    %% Compute model evidence of the fully-stochastic model
    %  ====================================================
    
    % The likelihood of a sequence considered as fully random is simply the
    % probability of observing any observation (here 1/2 since the sequence
    % is binary) powered by the number of observations in the sequence. It
    % is strictly equivalent to making coin tosses with an unbiased coin.
    % p(y|Mss) = (1/2)^K
    % <=> log(p(y|Mss)) = -K * log(2)
    if     islin, pYgMss(K) = Emergence_IO_Null(K, 'lin');
    elseif islog, pYgMss(K) = Emergence_IO_Null(K, 'log');
    end
    
    %% Loop over possible change point's positions
    %  ===========================================
    
    % Prepare outputs for change point's position:
    % p(y|Jk,Msp) & p(y|Jk,Msd)
    pYgJkMsp = NaN(1,K-1);
    pYgJkMsd = NaN(1,K-1);
    
    % Prepare outputs for predictions:
    % if Jk >= K, p(A|y,Jk,Msp) = (1/2)
    % if Jk <  K, see below...
    pAgJkYp2Mp = ones(1,N) ./ 2;
    pAgJkYp2Md = ones(1,N) ./ 2;
    
    % Prepare the output variable for the posterior distribution over
    % model's parameter(s):
    if ~isempty(postprec), pTgY = NaN(nTprec,K-1,nTdim); end % probabilistic observer
    pRgY = NaN(K-1,patlen); % deterministic observer
    
    % For each possible change point's position:
    % for all Jk E [1,K-1]
    for k = 1:K-1
        
        % Get the second part of the sequence:
        L   = K-k;                 % length of the second part
        Yp2 = Y(k+1:end);          % observations in the second part
        Wp  = decayP(end-L+1:end); % decaying weights for the second part
        Wd  = decayD(end-L+1:end); % decaying weights for the second part
        
        %% Compute the marginal likelihood of the 1st part under a fully-stochastic model
        %  ==============================================================================
        
        % The marginal likelihood of a fully stochastic model:
        % p(y1|Ms) = (1/2)^k
        % <=> log(p(y1|Ms)) = -k * log(2)
        if     islin, pYp1gMs = Emergence_IO_Null(k, 'lin');
        elseif islog, pYp1gMs = Emergence_IO_Null(k, 'log');
        end
        
        %% Compute the marginal likelihood of the 2nd part under a probabilistic model
        %  ===========================================================================
        
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
        % Compute the likelihood of the second part of the sequence under %
        % a probabilistic model that learns one of the following          %
        % statistic                                                       %
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
        
        % If the statistic to learn is the frequency of items:
        if probamod == 1
            
            % A Bernoulli process:
            % p(y2|Mp) ~ Beta(N(A),N(B))
            % p(y(K+1)=A|y2,Mp) = N(A) / (N(A) + N(B))
            [pYp2gMp, pAgJkYp2Mp(k), pTgY(:,k)] = Emergence_IO_Bernoulli(...
                Yp2, ...        % second part of the sequence
                scaleme, ...	% scale for p(y|Mp)
                false, ...      % use analytical solutions instead of numerical ones
                p_pTi, ...      % prior distribution over theta values
                Wp, ...         % vector of decaying weights
                postprec);      % precision of the grid for theta values
            
        % If the statistic to learn is the frequency of alternations:
        elseif probamod == 2
            
            % Recode the sequence with alternations/repetitions
            Yp2alt = abs(abs(diff(Yp2))-2);
            % N.B. Alternations are denoted by 1s and repetitions by 2s
            
            % A Bernoulli process:
            % p(y2|Mp) ~ Beta(N(alt.),N(rep.))
            % p(alt.|y2,Mp) = N(alt.) / (N(alt.) + N(rep.))
            [pYp2gMp, pAltgYp2JkMp, pTgY(:,k)] = Emergence_IO_Bernoulli(...
                Yp2alt, ...     % second part of the sequence
                scaleme, ...	% scale for p(y|Mp)
                false, ...      % use analytical solutions instead of numerical ones
                p_pTi, ...      % prior distribution over theta values
                Wp, ...         % vector of decaying weights
                postprec);      % precision of the grid for theta values
            
            % The prediction of the forthcoming observation depends on (i)
            % the identity of the last observation and (ii) the inferred
            % probability of alternation:
            % p(y(K+1)=A|y2,Mp) = 1 - p(alt.|y2,Mp) if y2(K) = A
            % p(y(K+1)=A|y2,Mp) =     p(alt.|y2,Mp) if y2(K) = B
            LastObs = Yp2(end);
            if     LastObs == 1, pAgJkYp2Mp(k) = 1 - pAltgYp2JkMp; % repetition  from A to A
            elseif LastObs == 2, pAgJkYp2Mp(k) =     pAltgYp2JkMp; % alternation from B to A
            end
            
        % If the statistics to learn are transition probabilities:
        elseif probamod == 3
            
            % A first-order Markov-chain:
            % p(y2|Mp) ~ Beta(N(A|A),N(B|A)) * Beta(N(A|B),N(B|B))
            %          = p(y2(1)) * p(y2(2:end)|t(A|B),t(B|A))
            % p(y(K+1)=A|y2,Mp) = N(A|X) / (N(A|X) + N(B|X))
            [pYp2gMp, pAgJkYp2Mp(k), pTgY(:,k,:)] = Emergence_IO_Markov(...
                Yp2, ...        % second part of the sequence
                scaleme, ...	% scale for p(y|Mp)
                false, ...      % use analytical solutions instead of numerical ones
                p_pTi, ...      % prior distribution over theta values
                Wp, ...         % vector of decaying weights
                postprec);      % precision of the grid for theta values
            
        % If the statistics to learn are higher-order transition probabilities:
        elseif probamod == 4
            
            % A n-order Markov-chain:
            % p(y2|Mp) ~ prod_i Beta(N(A|T_i), N(B|T_i))
            %          = p(y2(1)) * p(y2(2:end)|t(A|T_1),t(A|T_2),...,t(B|T_n))
            % p(y(K+1)=A|y2,Mp) = N(A|X) / (N(A|X) + N(B|X))
            [pYp2gMp, pAgJkYp2Mp(k), pTgY(:,k,:)] = Emergence_IO_Chain(...
                Yp2, ...        % second part of the sequence
                scaleme, ...	% scale for p(y|Mp)
                p_pTi, ...      % prior distribution over theta values
                nOrd, ...       % order of the Markov chain
                postprec);      % the function does not handle joint posterior distribution
        end
        
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
        % Compute the probability of observing the sequence as a function %
        % of a particular change point's position under a stochastic-to-  %
        % probabilistic model                                             %
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
        
        % For each change point's position, sum the 2 marginal likelihoods
        % p(y,Jk|Msp) = p(y1|Mr) * p(y2|Mp)
        % <=> log(p(y,Jk|Msp)) = log(p(y1|Mr)) + log(p(y2|Mp))
        if     islin, pYgJkMsp(k) = pYp1gMs * pYp2gMp;
        elseif islog, pYgJkMsp(k) = pYp1gMs + pYp2gMp;
        end
        
        %% Compute the marginal likelihood of the 2nd part under a deterministic model
        %  ===========================================================================
        
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
        % Compute the marginal likelihood for the second part of the  %
        % sequence under a deterministic model that supposes that the %
        % sequence is generated according to the repetition of a      %
        % particular pattern                                          %
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
        
        % Marginalize over all possible patterns
        % p(Md|y2) = sum_i p(y|Ri) * p(Ri)
        [pYp2gMd, pAgJkYp2Md(k), pRgY(k,:)] = Emergence_IO_Tree(...
            Yp2, ...     % second part of the sequence
            patlen, ...  % depth of the tree to explore
            scaleme, ... % scale for p(y|Md)
            false, ...   % use analytical solutions instead of numerical ones
            p_pRi, ...   % type of prior over patterns
            Wd, ...      % vector of decaying weights
            true);       % make sure that the output variables are not singular
        
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
        % Compute the probability of observing the sequence as a function %
        % of a particular change point's position under a stochastic-to-  %
        % deterministic model                                             %
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
        
        % For each change point's position, sum the 2 marginal likelihoods
        % p(y|Jk,Msd) = p(y1|Ms) * p(y2|Md)
        % <=> log(p(y|Jk,Msd)) = log(p(y1|Mr)) + log(p(y2|Md))
        if     islin, pYgJkMsd(k) = pYp1gMs * pYp2gMd;
        elseif islog, pYgJkMsd(k) = pYp1gMs + pYp2gMd;
        end
    end
    
    %% Combine model evidence of stochastic-to-regular models
    %  ======================================================
    %  Because there are 2 parts (one random and one non-random), we
    %  must combine beliefs over all positions the change point may have
    %  already taken:
    %  for all i E {P,D}, p(y|Msi) = U_{Jk=1:K} p(y|Jk,Msi)
    %  Because the different possible change point's positions are
    %  independent from one another, it is just required to sum beliefs
    %  over all observed change point's position weighted by the prior
    %  probability of observing a change point at that point in time:
    %  for all i E {P,D}: p(y|Msi) = sum_{Jk=1:K-1} p(y|Jk,Msi) * p(Jk)
    
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
    % Posterior distribution of already observed change point's positions %
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %    
    % For each position of the change point that have already been
    % observed, compute the posterior probability that the change point
    % would be at that position:
    % for all Jk E {1,2,...,K-1}: p(Jk|y,Msp) propto p(y|Jk,Msp) * p(Jk)
    % <=> log(p(Jk|y,Msp)) propto log(p(y|Jk,Msp)) + log(p(Jk))
    
    % Prior probability over observed and unobserved change point's
    % positions
    p_pobsJk   = p_pJk(1:K-1); %   observed positions
    p_punobsJk = p_pJk(K:N);   % unobserved positions
    
    % For the stochastic-to-probabilistic observer
    if     islin, ppobsJkgYMsp = pYgJkMsp .* p_pobsJk;
    elseif islog, ppobsJkgYMsp = pYgJkMsp  + p_pobsJk;
    end
    
    % For the stochastic-to-deterministic observer
    if     islin, ppobsJkgYMsd = pYgJkMsd .* p_pobsJk;
    elseif islog, ppobsJkgYMsd = pYgJkMsd  + p_pobsJk;
    end
    
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
    % Compute the model evidence of both models by marginalizing over %
    % change points' positions                                        %
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %    
    % To get the likelihood of the sequence under a particular model, we
    % just need to sum over all possible positions of the change point
    % p(y|Msp) = sum_{Jk=1:K-1} p(Jk|y,Msp)
    % <=> log(p(y|Msp)) = log(sum_{Jk=1:K-1} p(Jk|y,Msp))
    
    % For the stochastic-to-probabilistic observer
    if     islin, pYgMsp(K) =     sum(    ppobsJkgYMsp  );
    elseif islog, pYgMsp(K) = log(sum(exp(ppobsJkgYMsp)));
    end
    
    % For the stochastic-to-deterministic observer
    if     islin, pYgMsd(K) =     sum(    ppobsJkgYMsd  );
	elseif islog, pYgMsd(K) = log(sum(exp(ppobsJkgYMsd)));
    end
    
    %% Compute posterior beliefs regarding the change point's position
    %  ===============================================================
    
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
    % Posterior distribution of yet unobserved change point's positions %
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %    
    % The probability of the change point being at a yet unseen position
    % is the prior probability of that position times the likelihood of a
    % fully-stochastic sequence of the corresponding length. This has the
    % effect of increasing the probability of an alredy observed change
    % point has observations are received. Note that this is increase is 
    % striclty linear if the prior is uniform.
    % for all Jk E {K,K+1,...,N}: p(y|Jk) = (1/2)^K * p(Jk)
    % <=> log(p(y|Jk>K)) = -K * log(2) + p(Jk)
    if     islin, ppunobsJkgYKMsp = Emergence_IO_Null(K, 'lin') .* p_punobsJk;
    elseif islog, ppunobsJkgYKMsp = Emergence_IO_Null(K, 'log')  + p_punobsJk;
    end
    
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
    % Combine posterior beliefs over both observed and yet unobserved %
    % change point's positions to get the full posterior distribution %
    % over change point's positions                                   %
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
    
    % Combine posterior distributions for the stochastic-to-probabilistic observer
    if     islin, pJkgYKMsp =     cat(2, ppobsJkgYMsp, ppunobsJkgYKMsp);
    elseif islog, pJkgYKMsp = exp(cat(2, ppobsJkgYMsp, ppunobsJkgYKMsp));
    end
    
    % Combine posterior distributions for the stochastic-to-deterministic observer
    if     islin, pJkgYKMsd =     cat(2, ppobsJkgYMsd, ppunobsJkgYKMsp);
    elseif islog, pJkgYKMsd = exp(cat(2, ppobsJkgYMsd, ppunobsJkgYKMsp));
    end
    
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
    % Normalize the posterior distribution over change point's positions %
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
    % for all i E {P,D} and Jk E {1,2,...,N}:
    % p(Jk|y,Msi) = (p(y|Jk,Msi) * p(Jk)) / p(y|Msi)
    % where p(y|Msi) = sum_{Jk=1:N} p(y|Jk,Msi) * p(Jk)
    
    % Normalize the posterior distribution for the stochastic-to-
    % probabilistic observer
    pJkgYMsp(:,K) = pJkgYKMsp ./ sum(pJkgYKMsp);
    
    % Normalize the posterior distribution for the stochastic-to-
    % deterministic observer
    pJkgYMsd(:,K) = pJkgYKMsd ./ sum(pJkgYKMsd);
    
    %% Compute model-specific posterior distributions
    %  ==============================================
    %  Note that here the posterior distributions over models' parameters
    %  are marginalized over change point's position. Importantly, we
    %  consider only change point's positions that are already observed.
    %  Another possibility would be to marginalize over all (including not
    %  yet observed change point's positions) but it is less interesting
    %  in that case because the maginalized posteriors would be strongly
    %  driven by the (uniform) prior distribution over change point's
    %  positions thus resulting in uniform distribution over models'
    %  parameters as well.
    %  for all i E {P,D} and Jk E {1,2,...,K-1}:
    %  p(theta|y,Msi) = sum_{Jk=1:K-1} p(theta|y,Jk,Msi) * p(Jk|y,Msi)
    
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
    % Normalize posterior distributions over observed change point's positions %
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
    
    % Normalize the posterior distribution for the stochastic-to-
    % probabilistic observer
    if     islin, pobsJkgYKMsp =     ppobsJkgYMsp ./ pYgMsp(K);
    elseif islog, pobsJkgYKMsp = exp(ppobsJkgYMsp  - pYgMsp(K));
    end
    
    % Normalize the posterior distribution for the stochastic-to-
    % deterministic observer
    if     islin, pobsJkgYKMsd =     ppobsJkgYMsd ./ pYgMsd(K);
    elseif islog, pobsJkgYKMsd = exp(ppobsJkgYMsd  - pYgMsd(K));
    end
    
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
    % Compute the posterior beliefs over model's parameters for the %
    % the stochastic-to-probabilistic model                         %
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
    % Posterior over model's parameters (theta) marginalized over change
    % point's positions
    pTgYMsp(:,K,:) = sum(pobsJkgYKMsp .* pTgY, 2);
    
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
    % Compute the posterior beliefs over model's parameters for the %
    % the stochastic-to-deterministic model                         %
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %   
    % Posterior over patterns marginalized over change point's positions
    pTgYMsd(:,K) = pobsJkgYKMsd * pRgY;
    
    %% Compute predictions regarding the identity of the forthcoming observation
    %  =========================================================================
    %  The prediction is computed after receiving the observation K:
    %  for all i E {P,D}: p(A|y,Msi) = p(y_K+1=A|y_1:K,Msi) = 1 - p(B|y,Msi)
    %  
    %  It is obtained by deriving the expectation of a A given the
    %  observations following the change point times the posterior
    %  probability of the change point being at that point:
    %  for all i E {P,D}: p(y(K+1)=A|y,Jk,Msi) = p(y_K+1=A|y,Msi) * p(Jk|y)
    
    % for the stochastic-to-probabilistic observer
    pAgYMsp(K) = pAgJkYp2Mp * pJkgYMsp(:,K);
    
    % for the stochastic-to-deterministic observer
    pAgYMsd(K) = pAgJkYp2Md * pJkgYMsd(:,K);
end

%% Compute the posterior probability of each model
%  ===============================================

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Define models' prior probability %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

% We use a uniform uninformative Bayes-Laplace prior, which is in that case:
% p(Msi) = 1/{M} = 1/3
% with i E {D,P,S} and where {M} is the number of possible models (i.e. 3).
if     islin, p_pMi = 1/3;
elseif islog, p_pMi = -log(3);
end
p_pMss = p_pMi; % fully-stochastic model
p_pMsp = p_pMi; % stochastic-to-probabilistic model
p_pMsd = p_pMi; % stochastic-to-deterministic model
p_pMsi = [p_pMss; p_pMsp; p_pMsd];

% ~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Get Bayes' rule numerator %
% ~~~~~~~~~~~~~~~~~~~~~~~~~ %

% This term is proportional to the posterior over models
pYgMsi = [pYgMss; pYgMsp; pYgMsd];
if     islin, ppMsigY = pYgMsi .* p_pMsi;
elseif islog, ppMsigY = pYgMsi  + p_pMsi;
end

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Compute the normalization factor %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

% The normalization factor for models' posterior probabilities is the sum
% of models' evidences over all possible models times models' prior
% probabilities
% p(y) = sum(i E {S,P,D}) p(y|Msi) * p(Msi)
if     islin, pY =     sum(    ppMsigY , 1);
elseif islog, pY = log(sum(exp(ppMsigY), 1));
end

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Compute the models' posterior probability %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

% Bayes' rule states that the model's posterior probability is
% proportional to its model evidence times its prior probability. Dividing
% this product by the normalization factor computed just above gives
% probabilities that evolve from 0 to 1.
% for all i E {s,p,d}: p(Msi|y) = (p(y|Msi) * p(Msi)) / p(y)
if     islin, pMsigY =     ppMsigY ./ pY;
elseif islog, pMsigY = exp(ppMsigY  - pY);
end
pMssgY = pMsigY(1,:);
pMspgY = pMsigY(2,:);
pMsdgY = pMsigY(3,:);

% N.B. The sequence begins with a fully stochastic part, therefore beliefs 
% at the very beginnning are simply the following ones.
pMssgY(1) = 1; % p(Mss|y_1) = 1
pMspgY(1) = 0; % p(Msp|y_1) = 0
pMsdgY(1) = 0; % p(Msd|y_1) = 0

% ~~~~~~~~~~~~~~~~~~ %
% Get the best model %
% ~~~~~~~~~~~~~~~~~~ %

% Get the model with the highest posterior probability
[~, Mhat] = max(pMsigY, [], 1);

% Detect when one of the model reaches a significant threshold
pThr = 1/2; % the minimal probability at which a model is better than any other ones
dMss = pMssgY > pThr;
dMsp = pMspgY > pThr;
dMsd = pMsdgY > pThr;

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Entropy of the posterior distribution over models %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

% The entropy of discrete distributions is defined as:
% H(p(Mi|y)) = sum_{i} p(Mi|y) * log(1/p(Mi|y))
HpMgY = cellfun(@Emergence_IO_Entropy, mat2cell(pMsigY, 3, ones(1,N)));

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Update of beliefs regarding the models' posterior probability %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

% The iterative update is defined as:
% D_JS = (1/2) * (D_KL(p(Mi|y_1:K-1)||((p(Mi|y_1:K-1)+p(Mi|y_1:K))/2) +
%                 D_KL(p(Mi|y_1:K)  ||((p(Mi|y_1:K-1)+p(Mi|y_1:K))/2))
% Where D_KL(p(Mi|y_1:K-1)||p(Mi|y_1:K) = sum_{i} p(Mi|y_1:K-1) * log2(p(Mi|y_1:K-1) / p(Mi|y_1:K))
JSpMgY = Emergence_IO_DistDist(pMsigY');

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Compute the iterative covariance between current and previous beliefs %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

% Convert posterior beliefs from matrix to cells
pB = mat2cell(pMsigY(:,1:end-1), 3, ones(1,N-1)); % previous beliefs
cB = mat2cell(pMsigY(:,2:end),   3, ones(1,N-1)); %  current beliefs

% Compute iterative covariance
covpMsi = cellfun(@(Xi,Xj) cov(Xi, Xj), pB, cB, 'UniformOutput', 0);
covpMsi = [NaN, cellfun(@(x) x(2,1), covpMsi, 'UniformOutput', 1)];

%% Compute the posterior probability of change point's position using BMA
%  ======================================================================

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Bayesian model averaging regarding posterior beliefs about the %
% position of the change point                                   %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

% Since, in the fully-stochastic model, there is no assumption about
% any change point, we just need to combine the beliefs of the
% probabilistic-to-stochastic and stochastic-to-deterministic models
% regarding the position of the change point:
% p(Jk|y) propto p(Jk|y,Msp) * p(Msp) + p(Jk|y,Msd) * p(Msd)
pJkgY = (pJkgYMsp .* pMspgY + ... % p(Jk|y,Msp) * p(Msp|y)
         pJkgYMsd .* pMsdgY) ...  % p(Jk|y,Msd) * p(Msd|y)
        ./ (pMspgY + pMsdgY);     % normalized by p(Msp|y) + p(Msd|y)

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Entropy of the posterior distribution over change point's position %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

% Compute the entropy of the posterior distributions over change point's
% positions (see previous chunks for the mathematical definition)
HpJkgYMsp = cellfun(@Emergence_IO_Entropy, mat2cell(pJkgYMsp, N, ones(1,N)));
HpJkgYMsd = cellfun(@Emergence_IO_Entropy, mat2cell(pJkgYMsd, N, ones(1,N)));
HpJkgY    = cellfun(@Emergence_IO_Entropy, mat2cell(pJkgY,    N, ones(1,N)));

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Update of beliefs regarding the position of the change point %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

% Compute the iterative updates of the posterior distributions over change
% point's positions (see previous chunks for the mathematical definition)
JSpJkgYMsp = Emergence_IO_DistDist(pJkgYMsp');
JSpJkgYMsd = Emergence_IO_DistDist(pJkgYMsd');
JSpJkgY    = Emergence_IO_DistDist(pJkgY');

%% Compute expectancy-related quantities
%  =====================================

% ~~~~~~~~~~ %
% Prediction %
% ~~~~~~~~~~ %
% This the expected probability of receiving an A at position K given the
% beliefs derived from the sequence up to position K-1:
% p(y_K=A|y_1:K-1,Msi)

% The prediction from the fully-stochastic observer is simply chance level
pAgYMss = ones(1,N) ./ 2; % p(A|y_0,Mss) = 1/2

% Before seeing any observation, the probability of observing a A is equal
% to chance level, which is 1/2 for a binary sequence
% for all i E {S,P,D}, p(A|y_0,Mi) = 1/2
pAgYMsp = cat(2, ones(1,2) ./ 2, pAgYMsp(2:N-1)); % p(A|y_0,Msp) = 1/2
pAgYMsd = cat(2, ones(1,2) ./ 2, pAgYMsd(2:N-1)); % p(A|y_0,Msd) = 1/2

% Probability that the forthcoming stimulus will be a A computed using
% Bayesian Model Averaging
% p(A|y) = sum(i) p(A|y,Mi) * p(Mi|y)
pAgY = (pAgYMss .* pMssgY) + ...   % p(A|y,Mss) * p(Mss|y)
       (pAgYMsp .* pMspgY) + ...   % p(A|y,Msp) * p(Msp|y)
       (pAgYMsd .* pMsdgY) ;       % p(A|y,Msd) * p(Msd|y)

% Posterior Bernoulli distributions
pXgYMsp = [pAgYMsp; 1-pAgYMsp];
pXgYMsd = [pAgYMsd; 1-pAgYMsd];
pXgY    = [pAgY   ; 1-pAgY   ];

% ~~~~~~~~~~~~~~~~~~~~~~~~ %
% Equivalent learning rate %
% ~~~~~~~~~~~~~~~~~~~~~~~~ %

% Compute the prediction error (distance to prediction in lin scale)
pe = NaN(1,N);
pe(s == 1) = 1 -      pAgY(S == 1);
pe(s == 2) = 1 - (1 - pAgY(S == 2));

% Equivalent learning rate
% equivalent alpha = update / prediction error
eqalpha = JSpMgY ./ pe;
eqalpha(pe == 0) = 0;

% ~~~~~~~~ %
% Surprise %
% ~~~~~~~~ %
% Compute the surprise (distance to prediction in log2 scale):
% -log(p(X|y,Msi)) = -log(p(A|y,Msi))                      if y_K+1 = A
% -log(p(X|y,Msi)) = -log(p(B|y,Msi)) = -log(1-p(A|y,Msi)) if y_K+1 = B

% Find indices corresponding to A and B observations
A = find(S == 1);
B = find(S == 2);

% Remove first index
A = A(A > 1);
B = B(B > 1);

% Surprise evoked by the current observation (K) given previous beliefs
% (K-1) from the stochastic-to-probabilistic model
IpAgYMsp    = NaN(1,N);
IpAgYMsp(A) = -log2(  pAgYMsp(A));
IpAgYMsp(B) = -log2(1-pAgYMsp(B));

% Surprise evoked by the current observation (K) given previous beliefs
% (K-1) from the stochastic-to-deterministic model
IpAgYMsd    = NaN(1,N);
IpAgYMsd(A) = -log2(  pAgYMsd(A));
IpAgYMsd(B) = -log2(1-pAgYMsd(B));

% Surprise evoked by the current observation (K) given previous beliefs
% (K-1) marginalized (through Bayesian Model Averaging) over models
IpAgY    = NaN(1,N);
IpAgY(A) = -log2(  pAgY(A));
IpAgY(B) = -log2(1-pAgY(B));

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Entropy of the distribution of next observation being A/B %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

% Compute the entropy of the posterior distributions over the identity of
% the next observation (see previous chunks for the mathematical definition)
HpAgYMsp = cellfun(@Emergence_IO_Entropy, mat2cell(pXgYMsp, 2, ones(1,N)));
HpAgYMsd = cellfun(@Emergence_IO_Entropy, mat2cell(pXgYMsd, 2, ones(1,N)));
HpAgY    = cellfun(@Emergence_IO_Entropy, mat2cell(pXgY,    2, ones(1,N)));

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Update of beliefs regarding the likelihood of the next observation being A/B %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

% Compute the entropy of the posterior distributions over the identity of
% the next observation (see previous chunks for the mathematical definition)
JSpAgYMsp = Emergence_IO_DistDist(pXgYMsp');
JSpAgYMsd = Emergence_IO_DistDist(pXgYMsd');
JSpAgY    = Emergence_IO_DistDist(pXgY');

%% Wrap things up
%  ==============

% Prepare the structure output
out = [];

% ~~~~~~ %
% Inputs %
% ~~~~~~ %

% Sequence
out.in.S = S;

% Parameters
out.in.nu     = patlen;
out.in.pEp    = pEp;
out.in.pEd    = pEd;
out.in.decayP = decayP;
out.in.decayD = decayD;
out.in.stat   = stat;

% Prior knowledge
out.in.Jprior = Jprior;
out.in.p_pJ   = p_pJk;
out.in.p_pR   = p_pRi;
out.in.p_pT   = p_pTi;
out.in.p_pMss = p_pMss;
out.in.p_pMsp = p_pMsp;
out.in.p_pMsd = p_pMsd;

% Computation options
out.in.postgrid  = postprec;
out.in.scaleme   = scaleme;
out.in.computfor = computfor;
out.in.verbose   = verbose;

% ~~~~~~ %
% Models %
% ~~~~~~ %

% Sequence marginal likelihood
out.pY = pY;

% Export model evidences
out.pYgMss = pYgMss;
out.pYgMsp = pYgMsp;
out.pYgMsd = pYgMsd;

% Compute Bayes factors (in log or lin scale)
% BF(M1,M2) = p(y|M1) / p(y|M2)
if islin % Bayes factors
    out.BFspss = pYgMsp ./ pYgMss;
    out.BFsdss = pYgMsd ./ pYgMss;
    out.BFsdsp = pYgMsd ./ pYgMsp;
% <=> log(BF) = log(p(y|M1) / p(y|M2)) = log(p(y|M1)) - log(p(y|M2))
elseif islog % Log Bayes factors
    out.BFspss = pYgMsp - pYgMss;
    out.BFsdss = pYgMsd - pYgMss;
    out.BFsdsp = pYgMsd - pYgMsp;
end

% Return the posterior probability of each model
out.pMssgY = pMssgY;
out.pMspgY = pMspgY;
out.pMsdgY = pMsdgY;

% Return entropy of the posterior distribution over models
out.HpMgY = HpMgY;

% Return iterative updates of the posterior distribution over models
out.JSpMgY  = JSpMgY;
out.covpMsi = covpMsi;

% Return equivalent learning rate
out.PE      = pe;
out.eqAlpha = eqalpha;

% Return best models 
out.dMss = dMss;
out.dMsp = dMsp;
out.dMsd = dMsd;
out.Mhat = Mhat;

% ~~~~~~~~~~~~~~~~~~ %
% Models' parameters %
% ~~~~~~~~~~~~~~~~~~ %

% Return posterior distributions over models' parameters (marginalized over
% change points' positions)
out.pTgYMsp = pTgYMsp;
out.pTgYMsd = pTgYMsd;

% ~~~~~~~~~~~~ %
% Change point %
% ~~~~~~~~~~~~ %

% Return the posterior distributions over change point's positions
out.pJkgYMsp  = pJkgYMsp;
out.pJkgYMsd  = pJkgYMsd;
out.pJkgY     = pJkgY;

% Return entropy of the posterior distribution over change point's
% positions
out.HpJkgYMsp = HpJkgYMsp;
out.HpJkgYMsd = HpJkgYMsd;
out.HpJkgY    = HpJkgY;

% Return iterative updates of the posterior distributions over change
% point's positions
out.JSpJkgYMsp = JSpJkgYMsp;
out.JSpJkgYMsd = JSpJkgYMsd;
out.JSpJkgY    = JSpJkgY;

% ~~~~~~~~~~~~ %
% Expectations %
% ~~~~~~~~~~~~ %

% Return expectancies
out.pAgYMsp   = pAgYMsp;
out.pAgYMsd   = pAgYMsd;
out.pAgY      = pAgY;

% Return surprise
out.IpAgYMsp  = IpAgYMsp;
out.IpAgYMsd  = IpAgYMsd;
out.IpAgY     = IpAgY;

% Return entropy of the prediction
out.HpAgYMsp  = HpAgYMsp;
out.HpAgYMsd  = HpAgYMsd;
out.HpAgY     = HpAgY;

% Return iterative updates of the prediction
out.JSpAgYMsp = JSpAgYMsp;
out.JSpAgYMsd = JSpAgYMsd;
out.JSpAgY    = JSpAgY;

% ~~~~~~~~ %
% Conclude %
% ~~~~~~~~ %

% In the case of computations made after each observation
if verbose == 1 && strcmpi(computfor, 'all'), fprintf('. Done! '); end
if verbose == 2, fprintf('\n%s\n', repmat('-', 1, numel(cols))); end
    
% Get the names of all the output variables
f = fieldnames(out);
f = f(~strcmpi(f, 'in')); % except the "input" field

% In the case of computations made only after the last observation
if strcmpi(computfor, 'last')
    
    % Squeeze the second dimension that corresponds to the different
    % observations in the sequence
    for iVar = 1:numel(f)
        out.(f{iVar}) = squeeze(out.(f{iVar})(:,end,:));
    end
end

% Get the elapsed time since the beginning of the function
toc;
out.in.exectime = toc;

% Go back to line
if verbose == 1, fprintf('%s\n', repmat('=', 1, eqnum)); end

end
