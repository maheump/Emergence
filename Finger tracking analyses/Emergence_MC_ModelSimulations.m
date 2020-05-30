% This script runs simulations of the alternative ideal observed models. In
% particular we constrast the repeating pattern detector to various Markov
% chains of different orders. The simulations are run on sequences that
% were presented to the subjects. We considered several alternatives.
%
% NORMATIVE SINGLE-SYSTEM MODEL (i.e. direct comparison between statistics
% and rules, but same hypothesis space relying on Markov chains of more or
% less high order)
%   - "PseudoDeterministic": replaces the deterministic hypothesis with a
%       probabilistic hypothesis learning higher-order transition
%       probabilities.
%   - "BiasedPseudoDeterministic": same as the previous one but in addition
%       uses a prior biases for predictable cases.
%   - "DifferentPriors": probabilistic and deterministic hypotheses use the
%       same Markov chain (we vary the order) but with different prior
%       distribution (respectively biased vs. flat).
%
% NON-COMMENSURATE TWO-SYSTEM MODEL (i.e. distinct hypothesis spaces for
% statistics and rules, but no direct comparison)
%   - "IndependentDiffContinuous": uses a sigmoid on the difference between
%       independently-computed likelihoods in the regular hypotheses in
%       order to combine them (instead of the rules of probability).
%   - "IndependentDiffDiscrete": uses example extreme cases of the previous
%       one (i.e. a max, linear and example sigmoid versions).
%   - "IndependentRatioContinuous": uses a sigmoid on the log-ratio between
%       independently-computed likelihoods in the regular hypotheses in
%       order to combine them (instead of the rules of probability).
%   - "IndependentRatioDiscrete": uses example extreme cases of the
%       previous one (i.e. a max, linear and example sigmoid versions).
%
% OTHER ALTERNATIVES:
%   - "TreeDepth": uses a deterministic hypothesis that considers different
%       maximum pattern length than the one of the longest pattern (i.e.
%       10) used in the experiment.
%   - "Leak": uses an exponential leak when counting observations in the
%       case of the probabilistic hypothesis.
% 
% Copyright (c) 2020 Maxime Maheu

%% DEFINE OPTIONS OF THE IDEAL OBSERVERS TO SIMULATE
%  =================================================

% Define default options for the observer (this is the *ideal* observer
% given the structure of the task, the one that is analysed in the main
% text of the paper)
defo        = [];                % create empty structure
defo.pEd    = 0;                 % probability of making a memory error at each observation
defo.pEp    = 0;                 % probability of making a memory error at each observation
defo.patlen = 10;                % depth of the rules' tree to explore
defo.stat   = 'Transitions';     % the type of statistics to learn
defo.pR     = 'Size-principle';  % the prior probability of each rule depends on its length
defo.pT     = 'Bayes-Laplace';   % prior over transitions
defo.pJ     = 'Bayes-Laplace';   % prior over change point's position
defo.comp   = 'all';             % compute after each observation
defo.scale  = 'log';             % scale of the model evidence
defo.pgrid  = [];                % precision of the posterior over theta
defo.verb   = 0;                 % do not output messages in the command window

% Get the names of the options
lab = fieldnames(defo);

% Create a useful compact function for the follcwing command
strfun = @(x) contains(SimuType, x, 'IgnoreCase', true);

% SIMU #1: Alternative deterministic hypotheses that estimate higher-order
% transitions
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if strfun('PseudoDeterministic')
    
    % Define models to test
    orders = 2:9;
    nMod = numel(orders);
    
    % Replicate the default options for all the models to simulate
    options = struct2cell(defo);
    options = repmat(options, 1, nMod);
    
    % Change the statistics to learn for these models
    lidx = strcmpi(lab, 'stat');
    options(lidx,:) = arrayfun(@(x) sprintf('Chain%1.0f', x), orders, 'uni', 0);
    
% SIMU #2: Alternative deterministic hypotheses that estimate higher-order
% transitions and preferring predictive cases
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	if strfun('Biased')
        lidx = strcmpi(lab, 'pT');
        options(lidx,:) = arrayfun(@(x) repmat(999/1000, 2, 2^x), orders, 'uni', 0);
    end
    
% SIMU #3: Alternative probabilistic and deterministic hypotheses which
% both estimate (low to high-order) transition probabilities but using
% respectively a prior flat or biased for predictive cases
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
elseif strfun('DifferentPriors')
    
    % Define models to test
    orders = repmat(1:9, 1, 2);
    nMod = numel(orders);
    
    % Replicate the default options for all the models to simulate
    options = struct2cell(defo);
    options = repmat(options, 1, nMod);
    
    % Change the statistics to learn for these models
    lidx = strcmpi(lab, 'stat');
    options(lidx,:) = arrayfun(@(x) sprintf('Chain%1.0f', x), orders, 'uni', 0);
    
    % Half of the models use a flat prior and the remaining half use a
    % prior biased for predictability
    lidx = strcmpi(lab, 'pT');
    options(lidx,nMod/2+1:end) = arrayfun(@(x) repmat(999/1000, 2, 2^x), orders(nMod/2+1:end), 'uni', 0);
    
    
% SIMU #4: Alternative probabilistic hypothesis assuming a local integration
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
elseif strfun('Leak')
    
    % Define leak/substitution error parameter to test
    pe = [0.3158158158158160 0.1966966966966970 0.1416416416416420 0.1106106106106110 ...
          0.0905905905905906 0.0765765765765766 0.0665665665665666 0.0585585585585586 ...
          0.0525525525525526 0.0475475475475475 0.0435435435435435 0.0400400400400400 ...
          0.0370370370370370 0.0345345345345345 0.0320320320320320 0.0305305305305305 ...
          0.0285285285285285 0.0270270270270270 0.0255255255255255 0.0245245245245245 ...
          0.0230230230230230 0.0220220220220220 0.0215215215215215 0.0205205205205205 ...
          0.0195195195195195 0.0190190190190190 0.0180180180180180 0.0175175175175175 ...
          0.0170170170170170 0.0165165165165165 0.0160160160160160 0.0155155155155155 ...
          0.0150150150150150 0.0145145145145145 0.0140140140140140 0.0135135135135135 ...
          0.0130130130130130 0.0125125125125125 0.0120120120120120 0.0115115115115115 ...
          0.0110110110110110 0.0105105105105105 0.0100100100100100 0];  
    nMod = numel(pe);
    
    % Replicate the default options for all the models to simulate
    options = struct2cell(defo);
    options = repmat(options, 1, nMod);
    
    % Change the substitution error probability for these models
    lidx = strcmpi(lab, 'pEp');
    options(lidx,:) = arrayfun(@(x) x, pe, 'uni', 0);
    
    % Also specify the precision grid since in the case of leaky
    % integration, only grid based computations are available
    lidx = strcmpi(lab, 'pgrid');
    options(lidx,:) = {1e-1};
    
    % For the ultimate model (no error), use analytical solutions
    options(lidx,end) = {[]};
    
% SIMU #5: Alternative deterministic hypothesis with different tree depth
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
elseif strfun('TreeDepth')
    
    % Define depth of trees to explore
    nu = [4 6 8 20 50 N];
    nMod = numel(nu);
    
    % Replicate the default options for all the models to simulate
    options = struct2cell(defo);
    options = repmat(options, 1, nMod);
    
    % Change the depth of the tree for these models
    lidx = strcmpi(lab, 'patlen');
    options(lidx,:) = arrayfun(@(x) x, nu, 'uni', 0);
    
% SIMU #6: Alternative (non-normative) weighting of the regular hypotheses
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
elseif strfun('Independent')
    
    % Define models to test
    if     strfun('Continuous'), sigmparam = logspace(1/2,10-1/2,100);
    elseif strfun('Discrete'),	 sigmparam = [-2, 10^1, -1];
    end
    nMod = numel(sigmparam);
    
    % Replicate the default options for all the models to simulate
    options = struct2cell(defo);
    options = repmat(options, 1, nMod);
    
    % Specify the type of models that are simulated
    if strfun('Discrete')
        options(end+1,:) = {'Lin', 'Sigm', 'Max'};
    elseif strfun('Continuous')
        if     strfun('Ratio'), labsigm = 'logratio';
        elseif strfun('Diff'),  labsigm = 'difference';
        end
        options(end+1,:) = arrayfun(@(x) sprintf('sigm(%s,%1.2f)', ...
            labsigm, x), sigmparam, 'uni', 0);
    end
end

% Create a name for the file
fname = sprintf('Emergence_MC_%s', SimuType);

% Prepare the output variable
mIO = cell(nSeq,nSub,nMod); % conditions x subjects cell matrix with IO's inference
mIO = cellfun(@(x) NaN(N,3), mIO, 'uni', 0);

%% RUN SIMULATIONS
%  ===============

% Try to load existing simulations with the same options
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
try
    % Try to load data
    exdata = load(fname);
    
    % Compare options of the previous simulation and specified options 
    compopt = cellfun(@(x,y) all(x(:) == y(:)), options, exdata.options);
    
    % If simulations have already been run, skip the simulations, otherwise
    % run all required simulations
    if all(compopt(:))
        mIO = exdata.mIO;
        clear('exdata');
    elseif any(~compopt(:))
        error('No preexisting simulations'); % generate an error
    end

% Otherwise, run the (time-consuming) simulations
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
catch
    
    % Get previously computed simulations
    if strfun('Independent')
        mIO = repmat(cellfun(@(x) x.SeqLLH, IO, 'uni', 0), [1,1,nMod]);
        
    % Run simulations
    else
        % For each subject
        for iSub = 1:nSub
            
            % For each condition
            for iSeq = 1:nSeq
                fprintf('- Running models on sequence #%2.0f/%2.0f from subject #%2.0f/%2.0f.\n', ...
                    iSeq, nSeq, iSub, nSub);
                
                % Get the sequence
                seq = G{iSeq,iSub}.Seq;
                
                % For each order of the Markov chain
                for iMod = 1:nMod
                    fprintf('\t* Model %1.0f/%1.0f... ', iMod, nMod);
                    
                    % Run the model with its specific options
                    io = Emergence_IO_FullIO(seq, options{:,iMod});
                    mIO{iSeq,iSub,iMod} = cat(1, io.pYgMsp, io.pYgMsd, io.pYgMss)';
                end
                
                % Save temporary file
                save(sprintf('%s_tmp.mat', fname), 'mIO', 'options');
            end
        end
        
        % Save file containing all the simulations and delete temporary file
        save(sprintf('%s.mat', fname), 'mIO', 'options');
        delete(sprintf('%s_tmp.mat', fname));
    end
end

%% COMPUTE HYPOTHESIS LIKELIHOOD
%  =============================

% Get sequence likelihood for the ideal observer model
fullIO = cellfun(@(x) x.SeqLLH, IO, 'uni', 0);

% For simulations with a pseudo-deterministic hypothesis which uses higher-
% order Markov chains
if strfun('PseudoDeterministic')
    
    % Get sequence likelihood under different hypotheses
    pYgMsp = cellfun(@(x) x(:,1), fullIO, 'uni', 0); % 1st-order Markov chain
    pYgMsd = cellfun(@(x) x(:,2), fullIO, 'uni', 0); % pattern learner
    pYgMss = cellfun(@(x) x(:,3), fullIO, 'uni', 0); % fully-stochastic
    pYgMsc = cellfun(@(x) x(:,1), mIO,    'uni', 0); % higher-order Markov chains
    
    % Combine sequence likelihood under the pseudo-deterministic
    % hypothesies, with the probabilistic and fully-stochastic hypotheses
    mIO = cellfun(@(pH,pD,pS) cat(2,pH,pD,pS), ...
        repmat(pYgMsp, [1,1,nMod]), pYgMsc, repmat(pYgMss, [1,1,nMod]), ...
        'uni', 0);
    
    % Attribute a color to each model
    modc = flipud(winter(nMod));
    
% For simulations with same low- to high-order Markov chains for both
% probabilistic and deterministic hypotheses but with different priors
% (respectively flat and biased)
elseif strfun('DifferentPriors')
    
    % Get simulations that use either a flat or biased prior
    lidx = strcmpi(lab, 'pT');
    unifprior = strcmpi(options(lidx,:), 'Bayes-Laplace');
    biasprior = ~unifprior;
    nMod = sum(unifprior);
    
    % Get sequence likelihood
    pYgMsp = cellfun(@(x) x(:,1), mIO(:,:,unifprior), 'uni', 0); % pseudo-probabilistic
    pYgMsd = cellfun(@(x) x(:,1), mIO(:,:,biasprior), 'uni', 0); % pseudo-deterministic
    pYgMss = repmat(cellfun(@(x) x(:,3), fullIO, 'uni', 0), [1,1,nMod]); % fully-stochastic
    
    % Combine sequence likelihood under each hypothesis
    mIO = cellfun(@(pH,pD,pS) cat(2,pH,pD,pS), pYgMsp, pYgMsd, pYgMss, 'uni', 0);
    
    % Adapt options to reflect modified simulations
    lidx1 = strcmpi(lab, 'patlen');
    lidx2 = strcmpi(lab, 'stat');
    options(lidx1,unifprior) = options(lidx2,biasprior);
    lidx1 = strcmpi(lab, 'pR');
    lidx2 = strcmpi(lab, 'pT');
    options(lidx1,unifprior) = options(lidx2,biasprior);
    options = options(:,unifprior);
    
    % Attribute a color to each model
    modc = flipud(winter(nMod));

% For stimulations with different maximum pattern length allowed
elseif strfun('TreeDepth')
    modc = flipud(Emergence_Colormap('Reds', nMod));
    
% For stimulations with exponential leak
elseif strfun('Leak')
    modc = flipud(autumn(nMod));
    
% For simulations with independent (non rational) combination of hypothesis
% likelihood
elseif strfun('Independent')
	
    % For each subject, each sequence and each model
    for iSub = 1:nSub
        for iSeq = 1:nSeq
            for iMod = 1:nMod
                
                % Get sequence likelihood
                pYgHp = mIO{iSeq,iSub,iMod}(:,1);
                pYgHd = mIO{iSeq,iSub,iMod}(:,2);
                pYgHs = mIO{iSeq,iSub,iMod}(:,3);
                
                % Compute posterior likelihood of regular hypotheses
                % independently
                qHpgY = exp(pYgHp) ./ (exp(pYgHs) + exp(pYgHp));
                qHdgY = exp(pYgHd) ./ (exp(pYgHs) + exp(pYgHd));
                
                % Get the slope parameter of the sigmoid function to use
                slope = sigmparam(iMod);
                
                % Relative weighting
                % ~~~~~~~~~~~~~~~~~~
                if slope == -2
                    Wp = qHpgY ./ (qHpgY + qHdgY);
                    
                % Maximum a posteriori
                % ~~~~~~~~~~~~~~~~~~~~
                elseif slope == -1
                    dif = qHpgY - qHdgY;
                    Wp = double(dif > 0);
                    Wp(dif == 0) = 1/2;
                    
                % Sigmoid weighting
                % ~~~~~~~~~~~~~~~~~
                elseif slope > 0
                    if strfun('Ratio') % log ratio between regular hypothesis likelihoods
                        Wp = 1 ./ (1 + exp(-slope .* log(qHpgY ./ qHdgY)));
                    elseif strfun('Diff') % difference between regular hypothesis likelihoods
                        Wp = 1 ./ (1 + exp(-slope .* (qHpgY - qHdgY)));
                    end
                end
                
                % Compute the likelihood of the three hypotheses
                pHsgY = (1 - qHpgY) .* Wp + (1 - qHdgY) .* (1 - Wp);
                pHpgY = (1 - pHsgY) .* Wp;
                pHdgY = (1 - pHsgY) .* (1 - Wp);
                mIO{iSeq,iSub,iMod} = cat(2, pHpgY, pHdgY, pHsgY);
            end
        end
    end
    
    % Append the fully ideal observer model
    pMgY = cat(3, mIO, cellfun(@(x) x.BarycCoord, IO, 'uni', 0));
    options = cat(2, options, cat(1, struct2cell(defo), 'Rational'));
    
    % Attribute a color to each model
    modc = cool(nMod);
end

% If the posterior probabilities have not yet been computed
if ~strfun('Independent')
    
    % Append the fully ideal observer model
    mIO = cat(3, mIO, fullIO);
    options = cat(2, options, struct2cell(defo));
    
    % Compute posterior probability over hypotheses
    pMgY = cellfun(@(x) exp(x) ./ sum(exp(x), 2), mIO, 'uni', 0);
end

% Wrap things up
modc = [modc; zeros(1,3)];
nMod = size(pMgY, 3);

% Make all trajectories to start at the bottom tip of the triangle
for i = 1:numel(pMgY), pMgY{i}(1,:) = [0,0,1]; end

% Remove useless datasets (i.e. sequence likelihoods)
clear('mIO');
