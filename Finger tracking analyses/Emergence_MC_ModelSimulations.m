% This script runs simulations of the alternative ideal observed models. In
% particular we constrast the repeating pattern detector to various Markov
% chains of different orders. The simulations are run on sequences that
% were presented to the subjects.
% 
% Copyright (c) 2018 Maxime Maheu

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

% SIMU #1: Alternative deterministic hypotheses that estimate higher-order
% transitions
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if strcmpi(SimuType, 'PseudoDeterministic')
    
    % Define models to test
    orders = 1:9;
    nMod = numel(orders);
    
    % Replicate the default options for all the models to simulate
    options = struct2cell(defo);
    options = repmat(options, 1, nMod);
    
    % Change the statistics to learn for these models
    lidx = strcmpi(lab, 'stat');
    options(lidx,:) = arrayfun(@(x) sprintf('Chain%1.0f', x), orders, 'UniformOutput', 0);
    
% SIMU #2: Alternative deterministic hypotheses that estimate higher-order
% transitions and preferring predictive cases
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
elseif strcmpi(SimuType, 'BiasedPseudoDeterministic')
    
    % Define models to test
    orders = 1:9;
    nMod = numel(orders);
    
    % Replicate the default options for all the models to simulate
    options = struct2cell(defo);
    options = repmat(options, 1, nMod);
    
    % Change the statistics to learn and the prior over transitions for
    % these models
    lidx = strcmpi(lab, 'stat');
    options(lidx,:) = arrayfun(@(x) sprintf('Chain%1.0f', x), orders, 'UniformOutput', 0);
    lidx = strcmpi(lab, 'pT');
    options(lidx,:) = arrayfun(@(x) repmat(999/1000, 2, 2^x), orders, 'UniformOutput', 0);
    
% SIMU #3: Alternative probabilistic hypothesis assuming a local integration
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
elseif strcmpi(SimuType, 'Leak')
    
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
    options(lidx,:) = arrayfun(@(x) x, pe, 'UniformOutput', 0);
    
    % Also specify the precision grid since in the case of leaky
    % integration, only grid based computations are available
    lidx = strcmpi(lab, 'pgrid');
    options(lidx,:) = {1e-1};
    
    % For the ultimate model (no error), use analytical solutions
    options(lidx,end) = {[]};
    
% SIMU #4: Alternative deterministic hypothesis with different tree depth
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
elseif strcmpi(SimuType, 'TreeDepth')
    
    % Define depth of trees to explore
    nu = [unique(cellfun(@numel, dr)') 20 50 N];
    nMod = numel(nu);
    
    % Replicate the default options for all the models to simulate
    options = struct2cell(defo);
    options = repmat(options, 1, nMod);
    
    % Change the depth of the tree for these models
    lidx = strcmpi(lab, 'patlen');
    options(lidx,:) = arrayfun(@(x) x, nu, 'UniformOutput', 0);
end

% Create a name for the file
fname = sprintf('Emergence_MC_%s', SimuType);

% Prepare the output variable
mIO = cell(nSeq,nSub,nMod); % conditions x subjects cell matrix with IO's inference
mIO = cellfun(@(x) NaN(N,3), mIO, 'UniformOutput', 0);

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
    % For each subject
    for iSub = 1:nSub
        
        % For each condition
        for iSeq = 1:nSeq
            fprintf('- Running the IO on sequence #%2.0f/%2.0f from subject #%2.0f/%2.0f.\n', ...
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

%% COMPUTE HYPOTHESIS LIKELIHOOD
%  =============================

% For simulations with high-order Markov chains use them as
% pseudo-deterministic hypothesis, not as a higher-order probabilistic
% hypothesis
if contains(SimuType, 'PseudoDeterministic', 'IgnoreCase', true)
    
    % Get sequence likelihood under different hypotheses
    pYgMsp = cellfun(@(x) x(:,1), mIO(:,:,1), 'UniformOutput', 0); % 1st-order Markov chain
    pYgMsd = cellfun(@(x) x(:,2), mIO(:,:,1), 'UniformOutput', 0); % pattern learner
    pYgMss = cellfun(@(x) x(:,3), mIO(:,:,1), 'UniformOutput', 0); % fully-stochastic
    pYgMsc = cellfun(@(x) x(:,1), mIO,        'UniformOutput', 0); % higher-order Markov chains
    
    % Get sequence likelihood for the fully ideal observer modem
    fullIO = mIO(:,:,1);
    
    % Combine sequence likelihood under the pseudo-deterministic
    % hypothesies, with the probabilistic and fully-stochastic hypotheses
    mIO = cellfun(@(pH,pD,pS) cat(2,pH,pD,pS), ...
        repmat(pYgMsp, [1,1,nMod]), pYgMsc, repmat(pYgMss, [1,1,nMod]), ...
        'UniformOutput', 0);
    
    % Append the fully ideal observer model
    mIO = cat(3, mIO, fullIO);
    options = cat(2, options, struct2cell(defo));
    
    % Attribute a color to each model
    modc = [flipud(autumn(nMod)); tricol(2,:)];
elseif strcmpi(SimuType, 'Leak')
    modc = [flipud(winter(nMod-1)); zeros(1,3)];
end

% Compute posterior probability over hypotheses
pMgY = cellfun(@(x) exp(x) ./ sum(exp(x), 2), mIO, 'UniformOutput', 0);
nMod = size(pMgY, 3);

% Define starting point
for i = 1:numel(pMgY), pMgY{i}(1,:) = [0,0,1]; end

% Remove datasets 
clear('mIO');
