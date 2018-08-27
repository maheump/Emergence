% Compare the performance of the Bayesian deterministic pattern learner to
% a Bayesian observer using a n-order Markov chain at detecting a repeating
% pattern from noise.
% 
% Copyright (c) 2018 Maxime Maheu

%% INITIALIZATION
%  ==============

% Clear the place
clear;
close('all');

% Add ideal observer functions to the MATLAB path
scriptpath = mfilename('fullpath');
folderpath = scriptpath(1:strfind(scriptpath,'Emergence')+8);
addpath(genpath(folderpath));

% Set default figure properties
Emergence_DefaultFigureProperties;

%% CREATE A SEQUENCE
%  ================= 

% Define the repeating rule
rule = 'AABABABB';
nr = numel(rule);

% Define sequence properties
cp = 100;   % position of the change point
N = 200;    % length of the sequence
lp1 = cp;   % length of the first fully-stochastic part
lp2 = N-cp; % length of the second deterministic part

% Generate the sequence
seqp1 = GenRandSeq(lp1, 1/2);                     % the first part
seqp2 = repmat(str2pat(rule), [1, ceil(lp2/nr)]); % the second part
seq = [seqp1, seqp2(1:lp2)];                      % the full sequence

%% RUN THE BAYESIAN IDEAL OBSERVERS
%  ================================

% Order of the Markov chains to estimate
orders = 1:8;
nOrd = numel(orders);

% Define options for the observer
pEd   = 0; % probability of making a memory error at each observation
pEp   = 0; % probability of making a memory error at each observation
treed = N; % depth of the rules' tree to explore
p_pR  = 'Size-principle'; % the prior probability of each rule depends on its length
p_pT  = 'Bayes-Laplace'; % the prior over statistics to be learnt
p_pJ  = 'Uniform'; % prior over change point's position
comp  = 'all'; % update beliefs after each observation
scale = 'log'; % scale of the model evidence
pgrid = []; % precision of the posterior over theta
verb  = 0; % output some messages in the command window

% Prepare the output variable
io = cell(nOrd,1);

% Run the observer with these options
for iOrd = 1:nOrd
    
    % Define the statistics to be learned by the probabilistic model
    stat  = sprintf('Chain%1.0f', orders(iOrd));
    
    % Run the full ideal observer
    io{iOrd} = Emergence_IO_FullIO(seq, ... % binary sequence
        pEd, pEp, treed, stat, p_pR, p_pT, p_pJ, comp, scale, pgrid, verb); % options
end

% Get sequence likelihood under each model
pYgMss = cell2mat(cellfun(@(x) x.pYgMss, io, 'UniformOutput', 0));
pYgMsp = cell2mat(cellfun(@(x) x.pYgMsp, io, 'UniformOutput', 0));
pYgMsd = cell2mat(cellfun(@(x) x.pYgMsd, io, 'UniformOutput', 0));

% Compute posterior model probabilities against each other
if     strcmpi(scale, 'lin'), fun = @(x) x;
elseif strcmpi(scale, 'log'); fun = @exp;
end
pMpvssgY = fun(pYgMsp) ./ (fun(pYgMsp) + fun(pYgMss));
pMdvssgY = fun(pYgMsd) ./ (fun(pYgMsd) + fun(pYgMss));

%% DISPLAY THE RESULT OF THE INFERENCE
%  ===================================

% Prepare a new window
figure('Units', 'Normalized', 'Position', [0.4 0.2 0.2 0.8]);

% For each order of the Markov chain
for iOrd = 1:nOrd
    subplot(nOrd, 1, iOrd); lgd = NaN(1,3);
    
    % Display chance level
    plot([1,N], ones(1,2)./2, '--', 'Color', ones(1,3)./2); hold('on');
    
    % Display change point's position
    plot(repmat(cp,1,2), [0,1], 'k-');
    
    % Display posterior beliefs post change point
    plot(1:N, pMdvssgY(iOrd,:), '-', 'Color', [239 059 033]./255, 'LineWidth', 1);
    plot(1:N, pMpvssgY(iOrd,:), '-', 'Color', [066 146 198]./255, 'LineWidth', 3);

    % Customize the axis
    set(gca, 'TickLabelInterpreter', 'LaTeX');

    % Add some text labels
    ylabel('$p(\mathcal{M}_{i}|y)$', 'Interpreter', 'LaTeX');
    title(sprintf('%i-order Markov chain', orders(iOrd)), 'Interpreter', 'LaTeX');
end
xlabel('Observation ($K$)', 'Interpreter', 'LaTeX');
