% Script showing the effect of the type of deterministic hypothesis and the
% strength of prior beliefs over frequency of transitions 
% 
% Copyright (c) 2018 Maxime Maheu

%% INITIALIZATION
%  ==============

% Clear the place
clear;
close('all');

% Add ideal observer functions to the MATLAB path
scriptpath = mfilename('fullpath');
ind = strfind(scriptpath,'Emergence');
folderpath = scriptpath(1:ind(end-1)+8);
addpath(genpath(folderpath));

% Set default figure properties
Emergence_DefaultFigureProperties;

%% DISPLAY EXAMPLE PRIOR DISTRIBUTIONS
%  ===================================

% Define the grid over prior (pseudo-)counts
grid = [999/1000, 1];
ng = numel(grid);

% Define the colors to use
cmap = flipud(winter(ng));

% Prepare a new window
figure('Position', [119 678 560 420]);

% Display prior distributions
pgrid = 0:0.01:1;
for i = ng:-1:1
    plot(pgrid, Emergence_IO_BetaPDF(pgrid, grid(i), grid(i)), ...
        'Color', cmap(i,:)); hold('on');
end

% Customize the axes
set(gca, 'Box', 'Off', 'YTick', [], 'TickLabelInterpreter', 'LaTeX');

% Add some text labels
legend(cellfun(@num2str, num2cell(grid), 'UniformOutput', 0), 'Interpreter', 'LaTeX');
xlabel('Outcome probability', 'Interpreter', 'LaTeX');
ylabel('Density', 'Interpreter', 'LaTeX');

%% EXAMPLE INFERENCE ON A SIMPLE FULLY-STOCHASTIC SEQUENCE WITH DIFFERENT
%  PRIOR BELIEFS
%  ======================================================================

% Create a fully-stochastic sequence
N = 200;
y = GenRandSeq(N, 1/2);

% Define the order of the Markov chain to estimate
order = 1;

% Define pseudo-prior counts
pN = 1/2; % in number of observations
customprior = repmat(pN, 2, 2^order); % same for all transitions

% Run the different observers
[pYgHpf, pAgYHpf] = Emergence_IO_RunIO(@Emergence_IO_Chain, y, {'log', 'Bayes-Laplace', order});
[pYgHpe, pAgYHpe] = Emergence_IO_RunIO(@Emergence_IO_Chain, y, {'log', customprior,     order});
[pYgHs,  pAgYMs ] = Emergence_IO_Null(1:200, 'log');

% Compute models posterior probability
pHpfgY = exp(pYgHpf) ./ (exp(pYgHpf) + exp(pYgHs));
pHpegY = exp(pYgHpe) ./ (exp(pYgHpe) + exp(pYgHs));

% Prepare a new window
figure('Position', [679.5000 678 560 420]);

% Display models' predictions
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~
subplot(3,1,1); l = NaN(1,2);
plot([0,N+1], ones(1,2)./2, 'k--'); hold('on');
l(1) = plot(1:N, pAgYHpf);
l(2) = plot(1:N, pAgYHpe, '.');

% Customize the axes
set(gca, 'XLim', [1,N], 'Ylim', [0,1], 'TickLabelInterpreter', 'LaTeX');

% Add some text labels
legend(l, {'Bayes-Laplace', 'Jeffreys'}, 'Interpreter', 'LaTeX');
ylabel('$p(\mathrm{A})$', 'Interpreter', 'LaTeX');
title(sprintf('Markov chain of order %1.0f', order), 'Interpreter', 'LaTeX');

% Display sequence marginal likelihood
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subplot(3,1,2);
plot(1:N, pYgHpf); hold('on');
plot(1:N, pYgHpe);
plot(1:N, pYgHs);

% Customize the axes
set(gca, 'XLim', [1,N], 'TickLabelInterpreter', 'LaTeX');

% Add some text labels
legend({'Bayes-Laplace', 'Jeffreys', 'Fully-stochastic'}, 'Interpreter', 'LaTeX');
ylabel('$p(y|\mathcal{M})$', 'Interpreter', 'LaTeX');

% Display posterior probability of hypotheses
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subplot(3,1,3);
plot([0,N+1], ones(1,2)./2, 'k--'); hold('on');
plot(1:N, pHpfgY); 
plot(1:N, pHpegY);

% Customize the axes
set(gca, 'XLim', [1,N], 'TickLabelInterpreter', 'LaTeX');

% Add some text labels
xlabel('Observation \#', 'Interpreter', 'LaTeX');
ylabel('$p(\mathcal{M}|y)$ vs $p(\mathcal{M}_{0}|y)$', 'Interpreter', 'LaTeX');

%% EXAMPLE INFERENCE ON A SEQUENCE WITH CHANGE POINT WITH DIFFERENT PRIOR
%  BELIEFS AND LEARNING DIFFERENT TRANSITIONS OF DIFFERENT ORDERS
%  ======================================================================

% Prepare a sequence with a change point
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

seq = [GenRandSeq(100, 1/2), GenRandSeq(100, [3/4, 3/4])];

% Run the ideal observers
% ~~~~~~~~~~~~~~~~~~~~~~~

% Define general options for the observer
pEd    = 0;                 % probability of making a memory error at each observation
pEp    = 0;                 % probability of making a memory error at each observation
patlen = 10;                % depth of the rules' tree to explore
pR     = 'Size-principle';  % the prior probability of each rule depends on its length
pJ     = 'Bayes-Laplace';   % prior over change point's position
comp   = 'all';             % compute after each observation
scale  = 'log';             % scale of the model evidence
pgrid  = [];                % precision of the posterior over theta
verb   = 0;                 % do not output messages in the command window

% Define the orders of the chain to test
order = [1,5,9];
no = numel(order);

% Define the type of prior to tests
customprior = {1, 999/1000};
np = numel(customprior);

% For each pair of prior type and statistics
io = cell(np,no);
for p = 1:np
    for o = 1:no
        
        % Get the prior over transitions to use
        if isnumeric(customprior{p}), pT = repmat(customprior{p}, 2, 2^order(o));
        else, pT = customprior{p}; end
        
        % Get the type of statistics to estimate
        stat = sprintf('Chain%1.0f', order(o));
        
        % Run the ideal observer model
        io{p,o} = Emergence_IO_FullIO(seq, ...                   % input binary sequence
                                      pEd, pEp, ...              % memory errors
                                      patlen, stat, ...          % what to learn
                                      pR, pT, pJ, ...            % prior distributions
                                      comp, scale, pgrid, verb); % options
    end
end

% Display posterior probabilities
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Prepare a new figure
figure('Position', [1241 678 560 420]);

% For each prior distribution
for p = 1:np
    
    % For each order of the chain
    for o = 1:no
        subplot(numel(order), numel(customprior), p+numel(customprior)*(o-1));
        
        % Display the posterior probability over hypotheses
        pMgY = [io{p,o}.pMspgY; io{p,o}.pMsdgY; io{p,o}.pMssgY]';
        Emergence_PlotBarycTraj(pMgY);
        
        % Display the position of the change point
        plot([100,100], [0,1], 'k-');
        
        % Customize the axes
        set(gca , 'TickLabelInterpreter', 'LaTeX');
        
        % Add some text labels
        xlabel('Observation \#', 'Interpreter', 'LaTeX');
        ylabel('$p(\mathcal{H}|y)$', 'Interpreter', 'LaTeX');
        title([sprintf('Markov chain of order %1.0f', order(o)), ...
            sprintf(' with a [%g %g] prior', repmat(customprior{p},1,2))], ...
            'Interpreter', 'LaTeX');
    end
end
