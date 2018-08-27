% This script displays the inference of different observers learning
% different types of probabilistic regularities when faced to the same
% input sequence that is biased in terms of first-order transition
% probabilities.
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

% Create a sequence biased toward repetitions
nObs = 200;
pAgB = 1/3;
pBgA = 1/3;
y = GenRandSeq(nObs, [pAgB, pBgA]);

%% RUN THE BAYESIAN IDEAL OBSERVERS LEARNING DIFFERENT PROBABILISTIC REGULARITIES
%  ==============================================================================

% Define the models to use
% N.B. order 0 is Bernoulli, 1 is Markov and greater than 1 are
% higher-order Markov chains
nOrder = [0, 1:4];
nMod = numel(nOrder);

% Prepare outputs
pYgMp  = NaN(nMod,nObs);
pAgYMp = NaN(nMod,nObs);
IgYMp  = NaN(nMod,nObs);
pTgYMp = cell(1,nMod);

% Define the observers to 
iofun = cat(2, {@Emergence_IO_Bernoulli}, {@Emergence_IO_Markov}, ...
    repmat({@Emergence_IO_Chain}, 1, sum(nOrder > 1)));

% Define the same properties for these observers
dt      = 0.01; % precision of the posterior
scaleme = 'log'; % scale of the model evidence
prior   = 'Bayes-Laplace'; % prior over observers' parameter(s)
options = cat(2, {{scaleme, false, prior, [], dt}}, ...
                 {{scaleme, false, prior, [], dt}}, ...
     cellfun(@(x) {scaleme,        prior,  x, dt}, ...
     num2cell(nOrder(nOrder>1)), 'UniformOutput', 0));

% Present the same input sequence to all the different observers
for iMod = 1:nMod
    [pYgMp(iMod,:), pAgYMp(iMod,:), ~, ~, IgYMp(iMod,:), pTgYMp{iMod}] = ...
        Emergence_IO_RunIO(iofun{iMod}, y, options{iMod});
end

% (log-)Likelihood of the sequence under a null model
pYgMs = Emergence_IO_Null(1:nObs, scaleme);

% (log-)Likelihood ratio
if     strcmpi(scaleme, 'lin'), LR = pYgMp ./ pYgMs;
elseif strcmpi(scaleme, 'log'), LR = pYgMp  - pYgMs;
end

% Observers posterior probabilities
if     strcmpi(scaleme, 'lin'), pMpgY =     pYgMp  ./ (    pYgMp  +     pYgMs );
elseif strcmpi(scaleme, 'log'), pMpgY = exp(pYgMp) ./ (exp(pYgMp) + exp(pYgMs));
end

%% DISPLAY THE RESULT OF THE INFERENCE
%  ===================================

% Define the variables to look at
Vars = {'LR', 'pMpgY', 'pAgYMp', 'IgYMp'};
VarLab = {{'Likelihood ratio', '$\frac{p(y|\mathcal{M}_{i})}{p(y|\mathcal{M}_{\rm{S}})}$'}, ...
          {'Posterior probability', '$p(\mathcal{M}_{i}|y)$'}, ...
          {'Prediction', '$p(y_{k}=\mathrm{A}|y_{1:k-1},\mathcal{M}_{i})$'}, ...
          {'Surprise', '$-\log_{2} p(y_{k}|y_{1:k-1}\mathcal{M}_{i})$'}};
nVar = numel(VarLab);

% Prepare a new figure
figure('Units', 'Normalized', 'Position', [0.1 0.5 0.8 0.4]);
col = lines(nMod);

% For each type of observer
for iMod = 1:nMod
    
    % Variables
    % ~~~~~~~~~
    
    % For each variable from the ideal observer
    for iVar = 1:nVar
        subplot(nVar+2, nMod, iMod+(iVar-1)*nMod);
        
        % Define vertical limits of the plot
        tp = eval(Vars{iVar});
        limy1 = [min([0,min(tp(:))]), max([1,max(tp(:))])];
        margin = diff(limy1).*(1/4);
        limy2 = limy1 + [-1,1].*margin;
        
        % Display the sequence
        pos = (limy2 - limy1) / 2 + limy1;
        plot(find(y == 2), pos(1), 'k.', 'MarkerSize', 6); hold('on');
        plot(find(y == 1), pos(2), 'k.', 'MarkerSize', 6);
        
        % Display the beliefs of the ideal observer
        plot(1:nObs, tp(iMod,:), '-', 'Color', col(iMod,:), ...
            'LineWidth', 2, 'MarkerSize', 8); hold('on');
        
        % Customize the axes
        axis([1/2, nObs+1/2, limy2]);
        set(gca, 'XTick', [1, get(gca, 'XTick')], 'TickLabelInterpreter', 'LaTeX');
            
        
        % Add some text labels
        if iMod == 1
            ylabel(VarLab{iVar}, 'Interpreter', ...
                'LaTeX', 'Rotation', 0, 'HorizontalAlignment', 'Right', ...
                'VerticalAlignment', 'Middle');
        end
        if iVar == 1
            ttl = func2str(iofun{iMod});
            title({ttl(max(strfind(ttl,'_'))+1:end), ...
                sprintf('Order of the transition: %i', nOrder(iMod))}, ...
                'Interpreter', 'LaTeX');
        end
    end
    
    % Posterior distribution over observer's parameters
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    % Display the posterior beliefs over patterns
    subplot(nVar+2, nMod, iMod + nMod*(nVar) + [0,nMod]);
    tp = pTgYMp{iMod}; nTdim = size(tp,3);
    tp = reshape(permute(tp, [2,1,3]), [nObs, size(tp,1)*nTdim])';
    imagesc(1:nObs, [], tp); hold('on');
    
    % Display the limit between the different marginal distributions
    plot([1,nObs], repmat(cumsum(repmat(1/dt, 1, nTdim)),2,1), 'k-');
    
    % Customize the colormap
    colormap(parula); caxis([min(tp(:)), max(tp(:))]);
    
    % Customize the axes
    axis('xy');
    set(gca, 'XTick', [1, get(gca, 'XTick')], 'TickLabelInterpreter', 'LaTeX');
    
    % Add some text labels
    xlabel('Observation ($K$)', 'Interpreter', 'LaTeX');
    if iMod == 1
        ylabel({'Posterior', 'distribution', '$p(\theta|y,\mathcal{M}_{i})$'}, ...
            'Interpreter', 'LaTeX', 'Rotation', 0, ...
            'HorizontalAlignment', 'Right', 'VerticalAlignment', 'Middle');
    end
end
