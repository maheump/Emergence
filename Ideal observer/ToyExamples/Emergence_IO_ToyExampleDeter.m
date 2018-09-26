% This script shows the learning dynamic of a deterministic (i.e. assuming
% repetitions of the same pattern) ideal Bayesian observer presented with
% various sequences made of the repetition of different patterns.
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

%% PREPARE SEQUENCES
%  =================

% Define patterns to test
patterns = {'AABB', ...
            'AAABBB', 'AABABB', 'AAABAB'...
            'AAAABBBB', 'AABABABB', 'AAABBABB', ...
            'AAAAABBBBB', 'AABABABABB', 'AAABAABBAB'};

% Get size of sequences
L = cellfun(@numel, patterns); % length of each pattern
maxL = max(L);                 % maximum patterns' length
nRep = 3;                      % number of repetition of the longest pattern
nObs = nRep*maxL;              % number of observations

% Create sequences based on these patterns
patternsrecod = cellfun(@str2pat, patterns, 'UniformOutput', 0);
Seq = cellfun(@(x) repmat(x, 1, nObs), patternsrecod, 'UniformOutput', 0);
Seq = cellfun(@(x) x(1:nObs), Seq, 'UniformOutput', 0);
nSeq = numel(Seq);

%% RUN THE BAYESIAN IDEAL OBSERVER
%  ===============================

% Properties of the Bayesian ideal observer
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Define the probability that the observer will make a mistake at each
% received observation (0 => perfect memory)
pErr = 0;
leak = Emergence_IO_Leak(pErr, nObs);

% Define properties of the observer
nu      = nObs;             % the longest possible pattern considered
scaleme = 'log';            % scale for the model evidence
prior   = 'Size-principle'; % the prior distribution over patterns

% Prepare output variables
pRgY   = NaN(nu,nObs,nSeq);
pYgMd  = NaN(nSeq,nObs);
pAgYMd = NaN(nSeq,nObs);
IgYMd  = NaN(nSeq,nObs);
JSdiv  = NaN(nSeq,nObs);
HpRgY  = NaN(nSeq,nObs);

% Run the Bayesian ideal observer
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% For each sequence
for iSeq = 1:nSeq
    
    % Run the Bayesian ideal observer
    [pYgMd(iSeq,:), pAgYMd(iSeq,:), ~, IgYMd(iSeq,:), ~, ...
	 pRgY(:,:,iSeq), ~, HpRgY(iSeq,:), JSdiv(iSeq,:)] = ...
        Emergence_IO_RunIO(@Emergence_IO_Tree, ... % IO learning repeating patterns
        Seq{iSeq}, ...                             % current binary sequence
        {nu, scaleme, false, prior, leak, true});  % properties of the IO
end

% For the update, take the square root value of the Jensen-Shannon divergence
JSdiv = sqrt(JSdiv);

% Compare to a random toss scenario
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% (log-)Likelihood of any sequence under a null model
pYgMs = Emergence_IO_Null(1:nObs, scaleme);

% (log-)Likelihood ratio
if     strcmpi(scaleme, 'lin'), LR = pYgMd ./ pYgMs;
elseif strcmpi(scaleme, 'log'), LR = pYgMd  - pYgMs;
end

% Observers posterior probabilities
if     strcmpi(scaleme, 'lin'), pMdgY =     pYgMd  ./ (    pYgMd  +     pYgMs );
elseif strcmpi(scaleme, 'log'), pMdgY = exp(pYgMd) ./ (exp(pYgMd) + exp(pYgMs));
end

%% DISPLAY THE RESULT OF THE INFERENCE
%  ===================================

% Define the different variables from the IO to display
Vars = {'LR', 'pMdgY', 'pAgYMd', 'IgYMd', 'HpRgY'};
nVar = numel(Vars);
VarLab = {{'Likelihood ratio', '$\frac{p(y|\mathcal{M_{\mathrm{D}}})}{p(y|\mathcal{M}_{\rm{S}})}$'}, ...
          {'Posterior probability', '$\frac{p(\mathcal{M_{\mathrm{D}}}|y)}{p(\mathcal{M}_{\rm{S}}|y)}$'}, ...
          {'Prediction', '$p(y_{k}=\mathrm{A}|y_{1:k-1},\mathcal{M_{\mathrm{D}}})$'}, ...
          {'Surprise', '$-\log_{2} p(y_{k}|y_{1:k-1}\mathcal{M_{\mathrm{D}}})$'}, ...
          {'Entropy of the posterior', '$H(p(R|y,\mathcal{M_{\mathrm{D}}}))$'}, ...
          {'Model update', ['$\sqrt{D_\mathrm{JS}(p(R|y_{1:k},\mathcal{M_', ...
          '{\mathrm{D}}})||p(R|y_{1:k-1},\mathcal{M_{\mathrm{D}}}))}$']}};

% Create a new window
figure('Units', 'Normalized', 'Position', [0 1/2-0.2 1 0.4]);

% Choose different colors for each sequence (i.e. each pattern)
col = lines(nSeq);

% For each sequence (i.e. each pattern)
for iSeq = 1:nSeq
    
    % For each variable from the ideal observer
    for iVar = 1:nVar
        subplot(nVar+2, nSeq, iSeq+nSeq*(iVar-1));
        
        % Define vertical limits of the plot
        tp = eval(Vars{iVar});
        limy1 = [min([0,min(tp(:))]), max([1,max(tp(:))])];
        margin = diff(limy1).*(1/4);
        limy2 = limy1 + [-1,1].*margin;
        
        % Display each time a pattern is repeated
        rep = (L(iSeq):L(iSeq):nObs) + 1/2;
        plot(repmat(rep', [1,2])', repmat(limy2, [numel(rep),1])', ...
            '-', 'Color', ones(1,3)./2, 'LineWidth', 1/2); hold('on');
        
        % Display the sequence
        pos = (limy2 - limy1) / 2 + limy1;
        plot(find(Seq{iSeq} == 2), pos(1), 'k.', 'MarkerSize', 6);
        plot(find(Seq{iSeq} == 1), pos(2), 'k.', 'MarkerSize', 6);
        
        % Display the beliefs of the ideal observer
        plot(1:nObs, tp(iSeq,:), '.-', 'Color', col(iSeq,:), ...
            'LineWidth', 1, 'MarkerSize', 8); hold('on');
        
        % Customize the axes
        axis([1/2, nObs+1/2, limy2]);
        set(gca, 'XTick', [1, 5:5:nObs], 'TickLabelInterpreter', 'LaTeX');

        % Add some text labels
        if iVar == 1
            title(sprintf('[%s]$^n$', patterns{iSeq}), ...
                'Interpreter', 'LaTeX', 'Color', col(iSeq,:));
        end
        if iSeq == 1
            ylabel(VarLab{iVar}, 'Interpreter', ...
                'LaTeX', 'Rotation', 0, 'HorizontalAlignment', 'Right', ...
                'VerticalAlignment', 'Middle');
        end
    end
    
    % Display the posterior beliefs over patterns
    subplot(nVar+2, nSeq, iSeq + nSeq*nVar + [0,nSeq]);
    imagesc(1:nObs, 1:nu, pRgY(:,:,iSeq)); hold('on');
    
    % Display each time a pattern is repeated
    plot(repmat(rep', [1,2])', repmat([1/2, nObs+1/2], [numel(rep),1])', ...
        '-', 'Color', ones(1,3)./2, 'LineWidth', 1/2);
    
    % Customize the colormap
    colormap(parula); caxis([min(pRgY(:)), max(pRgY(:))]);
    
    % Customize the axes
    axis('xy');
    set(gca, 'XTick', [1, 5:5:nObs], 'TickLabelInterpreter', 'LaTeX');
    
    % Add some text labels
    xlabel('Observation ($K$)', 'Interpreter', 'LaTeX');
    if iSeq == 1
        ylabel({'Posterior', 'distribution', '$p(R_i|y,\mathcal{M_{\mathrm{D}}})$'}, ...
            'Interpreter', 'LaTeX', 'Rotation', 0, ...
            'HorizontalAlignment', 'Right', 'VerticalAlignment', 'Middle');
    end
end
