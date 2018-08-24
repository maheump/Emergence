% This script compares the learning dynamic of different patterns in a
% binary sequence detection by a deterministic (supposing repetitions of
% the same pattern) ideal (Bayesian) observer.
% 
% Copyright (c) 2018 Maxime Maheu

%% Initialization
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

% Define patterns to test
patterns = {'AABB', ...
            'AAABBB', 'AABABB', 'AAABAB'...
            'AAAABBBB', 'AABABABB', 'AAABBABB', ...
            'AAAAABBBBB', 'AABABABABB', 'AAABAABBAB'};

% Get size of sequences
L = cellfun(@numel, patterns); % length of each pattern
maxL = max(L); % maximum patterns' length
nRep = 3; % number of repetition of the longest pattern
nObs = nRep*maxL; % number of observations

% Create sequences based on these patterns
patternsrecod = cellfun(@str2pat, patterns, 'UniformOutput', 0);
Seq = cellfun(@(x) repmat(x, 1, nObs), patternsrecod, 'UniformOutput', 0);
Seq = cellfun(@(x) x(1:nObs), Seq, 'UniformOutput', 0);

% Add a fully-stochastic sequence
L(end+1) = NaN;
patterns{end+1} = 'Stochastic';
Seq{end+1} = GenRandSeq(nObs, 1/2); 
nSeq = numel(Seq);

%% Run the Bayesian ideal observer
%  ===============================

% Properties of the Bayesian ideal observer
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Define the probability that the observer will make a mistake at each
% received observation (0 => perfect memory)
pErr = 0;
leak = Emergence_IO_Leak(pErr, nObs);

% Define the depth of the tree to explore
nu = round(nObs/2);

% Prepare output variables
pRgY       = NaN(nu,nObs,nSeq);
log_pYgMd  = NaN(nSeq,nObs);
pAgYMd     = NaN(nSeq,nObs);
IgYMd      = NaN(nSeq,nObs);
JSdiv      = NaN(nSeq,nObs);
HpRgY      = NaN(nSeq,nObs);

% Run the Bayesian ideal observer
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% For each sequence
for iSeq = 1:nSeq
    
    % Run the Bayesian ideal observer
    [log_pYgMd(iSeq,:), pRgY(:,:,iSeq), pAgYMd(iSeq,:), ...
     IgYMd(iSeq,:), ~, JSdiv(iSeq,:), ~, HpRgY(iSeq,:)] = ...
        Emergence_IO_RunIO(@Emergence_IO_Tree, ... % IO learning repeating patterns
        Seq{iSeq}, ... % current binary sequence
        {nu, 'log', false, 'Size-principle', leak, true}); % properties of the IO
end

% Compare to a random toss scenario
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Sequence likelihood for a fully-stochastic observer
log_pYgMs = -log(2) .* (1:nObs);

% Log-likelihood ratio
LLR = log_pYgMd - log_pYgMs;

% Posterior probability
PP = exp(log_pYgMd) ./ (exp(log_pYgMd) + exp(log_pYgMs));

%% Display the inference process
%  =============================

% Define the different variables from the IO to display
Var = {PP, pAgYMd, IgYMd, HpRgY, sqrt(JSdiv)};
nVar = numel(Var);
VarLab = {{'Bayes factor', '$\frac{p(y|\mathcal{M_{\mathrm{D}}})}{p(y|\mathcal{M}_{\rm{S}})}$'}, ...
          {'Prediction', '$p(y_{k+1}=\mathrm{A}|y_{1:k},\mathcal{M_{\mathrm{D}}})$'}, ...
          {'Surprise', '$-\log_{2} p(y_{k+1}|\mathcal{M_{\mathrm{D}}})$'}, ...
          {'Entropy of the posterior', '$H(p(R|y,\mathcal{M_{\mathrm{D}}}))$'}, ...
          {'Model update', '$\sqrt{D_\mathrm{JS}(p(R|y_{1:k},\mathcal{M_{\mathrm{D}}})||p(R|y_{1:k-1},\mathcal{M_{\mathrm{D}}}))}$'}};

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
        y = Var{iVar}(~isinf(Var{iVar}));
        limy1 = [min(y), max(y)];
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
        plot(1:nObs, Var{iVar}(iSeq,:), '.-', 'Color', col(iSeq,:), ...
            'LineWidth', 1, 'MarkerSize', 8); hold('on');
        
        % Customize the axes
        axis([1/2, nObs+1/2, limy2]);
        set(gca, 'XTick', [1, 5:5:nObs]);
        set(gca, 'Color', 'None', 'LineWidth', 1, 'Layer', 'Top', ...
            'TickLabelInterpreter', 'LaTeX');

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
    sp = subplot(nVar+2, nSeq, iSeq + nSeq*nVar + [0,nSeq]);
    imagesc(1:nObs, 1:nu, pRgY(:,:,iSeq)); hold('on');
    
    % Display each time a pattern is repeated
    plot(repmat(rep', [1,2])', repmat([1/2, nObs+1/2], [numel(rep),1])', ...
        '-', 'Color', ones(1,3)./2, 'LineWidth', 1/2);
    
    % Customize the colormap
    colormap(viridis); caxis([min(pRgY(:)), max(pRgY(:))]);
    
    % Customize the axes
    set(gca, 'XTick', [1, 5:5:nObs]);
    set(gca, 'Color', 'None', 'LineWidth', 1, 'Layer', 'Top', ...
        'TickLabelInterpreter', 'LaTeX');
    
    % Add some text labels
    xlabel('Observation ($K$)', 'Interpreter', 'LaTeX');
    if iSeq == 1
        ylabel({'Posterior', 'distribution', '$p(R_i|y,\mathcal{M_{\mathrm{D}}})$'}, ...
            'Interpreter', 'LaTeX', 'Rotation', 0, ...
            'HorizontalAlignment', 'Right', 'VerticalAlignment', 'Middle');
    end
end
