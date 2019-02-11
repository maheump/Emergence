% This script displays various quantities from the different ideal
% observers when faced with the same small input sequence that we use in
% the introduction of the paper (based on the repetition of the pattern
% AAB). In particular, we show how an ultimate observation that either
% violates or confirms the pattern is processed by the deterministic versus
% probabilistic observers.
% 
% Copyright (c) 2018 Maxime Maheu

%% INITIALIZATION
%  ==============

% Clear the place
clear;
close('all');

% Add functions to the MATLAB path
scriptpath = mfilename('fullpath');
ind = strfind(scriptpath,'Emergence');
folderpath = scriptpath(1:ind(end-1)+8);
addpath(genpath(folderpath));

% Set default figure properties
Emergence_DefaultFigureProperties;

% Define the pattern over which to loop
pat = 'AAB';
nrepet = 3;

% Create a sequence based on a repetition of a pattern
Y = repmat(str2pat(pat), [1, nrepet]);

% Remove the last observation because we want to compare effect of a
% violating versus confirming observation
Y = Y(1:end-1);

% Get the length of that sequence
N = numel(Y);

%% SIMULATION LOOP
%  ===============

% Define the observers to use
IOfun = {@Emergence_IO_Bernoulli, @Emergence_IO_Tree};

% Define properties of these observers
nu      = 5;     % longest possible pattern allowed
dt      = 0.01;  % precision of the posterior over theta
scaleme = 'log'; % scale for the sequence marginal likelihood
usegrid = false; % whether to use grid-based or analytical computations
inputs  = {{    scaleme, usegrid, 'Bayes-Laplace', [], dt}, ...
           {nu, scaleme, usegrid, 'Size-principle'      }};

% Prepare output variables
pYgMr    = NaN(2,N+1,2); % likelihood of the sequence
pAgYMr   = NaN(2,N+1,2); % expectation
PEpAgYMr = NaN(2,N+1,2); % prediction error
pTgYMr   = cell(2,2);    % posterior over models' parameters

% Whether the last observation is an A or a B
for lastit = 1:2
    seq = [Y, abs(lastit-3)];
    
    % For each type of observer
    for iMod = 1:2
        
        % Iteratively estimate different quantities from the observer
        [pYgMr(iMod,:,lastit), pAgYMr(iMod,:,lastit), ~, ...
         PEpAgYMr(iMod,:,lastit), ~, pTgYMr{iMod,lastit}] = ...
            Emergence_IO_RunIO(IOfun{iMod}, seq, inputs{iMod});
    end
    
    % If the probabilistic observer learns transition probabilities,
    % combine marginal posterior distributions by conditioning them on the
    % previous observation of the sequence
    if contains(func2str(IOfun{1}), 'Markov', 'IgnoreCase', true)
        pAgB = pTgYMr{1,lastit}(:,:,1);
        pBgA = pTgYMr{1,lastit}(:,:,2);
        pTgYMr{1,lastit} = margpost2condpost(pAgB, pBgA, seq);
    end
end

% Estimate quantities from a fully-stochastic process
[pYgMs, pAgYMs] = Emergence_IO_Null(1:N+1, scaleme);
IpAgYMs = NaN(size(pAgYMs));
IpAgYMs(Y == 1) = abs(Y(Y == 1) -      pAgYMs(Y == 1));
IpAgYMs(Y == 2) = abs(Y(Y == 2) - (1 - pAgYMs(Y == 2)));

% Compute posterior probabilities
if strcmpi(scaleme, 'log')
    pY = sum(exp(pYgMr),1) + exp(pYgMs);
    pMgY = exp(cat(1, pYgMr, repmat(pYgMs, [1 1 2]))) ./ pY;
elseif strcmpi(scaleme, 'lin')
    pY = sum(pYgMr,1) + pYgMs;
    pMgY = cat(1, pYgMr, repmat(pYgMs, [1 1 2])) ./ pY;
end

%% DISPLAY RESULTS
%  ===============

% Prepare a new window
figure('Position', [50 570 150 500]);

% Define useful variables
tricol = [049, 130, 189; 222, 045, 038; 049, 163, 084] ./ 255;
w = 10; % width of the first column of subplots

% Sequence
% ~~~~~~~~

% Sequence
subplot(6, w, 1:(w-2));
plot(1:N, Y, 'k-', 'LineWidth', 1/2); hold('on');
for i = unique(Y)
    plot(find(Y == i), i, 'ko', 'MarkerFaceColor', 'w');
end
set(gca, 'YDir', 'Reverse', 'XColor', 'None', 'Box', 'Off');
set(gca, 'YTick', 1:2, 'YTickLabel', {'A','B'});
axis([1/2, N+1/2, 0, 3]);

% Ultimate observation
for i = 1:2
    subplot(5, w, w-2+i);
    plot(N+1, abs(i-3), 'ko', 'MarkerFaceColor', 'w');
    axis([N+1-1/2, N+1+1/2, 0, 3]); axis('off');
    set(gca, 'XTick', [], 'YTick', [], 'YDir', 'Reverse');
end

% Posterior model probabilities
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Sequence
subplot(6, w, (1:(w-2)) + w);
l = bar(1:N, pMgY(:,1:N)', 1, 'stacked');
for j = 1:3, set(l(j), 'FaceColor', tricol(j,:)); end
axis([1/2, N+1/2, 0, 1]);
set(gca, 'XTickLabel', {});
ylabel('$p\left(\mathcal{M}|y\right)$', 'Interpreter', 'LaTeX');

% Ultimate observation
for i = 1:2
    subplot(6, w, w-2+i + w);
    l = bar([0,N+1], [zeros(3,1), pMgY(:,N+1,i)]', 1, 'stacked');
    for j = 1:3, set(l(j), 'FaceColor', tricol(j,:)); end
    axis([N+1-1/2, N+1+1/2, 0, 1]);
    set(gca, 'XTickLabel', {}, 'YTickLabel', {});
end

% Expectation and surprise
% ~~~~~~~~~~~~~~~~~~~~~~~~

% Separately for expectation and surprise
txtlab = {'$p\left(\mathrm{A}|y_{1:K-1}\mathcal{M}\right)$', 'PE'};
lab = {'', 'PE'};
for h = 1:2
    var = [lab{h}, 'pAgYMr'];
    
    % Sequence
    subplot(6, w, (1:(w-2)) + (1+h)*w);
    plot([0,N+1], ones(1,2)./2, 'k', 'Color', tricol(3,:)); hold('on');
    l = plot(1:N, eval([var,'(:,1:N,1)']), '.-', 'MarkerSize', 15, 'LineWidth', 2);
    for j = 1:2, set(l(j), 'Color', tricol(j,:)); end
    axis([1/2, N+1/2, 0, 1]);
    set(gca, 'Box', 'Off', 'Layer', 'Bottom');
    set(gca, 'XTick', 1:N, 'XTickLabel', {});
    ylabel(txtlab{h}, 'Interpreter', 'LaTeX');
    
    % Ultimate observation
    for i = 1:2
        subplot(6, w, w-2+i + (1+h)*w);
        plot([N-1,N+2], ones(1,2)./2, 'k-', 'Color', tricol(3,:)); hold('on');
        plot(N+1, eval([var,'(1,N+1,i)']), '.', 'Color', tricol(1,:), 'MarkerSize', 15);
        plot(N+1, eval([var,'(2,N+1,i)']), '.', 'Color', tricol(2,:), 'MarkerSize', 15);
        axis([N+1-1/2, N+1+1/2, 0, 1]);
        set(gca, 'Box', 'Off', 'Layer', 'Bottom');
        set(gca, 'XTickLabel', {}, 'YTickLabel', {}, 'YColor', 'None');
    end
end

% Posterior over models
% ~~~~~~~~~~~~~~~~~~~~~

% For both types of regular observers
col = {'Blues', 'Reds'};
for iMod = 1:2
    
    % Sequence
    sp = subplot(6, w, (1:(w-2)) + w*(iMod+3));
    if     iMod == 1, yt = 0:dt:1;
    elseif iMod == 2, yt = 1:nu;
    end
    imagesc(1:N, yt, pTgYMr{iMod,1}(:,1:N)); hold('on');
    plot(repmat(1/2:1:(N+1/2), 2, 1), repmat(ylim', 1, N+1), 'k-');
    set(gca, 'XTick', 1:N);
    if iMod == 1
        set(gca, 'XTickLabel', {}); axis('xy');
    elseif iMod == 2
        set(gca, 'YTick', 1:nu);
        xlabel('Observation (K)');
    end
    caxis([0,max(cellfun(@(x) max(x(:)), pTgYMr(iMod,:)))]);
    colormap(sp, cbrewer2(col{iMod}, 100));
    
    % Ultimate observation
    for i = 1:2
        sp = subplot(6, w, w-2+i + w*(iMod+3));
        imagesc(N+1, yt, pTgYMr{iMod,i}(:,N+1));
        set(gca, 'YTickLabel', {});
        if iMod == 1, set(gca, 'XTickLabel', {}); end
        if     iMod == 1, axis('xy');
        elseif iMod == 2, set(gca, 'YTick', 1:nu);
        end
        caxis([0,max(cellfun(@(x) max(x(:)), pTgYMr(iMod,:)))]);
        colormap(sp, cbrewer2(col{iMod}, 100));
    end
end
colorbar('Position', [3/4,1/20,1/6,1/50], 'Orientation', 'Horizontal');
