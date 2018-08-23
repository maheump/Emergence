% This script displays the inference a n-order Markov chain on an example
% binary sequence.
% 
% Copyright (c) 2018 Maxime Maheu

%% Initialization
%  ==============

% Clear the place
clear;
close('all');

% Create a sequence
nObs = 200;
pAgB = 1/3;
pBgA = 1/3;
seq = GenRandSeq(nObs, [pAgB, pBgA]);

%% Run the n-order Markov chain
%  ============================

% Define properties of the n-order Markov chain
N = 3; % order of the chain
options = {'log', 'Bayes-Laplace', N, 0.01};

% Run the chain iteratively each time an observation is received
[pYgM, post, pred, surp, entr] = Emergence_IO_RunIO(@Emergence_IO_Chains, seq, options);

% Likelihood of the null model
pYgM0 = -(1:nObs)*log(2);

% Log-likelihood ratio
LLR = pYgM - pYgM0;

% Posterior belief
pMgY = exp(pYgM) ./ (exp(pYgM) + exp(pYgM0));

%% Display sequential inference
%  ============================

% Likelihood against a null model
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

figure('Units', 'Normalized', 'Position', [0.2 0.07 0.3 0.8]);

subplot(5,1,1);
plot([1,nObs], zeros(1,2), '-', 'Color', ones(1,3)/2); hold('on');
plot(1:nObs, LLR, 'k.-', 'MarkerSize', 10);
set(gca, 'Box', 'Off');
ylabel('Log likelihood ratio');

subplot(5,1,2);
plot([1,nObs], ones(1,2)/2, '-', 'Color', ones(1,3)/2); hold('on');
plot(1:nObs, pMgY, 'k.-', 'MarkerSize', 10);
set(gca, 'Box', 'Off', 'YLim', [0,1]);
ylabel('Posterior probability');

% Predictions
% ~~~~~~~~~~~

subplot(5,1,3);
plot([1,nObs], ones(1,2)/2, '-', 'Color', ones(1,3)/2); hold('on');
plot(find(seq == 1), 1.1, 'b.', 'MarkerSize', 10);
plot(find(seq == 2), -.1, 'r.', 'MarkerSize', 10); 
plot(1:nObs, pred, 'k.-', 'MarkerSize', 10); 
set(gca, 'Box', 'Off', 'YLim', [-0.2,1.2]);
ylabel('Prediction');

% Entropy
% ~~~~~~~

subplot(5,1,4);
plot(1:nObs, entr, 'k.-', 'MarkerSize', 10); 
set(gca, 'Box', 'Off', 'YLim', [0,1]);
xlabel('Observation #');
ylabel('Entropy');

% Surprise
% ~~~~~~~~

subplot(5,1,5);
plot([1,nObs], ones(1,2), '-', 'Color', ones(1,3)/2); hold('on');
plot(1:nObs, surp, 'k.-', 'MarkerSize', 10); 
set(gca, 'Box', 'Off');
xlabel('Observation #');
ylabel('Surprise');

% Posterior marginal distributions
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if N < 5
figure('Units', 'Normalized', 'Position', [0.5 0.07 0.3 0.8]);

pgrid = linspace(0, 1, size(post, 2));
trans = mat2cell(ff2n(N)+1, ones(2^N,1), N);
trans = cellfun(@(x) pat2str(x, {'A','B'}, 1:2), trans, 'UniformOutput', 0);
ntrans = 2^N;

for t = 1:ntrans
    subplot(ntrans, 1, t)
    imagesc(1:nObs, pgrid, squeeze(post(t,:,:)));
    caxis([0, max(post(:))]);
    axis('xy');
    if t ~= ntrans, set(gca, 'XTickLabel', {}); end
    set(gca, 'YTick', [0,1]);
    ylabel(['p(\theta_{A|', trans{t}, '}|y)'], 'Rotation', 0, ...
        'HorizontalAlignment', 'Right', 'VerticalAlignment', 'Middle');
end
ScaleAxis('c');
xlabel('Observation #');
end
