% A script showing the effect of the memory leak onto different variables
% reflecting the inference scenario of ideal observers learning either
% probabilistic or deterministic regularities.
%
% Copyright (c) 2018 Maxime Maheu

%% Initialization
%  ==============

% Clear the place
clear;
close('all');

% Array of forgetting parameters to test
pMemError = [0, 0.003, 0.009, 0.025, 0.07, 0.185, 0.5];
nParam = numel(pMemError);

% Labels for each leak parameter
leakval = cellfun(@(x) sprintf('%1.3f', x), ...
    num2cell(pMemError), 'UniformOutput', 0);

% Length of sequences to consider
nObs = 50;

%% Run ideal observers using a leaky integration
%  =============================================

% Prepare an output variable
leak = NaN(nParam,nObs);

% Compute the corresponding leak, for each value of the error parameter
for iParam = 1:nParam, leak(iParam,:) = Emergence_IO_Leak(pMemError(iParam), nObs); end

% Example sequences
nSeq = 2;
y = cell(1,nSeq);
y{1} = GenRandSeq(nObs, [1/3, 2/3]); % probabilistic regularity
y{2} = repmat([1 1 2 2 1], [1, 10]); % deterministic regularity

% Prepare outputs
seqLH  = NaN(nParam,nObs,2);
predA  = NaN(nParam,nObs,2);
surp   = NaN(nParam,nObs,2);
entr   = NaN(nParam,nObs,2);
update = NaN(nParam,nObs,2);
post   = cell(1,2);

% Define properties of the ideal observers' inference
nu = nObs; % depth of the tree
dt = 0.01; % precision of the grid

% For each value of the leak
for iParam = 1:nParam
    
    % For each Bayesian observer
    for iObs = 1:2
        
        % For an observer learning probabilistic regularities
        if iObs == 1, IOfun = @Emergence_IO_Markov;
            inputs = {'log', true,  'Bayes-Laplace', leak(iParam,:), dt};
        
        % For an observer learning deterministic regularities
        elseif iObs == 2, IOfun = @Emergence_IO_Tree;
            inputs = {nu, 'log', false, 'Size-principle', leak(iParam,:)};
        end
        
        % Run the observer iteratively after each observation of the
        % sequence
        [seqLH(iParam,:,iObs), post{iObs}(:,:,:,iParam), predA(iParam,:,iObs), ...
            surp(iParam,:,iObs), ~, update(iParam,:,iObs)] = ...
            Emergence_IO_RunMi(IOfun, y{iObs}, inputs);
    end
end

% Get rid of non-relevant dimensions
post = cellfun(@(x) squeeze(x), post, 'UniformOutput', 0);

% Compare inference to a dumb observer hypothesizing that the sequence is
% purely random
batch = log((1/2).^(1:nObs));

%% Display the different weighting functions
%  =========================================

% Prepare window
figure('Units', 'Normalized', 'Position', [0.05 0.5 0.3 0.4]);

% Display observations
g = repmat(2/3,1,3);
plot(repmat(1:nObs,2,1)', [0,1], '-', 'Color', g); hold('on');

% Display leaking functions implementing the imperfect memory
col = winter(nParam);
lgd = plot(1:nObs, leak, '.-', 'MarkerSize', 15);
for iParam = 1:nParam, set(lgd(iParam), 'Color', col(iParam,:)); end

% Display chance level
plot([1,nObs], ones(1,2)/2, 'k--', 'LineWidth', 1);

% Customize the axes
set(gca, 'XLim', [1,nObs], 'YLim', [0,1]);
set(gca, 'Box', 'Off', 'LineWidth', 1, 'TickLabelInterpreter', 'LaTeX');

% Add some text labels
xlabel('Observations ($k$)', 'Interpreter', 'LaTeX');
ylabel({'Weights $w_{k}$', '$p\left(\hat{y}_{k} = {y}_{k}\right)$'}, ...
    'Interpreter', 'LaTeX', 'Rotation', 0, 'HorizontalAlignment', 'Right');
legend(lgd, leakval, 'Interpreter', 'LaTeX', 'Orientation', ...
    'Horizontal', 'Location', 'SouthOutside');

%% Display the effect of the leaky integration of the dynamics of the inference
%  ============================================================================

% Properties of the figure
fs = 12;
col = lines(1);

% For each type of observer
for iObs = 1:2
    figure('Units', 'Normalized', 'Position', [0.35 0.5*(iObs-1) 0.6 0.4]);
    
    % For each memory parameter
    for iParam = 1:nParam

        % Likelihood ratio
        % ~~~~~~~~~~~~~~~~
        
        subplot(6,nParam,iParam); hold('on');
        LLR = seqLH(iParam,:,iObs) - batch;
        mLLR = max(max(abs(seqLH(:,:,iObs))));
        plot(repmat(1:nObs,2,1)', [-mLLR,mLLR], '-', 'Color', g);
        plot([1,nObs], zeros(1,2), 'k--');
        plot(1:nObs, LLR, '-', 'Color', col, 'LineWidth', 2);
        set(gca, 'XLim', [1,nObs], 'YLim', [-mLLR,mLLR]);
        set(gca, 'FontSize', fs, 'LineWidth', 1, 'TickLabelInterpreter', 'LaTeX');
        if iParam == 1
            ylabel({'Bayes factor', ['$\frac{p(y|\mathcal{M})}', ...
                '{p(y|\mathcal{M}_{0})}$']}, 'Interpreter', 'LaTeX', 'Rotation', 0, ...
                'HorizontalAlignment', 'Right', 'VerticalAlignment', 'Middle');
        end
        title(leakval{iParam}, 'Interpreter', 'LaTeX');

        % Prediction
        % ~~~~~~~~~~
        
        subplot(6,nParam,nParam+iParam); hold('on');
        plot(repmat(1:nObs,2,1)', [0,1], '-', 'Color', g);
        plot([1,nObs], ones(1,2)/2, 'k--');
        plot(1:nObs, predA(iParam,:,iObs), '-', 'Color', col, 'LineWidth', 2);
        set(gca, 'XLim', [1,nObs], 'YLim', [0,1]);
        set(gca, 'FontSize', fs, 'LineWidth', 1, 'TickLabelInterpreter', 'LaTeX');
        if iParam == 1
            ylabel({'Prediction', '$p(y_{k+1}=\rm{A}\it{|}y_{1:k})$'}, ...
                'Interpreter', 'LaTeX', 'Rotation', 0, ...
                'HorizontalAlignment', 'Right', 'VerticalAlignment', 'Middle');
        end
        
        % Surprise
        % ~~~~~~~~
        
        subplot(6,nParam,(nParam*2)+iParam); hold('on');
        plot(repmat(1:nObs,2,1)', [0,6], '-', 'Color', g);
        plot([1,nObs], ones(1,2), 'k--');
        plot(1:nObs, surp(iParam,:,iObs), '-', 'Color', col, 'LineWidth', 2);
        set(gca, 'XLim', [1,nObs], 'YLim', [0,6]);
        set(gca, 'FontSize', fs, 'LineWidth', 1, 'TickLabelInterpreter', 'LaTeX');
        if iParam == 1
            ylabel({'Surprise', '$-\log_{2} p(y_{k+1})$'}, ...
                'Interpreter', 'LaTeX', 'Rotation', 0, ...
                'HorizontalAlignment', 'Right', 'VerticalAlignment', 'Middle');
        end
        
        % Model update
        % ~~~~~~~~~~~~
        
        subplot(6,nParam,(nParam*3)+iParam); hold('on');
        plot(repmat(1:nObs,2,1)', [0,1], '-', 'Color', g);
        plot(1:nObs, sqrt(update(iParam,:,iObs)), '-', 'Color', col, 'LineWidth', 2);
        set(gca, 'XLim', [1,nObs], 'YLim', [0,1]);
        set(gca, 'FontSize', fs, 'LineWidth', 1, 'TickLabelInterpreter', 'LaTeX');
        if iParam == 1
            ylabel({'Model', 'update', '$\sqrt{\rm{JS}(\it{p}(\theta|y_{1:k})}$', ...
                '$\overline{||p(\theta|y_{1:k-1}))}$'}, 'Interpreter', 'LaTeX', 'Rotation', 0, ...
                'HorizontalAlignment', 'Right', 'VerticalAlignment', 'Middle');
        end
        
        % Posterior distribution
        % ~~~~~~~~~~~~~~~~~~~~~~
        subplot(6,nParam,(nParam*(4:5))+iParam); hold('on');
        
        % Customize the axes
        set(gca, 'FontSize', fs, 'LineWidth', 1, 'TickLabelInterpreter', 'LaTeX');
        
        % Model learning probabilistic regularities
        if iObs == 1
            if ndims(post{1}) == 4
                imagesc(0:dt:1, 0:dt:1, post{iObs}(:,:,end,iParam));
                plot([0,1], [0,1], 'k-');
                plot([0,1], [1,0], 'k-');
                plot([0,1], ones(1,2)/2, 'k--');
                plot(ones(1,2)/2, [0,1], 'k--');
                axis([0,1,0,1]);
                axis('xy'); axis('square'); 
                if iParam == 1
                    xlabel('$\theta_{A|B}$', 'Interpreter', 'latex');
                    ylabel('$\theta_{B|B}$', 'Interpreter', 'latex');
                end
            elseif ndims(post{1}) == 3
                imagesc(1:nObs, 0:dt:1, post{iObs}(:,:,iParam));
                plot([1,nObs], ones(1,2)./2, 'k--');
                axis([1,nObs,0,1]);
                if iParam == 1
                    xlabel('Observation \#', 'Interpreter', 'latex');
                    ylabel({'Posterior', 'distribution', '$\theta_{A}$'}, ...
                        'Interpreter', 'latex', 'Rotation', 0, ...
                        'HorizontalAlignment', 'Right', 'VerticalAlignment', 'Middle');
                end
            end
            
        % Model learning deterministic regularities
        elseif iObs == 2
            plot([0,nu+1], repmat(1/nu,1,2), 'k--');
            bar(1:nu, post{iObs}(:,end,iParam), 1, 'FaceColor', col);
            set(gca, 'XLim', [1/2,nu+1/2], 'YLim', [0,1]);
            xlabel('$R_{i}$', 'Interpreter', 'LaTeX');
            if iParam == 1
                ylabel({'Posterior', 'distribution', '$p(R_{i}|y_{1:N})$'}, ...
                    'Interpreter', 'LaTeX', 'Rotation', 0, ...
                    'HorizontalAlignment', 'Right', 'VerticalAlignment', 'Middle');
            end
        end
    end
end
