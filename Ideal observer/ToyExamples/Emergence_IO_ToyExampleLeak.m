% A script showing the effect of the memory leak onto different variables
% reflecting the inference scenario of ideal observers learning either
% probabilistic or deterministic regularities.
%
% Copyright (c) 2018 Maxime Maheu

%% INITIALIZATION
%  ==============

% Clear the place
clear;
close('all');
 
% Add Emergence functions to the MATLAB path
scriptpath = mfilename('fullpath');
folderpath = scriptpath(1:strfind(scriptpath,'Emergence')+8);
addpath(genpath(folderpath));

% Set default figure properties
Emergence_DefaultFigureProperties;

%% DISPLAY THE EFFECT OF THE LEAKY WEIGHTS ON INFERRED VALUES OF THETA
%  ===================================================================

% Define grids and compute remembered values of theta depending on (1) the
% true value of theta and (2) the weight parameter
theta  = 0:0.01:1;
weight = 1/2:0.005:1;
thetastar = (theta'*weight) + ((1-theta')*(1-weight));

% Prepare a new window
figure('Units', 'Normalized', 'Position', [0.4 1/3 1/5 1/3]);

% Display the heat map
imagesc(weight, theta, thetastar); hold('on');

% Add a colorbar and customize the axes 
cbr = colorbar('Location', 'SouthOutside'); caxis([0,1]);
axis('xy'); axis('square');

% Add some text labels
cbr.Label.String = 'Remembered value of $\theta$';
cbr.Label.Interpreter = 'LaTeX'; 
xlabel('Weights $w = p(\hat{y}_{k} = {y}_{k})$', 'Interpreter', 'LaTeX');
ylabel('True value of $\theta$', 'Interpreter', 'LaTeX');

%% DEFINE A SET OF LEAKY PARAMETERS
%  ================================

% Array of forgetting parameters to test
pMemError = [0, 0.003, 0.009, 0.015, 0.025, 0.07, 0.185];
nParam = numel(pMemError);
col = lines(nParam); % choose different colors for each leak parameter

% Labels for each leak parameter
leakval = cellfun(@(x) sprintf('%1.3f', x), num2cell(pMemError), 'UniformOutput', 0);

%% CREATE 2 SEQUENCES WITH DIFFERENT REGULARITIES
%  ==============================================

% Length of sequences to consider
nObs = 48;

% Prepare an output variable
leak = NaN(nParam,nObs);

% Compute the corresponding leak, for each value of the error parameter
for iParam = 1:nParam, leak(iParam,:) = Emergence_IO_Leak(pMemError(iParam), nObs); end

% Example sequences
y = cell(1,2);
y{1} = GenRandSeq(nObs, [1/4, 1/4]);        % probabilistic regularity
y{2} = repmat(str2pat('AAAABBBB'), [1, 6]); % deterministic regularity

%% RUN BAYESIAN IDEAL OBSERVES WITH A LEAKY INTEGRATION
%  ====================================================

% Prepare outputs
pYgM  = NaN( nParam,nObs,2);
pAgYM = NaN( nParam,nObs,2);
IgYM  = NaN( nParam,nObs,2);
JSdiv = NaN( nParam,nObs,2);
pTgY  = cell(nParam,     2);

% Define properties of the ideal observers' inference
nu = 10; % depth of the tree
dt = 0.01; % precision of the grid

% Scale to use for the model evidences
scaleme = 'log';

% For each value of the leak
for iParam = 1:nParam
    
    % For each Bayesian observer
    for iObs = 1:2
        
        % For an observer learning probabilistic regularities
        if iObs == 1, IOfun = @Emergence_IO_Markov;
            inputs = {scaleme, true, 'Bayes-Laplace', leak(iParam,:), dt};
        
        % For an observer learning deterministic regularities
        elseif iObs == 2, IOfun = @Emergence_IO_Tree;
            inputs = {nu, scaleme, true, 'Size-principle', leak(iParam,:)};
        end
        
        % Run the observer iteratively after each observation of the
        % sequence
        [pYgM(iParam,:,iObs), pAgYM(iParam,:,iObs), ~, ~, IgYM(iParam,:,iObs), ...
         pTgY{iParam,iObs}, ~, ~, JSdiv(iParam,:,iObs)] = ...
            Emergence_IO_RunIO(IOfun, y{iObs}, inputs);
    end
end

% For the update, take the square root value of the Jensen-Shannon divergence
JSdiv = sqrt(JSdiv);

% Compare inference to a dumb observer hypothesizing that the sequence is
% purely random
if     strcmpi(scaleme, 'lin'), pYgMr = 1./2.^(1:nObs);
elseif strcmpi(scaleme, 'log'), pYgMr = -(1:nObs) * log(2);
end

% Compute posterior probabilities
if     strcmpi(scaleme, 'lin'), PP = pYgM ./ (pYgM + pYgMr);
elseif strcmpi(scaleme, 'log'), PP = exp(pYgM) ./ (exp(pYgM) + exp(pYgMr));
end

%% DISPLAY THE RESULT OF THE INFERENCE
%  ===================================

% Define the variables to look at
vars = {'PP', 'pAgYM', 'IgYM', 'JSdiv'};
VarLab = {{'Posterior probability', '$\frac{p(y|\mathcal{M}_{i})}{p(y|\mathcal{M}_{\rm{S}})}$'}, ...
          {'Prediction', '$p(y_{k+1}=\mathrm{A}|y_{1:k},\mathcal{M}_{i})$'}, ...
          {'Surprise', '$-\log_{2} p(y_{k+1}|\mathcal{M}_{i})$'}, ...
          {'Model update', '$\sqrt{D_\mathrm{JS}(p(\theta|y_{1:k},\mathcal{M}_{i}}$', ...
          '$\overline{||p(\theta|y_{1:k-1},\mathcal{M}_{i}))}$'}};
nVar = numel(VarLab);

% For each type of observer
for iObs = 1:2
    figure('Units', 'Normalized', 'Position', [0.2 0.5*(iObs-1) 0.6 0.4]);

    % For each memory parameter
    for iParam = 1:nParam
        
        % Leak function
        % ~~~~~~~~~~~~~
        
        % Display leaking functions implementing the imperfect memory
        subplot(nVar+3, nParam, iParam);
        lgd = plot(1:nObs, leak(iParam,:), '.-', 'Color', ...
            col(iParam,:), 'MarkerSize', 8); hold('on');
        
        % Customize the axes
        axis([1,nObs,1/2,1]);
        set(gca, 'Color', 'None', 'LineWidth', 1, 'Layer', 'Top', ...
            'TickLabelInterpreter', 'LaTeX');
        
        % Add some text labels
        if iParam == 1
            ylabel({'Weights $w_{k}$', '$p(\hat{y}_{k} = {y}_{k})$'}, ...
                'Interpreter', 'LaTeX', 'Rotation', 0, ...
                'HorizontalAlignment', 'Right', 'VerticalAlignment', 'Middle');
        end
        title(['$p_{\mathrm{error}} = ', leakval{iParam}, '$'], ...
            'Interpreter', 'LaTeX', 'Color', col(iParam,:));
        
        % Variables
        % ~~~~~~~~~~
        
        % For each variable from the ideal observer
        for iVar = 1:nVar
            subplot(nVar+3, nParam, iParam+nParam*iVar);
            
            % Define vertical limits of the plot
            tp = eval([vars{iVar}, '(iParam,:,iObs)']);
            tp = tp(~isinf(tp));
            limy1 = [min([0,min(tp)]), max([1,max(tp)])];
            margin = diff(limy1).*(1/4);
            limy2 = limy1 + [-1,1].*margin;
            
            % Display the sequence
            pos = (limy2 - limy1) / 2 + limy1;
            plot(find(y{iObs} == 2), pos(1), 'k.', 'MarkerSize', 6);
            plot(find(y{iObs} == 1), pos(2), 'k.', 'MarkerSize', 6);
            
            % Display the beliefs of the ideal observer
            plot(1:nObs, tp, '.-', 'Color', col(iParam,:), ...
                'LineWidth', 1, 'MarkerSize', 8); hold('on');
            
            % Customize the axes
            axis([1/2, nObs+1/2, limy2]);
            set(gca, 'XTick', [1, 10:10:nObs]);
            set(gca, 'Color', 'None', 'LineWidth', 1, 'Layer', 'Top', ...
                'TickLabelInterpreter', 'LaTeX');
            
            % Add some text labels
            if iParam == 1
                ylabel(VarLab{iVar}, 'Interpreter', ...
                    'LaTeX', 'Rotation', 0, 'HorizontalAlignment', 'Right', ...
                    'VerticalAlignment', 'Middle');
            end
        end
        
        % Posterior distribution over observer's parameters
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        % Display the posterior beliefs over patterns
        sp = subplot(nVar+3, nParam, iParam + nParam*(nVar+1) + [0,nParam]);
        tp = pTgY{iParam,iObs}; nTdim = size(tp,3);
        tp = reshape(permute(tp, [2,1,3]), [nObs, size(tp,1)*nTdim])';
        imagesc(1:nObs, [], tp); hold('on');
        
        % Customize the colormap
        colormap(parula); caxis([min(tp(:)), max(tp(:))]);
        
        % Customize the axes
        set(gca, 'XTick', [1, 10:10:nObs]); axis('xy');
        set(gca, 'Color', 'None', 'LineWidth', 1, 'Layer', 'Top', ...
            'TickLabelInterpreter', 'LaTeX');
        
        % Add some text labels
        xlabel('Observation ($K$)', 'Interpreter', 'LaTeX');
        if iParam == 1
            ylabel({'Posterior', 'distribution', '$p(\theta|y,\mathcal{M}_{i})$'}, ...
                'Interpreter', 'LaTeX', 'Rotation', 0, ...
                'HorizontalAlignment', 'Right', 'VerticalAlignment', 'Middle');
        end
    end
end
