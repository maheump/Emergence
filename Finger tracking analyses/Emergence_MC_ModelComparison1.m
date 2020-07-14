% This script shows mean squared error between different alternaive models
% and the subjects measured either on the abruptness of the detection
% dynamics of deterministic rules (pattern-specific) or on the average
% probability of the deterministic rule hypothesis in the fully-random
% parts of sequences
% 
% Copyright (c) 2020 Maxime Maheu

% Initialization
% ~~~~~~~~~~~~~~

% On which measure to quantify the MSE
mes = 'Categorisation'; % 'Categorisation', 'Abruptness', or 'FalseAlarm'

% Get model simulations, i.e. choose among:
% - PredictabilityLin
% - PredictabilityMax
% - PredictabilityUshape
% - PredictabilityAshape
% - PseudoDeterministic
% - BiasedPseudoDeterministic
% - DifferentPriors
% - IndependentDiffDiscrete
% - IndependentDiffContinuous
% - IndependentRatioDiscrete
% - IndependentRatioContinuous
% - Leak
% - TreeDepth
SimuType = 'PredictabilityMax';

% Define properties of the plots based on the model to simulate
lab = {'Continuous', 'Shape', 'Leak'};
dispcont = any(contains(SimuType, lab, 'IgnoreCase', true));
useXlog = false;
useYlog = false;
lab = {'Ashape', 'Independent'};
if any(contains(SimuType, lab, 'IgnoreCase', true)), useXlog = true; end
lab = {'Predictability', 'PseudoDeterministic', 'DifferentPriors'};
if any(contains(SimuType, lab, 'IgnoreCase', true)), useYlog = true; end

% Run simulations
% ~~~~~~~~~~~~~~~
Emergence_MC_ModelSimulations;

% Compute abruptness of deterministic detection dynamics
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if strcmpi(mes, 'abruptness')
    
    % Get positions of the change point
    cp = cellfun(@(x) x.Jump, G(cidx{2},:), 'uni', 0);
    
    % Measure of abruptness between change and end points of deterministic
    % regularities
    abruptness = @(p) mean(abs(diff(p,2,1)));
    modbeh = cellfun(@(p,c) abruptness(p(c:end,2)), pMgY(cidx{2},:,:), repmat(cp, [1 1 nMod]));
    subbeh = cellfun(@(p,c) abruptness(p.BarycCoord(c:end,2)), G(cidx{2},:), cp);
    
    % Remove deterministic regularities that were not identified by
    % subjects
    detecmask = (filter{2} == 3);
    subbeh(~detecmask) = NaN;
    modbeh(~detecmask) = NaN;
    
    % Compute the error between models and subjects
    MSE = (subbeh - modbeh) .^ 2;
    
    % Average over sequences
    MSE = squeeze(mean(MSE, 1, 'OmitNaN'));
    
% Compute average false alarm position
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
elseif strcmpi(mes, 'falsealarm')
    
    % Select fully-stochastic parts of sequences
    idxtrimap = Emergence_SelectFullyStochSeq(G, filter, 1);
    
    % Get average positions in the (pseudo-)deterministic hypothesis from
    % the models
    modbeh = cellfun(@(x,i) x(i,:), pMgY, repmat(idxtrimap, [1,1,nMod]), 'uni', 0);
    modbeh = squeeze(mat2cell(modbeh, nSeq, ones(1,nSub,1), ones(1,1,nMod)));
    modbeh = cellfun(@(x) mean(cell2mat(x)), modbeh, 'uni', 0);
    
    % Get average positions in the (pseudo-)deterministic hypothesis from
    % the subjects
    subbeh = cellfun(@(x,i) x.BarycCoord(i,:), G, idxtrimap, 'uni', 0);
    subbeh = squeeze(mat2cell(subbeh, nSeq, ones(1,nSub)))';
    subbeh = cellfun(@(x) mean(cell2mat(x)), subbeh, 'uni', 0);
    
    % Tranform to cartesian coordinates to 
    subbeh = cellfun(@(x) x*tricc, subbeh, 'uni', 0);
    modbeh = cellfun(@(x) x*tricc, modbeh, 'uni', 0);
    
    % Compute the error between models and subjects
    MSE = cellfun(@(x,y) sum((x-y).^2), repmat(subbeh, [1,nMod]), modbeh);
    
% Compute categorisation profiles
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
elseif strcmpi(mes, 'categorisation')
    
    % Define error table wich relates subjects (columns) x model (rows)
    % responses
    %           P   D   S
    errtbl = [  0	1   1	;... % P
                1	0   1	;... % D
                1	1   0	;... % S
                1	1/2 1/2	;... % D/S
                1/2	1   1/2	;... % P/S
                1/2	1/2 1	;... % P/D
                1/3 1/3 1/3 ];   % P/D/S
    
    % Get models' categorisation profiles
    modbeh = cellfun(@(p) p(end,:), pMgY, 'uni', 0); % 1=S/2=P/3=D
    modbeh = cellfun(@(x) Emergence_GetModelCategorisation(x, 0.001), modbeh);
    
    % Get subjects' categorisation profiles
    subbeh = cellfun(@(x) x.Questions(2), G);
    subbeh(isnan(subbeh)) = 3; % 1=S/2=P/3=D
    subbeh = repmat(subbeh, [1,1,nMod]);
    
    % Get error for each sequence and each subject
    MSE = arrayfun(@(m,s) errtbl(m,s), modbeh, subbeh);
    
    % Average over sequences
    MSE = squeeze(mean(MSE, 1));
end

% Center the MSE from the alternative models relatively to the normative
% two-system model
MSE = MSE(:,1:end-1) - MSE(:,end);

% Perform paired t-tests
[h,pval,tci,stats] = ttest(MSE);
Emergence_PrintTstats(pval,tci,stats);

% In case of models with 2 dimensions
nCrv = 1;
if dispcont && contains(SimuType, 'Predictability', 'IgnoreCase', true)
   nCrv = numel(unique(options(4,1:end-1)));
   nMod = (nMod-1) / nCrv;
   MSE = reshape(MSE, [nSub,nCrv,nMod]);
end

% Display the mean squared error
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Average MSE over subjects
m = squeeze(mean(MSE));
s = squeeze(sem( MSE));

% Prepare a new window
figure('Position', [1 905 300 200]); hold('on');

% If continuous models, use lines
if dispcont
    x = contparam';
    for i = 1:nCrv
        plotMSEM(x, m(i,:), s(i,:), 0.15, modc(i,:), modc(i,:), 2);
    end
    xlim(sort(x([1,end]))); xlabel('Parameters');
    
% If discrete models, use dots instead
elseif ~dispcont
    x = 1:size(MSE,2);
    plot(x, m, 'k-', 'LineWidth', 2);
    plot(repmat(x, [2,1]), m+s.*[-1;1], 'k-');
    scatter(x, m, 100, modc(x,:), 'Filled', 'MarkerEdgeColor', 'k');
    xlim([0,nMod]); xlabel('Models');
end

% Overlap p-values on top of the plot
if nCrv == 1
    lima = [0.001 0.01 0.05];
    sz = [6,4,2];
    for i = 1:numel(lima)
        idx = pval < lima(i);
        limy = ylim;
        plot(x(idx), repmat(max(limy), [1,sum(idx)]), 'ks', ...
           'MarkerEdgeColor', 'None', ...
           'MarkerFaceColor', 'k', 'MarkerSize', sz(i));
       ylim(limy);
    end
end

% Customize the axes
set(gca, 'Box', 'Off', 'Layer', 'Top');
lab = {'Continuous', 'Predictability'};
if dispcont && useXlog, set(gca, 'XScale', 'log'); end
if useYlog, set(gca, 'YScale', 'log'); end

% Add some text labels
ylabel('Mean squared error');
title(sprintf('%s - %s', SimuType, mes));
