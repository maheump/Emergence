% To prevent averaging effects to mask effects or lead to incorrect
% conclusions, individual trials are inspected. Sigmoid functions are then
% fitted to those individual trials. Resulting fitted slope parameters are
% then averaged over the different sequences in order to obtain an average
% detection scenario for each type of regularity. This confimed that the
% dynamics of detection is much more ballistic in case of deterministic
% compared to probabilistic regularities.
% 
% Copyright (c) 2018 Maxime Maheu

%% DISPLAY INDIVIDUAL TRIALS LIKELIHOODS IN THE RELEVANT HYPOTHESIS
%  ===============================================================

% Prepare output variable
cp          = cell(1,2);
sortedcp    = cell(1,2);
idxcp       = cell(1,2);
belincorhyp = cell(1,2);

% For each type of regularity
for iHyp = 1:2
    
    % Get sequences that were correctly labelled
    detecmask = (filter{iHyp} == 1 | filter{iHyp} == 3);
    
    % Order according to change point's position and detected/undetected
    cp{iHyp} = cellfun(@(x) x.Jump, G(cidx{iHyp},:), 'UniformOutput', 1);
    [sortedcp{iHyp}, idxcp{iHyp}] = sortrows([cp{iHyp}(:), ...
        detecmask(:)], [2,1]);
    
    % Get the beliefs in the corresponding (correct) hypothesis ordered
    % according to the position the change point
    belincorhyp{iHyp} = cellfun(@(x) x.BarycCoord(:,iHyp)', ...
        D(cidx{iHyp},:), 'UniformOutput', 0);
end

% Prepare a new window
figure('Position', [702 705 230 400]);
cmapcol = {'Blues', 'Reds'};

% For sequences with a probabilistic/deterministic regularity
for iHyp = 1:2
    
    % Display the change in beliefs as a heatmap 
    sp = subplot(2,1,iHyp);
    bel = belincorhyp{iHyp}(idxcp{iHyp});
    imagesc(cell2mat(bel)); hold('on');
    
    % Customize the colormap
    colorbar('Location', 'EastOutside'); caxis([0,1]);
    colormap(sp, cbrewer2(cmapcol{iHyp}));
    
    % Display the position of the change points
    for d = [0,1]
        x = cp{iHyp}(idxcp{iHyp}(sortedcp{iHyp}(:,2) == d));
        y = find(sortedcp{iHyp}(:,2) == d);
        stairs([x(1);x], [y(1)-1;y], 'k-');
    end
	
    % Display limits between detected and undetected sequences
    lim = find(abs(diff(sortedcp{iHyp}(:,2))) == 1) + 1/2;
    plot([0,N+1], repmat(lim, 1, 2), 'k-');
    
    % Add some text labels
    axis('xy'); set(gca, 'XTick', [1, get(gca, 'XTick')], 'YTick', [1,50]);
    xlabel('Observation #');
    ylabel({'Sequence # (sorted by', 'change point''s position)'});
    title(sprintf('%s sequences', proclab{iHyp}));
end

% Save the figure
if isfield(D{1}, 'Seq'), save2pdf(fullfile(ftapath, 'figs', 'F_Dyn_MapS.pdf'));
else, save2pdf(fullfile(ftapath, 'figs', 'F_Dyn_MapIO.pdf'));
end

%% FIT SIGMOID FUNCTIONS TO INDIVIDUAL TRIAL LIKELIHOODS IN THE RELEVANT HYPOTHESIS
%  ================================================================================

% Fit a sigmoid function to beliefs recorded for each individual trial 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Prepare output variables
p_sub = cell(1,2); p_gp = cell(1,2); % posterior over parameters
o_sub = cell(1,2); o_gp = cell(1,2); % quality of fit

% For each sequence that entailed a regularity
for iHyp = 1:2
    
    % Get beliefs in the relevant hypothesis post-change point
	cp = cellfun(@(x) x.Jump+1/2, G(cidx{iHyp},:), 'UniformOutput', 0);
    Y = cellfun(@(x,c) x.BarycCoord(c:end,iHyp), D(cidx{iHyp},:), cp, 'UniformOutput', 0);
    
    % Select sequences that were correctly identified by subjects AND for
    % which we were able to find a detection point for subjects and the
    % ideal observer
    detecmask = (filter{iHyp} == 3);
    
    % Separately for each subject
    for iSub = 1:nSub
        
        % Get data over the different sequences
        y = Y(:,iSub);
        idx = detecmask(:,iSub);
        y = y(idx);
        nseq = numel(y);
        
        % Specify the dimensions of the problem
        dim         = []; % empty variable
        dim.n_phi   = 4;  % number of observation parameters
        dim.n       = 0;  % number of hidden states
        dim.n_theta = 0;  % number of evolution parameters
        dim.n_t     = 1;  % number of time samples
        
        % Specify options
        options             = [];
        options.DisplayWin  = 0;
        options.verbose     = 0;
        options             = repmat({options}, [nseq,1]);
        optiongp            = [];
        optiongp.DisplayWin = 0;
        optiongp.verbose    = 0;
        
        % Explaining variable: position of the observation post-change-point
        for iseq = 1:nseq, options{iseq}.inG.p = (1:numel(y{iseq})); end
        
        % Run the mixed-effect fitting scheme (it deduces priors through
        % empirical Bayes) using a variational procedure (it uses Laplace
        % approximation) in order to estimate the best parameters for each
        % model, each condition and each subject
        [p_sub{iHyp}(idx,iSub), o_sub{iHyp}(idx,iSub), ...
         p_gp{iHyp}(idx,iSub),  o_gp{iHyp}(idx,iSub)] = ...
            VBA_MFX(y, [], [], @g_SIGM, dim, options, [], optiongp);
    end
end

% Get back the fitted parameters of the sigmoid
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% For each type of regularity
sigparam = cell(1,2);
for iHyp = 1:2
    
    % For non-identified sequences, fill with NaNs
    idx = find(cellfun(@isempty, p_sub{iHyp}));
    for d = 1:numel(idx)
        p_sub{iHyp}{idx(d)}.muPhi = NaN(4,1);
    end
    
    % Get the (4) fitted parameters of the sigmoid function
    sigparam{iHyp} = cellfun(@(x) reshape(x.muPhi, [1,1,4]), ...
        p_sub{iHyp}, 'UniformOutput', 0);
    sigparam{iHyp} = cell2mat(sigparam{iHyp});
end

% Display average sigmoid function (based on averaged fitted parameters)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Prepare a new window
figure('Position', [933 905 220 200]); lgd = NaN(1,2); hold('on');
xcp = 1:80;

% Draw some help lines
plot(xcp([1,end]),    ones(1,2)./2, '-',  'Color', g, 'LineWidth', 1/2); 
plot(xcp([1,end]),    ones(1,2)./3, '--', 'Color', g, 'LineWidth', 1/2);
plot(xcp([1,end]), 2.*ones(1,2)./3, '--', 'Color', g, 'LineWidth', 1/2);

% For each type of regularity
for iHyp = 1:2
    
    % Average sigmoid parameters over sequences
    param = squeeze(mean(sigparam{iHyp}, 1, 'OmitNaN'));
    avgparam = mean(param, 1, 'OmitNaN');
    semparam = sem(param, 1);
    
    % Draw error bars
    upperlim = g_SIGM([], avgparam-semparam, [], xcp);
    lowerlim = g_SIGM([], avgparam+semparam, [], xcp);    
    fill([xcp, fliplr(xcp)], [upperlim', fliplr(lowerlim')], 'k', ...
        'FaceColor', tricol(iHyp,:), 'EdgeColor', 'None', 'FaceAlpha', 0.15);
    
    % Draw the average sigmoid curve
    avgpred = g_SIGM([], avgparam, [], xcp);
    lgd(iHyp) = plot(xcp, avgpred, 'Color', tricol(iHyp,:), 'LineWidth', 3); 
end

% Customize the axes
set(gca, 'Box', 'Off');
axis([xcp([1,end]), 0, 1.001]);

% Add some text labels
xlabel('# observation w.r.t. change point');
ylabel('Fitted posterior belief p(M_i|y)');
legend(lgd, proclab, 'Location', 'SouthEast');

% Save the figure
if isfield(D{1}, 'Seq'), save2pdf(fullfile(ftapath, 'figs', 'F_Dyn_SigmS.pdf'));
else, save2pdf(fullfile(ftapath, 'figs', 'F_Dyn_SigmIO.pdf'));
end

% Compare the slope of the sigmoid in the two types of regularity
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Average the slope parameter over sequences for each subject and each type
% of regularity (probabilistic and deterministic ones)
slope = cellfun(@(x) exp(x(:,:,1)), sigparam, 'UniformOutput', 0);
avgslope = cell2mat(cellfun(@(x) mean(x, 1, 'OmitNaN')', slope, 'UniformOutput', 0));

% Paired t-test on the slope parameter of sigmoids
[~,pval,tci,stats] = ttest(diff(avgslope, 1, 2));
Emergence_PrintTstats(pval,tci,stats);

% Prepare a new window
figure('Position', [1154 905 120 200]);

% Display the difference of the sigmoids' slope in the two types of
% regularity
Emergence_PlotSubGp(avgslope, tricol(1:2,:));

% Customize the axes
set(gca, 'XTick', [], 'XColor', 'None', 'Box', 'Off');
axis([0,3,0,1.8]);

% Display whether the difference is significant or not
Emergence_DispStatTest(avgslope);

% Add some text labels
ylabel('Slope parameter');

% Save the figure
if isfield(D{1}, 'Seq'), save2pdf(fullfile(ftapath, 'figs', 'F_Dyn_SlopeS.pdf'));
else, save2pdf(fullfile(ftapath, 'figs', 'F_Dyn_SlopeIO.pdf'));
end

% The sigmoid function
% ~~~~~~~~~~~~~~~~~~~~
function [g, dgdx, dgdP] = g_SIGM(x, P, u, in)

% Make sure that the values at which to evaluate the function are within a structure
if ~isstruct(in)
    p = in;
    in = [];
    in.p = p;
end

% Build up the relationship between x and y
slope = exp(P(1)); % E [0,+inf[
intcp = P(2); % E ]-inf,+inf[
y = slope .* in.p + intcp; % linear relationship

% Compute predictions of the sigmoid functions
g = 1 ./ (1 + exp(-y));

% If provided, allows the sigmoid to start above zero
if numel(P) > 2
    lowerlim = 1 / (1 + exp(-P(3))); % E [0,1]
    g = g .* lowerlim;
    
    % If provided, allows the sigmoid to end below one
    if numel(P) > 3
        upperlim = 1 / (1 + exp(-P(4))); % E [0,1]
        g = g + upperlim;
    end
end

% Make sure the output is a column vector
g = g(:);

% Rely on numerical derivatives
dgdx = [];
dgdP = [];

end
