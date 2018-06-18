% This script looks at the dynamics of finger trajectory aroung change and
% detection points. We show that depending if we lock on the position of
% the change point or the detection point (when the beliefs in the true
% generative process crosses some threshold), the dynamics does not look
% the same even though detection dynamics look steeper in the case of
% deterministic regularities. To prevent averaging effects to mask effects
% or lead to incorrect conclusions, sigmoid functions are then fitted to
% individual trials. Resulting fitted parameters are then averaged over the
% different sequences in order to obtain an average detection scenario for
% each type of regularity. This confimed that the dynamics of detection is
% much more ballistic in case of deterministic compared to probabilistic
% regularities.
%
% Copyright (c) 2018 Maxime Maheu

%% BARYCENTRIC COORDINATES LOCKED ON DETECTION/CHANGE POINT
%  ========================================================

% Define properties of windows to look into
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% The number of observation to consider
nSamp = 80;

% Define x-axes for trajectories locked on change point
xcp = 0:nSamp; % x vector

% Define x-axes for trajectories locked on detection point
xdp = -nSamp/2:nSamp/2; 

% Threshold for detection on the relevant dimension
detecthr = 1/2;

% Get the trajectories around the points of interest
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Prepare the output variable
cp = cell(1,2); dp = cell(1,2); lag = cell(1,2);
fingerwrtp = cell(2);

% For each type of regularity
for iHyp = 1:2
    
    % Sample (in the entire sequence) around which the change happened
    cp{iHyp} = cellfun(@(x) x.Jump, G(cidx{iHyp},:)) - 1/2;
    
    % Get the position of the subjective detection point
    lag{iHyp} = cell2mat(cellfun(@(x,y) nanmin([NaN, ...
        find(x.BarycCoord( x.Jump - 1/2 : end, iHyp) ...
        > detecthr, 1, 'first') - 1]), ... % 
        D(cidx{iHyp},:), G(cidx{iHyp},:), 'UniformOutput', 0));
    dp{iHyp} = cp{iHyp} + lag{iHyp}; % detection point = change point + lag
    
    % Consider only sequences for which regularities were correctly identified
    detecmask = cellfun(@(x) x.Questions(2) == iHyp, G(cidx{iHyp},:));
    lag{iHyp}(~detecmask) = NaN;
    cp{iHyp}(~detecmask) = NaN;
    dp{iHyp}(~detecmask) = NaN;
    
    % For change- and detection- points
    for lock = 1:2
    
        % Beginning and ending of the window to look in
        if lock == 1 % lock on change point
            begwin = cp{iHyp} + xcp(1); % in number of samples
            endwin = cp{iHyp} + xcp(end); % in number of samples
        elseif lock == 2 % lock on detection point
            begwin = dp{iHyp} + xdp(1); % in number of samples
            endwin = dp{iHyp} + xdp(end); % in number of samples
        end
        endwin(endwin > 200) = 200;
        
        % Get trajectory (i.e. beliefs in each hypothesis) in that window
        % of interest
        fing = cellfun(@(x,b,e) x.BarycCoord(nanmin([200,b]):nanmin([199,e]),:), ...
            D(cidx{iHyp},:), num2cell(begwin), num2cell(endwin), 'UniformOutput', 0);
        fing(cellfun(@isempty, fing)) = {NaN(nSamp+1,3)};
        fing = cellfun(@(x) [x; NaN(nSamp-size(x,1)+1,3)], fing, 'UniformOutput', 0);
        fingerwrtp{iHyp,lock} = cell2mat(reshape(fing, [1, 1, size(fing)]));
    end
end

% Average over sequences for each type of regularity
avgfingerwrtp = cellfun(@(x) squeeze(mean(x, 3, 'OmitNaN')), fingerwrtp, 'UniformOutput', 0);
avglag = cellfun(@(x) mean(x, 'OmitNaN'), lag, 'UniformOutput', 0);

% Average finger trajectory over subjects
avgsubtraj = cellfun(@(x) mean(x, 3), avgfingerwrtp, 'UniformOutput', 0);
semsubtraj = cellfun(@(x) sem(x ,3),  avgfingerwrtp, 'UniformOutput', 0);
avgsublag = cellfun(@mean, avglag);
semsublag = cellfun(@sem, avglag);

% Display beliefs in each hypothesis locked to both change- and detection- points
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Prepare a new window
figure('Position', [202 706 500 400]);

% For change- and detection- points
pt = {'change', 'detection'};
for lock = 1:2
    if     lock == 1, x = xcp; % change point
    elseif lock == 2, x = xdp; % detection point
    end
    
    % For each type of regularity
    for iHyp = 1:2
        subplot(2,2,iHyp+2*(lock-1)); hold('on');
        
        % Draw some help lines
        plot(x([1,end]), ones(1,2)./2,    '-',  'Color', g, 'LineWidth', 1/2);
        plot(x([1,end]), ones(1,2)./3,    '--', 'Color', g, 'LineWidth', 1/2);
        plot(x([1,end]), 2.*ones(1,2)./3, '--', 'Color', g, 'LineWidth', 1/2);
        if lock == 2, plot(zeros(1,2), [0,1], 'k'); end
        
        % Draw barycentric coordinates
        lt = repmat({'--'}, 1, 3);
        lt{iHyp} = '-';
        for iDim = 1:3
            plotMSEM(x, avgsubtraj{iHyp,lock}(:,iDim), ...
                        semsubtraj{iHyp,lock}(:,iDim), ...
                0.15, tricol(iDim,:), tricol(iDim,:), 1+1*(iDim == iHyp), 1, lt{iDim}, 'none');
        end
        
        % Customize the axes
        axis([x([1,end]),0,1]);
        set(gca, 'Box', 'Off');
        
        % Add some text labels
        xlabel(sprintf('Observation w.r.t. %s point', pt{lock}));
        ylabel('Posterior beliefs p(M_i|y)');
        if lock == 1, title(sprintf('%s regularities', proclab{iHyp})); end
        
        % Display the distribution of average detection points
        if lock == 1 % locked on change point
            ax = get(gca, 'Position');
            axes('Position', [ax(1) ax(2)+0.8*ax(4) ax(3) 0.2*ax(4)]);
            Emergence_PlotSubGp(avglag{iHyp}, tricol(iHyp,:));
            axis([0, 2, x([1,end])]); axis('off');
            view([90,90]);
            Emergence_DispStatTest(avglag{iHyp});
        end
    end
end

% Save the figure
if isfield(D{1}, 'Seq'), save2pdf('figs/F_Dyn_CoordS.pdf');
else, save2pdf('figs/F_Dyn_CoordIO.pdf');
end

% Display the same trajectories but in the triangular space
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Properties of the triangle
tricc  = [0, sqrt(3)/2; 1, sqrt(3)/2; 1/2, 0];
tricol = [066 146 198; 239 059 033; 065 171 093] ./ 255;

% Prepare a new window
figure('Position', [1 906 200 400]);

% For change- and detection- points
for lock = 1:2
    subplot(2,1,lock);
    
    % Display the triangle
    Emergence_PlotTriInfo(tricc, tricol);
    
    % Display a useful custom grid on the triangle
    for k = 1:3
        avglag = Emergence_PlotGridOnTri(2, k, tricol(k,:), tricc);
        set(avglag(1,k), 'EdgeAlpha', 3/4);
        avglag = Emergence_PlotGridOnTri(3, k, tricol(k,:), tricc);
        for kk = 1:3, set(avglag(kk,k), 'LineStyle', '--', 'EdgeAlpha', 3/4); end
    end
    
    % Display the trajectories
    for iHyp = 1:2
        cc = avgsubtraj{iHyp,lock}*tricc; % cartesian coordinates
        plot(cc(:,1), cc(:,2), '.-', 'Color', tricol(iHyp,:), ...
            'LineWidth', 1, 'MarkerSize', 10);
        
        % Make sure the triangle is equilateral and we don't see the axes
        axis('equal'); axis('off');
    end
end

% Save the figure
if isfield(D{1}, 'Seq'), save2pdf('figs/F_Dyn_TriS.pdf');
else, save2pdf('figs/F_Dyn_TriIO.pdf');
end

%% FIT INDIVIDUAL TRAJECTORIES USING SIGMOID FUNCTIONS
%  ===================================================

% Fit a sigmoid function to beliefs recorded for each individual trial 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Fit trajectories locked on detection point
lock = 1;

% Try to load results from the previous analysis
try
funname = @(x) fullfile(homedir, 'Finger tracking analyses', ...
    'ppdata', sprintf('Dyn1_MFX_%s.mat', x));
if isfield(D{1}, 'Seq'), load(funname('S'));
else, load(funname('IO'));
end
catch

% Prepare output variables
p_sub = cell(1,2); p_gp = cell(1,2); % posterior over parameters
o_sub = cell(1,2); o_gp = cell(1,2); % quality of fit

% For each sequence that entailed a regularity
for iHyp = 1:2
    for iSeq = 1:numel(cidx{iHyp})
        
        % Specify the dimensions of the problem
        dim         = []; % empty variable
        dim.n_phi   = 4;  % number of observation parameters
        dim.n       = 0;  % number of hidden states
        dim.n_theta = 0;  % number of evolution parameters
        dim.n_t     = 1;  % number of time samples
        
        % Define minimal options
        options             = [];
        options.DisplayWin  = 0;
        options.verbose     = 0;
        options             = repmat({options}, [nSub,1]);
        optiongp            = [];
        optiongp.DisplayWin = 0;
        optiongp.verbose    = 0;
        
        % Variable to explain: barycentric coordinates along the relevant
        % dimension
        y = mat2cell(squeeze(fingerwrtp{iHyp,lock}(:,iHyp,iSeq,:)), ...
            nSamp+1, ones(nSub,1))';
        
        % Remove NaNs
        nanloc = cellfun(@(x) isnan(x), y, 'UniformOutput', 0);
        y = cellfun(@(x,i) x(~i), y, nanloc, 'UniformOutput', 0);
        
        % Explaining variable: position of the observation post-change-point
        for iSub = 1:nSub, options{iSub}.inG.p = find(~nanloc{iSub})-1; end
        
        % Run the mixed-effect fitting scheme (it deduces priors through
        % empirical Bayes) using a variational procedure (it uses Laplace
        % approximation) in order to estimate the best parameters for each
        % model, each condition and each subject
        [p_sub{iHyp}(iSeq,:), o_sub{iHyp}(iSeq,:), ...
         p_gp{iHyp}(iSeq,:),  o_gp{iHyp}(iSeq,:)] = ...
            VBA_MFX(y, [], [], @g_SIGM, dim, options, [], optiongp);
    end
end

% Save the result of this analysis in a MATLAB file
if isfield(D{1}, 'Seq'), save(funname('S'), 'p_sub', 'o_sub', 'p_gp', 'o_gp');
else, save(funname('IO'), 'p_sub', 'o_sub', 'p_gp', 'o_gp');
end
end

% Get back the fitted parameters of the sigmoid
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% For each type of regularity
sigparam = cell(1,2);
for iHyp = 1:2
    
    % Get the (4) fitted parameters of the sigmoid function
    sigparam{iHyp} = cellfun(@(x) reshape(x.muPhi, [1,1,4]), ...
        p_sub{iHyp}, 'UniformOutput', 0);
    
    % Consider only sequences for which regularities were correctly identified
    detecmask = cellfun(@(x) x.Questions(2) == iHyp, G(cidx{iHyp},:));
    sigparam{iHyp}(~detecmask) = {NaN(1,1,4)};
    sigparam{iHyp} = cell2mat(sigparam{iHyp});
end

% Display average sigmoid function (based on averaged fitted parameters)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Prepare a new window
figure('Position', [703 905 220 200]); lgd = NaN(1,2); hold('on');

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
if isfield(D{1}, 'Seq'), save2pdf('figs/F_Dyn_SigmS.pdf');
else, save2pdf('figs/F_Dyn_SigmIO.pdf');
end

% Compare the slope of the sigmoid in the two types of regularity
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Average the slope parameter over sequences for each subject and each type
% of regularity (probabilistic and deterministic ones)
slope = cell2mat(cellfun(@(x) mean(exp(x(:,:,1)), 1, 'OmitNaN')', ...
    sigparam, 'UniformOutput', 0));

% Paired t-test
[~,pval,tci,stats] = ttest(diff(slope, 1, 2));
disptstats(pval,tci,stats);

% Prepare a new window
figure('Position', [924 905 100 200]);

% Display the difference of the sigmoids' slope in the two types of
% regularity
Emergence_PlotSubGp(slope, tricol(1:2,:));

% Customize the axes
set(gca, 'XTick', [], 'Box', 'Off');
axis([0,3,0,1.2]);

% Display whether the difference is significant or not
Emergence_DispStatTest(slope);

% Add some text labels
ylabel('Slope parameter');

% Save the figure
if isfield(D{1}, 'Seq'), save2pdf('figs/F_Dyn_SlopeS.pdf');
else, save2pdf('figs/F_Dyn_SlopeIO.pdf');
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
    lowerlim = sigmoid(P(3)); % E [0,1]
    g = g .* lowerlim;
    
    % If provided, allows the sigmoid to end below one
    if numel(P) > 3
        upperlim = sigmoid(P(4)); % E [0,1]
        g = g + upperlim;
    end
end

% Make sure the output is a column vector
g = g(:);

% Rely on numerical derivatives
dgdx = [];
dgdP = [];

end
