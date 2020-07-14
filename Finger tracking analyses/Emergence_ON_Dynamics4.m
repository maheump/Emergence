% To prevent averaging effects to mask effects or lead to incorrect
% conclusions, sigmoid functions are fitted to individual trial detection
% dynamics. Resulting fitted slope parameters are then averaged over the
% different sequences in order to obtain an average detection scenario for
% each type of regularity. This confimed that the dynamics of detection is
% much more ballistic in case of deterministic compared to probabilistic
% regularities.
% 
% Copyright (c) 2020 Maxime Maheu

%% SIMULATE DIFFERENT DETECTION DYNAMICS USING SIGMOID FUNCTIONS
%  =============================================================

% Define parameter grid
% ~~~~~~~~~~~~~~~~~~~~~

% Prepare output variable
pgrid = cell(1,4);

% Define grid over slope parameter
slope  = logspace(-3,0,31);
pgrid{1} = slope;

% Define grid over intercept parameter
intcp  = 0:N;
pgrid{2} = intcp;
intcp  = reshape(intcp, [1,1,numel(intcp)]);

% Define grid over offset parameter
offset = 0:0.05:0.5;
pgrid{3} = offset;
offset = reshape(offset, [1,1,1,numel(offset)]);

% Define grid over scaling
scaling = 0.5:0.05:1;
pgrid{4} = scaling;
scaling  = reshape(scaling, [1,1,1,1,numel(scaling)]);

% Define values to evaluate with the sigmoid functions, here, the
% observation number post change point
x = 0:145;

% Compute sigmoid functions over the grid of parameters
simu = 1 ./ (1 + exp(slope.*(-x'+intcp)));
simu = simu .* scaling + offset;

% Display range of parameters
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Prepare a new figure
figure('Position', [1 905 700 200]);

% Display grid of slope parameter
subplot(1,4,1);
plot(x, squeeze(simu(:,:,73,1,end)));

% Display grid of intercept parameter
subplot(1,4,2);
plot(x, squeeze(simu(:,30,:,1,end)));

% Display grid of offset parameter
subplot(1,4,3);
plot(x, squeeze(simu(:,30,73,:,end)));

% Display grid of scaling parameter
subplot(1,4,4);
plot(x, squeeze(simu(:,30,73,1,:)));

% Customize the axes
for iHyp = 1:4
    subplot(1,4,iHyp);
    axis([0,145,0,1]);
    xlabel('Observation #');
    ylabel('p(H|y)');
end

%% COMPUTE ERROR BETWEEN DETECTION DYNAMICS AND SIGMOID FUNCTIONS
%  ==============================================================

% Compute error
% ~~~~~~~~~~~~~

% Prepare output variables
MSE = cell(1,2);

% For each category of regularity
for iHyp = 1:2
	
    % Get post-change point posterior likelihood in the corresponding
    % (correct) hypothesis
    data = cellfun(@(p,c) Emergence_LockOnPoint(p.BarycCoord(:,iHyp), ...
        c.Jump, [0,145]), D(cidx{iHyp},:), G(cidx{iHyp},:), 'uni', 0);
    
    % Preallocate matrices
    MSE{iHyp} = cell(size(data));
    
    % For each subject and each regularity of this category
    for i = 1:numel(data)
        if any(i == round(linspace(1, numel(data), 20)))
            fprintf('%1.0f%%|', 100*i/numel(data));
        end
        
        % Compute the sum/mean squared error between all simulated sigmoid
        % functions and the observed beliefs dynamics
        SSE = squeeze(sum((simu - data{i}).^2, 1, 'OmitNaN'));
        n = sum(~isnan(data{i}));
        [i1,i2] = ind2sub(size(data), i);
        MSE{iHyp}{i1,i2} = SSE ./ n;
        % N.B. we save the MSE because it can be compared from one sequence
        % to the next (because it is normalized by the number of
        % observations post-change point, which can vary from one sequence
        % to the next)
    end
    fprintf('\n');
end

%% DISPLAY THE RESULTS OF THE GRID-SEARCH FIT PROCEDURE
%  ====================================================

% Find best parameters
% ~~~~~~~~~~~~~~~~~~~~

% Prepare the output variable
pmes = cell(1,2);

% For each type of regularity
for iHyp = 1:2

    % Find the set of sigmoid parameters that induces the smallest error
    % (indices along the grid)
    [~,ind] = cellfun(@(x) min(x(:)), MSE{iHyp});
    pind = cell(1,4);
    [pind{1},pind{2},pind{3},pind{4}] = ind2sub(size(MSE{1}{1}), ind);
    
    % Transform indices along the grid into real parameter values
    pmes{iHyp} = cellfun(@(x,y) x(y), pgrid, pind, 'uni', 0);
    
    % Restrict to regularities accurately labeled
    detecmask = (filter{iHyp} == 3);
    for iParam = 1:4, pmes{iHyp}{iParam}(~detecmask) = NaN; end
end

% Display best parameters
% ~~~~~~~~~~~~~~~~~~~~~~~

% Prepare a new window
figure('Position', [702 905 120*4 200]);
pname = {'Slope', 'Intercept', 'Offset', 'Scaling'};

% For each parameter
for iParam = 1:4
    subplot(1,4,iParam);
    
    % Get slope parameter for each group of regularities for each subject
    currparam = [mean(pmes{1}{iParam}, 1, 'OmitNaN')', ...
                 mean(pmes{2}{iParam}, 1, 'OmitNaN')'];
    
    % Perform a paired t-test between slope parameters for probabilistis vs.
    % deterministic regularities
    [h,pval,tci,stats] = ttest(diff(currparam, 1, 2));
    Emergence_PrintTstats(pval,tci,stats);
    
    % Display the paired difference in slope parameters
    Emergence_PlotSubGp(currparam, tricol(1:2,:));
    
    % Customize the axes
    set(gca, 'XLim', [0,3], 'XTick', [], 'XColor', 'None');
    ylabel(sprintf('%s parameter', pname{iParam}));
end

% Display average sigmoid function
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Define the sigmoid functions
sigm = @(x,p) (1 ./ (1 + exp(p(1).*(-x'+p(2))))) .* p(4) + p(3);

 % Prepare a new window
figure('Position', [1183 905 220 200]); lgd = NaN(1,2); hold('on');

 % Draw some help lines
plot(x([1,end]),    ones(1,2)./2, '-',  'Color', g, 'LineWidth', 1/2); 
plot(x([1,end]),    ones(1,2)./3, '--', 'Color', g, 'LineWidth', 1/2);
plot(x([1,end]), 2.*ones(1,2)./3, '--', 'Color', g, 'LineWidth', 1/2);

% For each type of regularity
for iHyp = 1:2
    
    % Average over sequences
    param = cell2mat(cellfun(@(x) mean(x, 1, 'OmitNaN')', pmes{iHyp}, 'uni', 0));
    
    % Average over subjects
    avgparam = mean(param, 1, 'OmitNaN');
    semparam = sem(param, 1);
    
    % Draw error bars
    paramcomb = mat2cell(avgparam + semparam .* (ff2n(4).*2-1), ones(2^4,1), 4);
    allsigm = cell2mat(cellfun(@(p) sigm(x,p)', paramcomb, 'uni', 0));
    lowerlim = max(allsigm, [], 1);
    upperlim = min(allsigm, [], 1);
    fill([x, fliplr(x)], [upperlim, fliplr(lowerlim)], 'k', 'FaceColor', ...
        tricol(iHyp,:), 'EdgeColor', 'None', 'FaceAlpha', 0.15); hold('on');

    % Draw the average sigmoid curve
    avgpred = sigm(x, avgparam);
    lgd(iHyp) = plot(x, avgpred, 'Color', tricol(iHyp,:), 'LineWidth', 3); 
end

 % Customize the axes
set(gca, 'Box', 'Off');
axis([0,60,0,1]);

 % Add some text labels
xlabel('# observation w.r.t. change point');
ylabel('Fitted posterior belief p(M_i|y)');
legend(lgd, proclab, 'Location', 'NorthWest');

% Display MSE over the grid for the slope parameter
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Prepare new figure
figure('Position', [1404 905 220 200]);

% For each regularity
for iHyp = 1:2
    
    % Minimize MSE over the other parameters
    paramcomb = cellfun(@(x) squeeze(min(x,[],2:4)), MSE{iHyp}, 'uni', 0);
    
    % Restrict to regularities accurately labeled
    detecmask = (filter{iHyp} == 3);
    paramcomb(~detecmask) = {NaN(size(paramcomb{1}))};
    
    % Average over sequences
    paramcomb = cell2mat(cellfun(@(x) reshape(x, [1,1,size(x,1)]), paramcomb ,'uni', 0));
    paramcomb = squeeze(mean(paramcomb, 1, 'omitnan'));
    
    % Overlap MSE for both types of regularity
    if     iHyp == 1, yyaxis('left');
    elseif iHyp == 2, yyaxis('right');
    end
    
    % Display MSE averaged over subjects
    plotMSEM(pgrid{1}, mean(paramcomb), sem(paramcomb), 0.15, tricol(iHyp,:));
    
    % Display the best slope parameter (i.e. with the smallest MSE)
    [m,s] = min(mean(paramcomb));
    plot(pgrid{1}(s), m, '.', 'MarkerEdgeColor', tricol(iHyp,:), 'MarkerSize', 20);
    
    % Customize the axes
    axis('tight'); xlim([0.1,0.9]);
    set(gca, 'YColor', tricol(iHyp,:), 'YScale', 'log', 'Box', 'Off');
end

% Add some text labels
xlabel('Slope parameter');
