% To prevent averaging effects to mask effects or lead to incorrect
% conclusions, sigmoid functions are fitted to individual trial detection
% dynamics. Resulting fitted slope parameters are then averaged over the
% different sequences in order to obtain an average detection scenario for
% each type of regularity. This confimed that the dynamics of detection is
% much more ballistic in case of deterministic compared to probabilistic
% regularities.
% 
% Copyright (c) 2018 Maxime Maheu

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
for i = 1:4
    subplot(1,4,i);
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
        c.Jump+1/2, [0,145]), D(cidx{iHyp},:), G(cidx{iHyp},:), 'uni', 0);
    
    % Preallocate matrices
    MSE{iHyp} = cell(size(data));
    
    % For each subject and each regularity of this category
    for i = 1:numel(data)
        disp(i/numel(data));
        
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
end

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
    
    % Transform indices along the grid into real paratemer values
    pmes{iHyp} = cellfun(@(x,y) x(y), pgrid, pind, 'uni', 0);
    
    % Restrict to regularities accurately labeled
    detecmask = (filter{iHyp} == 3);
    for iParam = 1:4, pmes{iHyp}{iParam}(~detecmask) = NaN; end
end

%% DISPLAY BEST PARAMETERS
%  =======================

% Get slope parameter for each group of regularities for each subject
slopeparam = [mean(pmes{1}{1}, 1, 'OmitNaN')', mean(pmes{2}{1}, 1, 'OmitNaN')'];

% Perform a paired t-test between slope parameters for probabilistis vs.
% deterministic regularities
[h,pval,ci,stats] = ttest(diff(slopeparam, 1, 2));
Emergence_PrintTstats(pval,ci,stats);

% Prepare a new window
figure('Position', [702 905 120 200]);

% Display the paired difference in slope parameters
Emergence_PlotSubGp(slopeparam, tricol(1:2,:));

% Customize the axes
set(gca, 'XLim', [0,3], 'XTick', [], 'XColor', 'None');

% Add text labels
ylabel('Slope parameter');
