% This script shows examples of the different link functions that are used
% to combine the independently-combined likelihood 
% 
% Copyright (c) 2020 Maxime Maheu

%% COMBINATION FUNCTIONS
%  =====================

% Prepare a new window
figure('Position', [1 895 200 450]);

% Scenario #1: Linear weighting between q(Hd|y) and q(Hp|y)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Display the link function
subplot(4,1,1); hold('on');
plot([-1,1], [0,1], 'Color', 'k', 'LineWidth', 1);

% Add some text labels
axis([-1,1,0,1]);
xlabel('p(Hd|y)-p(Hd|y)');
ylabel('p(Hd|y)');
title('Difference');

% Scenario #2: Winner-takes-all between q(Hd|y) and q(Hp|y)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Display the link function
subplot(4,1,2); hold('on');
plot([-1,0,0,1], [0,0,1,1], 'Color', 'k', 'LineWidth', 1);

% Add some text labels
axis([-1,1,0,1]);
xlabel('p(Hd|y)-p(Hp|y)');
ylabel('p(Hd|y)');
title('Maximum');

% Scenario #3: Sigmoid weighting on q(Hd|y)-q(Hp|y)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Compute sigmoid functions
sigmparam = 10.^[0;1/4;1/2;3/4;1];
nParam = numel(sigmparam);
differ = -1:0.001:1;
pHpgY = 1 ./ (1 + exp(-sigmparam .* differ));

% Display the link functions
subplot(4,1,3); hold('on');
l = plot(differ, pHpgY, 'LineWidth', 1);
col = cool(nParam);
for i = 1:nParam, set(l(i), 'Color', col(i,:)); end

% Add some text labels
axis([-1,1,0,1]);
xlabel('p(Hd|y)-p(Hp|y)');
ylabel('p(Hd|y)');
title('Softmax(Difference)');

% Scenario #3: Sigmoids weighting on log(q(Hd|y)/q(Hp|y))
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Compute sigmoid functions
sigmparam = 10.^[0;1/4;1/2;3/4;1];
nParam = numel(sigmparam);
ratio = 0:0.01:4;
pHpgY = 1 ./ (1 + exp(-sigmparam .* log(ratio)));

% Display the link functions
subplot(4,1,4); hold('on');
l = plot(ratio, pHpgY, 'LineWidth', 1);
col = cool(nParam);
for i = 1:nParam, set(l(i), 'Color', col(i,:)); end

% Add some text labels
xlabel('p(Hd|y)/p(Hp|y)');
ylabel('p(Hd|y)');
title('Sigmoid(log(Ratio))');

%% EXAMPLE PROJECTION IN THE TRIANGLE
%  ==================================

% Define independently computed pseudo posteriors
qHpgY = 3/4;
qHdgY = 1/2;

% Define value of the slope parameter for functions which use it
slope = 4.5;

% Prepare output variable
cc = NaN(4,2);

for iFun = 1:4
    
    % Get weight
    if     iFun == 1, Wp = qHpgY ./ (qHpgY + qHdgY);
    elseif iFun == 2, Wp = double((qHpgY - qHdgY) > 0);
    elseif iFun == 3, Wp = 1 ./ (1 + exp(-slope .* (qHpgY - qHdgY)));
    elseif iFun == 4, Wp = 1 ./ (1 + exp(-slope .* log(qHpgY ./ qHdgY)));
    end
    
    % Compute posterior probabilities
    pHsgY = (1 - qHpgY) .* Wp + (1 - qHdgY) .* (1 - Wp);
    pHpgY = (1 - pHsgY) .* Wp;
    pHdgY = (1 - pHsgY) .* (1 - Wp);
    bc = [pHpgY pHdgY pHsgY];
    
    % Transform to cartesian coordinates
    cc(iFun,:) = bc*tricc;
end

% Prepare a new window
figure;

% Display the pseudo posteriors
subplot(1,2,1); hold('on');
b = bar([zeros(1,2); ones(1,2)], [qHpgY, 1-qHpgY; zeros(1,2)], 'stacked');
set(b(1), 'FaceColor', tricol(1,:));
set(b(2), 'FaceColor', tricol(3,:));
b = bar([zeros(1,2); ones(1,2)], [zeros(1,2); qHdgY, 1-qHdgY], 'stacked');
set(b(1), 'FaceColor', tricol(2,:));
set(b(2), 'FaceColor', tricol(3,:));
set(gca, 'XColor', 'none'); axis('square');
ylabel('Pseudo posterior probabilities');

% Display the output posteriors
subplot(1,2,2);
Emergence_PlotTriInfo;
plot(cc(:,1), cc(:,2), 'ko');
