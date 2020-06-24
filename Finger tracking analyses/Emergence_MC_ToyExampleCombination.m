% This script shows examples of the different link functions that are used
% to combine the independently-combined likelihood 
% 
% Copyright (c) 2020 Maxime Maheu

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
