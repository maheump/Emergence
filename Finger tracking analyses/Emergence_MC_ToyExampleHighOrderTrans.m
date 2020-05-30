% This script is a toy example of (1) how a given pattern translate into
% biases of apparent transition probabilities of different orders and (2)
% how sequence likelihood changes with (a) the order of transition
% probabilities that is considered and (b) the prior that is used (flat or
% more or less biased).
% 
% Copyright (c) 2018 Maxime Maheu

%% APPARENT TRANSITION PROBABILITIES (OF DIFFERENT ORDERS) OF A GIVEN PATTERN
%  ==========================================================================

% Define the pattern from which transition probabilities will be estimated
pat = 'AAAAABBBBB';

% Almost infinite number of repetitions of that pattern
y = str2pat(pat);
y = repmat(y, [1, 1e5]);
K = numel(y);

% Define the orders of the Markov chains to test
orders = 1:5;
maxncond = 2^max(orders);

% Prepare a new window
figure('Position', [329 432 616 302], 'Name', pat);

% Define the order of the Markov chain
for o = 1:numel(orders)
    order = orders(o);
    cond = ff2n(order) + 1;
    ncond = 2^order;
    
    % Count the number of As and Bs after each transition
    o_nXgX = NaN(ncond,2);
    for c = 1:ncond
        
        % Locate positions of the current transition
        condpos = strfind(y, cond(c,:));
        
        % Get the position of the observations coming just after
        following = condpos + order;
        following = following(following <= K);
        
        % Count the number of such observations separately for As and Bs
        o_nXgX(c,1) = sum(y(following) == 1); % for As
        o_nXgX(c,2) = sum(y(following) == 2); % for Bs
    end
	
    % Add prior counts
    p_nXgX = ones(size(o_nXgX));
    nXgX = o_nXgX + p_nXgX;
    pXgX = o_nXgX ./ sum(o_nXgX, 2);
    
    % Get the name of each transition
    letters = {'A','B'};
    condlab = mat2cell(cond, ones(ncond,1), order);
    condlab = cellfun(@(x) pat2str(x, letters, [1,2]), condlab, 'uni', 0);
    
    % Compute bars' horizontal positions
    tmp = 2^(max(orders)-order);
    ypos = mean(reshape(1:maxncond, [tmp, maxncond/tmp]), 1);
    
    % Display probability of getting As/Bs after each transition
    subplot(1, numel(orders), o); hold('on');
    b = barh(ypos, pXgX, ncond/maxncond, 'Stacked', 'EdgeColor', 'k');
    set(b(1), 'FaceColor', 'k');
    set(b(2), 'FaceColor', 'w');
    
    % Customize the axes
    set(gca, 'YTick', ypos, 'YTickLabel', condlab, 'TickLen', zeros(1,2), ...
        'Box', 'Off', 'YDir', 'Reverse');
    axis([0,1,1/2,maxncond+1/2]);
    
    % Add some text labels
    xlabel('Probability');
    if o == 1, ylabel('Previous n-gram'); end
    title(sprintf('Order %1.0f', orders(o)));
end

%% HOW SEQUENCE LIKELIHOOD CHANGES WITH ORDER AND PRIOR BELIEFS
%  ============================================================

% Prior pseudo-counts to test
pN = 1/2:0.01:1;

% Orders of the Markov chain
orders = 1:9;

% Compute prior log-likelihood
out = (sum(gammaln([pN; pN]), 1) - gammaln(sum([pN; pN], 1))) .* orders';

% Disply prior likelihood as a function of prior pseudo-counts and order of
% the Markov chain
figure('Position', [946 432 354 302]);
l = plot(pN, out, 'LineWidth', 2);
chaincol = flipud(winter(numel(orders)));
col = [tricol(1,:); chaincol(2:end,:)];
for i = 1:numel(orders), set(l(i), 'Color', col(i,:)); end

% Customize the axes
xlim([1/2,1]); axis('tight');
axl = get(gca, 'XTickLabel');
axl([1,end]) = {'Jeffreys', 'Bayes-Laplace'};
set(gca, 'XTickLabel', axl, 'Box', 'Off');

% Add some text labels
xlabel('Prior pseudo-counts');
ylabel('Prior log-likelihood');
legend(arrayfun(@(x) sprintf('Order %1.0f', x), orders, 'uni', 0));
