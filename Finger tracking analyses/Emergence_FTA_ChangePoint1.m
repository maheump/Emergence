% This script derives predictions from the ideal observer and show that
% beliefs regarding the position of the change point are more precise in
% the case of the emergence of a deterministic compared to a probabilistic
% regularity.
%
% Copyright (c) 2018 Maxime Maheu

%% PREDICTIONS FROM THE IDEAL OBSERVER
%  ===================================

% Choose the size of the window to look into
winwidth = 50;

% Get the beliefs from the IO regarding the position of the change point
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Get a vector of position centered around the change point
x = -winwidth/2:winwidth/2;

% Prepare output variables
cppost = NaN(nSub, winwidth+1, 2);
varcppost = NaN(nSub, 2);

% For each type of regularity
for iHyp = 1:2

    % Get the position of the change point
    cp = cellfun(@(x) x.Jump-1/2, G(cidx{iHyp},:), 'UniformOutput', 0);

    % Get beliefs about the change point's position at the end of each
    % sequence
    lab = sprintf('pJkgYMs%s', lower(proclab{iHyp}(1)));
    bel = cellfun(@(x) x.(lab)(end,:), IO(cidx{iHyp},:), 'UniformOutput', 0);

    % Look around the change point's position in the window defined earlier
    belwrtcp = cellfun(@(b,c) b(c-winwidth/2:c+winwidth/2), bel, cp, 'UniformOutput', 0);

    % Average over sequences for each subject
    belwrtcp = mat2cell(belwrtcp, numel(cidx{iHyp}), ones(1,nSub));
    belwrtcp = cellfun(@(x) mean(cell2mat(x), 1), belwrtcp, 'UniformOutput', 0);
    cppost(:,:,iHyp) = cell2mat(belwrtcp');

    % Measure the precision of these posterior distributions
    varcppost(:,iHyp) = cellfun(@(y) 1/std(x,y), belwrtcp);
end

% Display precision of posterior distributions
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Get the length of the sequence
N = S{1}.Nstims;

% Prepare a new window
figure('Position', [1 906 340 200]);

% For each type of regularity
lgd = NaN(1,2);
for iHyp = fliplr(1:2)

    % Average posterior distributions
    d = cppost(:,:,iHyp);
    m = mean(d,1);

    % Display the posterior distribution of change point's positions
    % centered around the real position of the change point
    lgd(iHyp) = fill([x(1), x, x(end)], [0, m, 0], 'k', 'FaceColor', ...
        tricol(iHyp,:), 'EdgeColor', 'None', 'FaceAlpha', 1/3); hold('on');
    plot(x, m, '-', 'Color', tricol(iHyp,:), 'LineWidth', 2)

    % Perform a t-test against the uniform (prior) distribution at each
    % position around the true change point's position
    h = ttest(cppost(:,:,iHyp), 1/N, 'tail', 'right');
    plot(x(logical(h)), repmat(0.34+0.005*(iHyp-1),1,sum(h)), ...
        '-', 'Color', tricol(iHyp,:), 'LineWidth', 2);
end

% Display the position of the change point (i.e. 0)
plot(zeros(1,2), ylim, 'k--');

% Display the uniform scenario over change point's position
plot(x([1,end]), repmat(1/N,1,2), 'k:');

% Customize the axes
xlim(x([1,end]));
set(gca, 'Box', 'Off');

% Add some text labels
legend(lgd, proclab(1:2), 'Location', 'NorthWest');
xlabel('Position w.r.t. change point');
ylabel('p(j_k|H_i,y)');

% Save the figure
save2pdf('figs/F_CP_PredCorr.pdf');

% Display precision of posterior distributions
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Prepare a new window
figure('Position', [342 906 120 200]);

% Display difference in precision between the two types of regularity
Emergence_PlotSubGp(varcppost, tricol(1:2,:));

% Customize the axes
xlim([0,3]);
set(gca, 'XTick', [], 'XColor', 'None', 'Box', 'Off');

% Display whether the difference is significant or not
Emergence_DispStatTest(varcppost);

% Add some text labels
ylabel('1/std');

% Save the figure
save2pdf('figs/F_CP_PredConf.pdf');
