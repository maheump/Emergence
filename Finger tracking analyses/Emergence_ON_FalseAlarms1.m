% This script shows that the average position in the triangle during
% fully-stochastic sequences depends upon the position of the current
% observation in the sequence. More specifically, it shows that false
% alarms becomr more likely to happen as time has passed since the
% beginning of the sequence.
% 
% Copyright (c) 2020 Maxime Maheu

% Define options
% ~~~~~~~~~~~~~~

% Define the number of bins to use (must be a multiple of 200)
nBin = 10;

% Select the sequences to look at
% 1: all fully-stochastic parts
% 2: only fully-stochastic sequences
% 3: only fully-stochastic sequences that were correctly labeled
% 4: only first-part of stochastic-to-regular sequences
restopt = 1;

% Average hypotheses likelihood in different parts of the sequences
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Select observations
randidx = Emergence_SelectFullyStochSeq(G, filter, restopt);

% Get data from fully-stochastic sequences
data = cellfun(@(x,i) [x.BarycCoord(i,:); NaN(sum(~i),3)], D, randidx, 'uni', 0);

% Group positions in the triangle by consecutive bins of n trials
binavg = cell2mat(cellfun(@(x) reshape(x, [1 1 N 3]), data, 'uni', 0));
binavg = permute(binavg, [3 4 2 1]);
binavg = mat2cell(binavg, repmat(N/nBin,nBin,1,1,1), ...
    ones(1,3,1,1), ones(1,1,nSub,1), ones(1,1,1,nSeq));

% Average positions within the triangle
binavg = cellfun(@(x) mean(x, 1,'OmitNaN'), binavg); % average within each bin
binavg = mean(binavg, 4, 'OmitNaN'); % average over sequences
avgbc = mean(binavg, 3); % average over subjects
errbc = sem(binavg, 3); % measure dispersion over subjects

% Convert beliefs in each hypothesis to cartesian coordinates
avgcc = flipud(avgbc*tricc);
semcc = flipud(errbc*tricc);

% Display average position in the triangle according to time spend since
% the beginning of the sequence (in number of observations)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Define plot options
cmap = flipud(Emergence_Colormap('Viridis', nBin));
dotsize = 50;

% Prepare a new window
figure('Position', [1 805 200 300]);

% Display the triangular arena
Emergence_PlotTrajOnTri; alpha(0);
Emergence_PlotGridOnTri(10, 3, tricol(3,:));

% Display error bars
plot(repmat(avgcc(:,1), 1, 2)', avgcc(:,2)' + [-1;1] .* semcc(:,2)', 'k-', 'LineWidth', 1/2);
plot(avgcc(:,1)' + ([-1;1] .* semcc(:,1)'), repmat(avgcc(:,2), 1, 2)', 'k-', 'LineWidth', 1/2);

% Display averaged positions 
scatter(flipud(avgcc(:,1)), flipud(avgcc(:,2)), ...
    dotsize, cmap, 'Filled', 'MarkerEdgeColor', 'k');

% Add a colorbar
cbr = colorbar('Location', 'SouthOutside');
cbr.Label.String = 'Position within the sequence';
colormap(cmap); caxis([0,1]);

% Zoom in the bottom part of the triangle
axis('tight');
ylim([0,(sqrt(3)/2)/2]);

% Save the figure
if isfield(D{1}, 'Seq'), save2pdf(fullfile(ftapath, 'figs', 'F_FA_SeqGpS.pdf'));
else, save2pdf(fullfile(ftapath, 'figs', 'F_FA_SeqGpIO.pdf'));
end

% Display the distribution of correlation coefficients over subjects
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Concatenate data for each subject across the different sequences
inddata = cellfun(@(x) cell2mat(x), mat2cell(data, nSeq, ones(1,nSub)), 'uni', 0);

% Correlate likelihood of the fully-stochastic hypothesis 
nanidx = cellfun(@(x) ~isnan(x(:,3)), inddata, 'uni', 0);
obspos = repmat((1:N)', nSeq, 1);
coef = cellfun(@(toexplain,idx) ...
    Emergence_Regress(toexplain(idx,3), obspos(idx), 'CC', 'r'), inddata, nanidx)';

% Prepare a new window
figure('Position', [202 905 120 200]);

% Display chance level
plot([0,2], zeros(1,2), '-', 'Color', g); hold('on');

% Display distribution of correlation coefficients
Emergence_PlotSubGp(coef, zeros(1,3));

% Customize the axes
set(gca, 'Box', 'Off', 'XTick', [], 'XColor', 'None', 'XLim', [0,2]);

% Display whether the difference is significant or not
Emergence_DispStatTest(coef);

% Add some text label
ylabel('Correlation coefficient');

% Save the figure
if isfield(D{1}, 'Seq'), save2pdf(fullfile(ftapath, 'figs', 'F_FA_SeqCorrS.pdf'));
else, save2pdf(fullfile(ftapath, 'figs', 'F_FA_SeqCorrIO.pdf'));
end

% Perform a t-test (against chance) on correlation coefficients
[~,pval,tci,stats] = ttest(coef);
Emergence_PrintTstats(pval,tci,stats);
