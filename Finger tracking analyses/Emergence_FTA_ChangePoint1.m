% This script derives predictions from the ideal observer and show that
% beliefs regarding the position of the change point are more precise in
% the case of the emergence of a deterministic compared to a probabilistic
% regularity.
%
% Copyright (c) 2018 Maxime Maheu

%% GET PREDICTIONS FROM THE IDEAL OBSERVER
%  =======================================

% Get the length of the sequence
N = S{1}.Nstims;

% Choose the size of the window to look into
wwin = 40;
xwin = -wwin/2:wwin/2;
nwin = wwin+1;

% Prepare output variables
pcppos      = cell(2);
pcpposwrtcp = NaN(nSub,nwin,2);
precpcppos  = NaN(nSub,2);

% For each type of regularity
for iHyp = 1:2
    
    % Get the position of the change points
    cp = cellfun(@(x) x.Jump-1/2, G(cidx{iHyp},:));
    
    % Get beliefs about the change point's position at the end of each
    % sequence
    lab = sprintf('pJkgYMs%s', lower(proclab{iHyp}(1)));
    bel = cellfun(@(x) x.(lab)(:,end)', IO(cidx{iHyp},:), 'UniformOutput', 0);
    
    % Get sequences that were corretly labeled 
    detecmask = (filter{iHyp} == 1 | filter{iHyp} == 3);
    
    % Get the posterior distribution of change point's position separately
    % for (in)correctly labeled sequences
    for d = [0,1]
        cbel = bel(detecmask == d);
        [~,idx] = sort(cp(detecmask == d));
        pcppos{d+1,iHyp} = cell2mat(cbel(idx));
    end
    
    % Get a distibution of change point's position restricted around true
    % change point's position
    winbel = cellfun(@(b,c) b(c+xwin), bel, num2cell(cp), 'UniformOutput', 0);
    
    % Remove posterior beliefs about change point's positions from
    % sequences that were mislabeled
    winbel(~detecmask) = {NaN(1,nwin)};
    winbel = cellfun(@(x) reshape(x,[1,1,nwin]), winbel, 'UniformOutput', 0);
    
    % Save the averaged full posterior distribution over change point
    % position
    winbel = cell2mat(winbel);
    pcpposwrtcp(:,:,iHyp) = squeeze(mean(winbel, 1, 'OmitNaN'));
    
    % Measure the log-precision of the posterior distributions over change
    % point position
    % => log(precision) = log(1/variance) = -log(variance)
    obsidx = 1:N;
    mu = cellfun(@(p) sum(p.*obsidx), bel);
    precision = cellfun(@(p,m) -log(sum(p.*((obsidx-m).^2))), bel, num2cell(mu));
    precision(~detecmask) = NaN;
    precpcppos(:,iHyp) = mean(precision, 1, 'OmitNaN');
end

%% DISPLAY AVERAGED POSTERIOR DISTRIBUTIONS AROUND TRUE CHANGE POINT
%  =================================================================

% Prepare a new window
figure('Position', [1 906 340 200]);

% For each type of regularity
lgd = NaN(1,2);
for iHyp = [2,1]
    
    % Average posterior distributions
    d = pcpposwrtcp(:,:,iHyp);
    m = mean(d, 1);
    
    % Display the posterior distribution of change point's positions
    % centered around the real position of the change point
    lgd(iHyp) = fill([xwin(1), xwin, xwin(end)], [0, m, 0], 'k', 'FaceColor', ...
        tricol(iHyp,:), 'EdgeColor', 'None', 'FaceAlpha', 1/3); hold('on');
    plot(xwin, m, '-', 'Color', tricol(iHyp,:), 'LineWidth', 2)
    
    % Perform a t-test against the uniform (prior) distribution at each
    % position around the true change point's position
    h = ttest(d, 1/N, 'tail', 'right');
    h(isnan(h)) = 0;
    plot(xwin(logical(h)), repmat(0.34+0.005*(iHyp-1),1,sum(h)), ...
        '-', 'Color', tricol(iHyp,:), 'LineWidth', 2);
end

% Display the position of the change point (i.e. 0)
plot(zeros(1,2), ylim, 'k--');

% Display the uniform scenario over change point's position
plot(xwin([1,end]), repmat(1/N,1,2), '-', 'Color', g);

% Customize the axes
axis([xwin([1,end]),0,0.4]);
set(gca, 'Box', 'Off');

% Add some text labels
legend(lgd, proclab(1:2), 'Location', 'NorthWest');
xlabel('Position w.r.t. change point');
ylabel('p(j_k|H_i,y)');

% Save the figure
save2pdf(fullfile(ftapath, 'figs', 'F_CP_AvgCtrPostCP.pdf'));

%% DISPLAY SINGLE-TRIAL FULL POSTERIOR DISTRIBUTIONS
%  =================================================

% Prepare a new window
figure('Position', [1 431 230 400]);
cmapcol = {'Blues', 'Reds'};

% For sequences with a probabilistic/deterministic regularity
for iHyp = 1:2
    
    % Display the change in beliefs as a heatmap 
    sp = subplot(2,1,iHyp);
    imagesc(cell2mat(pcppos(:,iHyp))); hold('on');
    
    % Customize the colormap
    colorbar('Location', 'EastOutside');
    colormap(sp, cbrewer2(cmapcol{iHyp}));
    caxis([0,0.3+0.2*(iHyp==2)]);
    
    % Display limits between (in)correctly labeled sequences
    lim = size(pcppos{1,iHyp},1) + 1/2;
    plot([0,N+1], repmat(lim, 1, 2), 'k-');
    
    % Customize the axes
    axis('xy'); set(gca, 'XTick', [1, get(gca, 'XTick')], 'YTick', [1,50]);
    
    % Add some text labels
    xlabel('Observation #');
    ylabel({'Sequence # (sorted by', 'change point''s position)'});
    title(sprintf('%s sequences', proclab{iHyp}));
end

% Save the figure
save2pdf(fullfile(ftapath, 'figs', 'F_CP_SingTrlPostCP.pdf'));

%% DISPLAY AVERAGED PRECISION OF POSTERIOR DISTRIBUTIONS
%  =====================================================

% Prepare a new window
figure('Position', [342 906 120 200]);

% Display difference in precision between the two types of regularity
Emergence_PlotSubGp(precpcppos, tricol(1:2,:));

% Customize the axes
xlim([0,3]);
set(gca, 'XTick', [], 'XColor', 'None', 'Box', 'Off');

% Display whether the difference is significant or not
Emergence_DispStatTest(precpcppos);

% Add some text labels
ylabel('Log-precision');

% Save the figure
save2pdf(fullfile(ftapath, 'figs', 'F_CP_PrecPostCP.pdf'));

% Compute a paired t-test on the precision of posterior distributions
% between the 2 types of regularities
[~,pval,tci,stats] = ttest(diff(precpcppos, 1, 2));
Emergence_PrintTstats(pval, tci, stats);
