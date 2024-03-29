% This script derives predictions from the ideal observer and show that
% beliefs regarding the position of the change point are more precise in
% the case of the emergence of a deterministic compared to a probabilistic
% regularity.
%
% Copyright (c) 2020 Maxime Maheu

%% GET PREDICTIONS FROM THE IDEAL OBSERVER
%  =======================================

% Get the length of the sequence
N = S{1}.Nstims;

% Choose the size of the window to look into
wwin = 40;
xwin = -wwin/2:wwin/2;
nwin = wwin+1;

% Prepare output variables
cp          = cell(1,2);
cppos       = cell(2);
pcppos      = cell(2);
pcpposwrtcp = NaN(nSub,nwin,2);
precpcppos  = NaN(nSub,2);

% For each type of regularity
for iHyp = 1:2
    
    % Get the position of the change points
    cp{iHyp} = cellfun(@(x) x.Jump, G(cidx{iHyp},:));
    
    % Get beliefs about the change point's position at the end of each
    % sequence
    bel = cellfun(@(x) x.CPbelief(:,end)', IO(cidx{iHyp},:), 'uni', 0);
    
    % Get inferred change position from the model and the subjects
    cpposmod = cellfun(@(x) x.Questions(3), IO(cidx{iHyp},:));
    cppossub = cellfun(@(x) x.Questions(3), G( cidx{iHyp},:));
    
    % Get sequences that were correctly labeled 
    detecmask = (filter{iHyp} == 1 | filter{iHyp} == 3);
    
    % Get the posterior distribution of change point's position separately
    % for (in)correctly labeled sequences
    for d = [0,1]
        
        % Restric analyses to particular sequence (either incorrectly or
        % correctly labeled)
        cbel = bel(detecmask == d);
        
        % Sort sequences by change point position
        [~,idx] = sort(cp{iHyp}(detecmask == d));
        pcppos{d+1,iHyp} = cell2mat(cbel(idx));
        
        % Same for most likely position
        mlpos = cpposmod(detecmask == d);
        cppos{1,iHyp} = mlpos(idx);
        mlpos = cppossub(detecmask == d);
        cppos{2,iHyp} = mlpos(idx);
    end
    
    % Get a distibution of change point's position restricted around true
    % change point's position
    winbel = cellfun(@(b,c) b(c+xwin), bel, num2cell(cp{iHyp}), 'uni', 0);
    
    % Remove posterior beliefs about change point's positions from
    % sequences that were mislabeled
    winbel(~detecmask) = {NaN(1,nwin)};
    winbel = cellfun(@(x) reshape(x,[1,1,nwin]), winbel, 'uni', 0);
    
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
for iHyp = 1:2
    subplot(1,2,iHyp); hold('on');
    
    % Average posterior distributions
    d = pcpposwrtcp(:,:,iHyp);
    m = mean(d, 1);
    s = sem(d, 1);
    
    % Display the posterior distribution of change point's positions
    % centered around the real position of the change point
    lgd(iHyp) = fill([xwin(1), xwin, xwin(end)], [0, m, 0], 'k', 'FaceColor', ...
        tricol(iHyp,:), 'EdgeColor', 'None', 'FaceAlpha', 1/7);
    plotMSEM(xwin, m, s, 1/2, tricol(iHyp,:), tricol(iHyp,:), 2, 1);
    
    % Display the position of the change point (i.e. 0)
    plot(zeros(1,2), ylim, 'k--');
    
    % Display the uniform scenario over change point's position
    plot(xwin([1,end]), repmat(1/N,1,2), '-', 'Color', g);
    
    % Customize the axes
    set(gca, 'XLim', xwin([1,end]), 'Box', 'Off');
    
    % Add some text labels
    xlabel('Position w.r.t. change point');
    ylabel('p(j_k|H_i,y)');
end

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
    colormap(sp, Emergence_Colormap(cmapcol{iHyp}));
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

% Compute a paired t-test on the precision of posterior distributions
% between the 2 types of regularities
[~,pval,tci,stats] = ttest(diff(precpcppos, 1, 2));
Emergence_PrintTstats(pval, tci, stats);
