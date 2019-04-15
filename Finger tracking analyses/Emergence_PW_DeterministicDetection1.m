% The analysis implemented in this script focuses on the trajectories
% produced just before the detection of (effectively detected)
% deterministic regularities. In particular, we study the transient
% increase of beliefs in the probabilistic hypothesis that tend to occur in
% these periods. We show that the amount of such an increase in
% probabilistic beliefs is determined by the entropy of transition
% probabilities computed on the deterministic rule that is about to be
% detected. This suggests that subjects perform genuine hypothesis
% weighting (between probabilistic and deterministic regularities) and do
% not produce a trivial discrimination behviour betwen probabilistic and
% deterministic regularities.
%
% Copyright (c) 2018 Maxime Maheu

%% SHOW THAT DETERMINISTIC REGULARITIES ARE (ALSO) CHARACTERIZED BY DIFFERENT ENTROPY LEVELS
%  =========================================================================================

% Create a 2D entropy map
% ~~~~~~~~~~~~~~~~~~~~~~~

% Create colormap for the entropy
minH = 1.4;
maxH = Emergence_MarkovEntropy(1/2, 1/2);
prec = 101;
offset = round(prec * (maxH - minH));
EntCMap = flipud([flipud(Emergence_Colormap('Purples', offset)); Emergence_Colormap('Greys', prec)]);
prec = size(EntCMap,1);

% Compute 2D entropy map
ProbaGrid = linspace(0, 1, prec);
EntMap = arrayfun(@(x,y) Emergence_MarkovEntropy(x,y), ...
    repmat(ProbaGrid, [prec,1]), repmat(ProbaGrid', [1,prec]));

% Compute entropy of transition probabilities from the patterns
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Get the number of deterministic regularities
nR = numel(cidx{2});

% Get the list of patterns used in the experiment
numpat = cellfun(@(xval) xval.Rule, G(cidx{2},1), 'uni', 0); % numbers
strpat = cellfun(@(xval) pat2str(xval), numpat, 'uni', 0);   % strings

% Get theoretical probabilities associated with those patterns
[pA, pAlt, pAgB, pBgA] = cellfun(@(xval) pat2proba(xval, [1 2], true), numpat);

% Get corresponding entropy levels
TPent = arrayfun(@(x,y) Emergence_MarkovEntropy(x, y), pAgB, pBgA);

% Get colors corresponding to each entropy level
EntGrid = linspace(0, max(EntMap(:)), prec);
[~,colidx] = min(abs(TPent - EntGrid), [], 2);
EntCol = EntCMap(colidx,:);

% Define bins
% ~~~~~~~~~~~

% Group deterministic regularities together according to the entropy bins
EntBin = [1.4, 1.7665, 1.853, 2]';
EntLab = arrayfun(@(x) sprintf('%1.0f', x), 1:numel(EntBin), 'uni', 0);
nEnt = numel(EntBin)-1;
[~,~,EntIdx] = histcounts(TPent, EntBin);

% Deduce colors corresponding to the different entropy bins
EntVal = arrayfun(@(i) mean(TPent(EntIdx == i)), 1:nEnt)';
[~,colidx] = min(abs(EntVal - EntGrid), [], 2);
EntGpCol = EntCMap(colidx,:);

% Display the entropy map and the position of rules on that map
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Prepare a new window
figure('Position', [1 905 290 200]);

% Display the binned entropy map
imagesc(ProbaGrid, ProbaGrid, EntMap); hold('on');

% Display diagonals
plot([0,1], [0,1], 'k-', 'LineWidth', 1/2); hold('on');
plot([0,1], [1,0], 'k-', 'LineWidth', 1/2);

% Display limits of entropy bins
[~,c] = contour(ProbaGrid, ProbaGrid, EntMap, unique(TPent), ...
    'k-','ShowText', 'Off', 'LabelSpacing', 500, 'LineWidth', 1/4);

% Display each deterministic regularity
scatter(pAgB, pBgA, 100, EntCol, 'filled', 'MarkerEdgeColor', 'k');
text(pAgB, pBgA, cellfun(@(x) sprintf('  %s', x), strpat, 'uni', 0), ...
   'HorizontalAlignment', 'Left', 'VerticalAlignment', 'Top');

% Customize the colormap
colormap(EntCMap);
colorbar;
caxis([0,maxH]);

% Customize the axes
axis(repmat([0,1],1,2));
axis('square'); axis('xy');
set(gca, 'XTick', 0:0.2:1, 'YTick', 0:0.2:1);

% Add some text labels
xlabel('p(A|B)'); ylabel('p(B|A)');

% Save the figure
save2pdf(fullfile(ftapath, 'figs', 'F_HW_DeterRegTP.pdf'));

%% SHOW THE TRANSIENT TENSION BETWEEN REGULAR HYPOTHESES
%  =====================================================

% Find important points
% ~~~~~~~~~~~~~~~~~~~~~

% Find positions of change and detection points
cp = cellfun(@(x) x.Jump+1/2, G(cidx{2},:));
lag = cellfun(@(x,c) Emergence_FindDetecPoint(x.BarycCoord(c:end,2)), D(cidx{2},:), num2cell(cp));
dp = cp + lag;

% Restrict to sequences that were accurately classified by subjects and
% with a regular detection point for both the subjects and the ideal
% observer model (it is true in all cases for the ideal observer model in
% the case of deterministic regularities)
detecmask = (filter{2} == 3);
cp (~detecmask) = NaN;
dp (~detecmask) = NaN;
lag(~detecmask) = NaN;

% Align observers' trajectories in sequences with deterministic
% regularities to the detection points
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Window in which to look in (in # of observations)
xval = -130:60;
nsp = numel(xval);

% Beginning and ending of the window to look into
idx = arrayfun(@(x) x+xval, dp, 'uni', 0);

% Create a useful logical indexing that makes sure that the window around
% the detection point does not go beyond the observations' indices of the
% sequence
ok = cellfun(@(x) x >= 1 & x <= N, idx, 'uni', 0);

% Get the posterior beliefs over the 3 different hypotheses
beliefs = cellfun(@(x) x.BarycCoord, D(cidx{2},:), 'uni', 0);

% Select posterior beliefs around the detection point
subtraj = cellfun(@(x) NaN(nsp,3), cell(size(idx)), 'uni', 0);
for iSeq = 1:numel(idx)
    subtraj{iSeq}(ok{iSeq},:) = beliefs{iSeq}(idx{iSeq}(ok{iSeq}),:);
end

% Average trajectories over sequences and subjects
avgtraj = squeeze(mean(cell2mat(cellfun(@(x) reshape(x, [1 1 size(x)]), ...
    subtraj, 'uni', 0)), [1,2], 'OmitNaN'));
cc = avgtraj*tricc;

% Show the trajectory in the triangle locked on the detection of
% deterministic regularities
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% For each predictor
figure('Position', [292 905 300 200]);

% Display information about the triangle
Emergence_PlotTriInfo; hold('on');
l = Emergence_PlotGridOnTri(2, 2); set(l(1,2), 'LineWidth', 2);

% Display null ratio
plot(ones(1,2)./2, [0, sqrt(3)/2], 'k--');

% Display the average trajectory
plot(cc(:,1), cc(:,2), 'k-', 'LineWidth', 2);

% Customize the axes
axis([0, 1, 0, sqrt(3)/2]);
axis('xy'); axis('off'); axis('equal');

% Save the figure
if isfield(D{1}, 'Seq'), save2pdf(fullfile(ftapath, 'figs', 'F_HW_AvgTriS.pdf'));
else, save2pdf(fullfile(ftapath, 'figs', 'F_HW_AvgTriIO.pdf'));
end

% Display the P/D ratio locked on the detection of deterministic regularities
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Copute the P/D ratio as (p(Hd|y)-p(Hp|y)) / (p(Hd|y)+p(Hp|y))
ratioPD = cell2mat(cellfun(@(pHgY) reshape(...
    (pHgY(:,2) - pHgY(:,1)) ./ (pHgY(:,2) + pHgY(:,1)), ...
    [1 1 nsp]), subtraj, 'uni', 0));

% Average over sequences for each subject
ratioPD = squeeze(mean(ratioPD, 1, 'OmitNaN'));

% Average over subjects
avg = mean(ratioPD, 1);
err = sem( ratioPD, 1);

% Create a new window
figure('Position', [593 905 190 200]);

% Display origin
plot(xval([1,end]), zeros(1,2), '-', 'Color', g); hold('on');
plot(zeros(1,2), [-1,1], '-', 'Color', g);

% Display the P/D ratio centered on the detection point
plotMSEM(xval, avg, err, 1/5, 'k', 'k', 2);

% Display distribution of change point position
fout0 = ksdensity(-lag(:), ...              % which distribution to plot
    xval, ...                               % grid of positions
    'BandWidth',            8, ...          % bandwidth of the kernel smoothing window 
    'Support',              [-Inf,1], ...   % restrict the kernel to a certain range of values
    'BoundaryCorrection',   'Reflection');	% type of correction for the boundaries
fill([xval(1),xval,xval(end)], [0,fout0.*10,0]-1, 'k', 'FaceColor', g);

% Customize the axes
caxis([-1,1]);
axis([xval(1), xval(end), -1, 1]);
set(gca, 'Box', 'Off');

% Add text labels
xlabel('Observation w.r.t. detection point');
ylabel('Ratio P versus D');

% Save the figure
if isfield(D{1}, 'Seq'), save2pdf(fullfile(ftapath, 'figs', 'F_HW_PDratioS.pdf'));
else, save2pdf(fullfile(ftapath, 'figs', 'F_HW_PDratioIO.pdf'));
end

%% SHOW THAT THE TRANSIENT TENSION DEPENDS UPON THE ENTROPY OF THE RULES
%  =====================================================================

% Average trajectories around detection point for each entropy bin
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Prepare an output variable
subavgtrajent = NaN(nsp,3,nEnt,nSub);

% For each entropy bin, average the trajectories across subjects
for iEnt = 1:nEnt
    n = sum(EntIdx == iEnt);
    d = mat2cell(subtraj(EntIdx == iEnt,:), n, ones(nSub,1));
    d = cellfun(@(x) mean(cell2mat(reshape(x, [1,1,n])), 3, 'OmitNaN'), ...
        d, 'uni', 0);
    subavgtrajent(:,:,iEnt,:) = cell2mat(reshape(d, [1,1,nSub]));
end

% Average over subjects
gdavgtrajent = mean(subavgtrajent, 4, 'OmitNaN');
gdsemtrajent = sem(subavgtrajent, 4);

% Display group-averaged trajectories
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Prepare a new window
figure('Position', [784 905 240 200]);

% Display the triangle
Emergence_PlotTriInfo(tricc, tricol);

% For each entropy bin
for iEnt = 1:nEnt
    
    % Average over subjects and convert barycentric coordinates to cartesian
    ccM = gdavgtrajent(:,:,iEnt) * tricc;
    
    % Display the average trajectory
    col = EntCMap(colidx(iEnt),:);
    plot(ccM(:,1), ccM(:,2), '.-', 'Color', col, ...
        'MarkerSize', 10); hold('on');
end

% Display relevant triangular grid
Emergence_PlotGridOnTri(10, 1, tricol(1,:), tricc);
l = Emergence_PlotGridOnTri(2, 2); set(l(1,2), 'LineWidth', 2);

% Make sure the triangle is equilateral and the axes invisible
axis('equal'); axis('off');

% Save the figure
if isfield(D{1}, 'Seq'), save2pdf(fullfile(ftapath, 'figs', 'F_HW_TriS.pdf'));
else, save2pdf(fullfile(ftapath, 'figs', 'F_HW_TriIO.pdf'));
end

% Display corresponding barycentric coordinates
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Prepare a new window
figure('Position', [1025 905 240 200]);

% Draw some help lines
plot(xval([1,end]), ones(1,2)./2, '-',  'Color', g, 'LineWidth', 1/2); hold('on');
plot(zeros(1,2), [0,1], 'k');

% For each entropy bin
for iEnt = 1:nEnt
    col = EntCMap(colidx(iEnt),:);
    plotMSEM(xval, gdavgtrajent(:,1,iEnt), gdsemtrajent(:,1,iEnt), ...
        1/10, col, col, 3, 1, '-', 'none');
    plotMSEM(xval, gdavgtrajent(:,2,iEnt), gdsemtrajent(:,1,iEnt), ...
        1/10, col, col, 1, 1, '--', 'none');
end 

% Customize the axes
axis([xval([1,end]),0,1]);
set(gca, 'Box', 'Off');

% Add some text labels
xlabel('Observation w.r.t. detection point');
ylabel('Posterior beliefs p(M_i|y)');

% Save the figure
if isfield(D{1}, 'Seq'), save2pdf(fullfile(ftapath, 'figs', 'F_HW_DynS.pdf'));
else, save2pdf(fullfile(ftapath, 'figs', 'F_HW_DynIO.pdf'));
end

% Average beliefs in the probabilistic hypothesis between the change and
% the detection points
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Get the average beliefs in the probabilistic hypothesis between the
% change and the detection points
celldp = num2cell(dp);
celldp(cellfun(@isnan, celldp)) = {[]};
cellcp = num2cell(cp);
cellcp(cellfun(@isnan, cellcp)) = {[]};
transpbel = cellfun(@(p,c,d) ...
    (p.BarycCoord(c:d,2) - p.BarycCoord(c:d,1)) ./ (p.BarycCoord(c:d,2) + p.BarycCoord(c:d,1)), ...
    D(cidx{2},:), cellcp, celldp, 'uni', 0);
transpbel = cellfun(@mean, transpbel);

% For each subject, average the trajectories
data = NaN(nSub,nEnt);
for iEnt = 1:nEnt
    data(:,iEnt) = squeeze(mean(transpbel(EntIdx == iEnt,:), 1, 'OmitNaN'));
end

% Run an ANOVA
RMtbl = rmANOVA(data', 'SeqType');
Emergence_PrintFstats(RMtbl);

% Test whether there is significant correlations between beliefs in the
% probabilistic hypothesis and Shannon entropy
corcoef = cellfun(@(x) Emergence_Regress(x, TPent, 'CC', 'r'), ...
    mat2cell(transpbel, numel(cidx{2}), ones(1,nSub)));
[~,pval,tci,stats] = ttest(corcoef');
Emergence_PrintTstats(pval,tci,stats);

% Display group-averaged beliefs in the probabilistic hypothesis
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Prepare a new window
figure('Position', [1266 905 120 200]);

% Display the dispersion of of beliefs across subjects
Emergence_PlotSubGp(data, EntCMap(colidx,:));

% Customize the axes
xlim([0,nEnt+1]);
set(gca, 'XTick', 1:nEnt, 'XTickLabel', EntLab, 'Box', 'Off');

% Add some text labels
ylabel('Posterior beliefs p(M_P|y)');

% Save the figure
if isfield(D{1}, 'Seq'), save2pdf(fullfile(ftapath, 'figs', 'F_HW_EntS.pdf'));
else, save2pdf(fullfile(ftapath, 'figs', 'F_HW_EntIO.pdf'));
end
