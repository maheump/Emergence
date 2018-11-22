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

% Compute entropy of the transition probabilities from the rules
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Get the number of deterministic regularities
nR = numel(cidx{2});

% Get the list of patterns used in the experiment
numpat = cellfun(@(xval) xval.Rule, G(cidx{2},1), 'UniformOutput', 0); % numbers
strpat = cellfun(@(xval) pat2str(xval), numpat, 'UniformOutput', 0);   % strings

% Get theoretical probabilities associated with those patterns
[pA, pAlt, pAgB, pBgA] = cellfun(@(xval) pat2proba(xval, [1 2], true), numpat);

% Get corresponding entropy levels
TPent = arrayfun(@Emergence_IO_Entropy, pAgB) + ...
        arrayfun(@Emergence_IO_Entropy, pBgA);

% Define limits of the entropy bins
EntLab = {'Low', 'Med', 'High'};
EntBin = [1, 1.53, 1.73, 2]';
nEnt = numel(EntBin)-1;

% Group deterministic regularities together according to the entropy bins
EntIdx = cell(1,nEnt);
for iEnt = 1:nEnt
    EntIdx{iEnt} = find(TPent >= EntBin(iEnt) & TPent <= EntBin(iEnt+1));
end

% Get colormap for the entropy
EntGrid = linspace(0, 1, 111);
EntCMap = flipud(cbrewer2('RdPu', numel(EntGrid)));
EntCMap = EntCMap(1:101,:);
EntGrid = linspace(0, 1, 101);

% Compute 2D entropy map
EntMap = arrayfun(@Emergence_IO_Entropy, EntGrid);
EntMap = EntMap' + EntMap;
EntMap(EntMap < 1) = NaN;

% Deduce colors corresponding to each entropy level
[~,colidx] = min(abs(TPent - 1 - EntGrid), [], 2);
EntCol = EntCMap(colidx,:);

% Deduce colors corresponding to the different entropy bins
EntVal = mean([EntBin(1:end-1), EntBin(2:end)], 2);
[~,colidx] = min(abs(EntVal - 1 - EntGrid), [], 2);
EntGpCol = EntCMap(colidx,:);

% Colorbar for group of entropies
EntGpColBar = NaN(numel(EntGrid),3);
for iEnt = 1:nEnt
    idx = EntGrid > EntBin(iEnt)-1 & EntGrid <= EntBin(iEnt+1)-1;
    EntGpColBar(idx,:) = repmat(EntGpCol(iEnt,:), [sum(idx), 1]);
end

% Display the entropy map and the position of rules on that map
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Prepare a new window
figure('Position', [1 905 290 200]);

% Display bi-column colorbar
subplot(1,3,1);
image(1:2, 1+EntGrid, [reshape(EntCMap, [numel(EntGrid), 1, 3]), ...
                       reshape(EntGpColBar, [numel(EntGrid), 1, 3])]); hold('on');
plot(1/2+ones(1,2), [1,2], 'k-');
plot([0,1.5], repmat(EntBin(2:end-1), 1, 2)', 'k:');
plot([1.5,3], repmat(EntBin(2:end-1), 1, 2)', 'k-');
text(1+ones(1,nEnt), EntVal, EntLab, 'Rotation', 90, 'Color', 'w');
axis('xy'); set(gca, 'XTick', []); ylabel('Entropy (bits)');

% Display the binned 2D entropy map
subplot(1,3,2:3);
imagesc(EntGrid, EntGrid, EntMap, 'AlphaData', ~isnan(EntMap)); hold('on');
caxis([1,2]);

% Display diagonals
plot([0,1], [0,1], 'k-', 'LineWidth', 1/2);
plot([0,1], [1,0], 'k-', 'LineWidth', 1/2);

% Display limits of entropy bins
for iEnt = 1:nEnt
    contour(EntGrid, EntGrid, EntMap >= EntBin(iEnt), 1, 'k:', 'LineWidth', 1);
end

% Display each deterministic regularity
for iR = 1:nR
    plot(pAgB(iR), pBgA(iR), 'ko', 'MarkerFaceColor', EntCol(iR,:), 'MarkerSize', 8);
    text(pAgB(iR), pBgA(iR), sprintf('  %s', strpat{iR}), ...
        'HorizontalAlignment', 'Left', 'VerticalAlignment', 'Top');
end

% Customize the colormap
colormap(EntCMap);

% Customize the axes
axis(repmat([0,1],1,2));
axis('square'); axis('xy');
set(gca, 'Layer', 'Bottom')
set(gca, 'XTick', 0:0.2:1, 'YTick', 0:0.2:1);

% Add some text labels
xlabel('p(A|B)'); ylabel('p(B|A)');

% Save the figure
save2pdf(fullfile(ftapath, 'figs', 'F_HW_DeterRegTP.pdf'));

%% SHOW THE TRANSIENT TENSION BETWEEN REGULAR HYPOTHESES
%  =====================================================

% Align observers' trajectories in sequences with deterministic
% regularities to the detection points
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Theshold in beliefs that define the detection points
detecthr = 1/2;

% Window in which to look in (in # of observations)
xval = -130:60;
nsp = numel(xval);

% Get the detection point's position
% (note that we go back in time from the end of the sequence in order
% to avoid as much as possible the false alarms that might exist at the
% time of the change point)
cdp = cellfun(@(x) min([NaN, ...
    N - find(flipud(x.BarycCoord( x.Jump + 1/2: end, 2)) ...
    >= detecthr, 1, 'last')], [], 'OmitNaN'), ... % in # of observations
    D(cidx{2},:), 'UniformOutput', 0); % cells
    
% Only detected regularities
detecmask = cellfun(@(x) x.Questions(2) == 2, G(cidx{2},:));
cdp(~detecmask) = {NaN};

% Beginning and ending of the window to look into
idx = cellfun(@(x) x+xval, cdp, 'UniformOutput', 0);

% Create a useful logical indexing that makes sure that the window around
% the detection point does not go beyond the observations' indices of the
% sequence
ok = cellfun(@(x) x >= 1 & x <= N, idx, 'UniformOutput', 0);

% Get the posterior beliefs over the 3 different hypotheses
beliefs = cellfun(@(x) x.BarycCoord, D(cidx{2},:), 'UniformOutput', 0);

% Select posterior beliefs around the detection point
subtraj = cellfun(@(x) NaN(nsp,3), cell(size(idx)), 'UniformOutput', 0);
for iSeq = 1:numel(idx)
    subtraj{iSeq}(ok{iSeq},:) = beliefs{iSeq}(idx{iSeq}(ok{iSeq}),:);
end

% Average trajectories over sequences and subjects
avgtraj = squeeze(mean(cell2mat(cellfun(@(x) reshape(x, [1 1 size(x)]), ...
    subtraj, 'UniformOutput', 0)), [1,2], 'OmitNaN'));
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

% Copute the P/D ratio
ratioPD = cell2mat(cellfun(@(x) reshape(x(:,2) ./ sum(x(:,1:2), 2), ...
    [1 1 nsp]), subtraj, 'UniformOutput', 0));

% Transform it such that it evolves between -1 (Proba) and 1 (Deter)
ratioPD = 2 .* (ratioPD - 1/2);

% Average over sequences for each subject
ratioPD = squeeze(mean(ratioPD, 1, 'OmitNaN'));

% Average over subjects
avg = mean(ratioPD, 1);
err = sem( ratioPD, 1);

% Create a new window
figure('Position', [593 905 190 200]);

% Display origin
plot(xval([1,end]), zeros(1,2), 'k-'); hold('on');
plot(zeros(1,2), [-1,1], 'k-');

% Display the P/D ratio centered on the detection point
plotMSEM(xval, avg, err, 1/5, 'k', 'k', 2);

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
    n = numel(EntIdx{iEnt});
    d = mat2cell(subtraj(EntIdx{iEnt},:), n, ones(nSub,1));
    d = cellfun(@(x) mean(cell2mat(reshape(x, [1,1,n])), 3, 'OmitNaN'), ...
        d, 'UniformOutput', 0);
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
cdp(cellfun(@isnan, cdp)) = {[]};
transpbel = cellfun(@(x,d) mean(x.BarycCoord((x.Jump+1/2):d,1)), ...
    D(cidx{2},:), cdp);

% Remove sequences with deterministic regularities that were not detected
transpbel(~detecmask) = NaN;

% For each subject, average the trajectories
data = NaN(nSub,nEnt);
for iEnt = 1:nEnt
    data(:,iEnt) = squeeze(mean(transpbel(EntIdx{iEnt},:), 1, 'OmitNaN'));
end

% Run an ANOVA
RMtbl = rmANOVA(data', 'SeqType');
disp(RMtbl);

% Display group-averaged beliefs in the probabilistic hypothesis
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Prepare a new window
figure('Position', [1266 905 120 200]);

% Display the dispersion of of beliefs across subjects
Emergence_PlotSubGp(data, EntCMap(colidx,:));

% Customize the axes
xlim([0,nEnt+1]); ylim([0, max(ylim)]);
set(gca, 'Box', 'Off');
set(gca, 'XTick', 1:nEnt, 'XTickLabel', EntLab);

% Display whether the difference is significant or not
Emergence_DispStatTest(data);

% Add some text labels
ylabel('Posterior beliefs p(M_P|y)');

% Save the figure
if isfield(D{1}, 'Seq'), save2pdf(fullfile(ftapath, 'figs', 'F_HW_EntS.pdf'));
else, save2pdf(fullfile(ftapath, 'figs', 'F_HW_EntIO.pdf'));
end
