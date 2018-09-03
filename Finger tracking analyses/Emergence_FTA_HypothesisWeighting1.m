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

% Align observers' trajectories in sequences with deterministic
% regularities to the detection points
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Theshold in beliefs that define the detection points
detecthr = 1/2;

% Window in which to look in (in # of observations)
winwidth = 80;
xval = -winwidth/2:winwidth/2;

% Get the position of the change point
cp = cellfun(@(xval) xval.Jump-1/2, G(cidx{2},:));

% Get the position of the detection point
cdp = cellfun(@(xval,c) c + find(xval.BarycCoord(c:end,2) > detecthr, 1, ...
    'first') - 1, D(cidx{2},:), num2cell(cp), 'UniformOutput', 0); % cells

% Remove sequences with deterministic regularities that were not detected
detecmask = cellfun(@(xval) xval.Questions(2) == 2, G(cidx{2},:));
cdp(~detecmask) = {[]};

% Beginning and ending of the window to look into
idx = cellfun(@(xval) xval-winwidth/2:xval+winwidth/2, cdp, 'UniformOutput', 0);

% Create a useful logical indexing that makes sure that the window around
% the detection point does not go beyond the observations' indices of the
% sequence
ok = cellfun(@(xval) xval >= 1 & xval <= N, idx, 'UniformOutput', 0);

% Get the posterior beliefs over the 3 different hypotheses
beliefs = cellfun(@(xval) xval.BarycCoord, D(cidx{2},:), 'UniformOutput', 0);

% Select posterior beliefs around the detection point
subtraj = cellfun(@(xval) NaN(winwidth+1,3), cell(size(idx)), 'UniformOutput', 0);
for iSeq = 1:numel(cdp)
    subtraj{iSeq}(ok{iSeq},:) = beliefs{iSeq}(idx{iSeq}(ok{iSeq}),:);
end

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
EntGrid = linspace(0, 1, 1001);
EntCMap = hot(numel(EntGrid));

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

%% SHOW THAT THE TRANSIENT BELIEFS DEPENDS UPON THE ENTROPY OF THE RULES
%  =====================================================================

% Average trajectories around detection point for each entropy bin
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Prepare an output variable
subavgtrajent = NaN(winwidth+1,3,nEnt,nSub);

% For each entropy bin, average the trajectories across subjects
for iEnt = 1:nEnt
    n = numel(EntIdx{iEnt});
    d = mat2cell(subtraj(EntIdx{iEnt},:), n, ones(nSub,1));
    d = cellfun(@(xval) mean(cell2mat(reshape(xval, [1,1,n])), 3, 'OmitNaN'), d, 'UniformOutput', 0);
    subavgtrajent(:,:,iEnt,:) = cell2mat(reshape(d, [1,1,nSub]));
end

% Average over subjects
gdavgtrajent = mean(subavgtrajent, 4, 'OmitNaN');
gdsemtrajent = sem(subavgtrajent, 4);

% Display group-averaged trajectories
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Prepare a new window
figure('Position', [292 905 240 200]);

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
Emergence_PlotGridOnTri(2,  2, tricol(2,:), tricc);

% Make sure the triangle is equilateral and the axes invisible
axis('equal'); axis('off');

% Save the figure
if isfield(D{1}, 'Seq'), save2pdf(fullfile(ftapath, 'figs', 'F_HW_TriS.pdf'));
else, save2pdf(fullfile(ftapath, 'figs', 'F_HW_TriIO.pdf'));
end

% Display corresponding barycentric coordinates
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Prepare a new window
figure('Position', [533 905 240 200]);

% Draw some help lines
plot(xval([1,end]),    ones(1,2)./2, '-',  'Color', g, 'LineWidth', 1/2); hold('on');
plot(xval([1,end]),    ones(1,2)./3, '--', 'Color', g, 'LineWidth', 1/2);
plot(xval([1,end]), 2.*ones(1,2)./3, '--', 'Color', g, 'LineWidth', 1/2);
plot(zeros(1,2), [0,1], 'k');

% For each entropy bin
for iEnt = 1:nEnt
    col = EntCMap(colidx(iEnt),:);
    plotMSEM(xval, gdavgtrajent(:,1,iEnt), gdsemtrajent(:,1,iEnt), 1/10, col, col, 3, 1, '-', 'none');
    plotMSEM(xval, gdavgtrajent(:,2,iEnt), gdsemtrajent(:,1,iEnt), 1/10, col, col, 1, 1, '--', 'none');
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
transpbel = cellfun(@(b,c,d) mean(b.BarycCoord(c:d,1)), ...
    D(cidx{2},:), num2cell(cp), cdp);

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
figure('Position', [774 905 120 200]);

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
