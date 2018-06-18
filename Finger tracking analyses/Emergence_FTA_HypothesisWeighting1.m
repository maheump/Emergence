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
% regularities to the detection point
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Theshold in beliefs that define the detection points
detecthr = 1/2;

% Window in which to look in (in # of observations)
winwidth = 80;
xval = -winwidth/2:winwidth/2;

% Get the position of the change point
cp = cellfun(@(x) x.Jump-1/2, G(cidx{2},:));

% Get the position of the detection point
cdp = cellfun(@(x,c) c + find(x.BarycCoord(c:end,2) > detecthr, 1, ...
    'first') - 1, D(cidx{2},:), num2cell(cp), 'UniformOutput', 0); % cells

% Remove sequences with deterministic regularities that were not detected
detecmask = cellfun(@(x) x.Questions(2) == 2, G(cidx{2},:));
cdp(~detecmask) = {[]};

% Beginning and ending of the window to look in
idx = cellfun(@(x) x-winwidth/2:x+winwidth/2, cdp, 'UniformOutput', 0);

% Create a useful logical indexing that makes sure that the window around
% the detection point does not go beyond the observations' indices of the
% sequence
ok = cellfun(@(x) x >= 1 & x <= N, idx, 'UniformOutput', 0);

% Get the posterior beliefs over the 3 different hypotheses
beliefs = cellfun(@(x) x.BarycCoord, D(cidx{2},:), 'UniformOutput', 0);

% Select posterior beliefs around the detection point
subtraj = cellfun(@(x) NaN(winwidth+1,3), cell(size(idx)), 'UniformOutput', 0);
for iSeq = 1:numel(cdp)
    subtraj{iSeq}(ok{iSeq},:) = beliefs{iSeq}(idx{iSeq}(ok{iSeq}),:);
end

% Compute entropy of the transition probabilities from the rules
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Get the number of deterministic regularities
nR = numel(cidx{2});

% Get the list of patterns used in the experiment
numpat = cellfun(@(x) x.Rule, G(cidx{2},1), 'UniformOutput', 0); % numbers
strpat = cellfun(@(x) pat2str(x), numpat, 'UniformOutput', 0); % strings

% Get theoretical probabilities associated with those patterns
[pA, pAlt, pAgB, pBgA] = cellfun(@(x) pat2proba(x, [1 2], true), numpat);

% Get corresponding entropy levels
H = @(p) -(p .* log2(p) + (1-p) .* log2(1-p)); % entropy function
TPent = H(pAgB) .* H(pBgA); % entropy of transition probabilities
TPent(isnan(TPent)) = 0; % if not computable it means that one of p E {0,1}

% Define limits of the entropy bins
EntBin = [0, 0.5897, 0.7508, 1]';
nEnt = numel(EntBin)-1;

% Group deterministic regularities together according to the entropy bins
entidx = cell(1,nEnt);
for iEnt = 1:nEnt
    entidx{iEnt} = find(TPent >= EntBin(iEnt) & TPent < EntBin(iEnt+1));
end

% Get colormap for the entropy
entgrid = linspace(0, 1, 1001);
entcmap = inferno(numel(entgrid));

% Deduce colors corresponding to each entropy level
entcol = NaN(nR,3);
for iR = 1:nR
    [~,i] = min(abs(TPent(iR)-entgrid));
    entcol(iR,:) = entcmap(i,:);
end

% Deduce colors corresponding to the different entropy bins
entval = (EntBin(1:end-1)+EntBin(2:end)) ./ 2;
[~,colidx] = min(abs(entval'-entgrid'));
entgpcol = NaN(nEnt,3);
for iEnt = 1:nEnt
    entgpcol(iEnt,:) = entcmap(colidx(iEnt),:);
end

% Display the entropy map
% ~~~~~~~~~~~~~~~~~~~~~~~

% Prepare a new window
figure('Position', [1 905 290 200]);

% Compute 2D entropy map
entmap = H(entgrid') * H(entgrid);
entmap(isnan(entmap)) = 0;

% Get another 2D entropy map transformed according to the entropy bins
entmap2 = repmat({zeros(size(entmap))}, 1, 3);
for iEnt = 1:nEnt
    idx = entmap >= EntBin(iEnt) & entmap < EntBin(iEnt+1);
    for rgb = 1:3
        entmap2{rgb}(idx) = entgpcol(iEnt,rgb);
    end
end
entmap2 = cell2mat(reshape(entmap2, [1,1,3]));

% Display the binned 2D entropy map
image(entgrid, entgrid, entmap2); alpha(1/2); hold('on');

% Display diagonals
plot([0,1], [0,1], 'k-');
plot([0,1], [1,0], 'k-');

% Display each deterministic regularity
for iR = 1:nR
    plot(pAgB(iR), pBgA(iR), 'ko', 'MarkerFaceColor', entcol(iR,:), 'MarkerSize', 15);
    text(pAgB(iR), pBgA(iR), sprintf('   %s', strpat{iR}), ...
        'HorizontalAlignment', 'Left', 'VerticalAlignment', 'Top');
end

% Customize the colormap
colorbar('Location', 'WestOutside')
colormap(entcmap);

% Customize the axes
axis(repmat([0,1],1,2));
axis('square'); axis('xy');
set(gca, 'XTick', 0:0.2:1, 'YTick', 0:0.2:1);

% Add some text labels
xlabel('p(A|B)'); ylabel('p(B|A)');

% Save the figure
save2pdf('figs/F_HW_DeterRegTP.pdf');

%% SHOW THAT THE TRANSIENT BELIEFS DEPENDS UPON THE ENTROPY OF THE RULES
%  =====================================================================

% Average trajectories around detection point for each entropy bin
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Prepare an output variable
subavgtrajent = NaN(winwidth+1,3,nEnt,nSub);

% For each entropy bin, average the trajectories across subjects
for iEnt = 1:nEnt
    n = numel(entidx{iEnt});
    d = mat2cell(subtraj(entidx{iEnt},:), n, ones(nSub,1));
    d = cellfun(@(x) mean(cell2mat(reshape(x, [1,1,n])), 3, 'OmitNaN'), d, 'UniformOutput', 0);
    subavgtrajent(:,:,iEnt,:) = cell2mat(reshape(d, [1,1,nSub]));
end

% Display group-averaged trajectories
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Prepare a new window
figure('Position', [292 905 240 200]);

% Display the triangle
Emergence_PlotTriInfo(tricc, tricol);

% For each entropy bin
for iEnt = 1:nEnt

    % Average over subjects and convert barycentric coordinates to cartesian
    d = subavgtrajent(:,:,iEnt,:);
    m = mean(d, 4, 'OmitNaN');
    ccM = m * tricc;

    % Display the average trajectory
    col = entcmap(colidx(iEnt),:);
    plot(ccM(:,1), ccM(:,2), '.-', 'Color', col, ...
        'MarkerSize', 10); hold('on');

    % Display the detection point
    i = xval == 0;
    plot(ccM(i,1), ccM(i,2), 'o', 'Color', col, 'MarkerSize', 10);
end

% Display relevant triangular grid
Emergence_PlotGridOnTri(10, 1, tricol(1,:), tricc);
Emergence_PlotGridOnTri(2,  2, tricol(2,:), tricc);

% Make sure the triangle is equilateral and the axes invisible
axis('equal'); axis('off');

% Save the figure
if isfield(D{1}, 'Seq'), save2pdf('figs/F_HW_TriS.pdf');
else, save2pdf('figs/F_HW_TriIO.pdf');
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
    data(:,iEnt) = squeeze(mean(transpbel(entidx{iEnt},:), 1, 'OmitNaN'));
end

% Run an ANOVA
RMtbl = rmANOVA(data', 'SeqType');
disp(RMtbl);

% Display group-averaged beliefs in the probabilistic hypothesis
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Prepare a new window
figure('Position', [533 905 120 200]);

% Display the dispersion of of beliefs across subjects
Emergence_PlotSubGp(data, entcmap(colidx,:));

% Customize the axes
xlim([0,nEnt+1]); ylim([0, max(ylim)]);
set(gca, 'Box', 'Off');
set(gca, 'XTick', 1:nEnt, 'XTickLabel', {'Low', 'Med', 'High'});

% Display whether the difference is significant or not
Emergence_DispStatTest(data);

% Save the figure
if isfield(D{1}, 'Seq'), save2pdf('figs/F_HW_EntS.pdf');
else, save2pdf('figs/F_HW_EntIO.pdf');
end
