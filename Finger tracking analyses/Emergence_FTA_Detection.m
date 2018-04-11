% This script implements analyses looking at the inferred generative
% process of the sequence once the sequence has been entirely presented
% (post-sequence questions). First, we show that the correct labeling of
% sequences is better for those that entail a deterministic compared to
% a probabilistic regularity. This is then confirmed using a signal
% detection theory approach in which possible (and existing) response
% biases are taken into account. Second, we show that probabilistic
% regularity that are missed by subjects are associated with lower beliefs
% from the ideal observer. Finally, we investigate the determinants of
% missed and accurately detected regularities. In particuler, we show that,
% for probabilistic regularities, it depends on the entropy level and the
% type of bias (alternation are more often biased) while, for deterministic
% regularities, it depends on the length of the repeating rule and its
% inner organisation.
% 
% Copyright (c) 2018 Maxime Maheu

%% LABELING THE 3 TYPES OF SEQUENCES
%  =================================

% Get the proportions of each label that has been atrributed to each
% sequence type
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Get the inferred generative process once the sequence is over
estgenproc = cellfun(@(x) x.Questions(2), D);
estgenproc(isnan(estgenproc)) = 3;

% Average over sequences within each condition
avgcordet = NaN(3,3,nSub);
for iHyp = 1:3 % generative process
    data = mat2cell(estgenproc(cidx{iHyp},:), numel(cidx{iHyp}), ones(nSub,1));
    avgcordet(iHyp,:,:) = cell2mat(cellfun(@(x) histc(x, 1:3) ...
        ./ numel(x), data, 'UniformOutput', 0));
end

% Run a 1-way ANOVA on correct labeling between conditions
data = squeeze(cat(1, avgcordet(1,1,:), avgcordet(2,2,:), avgcordet(3,3,:)))';
RMtbl = rmANOVA(data, 'SeqType');
disp(RMtbl);

% Display results using a home-made stacked barplot
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Prepare the window
figure('Position', [1 806 120 300]); lgd = NaN(1,3);

% Display chance level
xpos = [2,3,1];
plot([min(xpos)-1,max(xpos)+1], repmat(1/3,1,2), '-', ...
    'Color', g, 'LineWidth', 1/2); hold('on');

% For each type of sequence
for iHyp = 1:3
    
    % Average over subjects
    m = mean(avgcordet(iHyp,:,:), 3);
    s = sem(avgcordet(iHyp,:,:), 3);
    
    % Display correct labeling
    lgd(iHyp) = bar(xpos(iHyp), m(iHyp), 0.8, 'FaceColor', tricol(iHyp,:));
    
    % Display error bar for correct reports
    plot(repmat(xpos(iHyp),1,2), m(iHyp)+s(iHyp).*[-1,1], 'k-');
    
    % Display incorrect labeling
    idx = setdiff(1:3, iHyp); % other responses than the real one
    [~,nd] = sort(m(idx), 2, 'descend');
    idx = idx(nd);
    h = bar([repmat(xpos(iHyp),1,2); repmat(xpos(iHyp)-1,1,2)], ...
        [-m(idx); NaN(1,2)], 0.8, 'Stacked');
    for i = 1:2, set(h(i), 'FaceColor', tricol(idx(i),:)); end
    
    % Display error bars for wrong reports
    cm = cumsum(-m(idx));
    for j = 1:2
        plot(repmat(xpos(iHyp),1,2), cm(j)+s(j).*[-1,1], 'k-');
    end
end

% Customize the axes
axis([min(xpos)-1,max(xpos)+1,-0.5,1]);
set(gca, 'XTick', 1:3, 'XTickLabel', proclab, 'XTickLabelRotation', 30);
set(gca, 'YTick', -0.4:0.2:1, 'Box', 'Off');

% Add some text labels
xlabel('Type of sequence');
ylabel('Average report rates');
legend(lgd, proclab, 'Location', 'SouthOutside');

% Save the figure
if isfield(D{1}, 'Seq'), save2pdf('figs/F_D_DetS.pdf');
else, save2pdf('figs/F_D_DetIO.pdf');
end

% Display labeling of each sequence and each subject
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Prepare the figure
figure('Position', [122 806 230 300]);

% For each category of regularity
for iHyp = 1:3
    subplot(3,1,iHyp);
    
    % Number of sequences with this type of regularity
    nR = numel(cidx{iHyp});
    
    % Display response
    pcolor(1:nSub, 1:nR, estgenproc(cidx{iHyp},:));
	
    % Customize the axes
    colormap(tricol);
    axis('equal'); axis('off');
    axis([1, nSub, 1, nR]);
    set(gca, 'YDir', 'Reverse', 'Box', 'Off');

    % Add some text labels
    text(mean(1:nSub), nR+1, 'Subjects');
    text(0, nR/2, 'Sequences', 'Rotation', 90);
    title(sprintf('%s sequences', proclab{iHyp}));
end

% Save the figure
if isfield(D{1}, 'Seq'), save2pdf('figs/F_D_SeqS.pdf');
else, save2pdf('figs/F_D_SeqIO.pdf');
end

%% MULTI-DIMENSIONAL SIGNAL DETECTION THEORY
%  =========================================

% Computed paired signal detection theory quantities
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Subjects' jugment of whether a regularity was observed in each sequence
regornot = cellfun(@(x) x.Questions(1), D);
regornot(regornot == 2) = 0;

% Prepare output variables
dprime = NaN(nSub,3);
crit   = NaN(nSub,3);

% Define the padding value for non-computable d' values
pad = N/(2*N);

% For each subject
for sub = 1:nSub
    
    % 1) Fully-stochastic versus one of the regularity
    
    % For each type of regularity
    for iHyp = 1:2
    
        % Get data
        signal = regornot(cidx{iHyp},sub); % sequences with a regularity
        noise  = regornot(cidx{3},sub);    % fully-stochastic sequences

        % Get frequencies for each possible answer
        H  = sum(signal == 1); % hit
        M  = sum(signal == 0); % miss
        CR = sum(noise  == 0); % correct rejection
        FA = sum(noise  == 1); % false alarm
        
        % Create a contingency table and pad unobserved values
        t = [H, M; FA, CR];
        t(t == 0) = pad;

        % Compute the proportions of hits and false alarms
        pH  = t(1,1) / sum(t(1,:));
        pFA = t(2,1) / sum(t(2,:));

        % Compute sensitivity and bias measures
        dprime(sub,iHyp) = (norminv(pH) - norminv(pFA));
        crit(sub,iHyp) = (1/2) * (norminv(pH) + norminv(pFA));
        % N.B. In that case we compare presence versus absence of
        % regularities. It is therefore equivalent to a yes/no type of
        % paradigme: no 1/sqrt(2) downward scaling factor is thus required
        % in the computation of the d' sensitivity measure.
    end 
    
    % 2) Probabilistic versus deterministic regularity
    
    % Get data
    proba = estgenproc(cidx{1},sub); % proba
    deter = estgenproc(cidx{2},sub); % deter
    
    % Get frequencies for each possible answer
    H  = sum(proba == 1); % hit
    M  = sum(proba == 2); % miss
    CR = sum(deter == 2); % correct rejection
    FA = sum(deter == 1); % false alarm
	
    % Create a contingency table and pad unobserved values 
    t = [H, M; FA, CR];
    t(t == 0) = pad;
    
    % Compute the proportions of hits and false alarms
    pH  = t(1,1) / sum(t(1,:));
    pFA = t(2,1) / sum(t(2,:));
    
    % Compute sensitivity and bias measures
    dprime(sub,end) = (1/sqrt(2)) * (norminv(pH) - norminv(pFA));
    crit(sub,end) = (1/2) * (norminv(pH) + norminv(pFA));
    % N.B. In that case we compare probabilistic versus deterministic
    % regularities. It is therefore equivalent to a 2AFC type of
    % paradigme: a 1/sqrt(2) downward scaling factor is thus required
    % in the computation of the d' sensitivity measure.
end

% Average over subjects
gavgdprime = mean(dprime,1); % average d'
gavgcrit = mean(crit,1); % average criterion

% Compute group-level multidimensional SDT variables
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Get average values 
Dsp = gavgdprime(1); % d' between fully-stochastic and probabilistic sequences
Dsd = gavgdprime(2); % d' between fully-stochastic and deterministic sequences
Dpd = gavgdprime(3); % d' between probabilistic and deterministic sequences
Csp = gavgcrit(1); % criterion between fully-stochastic and probabilistic sequences
Csd = gavgcrit(2); % criterion between fully-stochastic and deterministic sequences
Cpd = gavgcrit(3); % criterion between probabilistic and deterministic sequences

% Compute the angle between dimensions
Apsd = acosd((Dsd^2 + Dsp^2 - Dpd^2) / (2*Dsd*Dsp)); % "fully-stochastic" angle
Aspd = acosd((Dpd^2 + Dsp^2 - Dsd^2) / (2*Dpd*Dsp)); % "probabilistic" angle
Apds = acosd((Dpd^2 + Dsd^2 - Dsp^2) / (2*Dpd*Dsd)); % "deterministic" angle
theta = Apsd; % 

% A useful angle (for trigonometry below) that specify the rotation of the
% triangle from its horizontal position (in which the S/D axis is along the
% x-axis)
rotang = 90-theta/2;

% Compute coordinates in the 2D space of the 3 types of sequences given
% pairwise d'
Scoord = [0,0];
Dcoord = [ cosd(rotang)*Dsd, sind(rotang)*Dsd];
Pcoord = [-cosd(rotang)*Dsp, sind(rotang)*Dsp];
XYcoord = [Pcoord; Dcoord; Scoord]; % combine them in a single matrix

% Display the results
% ~~~~~~~~~~~~~~~~~~~

% Create colormaps based on used colors
cmaprec = 1000;
vec = linspace(0, 1, cmaprec);
cmap = repmat({NaN(cmaprec,3)}, 1, 3);
for iHyp = 1:3
    for rgb = 1:3
        cmap{iHyp}(:,rgb) = pchip(0:1, [1, tricol(iHyp,rgb)], vec);
    end
end

% Prepare the window
figure('Position', [353 806 360 300]);

% Display underlying distributions
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Prepare gaussian distributions
X = -3:0.01:3; % grid for the gaussian distribution
Y = normpdf(X, 0, 1); % 1D gaussian distribution
Y = setminmax(Y, [0,1]); % easier to scale below
Y = Y'*Y; % 2D gaussian distribution

% Overlap underlying gaussian distributions that gave rise to the observed
% d' sensitivity values
for iHyp = 1:3
    image(X+XYcoord(iHyp,1), X+XYcoord(iHyp,2), ...
        imagesc2image(Y, cmap{iHyp}), 'AlphaData', 3*Y/4); hold('on');
end

% Display standard deviations of these distributions
for iHyp = 1:3
    for iRad = 1/2:1/2:2
        l = circle(XYcoord(iHyp,1), XYcoord(iHyp,2), iRad, tricol(iHyp,:));
        set(l, 'LineWidth', 1/2, 'LineStyle', ':');
    end
end

% Display d'
% ~~~~~~~~~~

% Display the triangle whose sides are determined based on pairwise d'
fill(XYcoord(:,1), XYcoord(:,2), 'k-', 'FaceColor', 'None', 'LineWidth', 3);

% Display summits of the triangle
for iHyp = 1:3
    plot(XYcoord(iHyp,1), XYcoord(iHyp,2), 'o', 'MarkerEdgeColor', 'k', ...
        'MarkerFaceColor', tricol(iHyp,:), 'MarkerSize', 20);
    text(XYcoord(iHyp,1), XYcoord(iHyp,2), upper(proclab{iHyp}(1)), ...
        'Color', 'w', 'FontWeight', 'Bold', 'FontSize', 12);
end

% Display values of d'
text(mean([Scoord(1), Dcoord(1)]), mean([Scoord(2), Dcoord(2)]), ...
    sprintf('d'' = %1.2f', Dsd), 'Rotation', rotang, 'VerticalAlignment', 'Top');
text(mean([Scoord(1), Pcoord(1)]), mean([Scoord(2), Pcoord(2)]), ...
    sprintf('d'' = %1.2f', Dsp), 'Rotation', -rotang, 'VerticalAlignment', 'Top');
text(mean([Pcoord(1), Dcoord(1)]), mean([Pcoord(2), Dcoord(2)]), ...
    sprintf('d'' = %1.2f', Dpd), 'Rotation', Aspd-rotang, 'VerticalAlignment', 'Bottom');

% Display the value of the (non-)indenpendancy angle determined by d' values
text(0, -0.2, ['\theta', sprintf(' = %1.0f°', theta)], 'FontSize', 12, ...
    'VerticalAlignment', 'Top', 'HorizontalAlignment', 'Center');

% Display criterion
% ~~~~~~~~~~~~~~~~~

% Length of the criterion lines
cll = 1/2; % valid only when axis("equal") is called

% Criterion for S/D axis
Xmid = (Dcoord(1) - Scoord(1)) / 2;
Ymid = (Dcoord(2) - Scoord(2)) / 2;
plot(Xmid, Ymid, 'k.', 'MarkerSize', 15); % null criterion
if Csd > 0
    Xmid = Xmid + Csd * cosd(rotang);
    Ymid = Ymid + Csd * sind(rotang);
elseif Csd < 0
    Xmid = Xmid - Csd * cosd(rotang);
    Ymid = Ymid - Csd * sind(rotang);
end
a = cll/2 * sind(rotang);
b = cll/2 * cosd(rotang);
plot([Xmid+a, Xmid-a], [Ymid-b, Ymid+b], 'k-', 'LineWidth', 1/2);

% % Criterion for S/P axis
Xmid = (Pcoord(1) - Scoord(1)) / 2;
Ymid = (Pcoord(2) - Scoord(2)) / 2;
plot(Xmid, Ymid, 'k.', 'MarkerSize', 15); % null criterion
if Csp > 0
    Xmid = Xmid - Csp * cosd(rotang);
    Ymid = Ymid + Csp * sind(rotang);
elseif Csp < 0
    Xmid = Xmid + Csp * cosd(rotang);
    Ymid = Ymid - Csp * sind(rotang);    
end
a = cll/2 * sind(rotang);
b = cll/2 * cosd(rotang);
plot([Xmid-a, Xmid+a], [Ymid-b, Ymid+b], 'k-', 'LineWidth', 1/2);

% Customize the plot
% ~~~~~~~~~~~~~~~~~~

% Customize the axes
axis('xy'); axis('equal'); axis('off');
axis([-4,4,-2,5]);

% Add a scale at the bottom
ax = get(gca, 'Position');
axes('Position', [ax(1) ax(2) ax(3)/diff(xlim) ax(4)/diff(ylim)]);
axis(repmat([0,1],1,2));
set(gca, 'XTick', 0:1, 'YTick', 0:1, 'TickLen', zeros(1,2), 'Box', 'Off');
ylabel('z-units');

% Save the figure
if isfield(D{1}, 'Seq'), save2pdf('figs/F_D_MsdtS.pdf');
else, save2pdf('figs/F_D_MsdtIO.pdf');
end

% Compare d' between the detection of probabilistic versus deterministic
% regularities
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Compute a paired t-test between values of d' between fully-stochastic
% sequences and, respetiveley, probabilistic versus deterministic regularities
[~,pval,tci,stats] = ttest(diff(dprime(:,1:2), 1, 2));
disptstats(pval, tci, stats);

% Prepare a new window
figure('Position', [714 806 120 200]);

% Display difference in d' between the two types of regularity
Emergence_PlotSubGp(dprime(:,1:2), tricol(1:2,:));

% Customize the axes
xlim([0,3]);
set(gca, 'XTick', [], 'XColor', 'None', 'Box', 'Off');

% Add some text labels
ylabel('Detection d''');

% Save the figure
if isfield(D{1}, 'Seq'), save2pdf('figs/F_D_DpS.pdf');
else, save2pdf('figs/F_D_DpIO.pdf');
end

%% GET THE IO'S POTERIOR BELIEFS IN DETECTED AN UNDETECTED PROBABILISTIC REGULARITIES
%  ==================================================================================

% Focus on probabilistic regularities because these are those that are most
% often missed by the ideal observer
iHyp = 1;

% Get the average beliefs of the ideal observer for probabilistic
% regularities that have been detected versus missed by the subjects
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Get the hidden generative process
truegenproc = cellfun(@(x) find(strncmpi(x.Cond(1), proclab, 1)), G);

% Deduce wether the estimation was correct
corest = (truegenproc == estgenproc);
corest = mat2cell(corest(cidx{iHyp},:), numel(cidx{iHyp}), ones(nSub,1));

% Get the ultimate value of belief in the probabilistic hypothesis from the
% ideal observer
ultibel = cellfun(@(x,y) x.BarycCoord(end,y), IO, num2cell(truegenproc));
ultibel = mat2cell(ultibel(cidx{iHyp},:), numel(cidx{iHyp}), ones(nSub,1));

% Average those values separately in case of detected versus missed
% regularities (based on subjects' detection)
avgbel = cat(2, cellfun(@(x,y) mean(x(y == 0)), ultibel, corest)', ...
                cellfun(@(x,y) mean(x(y == 1)), ultibel, corest)');

% Perform a paired t-test
[~,pval,tci,stats] = ttest(diff(avgbel, 1, 2));
disptstats(pval, tci, stats);

% Display average IO's beliefs following sequences entailing a
% probabilistic regularity that was either detected or missed by subjects
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Prepare a new window
figure('Position', [835 806 120 200]);

% Display chance level
plot([0,3], ones(1,2)./3, '-', 'Color', g, 'LineWidth', 1/2); hold('on');

% Display difference in IO's beliefs between missed versus detected
% regulartieis
Emergence_PlotSubGp(avgbel, tricol(1,:));

% Customize the axes
axis([0,3,0,1]);
set(gca, 'XTick', 1:2, 'XTickLabel', {'Undetected', 'Detected'}, ...
    'XTickLabelRotation', 30);
set(gca, 'Box', 'Off');

% Add some text labels
xlabel({'probabilistic', 'regularities'});
ylabel('Belief p(M_P|y) from the IO');

% Save the figure
save2pdf('figs/F_D_BelUndet.pdf');

%% DETERMINANTS OF DETECTION FOR PROBABILISTIC REGULARITIES
%  ========================================================

% Depending on the entropy
% ~~~~~~~~~~~~~~~~~~~~~~~~

% Entropy function
H = @(p) -(p .* log2(p) + (1-p) .* log2(1-p));

% Compute entropy levels
h = cellfun(@(x) H(x(1))*H(x(2)), prob);

% Define bins
nEnt = 2;
[~,edges] = histcounts(h, nEnt);
entgpidx = {find(h < edges(2)), find(h > edges(2))};
entgpttl = {sprintf('Low (< %1.1f)', edges(2)), sprintf('High (> %1.1f)', edges(2))};

% Get the average detection rate in each entropy bins
entdetecrate = cell2mat(cellfun(@(x) mean(estgenproc(cidx{1}(x),:) == 1, 1)', ...
    entgpidx, 'UniformOutput', 0));

% Prepare a window
figure('Position', [956 806 240 200]);

% Display difference in detection rate according to entropy levels
subplot(1,2,1);
Emergence_PlotSubGp(entdetecrate, tricol(1,:));

% Customize the axes
axis([0,3,0,1]);
set(gca, 'XTick', 1:2, 'XTickLabel', entgpttl, 'XTickLabelRotation', 30);
set(gca, 'Box', 'Off');

% Add some labels
xlabel('Entropy levels');
ylabel('Average detection rate');

% Depending on the type of bias
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Create groups
pAlt = cellfun(@(x) (2.*x(1).*x(2)) ./ (x(1) + x(2)), prob);
biasgpidx = {find(pAlt < 1/2), find(pAlt > 1/2), 7:9};
biasgpttl = {'Repetition', 'Alternation', 'Frequency',};
nGp = numel(biasgpidx);

% Get the average detection rate in each group of regularities
biasdetecrate = cell2mat(cellfun(@(x) mean(estgenproc(cidx{1}(x),:) == 1, 1)', ...
    biasgpidx, 'UniformOutput', 0));

% Display difference in detection rate according to biased dimensions
subplot(1,2,2);
Emergence_PlotSubGp(biasdetecrate, tricol(1,:));

% Customize the axes
axis([0,4,0,1]);
set(gca, 'XTick', 1:nGp, 'XTickLabel', biasgpttl, 'XTickLabelRotation', 30);
set(gca, 'Box', 'Off');

% Add some labels
xlabel('Biased dimension');
ylabel('Average detection rate');

% Save the figure
if isfield(D{1}, 'Seq'), save2pdf('figs/F_D_ProbRegS.pdf');
else, save2pdf('figs/F_D_ProbRegIO.pdf');
end

%% DETERMINANTS OF DETECTION FOR DETERMINISTIC REGULARITIES
%  ========================================================

% Depending on the rules' length
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Create groups
len = cellfun(@numel, dr(1:end));
lengpidx = cellfun(@(x) find(x == len), num2cell(unique(len)), 'UniformOutput', 0)';
nGp = numel(lengpidx);

% Get the average detection rate in each group of regularities
lendetecrate = cell2mat(cellfun(@(x) mean(estgenproc(cidx{2}(x),:) == 2, 1)', ...
    lengpidx, 'UniformOutput', 0));

% Prepare a new window
figure('Position', [1196 806 240 200]);

% Display difference in detection rate according to rules' length
subplot(1,2,1);
Emergence_PlotSubGp(lendetecrate, tricol(2,:));

% Customize the axes
axis([0,nGp+1,0,1]);
set(gca, 'XTick', 1:nGp, 'XTickLabel', num2cell(unique(len)));
set(gca, 'Box', 'Off');

% Add some labels
xlabel('Rules'' length');
ylabel('Average detection rate');

% Depending on the rules' difficulty
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Create groups
cplxidx = {[1,2,5,8], [3,6,9], [4,7,10]};
cplxttl = {'Easy', 'Moderate', 'Difficult'};
nGp = numel(cplxidx);

% Get the average detection rate in each group of regularities
cplxetecrate = cell2mat(cellfun(@(x) mean(estgenproc(cidx{2}(x),:) == 2, 1)', ...
    cplxidx, 'UniformOutput', 0));

% Display difference in detection rate according to rules' difficulty
subplot(1,2,2);
Emergence_PlotSubGp(cplxetecrate, tricol(2,:));

% Customize the axes
axis([0,4,0,1]);
set(gca, 'XTick', 1:numel(cplxidx), 'XTickLabel', cplxttl, 'XTickLabelRotation', 30);
set(gca, 'Box', 'Off');

% Add some labels
xlabel('Rule''s difficulty');
ylabel('Average detection rate');

% Save the figure
if isfield(D{1}, 'Seq'), save2pdf('figs/F_D_DetRegS.pdf');
else, save2pdf('figs/F_D_DetRegIO.pdf');
end
