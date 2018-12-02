% This script inspects whether there is an effect of manipulated factors
% of probabilistic and deterministic regularities on  the detection speed
% of probabilistic regularities and the detection lag of deterministic
% regularities. These measures were chosen because they reflect how the
% beleifs are updated in both cases: progressively for probabilistic
% regularities (the speed can thus vary), abruptly for deterministic
% regularities (the lag can thus vary).
% 
% Copyright (c) 2018 Maxime Maheu

%% INITIALIZATION
%  ==============

% The following script must be run beforehand
Emergence_FTA_Dynamics1;

% Prepare a useful output variable
condmap = cell(1,2);
condmkr = cell(1,2);

% Define markers to use
mrkr = {'s', 'o', '^', 'd'};

%% FACTORS MODULATING DETECTION DELAY OF PROBABILISTIC REGULARITIES
%  ================================================================

% Compute entropy levels of theoretical transition probabilities
PpXgY = cell2mat(prob')';
PpAgB = PpXgY(1,:);
PpBgA = PpXgY(2,:);
TPent = arrayfun(@(x,y) Emergence_MarkovEntropy(x, y), PpAgB, PpBgA)';
    
% Define the type of biases that have been used
biases = NaN(1, numel(prob));
biases(PpXgY(1,:) < 1/2 & PpXgY(2,:) > 1/2) = 1; % frequency biases
biases(sum(PpXgY < 1/2) == 2) = 2;               % repetition biases
biases(sum(PpXgY > 1/2) == 2) = 3;               % alternation biases
biases(isnan(biases)) = 4;                       % outside diagonals

% Test whether there is an effect of bias strength (measured as entropy)
% on detection dpeed
[~,~,idx] = unique(TPent);
avggperstrength = cell2mat(cellfun(@(i) mean(speed{1}(idx == i,:), 1, ...
    'OmitNaN')', num2cell(1:numel(unique(TPent))), 'UniformOutput', 0));
avggperstrength = mat2cell(avggperstrength, ones(nSub,1), numel(unique(TPent)));
slope = cellfun(@(x) Emergence_Regress(x, unique(TPent)', 'OLS', 'beta1'), avggperstrength);
[~,pval,tci,stats] = ttest(slope);
disptstats(pval,tci,stats);

% Test whether there is an effect of bias type on detection dpeed
[~,~,idx] = unique(biases);
avggpertype = cell2mat(cellfun(@(i) mean(speed{1}(idx == i,:), 1, ...
    'OmitNaN')', num2cell(1:numel(unique(TPent))), 'UniformOutput', 0));
T = rmANOVA(avggpertype, 'BiasType');
disp(T);

% Get average detection speed over subjects
avg = mean(speed{1}, 2, 'OmitNaN');
err = sem(speed{1}, 2);

% Create the colormap
prec = 1001;
decal = round(prec*(max(TPent) - 1)); %500;
pmap = flipud([flipud(cbrewer2('Blues', decal)); cbrewer2('Greys', prec)]);
prec = size(pmap,1);

% Get dedicated color for each condition that is indexed on the entropy of
% the rule
[~,idx] = min(abs(linspace(0, Emergence_MarkovEntropy(1/2,1/2), prec) - TPent), [], 2);
condmap{1} = pmap(idx,:);

% Prepare a new window
figure('Position', [702 780 200 325]);
lgd = NaN(1,3); hold('on');

% Display the regression line between Shannon entropy and detection speed
confint  = Emergence_Regress(avg, TPent, 'OLS', 'confint');
confintx = Emergence_Regress(avg, TPent, 'OLS', 'confintx');
fill([confintx, fliplr(confintx)], [confint(1,:), fliplr(confint(2,:))], ...
    tricol(1,:), 'EdgeColor', 'none', 'FaceColor', 'k', 'FaceAlpha', 0.15);
b = Emergence_Regress(avg, TPent, 'OLS', {'beta0', 'beta1'});
plot(confintx([1,end]), confintx([1,end]).*b(2)+b(1), '-', ...
    'Color', tricol(1,:), 'LineWidth', 3);

% For each type of bias
for iGp = 1:4
    idx = find(biases == iGp);
    condmkr{1}(idx) = mrkr{iGp};
	
    % Conenct individual dots with a line
    if iGp ~= 4
        plot(TPent(idx), avg(idx), '-', 'Color', 'k', 'LineWidth', 1);
    end
    
    % For each level (strength) of that bias
    for j = 1:numel(idx)
        k = idx(j);
        
        % Display the SEM over subjects
        plot(TPent(k) + zeros(1,2), avg(k) + err(k) * [-1,1], '-', 'Color', 'k');
        
        % Display the average detection speed
        lgd(iGp) = plot(TPent(k), avg(k), 'k-', 'MarkerEdgeColor', 'k', ...
            'MarkerFaceColor', condmap{1}(k,:), 'Marker', mrkr{iGp}, 'MarkerSize', 10);
    end
end

% Display the colormap
caxis([0,1]);
colormap(pmap);

% Customize the axes
set(gca, 'Box', 'Off', 'XLim', [0.95, Emergence_MarkovEntropy(1/2,1/2)], 'YLim', [0,0.012]);

% Add some text labels
legend(lgd, cellfun(@(x) sprintf('%s diag.', x), {'Repetition', ...
    'Alternation', 'Frequency', 'Outside'}, 'UniformOutput', 0), ...
    'Location', 'South', 'Box', 'Off');
xlabel('Entropy (bits)');
ylabel('Speed');

% Save the figure
if isfield(D{1}, 'Seq'), save2pdf(fullfile(ftapath, 'figs', 'F_Dyn_DetProbaS.pdf'));
else, save2pdf(fullfile(ftapath, 'figs', 'F_Dyn_DetProbaIO.pdf'));
end

%% FACTORS MODULATING DETECTION DELAY OF DETERMINISTIC REGULARITIES
%  ================================================================

% Get variables that are manipulated
len = cellfun(@numel, dr);
group = [1, repmat(1:3, 1, 3)];

% Perform an ANOVA on the effet of pattern length and pattern type on 
% detection lag (in terms of number of repetitions one hand and number of
% observations on the other hand) 
Traw   = rmANOVA(reshape(lag{2}( 2:end,:)', [nSub,3,3]), {'Type', 'Length'}); % exclude AABB
Trelat = rmANOVA(reshape(lag2(   2:end,:)', [nSub,3,3]), {'Type', 'Length'}); % exclude AABB
Trelatall = rmANOVA(lag2'); % keep AABB
disp(Traw);
disp(Trelat);
disp(Trelatall);

% Compute the average detection lag across subjects
avg = mean(lag{2}, 2, 'OmitNaN');
err = sem(lag{2}, 2);

% Get dedicated color for each condition that is indexed on the length of
% the corresponding pattern
condmap{2} = flipud(cbrewer2('Reds', prec));
idx = round(1 + (len - 3) ./ (11 - 3) .* (prec-1));
condmap{2} = condmap{2}(idx,:);

% Prepare a new window
figure('Position', [702 381 200 325]);
lgd = NaN(1,3);

% Display a grid that expresses the detection delays in terms of the number
% of repetition of pattern length, and not simply the number of
% observations
x = min(len)-1:max(len)+1;
for iRep = 1:30
    plot(x, x*iRep, '-', 'Color', g, 'LineWidth', 1/2); hold('on');
end

% Display the regression line between Shannon entropy and detection speed
confint  = Emergence_Regress(avg, len, 'OLS', 'confint');
confintx = Emergence_Regress(avg, len, 'OLS', 'confintx');
fill([confintx, fliplr(confintx)], [confint(1,:), fliplr(confint(2,:))], ...
    tricol(2,:), 'EdgeColor', 'none', 'FaceColor', 'k', 'FaceAlpha', 0.15);
b = Emergence_Regress(avg, len, 'OLS', {'beta0', 'beta1'});
plot(confintx([1,end]), confintx([1,end]).*b(2)+b(1), '-', ...
    'Color', tricol(2,:), 'LineWidth', 3);

% For each group of conditions (based on pattern length)
for iGp = 1:max(group)
    idx = find(iGp == group);
    condmkr{2}(idx) = mrkr{iGp};
    
    % Display the change in detection delay as a function of pattern length
    plot(len(idx), avg(idx), '-', 'Color', 'k', 'LineWidth', 1);
    
    % For each condition, i.e. each repeating pattern
    for iCond = 1:numel(idx)
        k = idx(iCond);
        
        % Display the dispersion of detection delay across subjets
        plot(len(k) + zeros(1,2), avg(k) + err(k) * [-1,1], '-', 'Color', 'k');
        
        % Display the average detection delay
        lgd(iGp) = plot(len(k), avg(k), mrkr{iGp}, 'MarkerEdgeColor', 'k', ...
            'MarkerFaceColor', condmap{2}(k,:), 'MarkerSize', 10);
    end
end

% Customize the axes
set(gca, 'Box', 'Off', 'XLim', [min(len)-1, max(len)+1], 'YLim', [0,80]);

% Add some text labels
xlabel('Length (# obs.)');
ylabel('Detection delay (# obs.)');

% Save the figure
if isfield(D{1}, 'Seq'), save2pdf(fullfile(ftapath, 'figs', 'F_Dyn_DetDeterS.pdf'));
else, save2pdf(fullfile(ftapath, 'figs', 'F_Dyn_DetDeterIO.pdf'));
end

%% DISPLAY DETECTION SPEED OF MATCHED PROBABILISTIC/DETERMINISTIC REGULARITIES
%  ===========================================================================

% Get transition probabilities from the deterministic regularities
[~, DpAlt, DpAgB, DpBgA] = cellfun(@(x) pat2proba(x, 1:2, true), det');

% Restrict to regularities on the repetition/alternation diagonal
Pidx = find(PpAgB == PpBgA);
nidx = numel(Pidx);

% Find corresponding probabilistic biases 
Didx = NaN(1,nidx);
for i = 1:nidx
    Didx(i) = find(PpAgB(Pidx(i)) == DpAgB & PpBgA(Pidx(i)) == DpBgA);
end

% Sort according to the frequency of alternations
DpAlt = DpAlt(Didx);
[DpAlt, sortidx] = sort(DpAlt);
IDX = {Pidx(sortidx), Didx(sortidx)};

% Get speed measures (put NaN for negative speeds)
spd = cellfun(@(x,i) x(i,:), speed, IDX, 'UniformOutput', 0);
for i = 1:2, spd{i}(spd{i} < 0) = NaN; end
spd = cellfun(@(x) log(x), spd, 'UniformOutput', 0);

% Perform a t-test on average speed between probabilistic and deterministic
% versions
avgverspd = cell2mat(cellfun(@(x) mean(x, 1, 'OmitNaN')', spd, 'UniformOutput', 0));
[~,pval,tci,stats] = ttest(diff(avgverspd, 1, 2));
disptstats(pval,tci,stats);

% Perform a t-test on the effect of entropy on detection speed in
% probabilistic versus deterministic versions of regularities
h = TPent(IDX{1}(1:3)); % 3 entropy levels
probaver = fliplr(nanmean(cat(3, spd{1}(1:3,:)', flipud(spd{1}(4:6,:))'), 3));
deterver = fliplr(nanmean(cat(3, spd{2}(1:3,:)', flipud(spd{2}(4:6,:))'), 3));
probaslope = cellfun(@(x) Emergence_Regress(x', h, 'OLS', 'beta1'), mat2cell(probaver, ones(nSub,1), 3));
deterslope = cellfun(@(x) Emergence_Regress(x', h, 'OLS', 'beta1'), mat2cell(deterver, ones(nSub,1), 3));
[~,pval,tci,stats] = ttest(probaslope-deterslope);
disptstats(pval,tci,stats);

% Compute group average and error
avgregspd = cell2mat(cellfun(@(x) mean(x, 2, 'OmitNaN'), spd, 'UniformOutput', 0));
semregspd = cell2mat(cellfun(@(x) sem( x, 2           ), spd, 'UniformOutput', 0));

% Prepare a new figure
figure('Position', [903 858 400 300]); hold('on');

% Define colors indexed on entropy levels and that are different for
% repetition/alternation biases
prec = 200;
cmap = cbrewer2('PuOr', prec+1);
hgrid = [linspace(-2, -1, prec/2), linspace(1, 2, prec/2)];
[~,idx] = min(abs(hgrid - [-sort(h, 'Descend'); sort(h, 'Ascend')]), [], 2);
cmap = cmap(idx,:);


prec = 1001;
cmap = cbrewer2('BuPu', prec);
hgrid = linspace(1, Emergence_MarkovEntropy(1/2,1/2), prec);
[~,idx] = min(abs(hgrid - [sort(h, 'Descend'); sort(h, 'Ascend')]), [], 2);
cmap = cmap(idx,:);

% Display a distance map
lim = [-6.5, -4.5];
idx = linspace(lim(1), lim(2), 100);
bgmap = imagesc2image(abs(idx - idx'), flipud(gray));
image(idx, idx, bgmap); alpha(1/2);

% Display error bars in both dimensions
for i = 1:nidx
    plot(avgregspd(i,2)+semregspd(i,2).*[-1,1], repmat(avgregspd(i,1),1,2), ...
        '-', 'Color', 'k', 'LineWidth', 1/2);
    plot(repmat(avgregspd(i,2),1,2), avgregspd(i,1)+semregspd(i,1).*[-1,1], ...
        '-', 'Color', 'k', 'LineWidth', 1/2);
end

% Connect dots corresponding to repetition/alternation biases
for i = [1,4]
    idx = i:i+2;
    plot(avgregspd(idx,2), avgregspd(idx,1), '-', ...
        'LineWidth', 3, 'Color', mean(cmap(idx,:)));
end

% Display average speed
mkr = {'s', 'o'};
for i = 1:nidx
    plot(avgregspd(i,2), avgregspd(i,1), mkr{1 + (DpAlt(i) > 1/2)}, ...
        'MarkerSize', 8, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', cmap(i,:));
end

% Customize the axes
axis(repmat(lim, 1, 2)); axis('square');
set(gca, 'Box', 'On', 'XColor', tricol(2,:), 'YColor', tricol(1,:));

% Add an identity line
plot(lim, lim, '-', 'Color', g);
text(lim(1), lim(1), '   Identity', 'Rotation', 45, 'Color', g, ...
    'HorizontalAlignment', 'Left', 'VerticalAlignment', 'Bottom');

% Add text labels
xlabel({'log average speed along p(H_D|y)', 'for the deterministic versions'}, ...
    'Color', tricol(2,:));
ylabel({'log average speed along p(H_P|y)', 'for the probabilistic versions'}, ...
    'Color', tricol(1,:));

% Add a colorbar
colormap(cmap);
cbr = colorbar;
clab1 = dr(IDX{2})';
clab2 = arrayfun(@(p1,p2) sprintf('p(A|B) = %1.2f & p(A|B) = %1.2f', ...
    p1, p2), PpAgB(IDX{1}), PpBgA(IDX{1}), 'UniformOutput', 0);
set(cbr, 'Ticks', linspace(1/6/2, 1-1/6/2, 6), 'TickLabels', ...
    cellfun(@(x,y) sprintf('%s <=> %s', x, y), clab1, clab2, 'UniformOutput', 0));

% Save the figure
if isfield(D{1}, 'Seq'), save2pdf(fullfile(ftapath, 'figs', 'F_Dyn_DetVsProbaSpeedS.pdf'));
else, save2pdf(fullfile(ftapath, 'figs', 'F_Dyn_DetVsProbaSpeedIO.pdf'));
end
