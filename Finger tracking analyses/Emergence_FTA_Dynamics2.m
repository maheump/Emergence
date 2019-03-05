% This script inspects whether there is an effect of manipulated factors
% of probabilistic and deterministic regularities on the detection update
% of probabilistic regularities and the detection lag of deterministic
% regularities. These measures were chosen because they reflect how the
% beleifs are updated in both cases: progressively for probabilistic
% regularities (the update can thus vary), abruptly for deterministic
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

% Define usefule variables
mrkr = {'s', 'o', '^', 'd', 'p'};
lgdlab = {'Frequency', 'Alternation', 'Repetition', 'Outside', 'No'};
cmap = {'Blues', 'Purples'};

%% EFFECT OF SHANNON ENTROPY ON BELIEF UPDATE
%  ==========================================

% For each type of regularity
for iHyp = [2,1]
    
    % Get factors
    % ~~~~~~~~~~~
    
    % Compute entropy levels of theoretical transition probabilities
    if iHyp  == 1
        pXgY = cell2mat(prob');
        
    % Get transition probabilities from the deterministic regularities
    elseif iHyp == 2
        [~, ~, pAgB, pBgA] = cellfun(@(x) pat2proba(x, 1:2, true), det');
        pAgB(pAgB == 1) = 1-eps;
        pXgY = [pAgB, pBgA];
    end
    
    % Define the type of biases that have been used
    biases = NaN(numel(cidx{iHyp}), 1);
    biases(pXgY(:,1) == 1 - pXgY(:,2)) = 1;               % frequency biases
    biases(pXgY(:,1) > 1/2 & pXgY(:,1) == pXgY(:,2)) = 2; % repetition biases
    biases(pXgY(:,1) < 1/2 & pXgY(:,1) == pXgY(:,2)) = 3; % alternation biases
    biases(pXgY(:,1) == 1/2 & pXgY(:,1) == 1/2) = 5;      % no biases
    biases(isnan(biases)) = 4;                            % outside diagonals
    
    % Create colormap for the entropy
    minH = 1.4;
    maxH = Emergence_MarkovEntropy(1/2, 1/2);
    prec = 1001;
    offset = round(prec * (maxH - minH));
    EntCMap = flipud([flipud(cbrewer2(cmap{iHyp}, offset)); cbrewer2('Greys', prec)]);
    prec = size(EntCMap,1);
    
    % Compute Shannon entropy on transition probabilities
    TPent = arrayfun(@(x,y) Emergence_MarkovEntropy(x, y), pXgY(:,1), pXgY(:,2));
    
    % Get dedicated color for each condition that is indexed on the entropy of
    % the rule
    [~,idx] = min(abs(linspace(0, maxH, prec) - TPent), [], 2);
    condmap{1} = EntCMap(idx,:);
    
    % Get detection velocity
    % ~~~~~~~~~~~~~~~~~~~~~~
    
    % Regress entropy against detection velocity for each subject
    subvar = mat2cell(update{iHyp}', ones(nSub,1), numel(cidx{iHyp}));
    subvar = cellfun(@(x) x - mean(x, 2, 'OmitNaN'), subvar, 'UniformOutput', 0);
    offset = cellfun(@(x) Emergence_Regress(x, TPent', 'OLS', 'beta0'), subvar);
    slope  = cellfun(@(x) Emergence_Regress(x, TPent', 'OLS', 'beta1'), subvar);
    
    % Frequentist t-test on the distribution of statistics over subjects
    [~,pval,tci,stats] = ttest(slope);
    Emergence_PrintTstats(pval,tci,stats);
    
    % Frequentist t-test on the distribution of statistics over subjects
    [~,bf01] = BF_ttest(slope);
    fprintf('BF in favour of the null: %1.2f\n', bf01);
    
    % Get average detection update over subjects
    avg = mean(update{iHyp}, 2, 'OmitNaN');
    err = sem(update{iHyp}, 2);
    
    % Display detection velocity as a function of entropy
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    % Prepare a new window
    figure('Position', [501+201*iHyp 780 200 325]); hold('on');
    
    % Display the regression line between Shannon entropy and detection update
    confint  = Emergence_Regress(avg, TPent, 'OLS', 'confint');
    confintx = Emergence_Regress(avg, TPent, 'OLS', 'confintx');
    col = cbrewer2(cmap{iHyp}, 1);
    fill([confintx, fliplr(confintx)], [confint(1,:), fliplr(confint(2,:))], ...
        col, 'EdgeColor', 'none', 'FaceAlpha', 0.15);
    b = Emergence_Regress(avg, TPent, 'OLS', {'beta0', 'beta1'});
    plot(confintx([1,end]), confintx([1,end]).*b(2)+b(1), '-', ...
        'Color', col, 'LineWidth', 3);
    
    % Plot the distribution of 
    fill(minH + fout0./5e3, xout, 'k');
    
    % For each type of bias
    lgd = NaN(1,numel(lgdlab));
    for iGp = 1:numel(lgdlab)
        idx = find(biases == iGp);
        condmkr{1}(idx) = mrkr{iGp};
        
        % Connect individual dots with a line
        if iGp ~= 4
            plot(TPent(idx), avg(idx), '-', 'Color', 'k', 'LineWidth', 1);
        end
        
        % For each level (strength) of that bias
        for j = 1:numel(idx)
            k = idx(j);
            
            % Display the SEM over subjects
            plot(TPent(k) + zeros(1,2), avg(k) + err(k) * [-1,1], '-', 'Color', 'k');
            
            % Display the average detection update
            lgd(iGp) = plot(TPent(k), avg(k), 'k-', 'MarkerEdgeColor', 'k', ...
                'MarkerFaceColor', condmap{1}(k,:), 'Marker', mrkr{iGp}, 'MarkerSize', 10);
        end
    end
    
    % Display the colormap
    caxis([0,1]);
    colormap(EntCMap);
    
    % Customize the axes
    set(gca, 'Box', 'Off', 'XLim', [minH, maxH], 'YLim', [0.002,0.013]);
    
    % Add some text labels
    legend(lgd(unique(biases)), lgdlab(unique(biases)), 'Location', 'South', 'Box', 'Off');
    xlabel('Entropy (bits)');
    ylabel('Update');
    
    % Save the figure
    fname = sprintf('F_Dyn_Ent%sS.pdf', proclab{iHyp}(1));
    if isfield(D{1}, 'Seq'), save2pdf(fullfile(ftapath, 'figs', fname));
    else, save2pdf(fullfile(ftapath, 'figs', fname));
    end
end

%% EFFECT OF PATTERN LENGTH ON THE DETECTION OF DETERMINISTIC REGULARITIES
%  =======================================================================

% Get factors
% ~~~~~~~~~~~

% Get variables that are manipulated
len = cellfun(@numel, dr);
group = [4 1 2 3 1 2 4 1 2 3];

% Get dedicated colors for each condition that is indexed on the length of
% the corresponding pattern
condmap{2} = flipud(cbrewer2('Reds', prec));
idx = round(1 + (len - 3) ./ (11 - 3) .* (prec-1));
condmap{2} = condmap{2}(idx,:);

% Get detection delays
% ~~~~~~~~~~~~~~~~~~~~

% Regress pattern length against detection delay for each subject
delay = lag{2}; % lag2 or lag{2}
subvar = mat2cell(delay', ones(nSub,1), numel(cidx{2}));
subvar = cellfun(@(x) x - mean(x, 2, 'OmitNaN'), subvar, 'UniformOutput', 0);
offset = cellfun(@(x) Emergence_Regress(x, len', 'OLS', 'beta0'), subvar);
slope  = cellfun(@(x) Emergence_Regress(x, len', 'OLS', 'beta1'), subvar);
[~,pval,tci,stats] = ttest(slope);
Emergence_PrintTstats(pval,tci,stats);

% Regress out the effect of pattern length (our experimental design is not
% fully factorial)
yhat = cell2mat(cellfun(@(y,b0,b1) y - (b0 + len'.*b1), subvar, ...
   num2cell(offset), num2cell(slope), 'UniformOutput', 0))';

% Average detection lag for each group of patterns
avggpergp = cell2mat(arrayfun(@(i) mean(yhat(group == i,:), 1, ...
    'OmitNaN')', unique(group), 'UniformOutput', 0));

% Perform a repeated-measures 1-way ANOVA on detection lag depending on the
% group of pattern
T = rmANOVA(avggpergp);
disp(T);

% Compute the average detection lag across subjects
avg = mean(lag{2}, 2, 'OmitNaN');
err = sem(lag{2}, 2);

% Display detection delay as a function of pattern length
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Prepare a new window
figure('Position', [1103 780 200 325]);
lgd = NaN(1,3);

% Display a grid that expresses the detection delays in terms of the number
% of repetition of pattern length, and not simply the number of
% observations
x = min(len)-1:max(len)+1;
for iRep = 1:30
    plot(x, x*iRep, '-', 'Color', g, 'LineWidth', 1/2); hold('on');
end

% Display the regression line between pattern length and detection delay
confint  = Emergence_Regress(avg, len, 'OLS', 'confint');
confintx = Emergence_Regress(avg, len, 'OLS', 'confintx');
fill([confintx, fliplr(confintx)], [confint(1,:), fliplr(confint(2,:))], ...
    tricol(2,:), 'EdgeColor', 'none', 'FaceColor', tricol(2,:), 'FaceAlpha', 0.15);
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
set(gca, 'Box', 'Off', 'XLim', [min(len)-1, max(len)+1], 'YLim', [0,ceil(max(avg+err)/10)*10]);

% Add some text labels
xlabel('Length (# obs.)');
ylabel('Detection delay (# obs.)');

% Save the figure
if isfield(D{1}, 'Seq'), save2pdf(fullfile(ftapath, 'figs', 'F_Dyn_LagDS.pdf'));
else, save2pdf(fullfile(ftapath, 'figs', 'F_Dyn_LagDIO.pdf'));
end
