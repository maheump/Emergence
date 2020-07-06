% The script relates subjects' reported positions of the change point
% during post-sequence questions to the model's posterior probabilities
% over change point positions.
% 
% Copyright (c) 2020 Maxime Maheu

% Initialization
% ~~~~~~~~~~~~~~

% Define the window around most likely change pint position 
win = -50:50;
nObs = numel(win);

% Prepare output variables
repCP   = cell(1,2);
kernRep = NaN(nObs,nSub,2);
trueCP  = cell(1,2);
kernCP  = NaN(nObs,nSub,2);
postCP  = NaN(nObs,nSub,2);

% For each type of regularity
for iReg = 1:2
    
    % Get model's most likely change point position
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    % Get posterior probabilities over change point position
    fullpostCP = cellfun(@(p) p.CPbelief(:,end)', IO(cidx{iReg},:), 'uni', 0);
    
    % Find the most likely position of the change point
    [~,mlCPpos] = cellfun(@max, fullpostCP);
    
    % Restrict to sequences that were correctly labeled
    detecmask = (filter{iHyp} == 1 | filter{iHyp} == 3);
    mlCPpos(~detecmask) = NaN;
    
    % Get posterior probability of change positions
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    % Get indices of observations around most likely change point position
    idx = cellfun(@(i) win + i, num2cell(mlCPpos), 'uni', 0);
    
    % Get corresponding values of the posterior probabilities around the
    % most likely change point position (pad with NaNs)
    curpostCP = cellfun(@(i,p) [NaN(1, sum(i < 1)), ...    % pad before
                                p(i(i >= 1 & i <= N)), ... % posterior
                                NaN(1, sum(i > N))]', ...  % pad after
                                idx, fullpostCP, 'uni', 0);
    
	% Restrict to sequences that were correctly labeled
    curpostCP(~detecmask) = {NaN(nObs,1)};
    curpostCP(cellfun(@isempty, curpostCP)) = {NaN(nObs,1)};
    
    % Average over sequences
    curpostCP = cell2mat(reshape(curpostCP, [1,size(curpostCP)]));
    postCP(:,:,iReg) = squeeze(mean(curpostCP, 2, 'OmitNaN'));
	
    % Get reported change point positions from the subjects
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    % Get subjects' reports about most likely change point position
    repCP{iReg} = cellfun(@(x) x.Questions(3), G(cidx{iReg},:));
    
    % Center reports relative to the model's estimation of the most likely
    % change point position
    repCP{iReg} = repCP{iReg} - mlCPpos;
    
    % Restrict to sequences that were correctly labeled
    repCP{iReg}(~detecmask) = NaN;
    
    % For each subject
    for iSub = 1:nSub
        
        % Compute subject-specific kernel densities over reported change point
        % positions
        d = repCP{iReg}(:,iSub);
        kernRep(:,iSub,iReg) = ksdensity(d(~isnan(d)), ...
            win, ...                             % grid of positions
            'BandWidth',          12, ...        % bandwidth of the kernel smoothing window 
            'BoundaryCorrection', 'Reflection'); % type of correction for the boundaries
    end
    
    % Get true position of the change point
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    % True chance point position
    trueCP{iReg} = cellfun(@(x) x.Jump, G(cidx{iReg},:));
    trueCP{iReg} = trueCP{iReg} - mlCPpos;
    
    % Restrict to sequences that were correctly labeled
    trueCP{iReg}(~detecmask) = NaN;
    
    % For each subjects
    for iSub = 1:nSub
        
        % Compute kernel densities over true change point positions
        d = trueCP{iReg}(:,iSub);
        kernCP(:,iSub,iReg) = ksdensity(d(~isnan(d)), ...
            win, ...                             % grid of positions
            'BandWidth',          7, ...         % bandwidth of the kernel smoothing window 
            'BoundaryCorrection', 'Reflection'); % type of correction for the boundaries
    end
end

% Display the different distributions in a window around the most likely
% change point position given the model
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Label of distribution to plot
varlab = {'postCP', 'kernCP', 'kernRep'};
ylab = {'Posterior over CP positions', 'True CP positions (kernel)', ...
    'Reported positions (kernel)'};
nSp = numel(varlab);

% Prepare a new window
figure('Position', [463 631 660 200]);

% For eadch type of regularity
for iReg = 1:2
    
    % For each distribution
    for iDist = 1:nSp
        subplot(1, nSp, iDist);
        
        % Display the distribution (with shaded error bars)
        d = eval(cat(2, varlab{iDist}, '(:,:,iReg)'));
        m = mean(d, 2, 'OmitNaN');
        s = sem(d, 2);
        plotMSEM(win, m, s, 0.15, tricol(iReg,:), tricol(iReg,:), 2);
        
        % Display position of most likely change point position given the
        % model
        plot(zeros(1,2), [0, max(m)], 'k:'); hold('on');
        
        % Display true positions of change points
        nCP = histcounts(trueCP{iReg}(:), [-N,win(2:end),N]);
        nCP(nCP == 0) = NaN;
        scatter(win, ones(1,nObs) .* (0 + iReg*0.05*diff(ylim)), nCP*3, ...
            tricol(iReg,:), 'filled', 'MarkerEdgeColor', 'k');
        
        % Customize the axes
        set(gca, 'Box', 'Off');
        axis('tight'); xlim(win([1,end]));
        
        % Add some text labels
        xlabel({'Obs. # w.r.t. to most likely', 'change point given the model'});
        ylabel(ylab{iDist});
    end
end

% Display the difference in the variance of reported change point positions
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Get variance of reported change point position
varcp = cell2mat(cellfun(@(x) var(x, [], 1, 'OmitNaN')', repCP, 'uni', 0));

% Compare between types of regularity
[h,pval,tci,stats] = ttest(diff(varcp,1,2));
Emergence_PrintTstats(pval,tci,stats);

% Display difference in variance
figure('Position', [1124 631 120 200]);
Emergence_PlotSubGp(varcp, tricol(1:2,:));
set(gca, 'Yscale', 'log', 'Xcolor', 'None', 'Xlim', [0,3], 'Box', 'Off');
ylabel({'Variance of the reported', 'change point positions'});

% Display the delay between the model and the subjects' reports
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Get delay between most likely and reported change points positions
delay = cell2mat(cellfun(@(x) mean(x, 1, 'OmitNaN')', repCP, 'uni', 0));

% Compare delay to zero for both types of regularities
[h,pval,tci,stats] = ttest(delay);
Emergence_PrintTstats(pval,tci,stats);

% Compare delay between types of regularity
[h,pval,tci,stats] = ttest(diff(delay,1,2));
Emergence_PrintTstats(pval,tci,stats);

% Display distribution of 
figure('Position', [1245 631 120 200]);
plot([0,3], zeros(1,2), 'k-'); hold('on');
Emergence_PlotSubGp(delay, tricol(1:2,:));
set(gca, 'Xcolor', 'None', 'Xlim', [0,3], 'Box', 'Off');
ylim(max(abs(ylim)).*[-1,1]);
ylabel({'Delay between most likely and', 'reported change point positions'});
