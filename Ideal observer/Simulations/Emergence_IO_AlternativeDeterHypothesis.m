% This script runs simulations of the ideal observer model using different
% types of deterministic hypotheses. In particular we constrast the
% repeating pattern detector to various Markov chains of different orders.
% The simulations are run on sequences that were presented to the subjects.
% Post change point detection dynamics to probabilistic and deterministic
% regularities are then inspected. We conclude that the fully-deterministic
% version of the deterministic hypothesis (i.e. the repeating pattern
% detector) provides the best qualitative fit to the subjects' data.
% 
% Copyright (c) 2018 Maxime Maheu

%% INITIALIZATION
%  ==============

% Clear the place
clear; close('all');

% Load data
Emergence_FTA_LoadData;

% Define the order of the Markov chains to use
orders = 1:9;
nOrd = numel(orders);

% Restrict to sequences with deterministic regularities only
nCond = size(G,1);
conds = 1:nCond;

% Define options for the observer
pEd    = 0;                 % probability of making a memory error at each observation
pEp    = 0;                 % probability of making a memory error at each observation
patlen = 10;                % depth of the rules' tree to explore
pR     = 'Size-principle';  % the prior probability of each rule depends on its length
pT     = 'Bayes-Laplace';   % the prior over statistics to be learnt
pJ     = 'Uniform';         % prior over change point's position
comp   = 'all';             % compute after each observation
scale  = 'log';             % scale of the model evidence
pgrid  = [];                % precision of the posterior over theta
verb   = 0;                 % do not output messages in the command window

% Prepare the output variable
mIO = cell(nOrd,nCond,nSub); % conditions x subjects cell matrix with IO's inference

%% SIMULATIONS
%  ===========

% For each subject
for iSub = 1:nSub

    % For each condition
    for iCond = 1:nCond
        fprintf('- Running the IO on sequence #%2.0f/%2.0f from subject #%2.0f/%2.0f...\n', ...
            iCond, nCond, iSub, nSub);

        % Get the sequence
        seq = G{conds(iCond),iSub}.Seq;

        % For each order of the Markov chain
        for iOrd = 1:nOrd

            % Define the statistics to be learned by the probabilistic model
            stat  = sprintf('Chain%1.0f', orders(iOrd));

            % Run the observer with these options
            io = Emergence_IO_FullIO(seq, pEd, pEp, patlen, stat, pR, pT, pJ, comp, scale, pgrid, verb);
            mIO{iOrd,iCond,iSub} = cat(1, io.pYgMsp, io.pYgMsd, io.pYgMss)';
        end

        % Save temporary file
        save('Emergence_IO_AlternativeDeterHypothesis_tmp.mat', 'mIO');
    end
end

% Save file containing all the simulations and delete temporary file
save('Emergence_IO_AlternativeDeterHypothesis.mat', 'mIO');
delete('Emergence_IO_AlternativeDeterHypothesis_tmp.mat');

%% COMPUTE HYPOTHESES LIKELIHOOD
%  =============================

% Get sequence likelihood under different hypotheses
pYgMsp = cell2mat(cellfun(@(x) x(:,1), mIO(1,:,:), 'UniformOutput', 0)); % 1st-order Markov chain
pYgMsd = cell2mat(cellfun(@(x) x(:,2), mIO(1,:,:), 'UniformOutput', 0)); % pattern learner
pYgMss = cell2mat(cellfun(@(x) x(:,3), mIO(1,:,:), 'UniformOutput', 0)); % fully-stochastic
pYgMsc = cell2mat(cellfun(@(x) x(:,1), reshape(permute(mIO, ...          % higher-order Markov chains
    [2,3,1]), [1,nCond,nSub,nOrd]), 'UniformOutput', 0));

% Combine the together different deterministic hypothesis
pYgMsd = cat(4, pYgMsc, pYgMsd);

% Deduce posterior probabilities of the different hypotheses
pY = exp(pYgMss) + exp(pYgMsp) + exp(pYgMsd); % marginal likelihood
pMspgY = exp(pYgMsp) ./ pY; % posterior probability of the probabilistic hypothesis
pMsdgY = exp(pYgMsd) ./ pY; % posterior probability of the deterministic hypotheses
pMssgY = exp(pYgMss) ./ pY; % posterior probability of the fully-stochastic hypothesis
pMgY = cat(5, pMspgY, pMsdgY, pMssgY); % concatenate posterior probabilities

%% DISPLAY DYNAMICS AROUND CHANGE POINT
%  ====================================

% Window around change point to look into (in number of observations)
ObsWin = -50:50;

% Create an hybrid colormap
ColMap = [cbrewer2('Blues', nOrd+2); tricol(2,:)];
ColMap = ColMap(3:end,:);

% Define positions of figures' subplots
SpPos = {[6,7,8,10,11,12,14,15,16,2,1], ...
         [1,6,10,14,7,11,3,8,12,16]};

% For each type of regularity
for iReg = 1:2
    nRegCond = numel(cidx{iReg});
    
    % Select data correspinding to the current regularity type
    WinIO = mat2cell(pMgY(:,cidx{iReg},:,:,iReg), ...
        N, ones(nRegCond,1), ones(1,nSub), nOrd+1);
    
    % Select sequence likelihood around change point's position
    CP = cellfun(@(x) x.Jump - 1/2, G(cidx{iReg},:), 'UniformOutput', 0); 
    WinIO = cell2mat(cellfun(@(x,i) x(i+ObsWin,:,:,:), ...
        WinIO, reshape(CP, [1,size(CP)]), 'UniformOutput', 0));
    
    % Average over (pseudo-)subjects
    Avg = squeeze(mean(WinIO, 3, 'OmitNaN'));
    Err = squeeze(sem(WinIO, 3));
    
    % Prepare a new window
    figure('Position', [200+760*(iReg-1) 565 760 540]);
    
    % For each sequence with a deterministic regularity
    for iCond = 1:nRegCond
        subplot(4,4,SpPos{iReg}(iCond));
        
        % For each type of deterministic hypothesis
        for iOrd = 1:nOrd+1
            
            % Display averaged hypothesis likelihood
            plotMSEM(ObsWin, Avg(:,iCond,iOrd), Err(:,iCond,iOrd), ...
                0.15, ColMap(iOrd,:), ColMap(iOrd,:), 2);
        end
        
        % Display the position of the change point
        plot(zeros(1,2), [0,1], 'k-');
        
        % Display the detection threshold
        plot(ObsWin([1,end]), ones(1,2)./2, '-', 'Color', g);
        
        % Customize the axes
        axis([ObsWin([1,end]),0,1]);
        set(gca, 'Box', 'Off');
        
        % Add some text labels
        if SpPos{iReg}(iCond) == 1, xlabel('Position w.r.t. change point'); end
        if     iReg == 1, ylabel('p(Hp|y)'); title(pr{iCond});
        elseif iReg == 2, ylabel('p(Hd|y)'); title(dr{iCond});
        end
    end
end

%% QUANTITATIVELY COMPARE THE DIFFERENT TYPE OF DETERMINISTIC HYPOTHESES
%  =====================================================================

% Prepare the output variable
rho = NaN(nSub,nOrd+1,nCond);

% For each type of deterministic hypothesis
for iOrd = 1:nOrd+1
    
    % For each sequence
    for iSeq = 1:nCond
        
        % For each subject
        for iSub = 1:nSub
            
            % Get predictions of the ideal observer model
            Xi = squeeze(pMgY(:,iSeq,iSub,iOrd,:));
            Xi(1,:) = [0,0,1];
            Xi = Xi(:);
            
            % Get reported probabilities by the subject
            Yi = G{iSeq,iSub}.BarycCoord;
            Yi = Yi(:);
            
            % Measure the Pearson correlation between the two
            rho(iSub,iOrd,iSeq) = corr(Xi, Yi);
        end
    end
end

% Average over sequences
avgcoef = mean(rho, 3);

% Compare correlation coefficients obtained for the fully-deterministic
% hypothesis with correlation coefficients obtained with
% pseudo-deterministic hypotheses (high-order Markov chains)
[~,pval,ci,stats] = ttest(avgcoef(:,end) - avgcoef(:,1:end-1));
Emergence_PrintTstats(pval,ci,stats);

% Create a new window
figure('Position', [682 71 560 420]); hold('on');

% Display group-averaged correlation coefficients
for iOrd = 1:nOrd+1
    m = mean(avgcoef(:,iOrd));
    s = sem(avgcoef(:,iOrd));
    bar(iOrd, m, 'FaceColor', ColMap(iOrd,:));
    plot(repmat(iOrd,1,2), m+[-1,1].*s, 'k-');
end

% Customize the axes
set(gca, 'Box', 'Off', 'XTick', [], 'Ylim', [0.6, 0.8]);

% Add some text labels
xlabel('Deterministic hypotheses');
ylabel('Correlation coefficient');
