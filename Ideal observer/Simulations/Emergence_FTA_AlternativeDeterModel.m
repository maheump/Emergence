
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
    end
end

%
save('Emergence_FTA_AlternativeDeterModel.mat', 'mIO');

%% COMPUTE HYPOTHESES LIKELIHOOD
%  =============================

%
pYgMsp = cell2mat(cellfun(@(x) x(:,1), mIO(1,:,:), 'UniformOutput', 0));
pYgMss = cell2mat(cellfun(@(x) x(:,3), mIO(1,:,:), 'UniformOutput', 0));

% Get sequence likelihood under different deterministic hypotheses
pYgMsdD = cell2mat(cellfun(@(x) x(:,2), mIO(1,:,:), 'UniformOutput', 0)); % pattern learner
pYgMsdP = cell2mat(cellfun(@(x) x(:,1), reshape(permute(mIO, ...          % high-order Markov chain
    [2,3,1]), [1,nCond,nSub,nOrd]), 'UniformOutput', 0));
pYgMsd = cat(4, pYgMsdP, pYgMsdD); % combine them altogether

% Compute the likelihood of the deterministic hypothesis under different 
% type of deterministic hypotheses using Bayes' rule
pMsdgY = exp(pYgMsd) ./ (exp(pYgMss) + exp(pYgMsp) + exp(pYgMsd));

%% DISPLAY DYNAMICS AROUND DETECTION
%  =================================

%%%% for proba and deter regularities

%
obswin = -40:40;
j = cellfun(@(x) x.Jump - 1/2, G(cidx{2},:), 'UniformOutput', 0);

% Select sequence likelihood around change point's position
toto = mat2cell(pMsdgY, N, ones(nCond,1), ones(1,nSub), nOrd+1);
winIO = cellfun(@(x,i) x(i+obswin,:,:,:), toto, reshape(j, [1,size(j)]), 'UniformOutput', 0);
winIO = cell2mat(winIO);

% Average over subjects
avg = squeeze(mean(winIO, 3));
err = squeeze(sem(winIO, 3));

% Create an hybrid colormap
cmap = [cbrewer2('Blues', nOrd+2); tricol(2,:)];
cmap = cmap(3:end,:);

% Prepare a new window
figure('Position', [434 363 1000 430]);
pos = [1,6,10,14,7,11,3,8,12,16];

% For each sequence with a deterministic regularity
for iCond = 1:nCond
    subplot(4,4,pos(iCond));
    
    % For each type of deterministic hypothesis
    for iOrd = 6%1:nOrd+1
        plotMSEM(obswin, avg(:,iCond,iOrd), err(:,iCond,iOrd), ...
            1/4, cmap(iOrd,:), cmap(iOrd,:), 2);
    end
    
    % Display the position of the change point
    plot(zeros(1,2), [0,1], 'k-');
    
    % Display color coding
    colormap(cmap);
    colorbar;
    
    % Customize the axes
    axis([obswin([1,end]),0,1]);
    set(gca, 'Box', 'Off');
    
    % Add some text labels
    if iCond == 1 && iOrd == 1
        xlabel('Position w.r.t. change point');
        ylabel('p(Hd|y)');
    end
    title(dr{iCond});
end

%% COMPARE BASELINE AND DETECTION VELOCITY
%  =======================================

% Pre-change point baseline
% ~~~~~~~~~~~~~~~~~~~~~~~~~

%
fun = {@(p,i) mean(p, 1, 'OmitNaN'), ...
       @(p,i) mean(diff(p, 1, 1), 1, 'OmitNaN')};
idx = {'1:i-1', 'i:N'};
detecmask = cellfun(@(x) x.Questions(2) == 2, G(cidx{2},:));

figure;
for mes = 1:2

    io_mes = squeeze(cell2mat(cellfun(@(x,i) fun{mes}(x(i:N,:,:,:)), ...
        toto, reshape(j, [1,size(j)]), 'UniformOutput', 0)));
    sub_mes = cellfun(@(x,i) fun{mes}(x.BarycCoord(i:N,2)), ...
        G(cidx{2},:), j, 'UniformOutput', 1);
    
    % Remove data from undetected regularities
    io_mes(repmat(~detecmask, [1,1,nOrd+1])) = NaN;
    sub_mes(~detecmask) = NaN;
    
    % Average over sequences
    io_mes  = mean(io_mes, 1, 'OmitNaN');
    sub_mes = mean(sub_mes, 1, 'OmitNaN');

    % Average over subjects
    avg_mes = [squeeze(mean(io_mes, 2)); mean(sub_mes)];
    sem_mes = [squeeze(sem(io_mes, 2)); sem(sub_mes)];

    cmap = [cbrewer2('Blues', nOrd+2); tricol(2,:); g];
    cmap = cmap(3:end,:);

    %
    subplot(1,2,mes)
    for iOrd = 1:nOrd+2
        plot(repmat(iOrd, 1, 2), avg_mes(iOrd) + sem_mes(iOrd) .* [-1,1], 'k-'); hold('on');
        plot(iOrd, avg_mes(iOrd), 'ko', 'MarkerFaceColor', cmap(iOrd,:), 'MarkerSize', 10);
    end
    
    % Customize the axes
    xlim([0,nOrd+3]);
    if     mes == 1, ylabel('Baseline');
    elseif mes == 2, ylabel('Velocity');
    end
end
