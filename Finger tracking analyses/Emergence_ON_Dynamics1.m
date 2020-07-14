% This script looks at the dynamics of finger trajectory aroung change and
% detection points. We show that depending if we lock on the position of
% the change point or the detection point (when the beliefs in the true
% generative process crosses some threshold), the dynamics does not look
% the same even though detection dynamics look steeper in the case of
% deterministic regularities.
% 
% Copyright (c) 2020 Maxime Maheu

%% GET TRAJECTORIES LOCKED ON PARTICULAR POINTS
%  ============================================

% Define properties of windows to look into
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% The number of observation to consider
nSamp = 30;

% Define x-axes for trajectories...
xcp =  0:(nSamp*2); % ... locked on change    point
xdp = -nSamp:nSamp; % ... locked on detection point
xep = -(nSamp/2):0; % ... locked on end       point
xXp = {xcp,xdp,xep};

% Prepare the output variables
cp = cell(1,2); dp = cell(1,2); ep = cell(1,2); % points
lag = cell(1,2); update = cell(1,2); % measures of the trajectories
fingerposwrtp = cell(2,3); % trajectories locked on different point

% For each type of regularity
for iHyp = 1:2
    
    % Get end point positions
    ep{iHyp} = repmat(N, numel(cidx{iHyp}), nSub);
    
    % Get change point positions
    cp{iHyp} = cellfun(@(x) x.Jump, G(cidx{iHyp},:));
    
    % Measure important quantities
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    % Lock trajectories on change point
    lockbel = cellfun(@(x,c) Emergence_LockOnPoint(x.BarycCoord(:,iHyp),c,[0,N]), ...
            D(cidx{iHyp},:), num2cell(cp{iHyp}), 'uni', false);
    
    % Get detection point position
    lag{iHyp} = cellfun(@(p) Emergence_FindDetecPoint(p), lockbel);
    dp{iHyp} = cp{iHyp} + lag{iHyp};
    
    % Measure the build-up of beliefs along the relevant dimension
    update{iHyp} = cellfun(@(p) mean(diff(p), 'OmitNaN'), lockbel);
    
    % Remove sequences entailing undetected regularities
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % N.B. this is computed on subject such that we keep the same sequences
    % for subjects and the ideal observer
    
    % Apply those 2 criterions to select the sequences
    detecmask = (filter{iHyp} == 3);
    cp    {iHyp}(~detecmask) = NaN;
    dp    {iHyp}(~detecmask) = NaN;
    ep    {iHyp}(~detecmask) = NaN;
    lag   {iHyp}(~detecmask) = NaN;
    update{iHyp}(~detecmask) = NaN;
    
    % Get the trajectories around the points of interest
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    % For the important points we want to lock onto
    for lock = 1:3
        
        % Get observations to lock onto
        if     lock == 1, lockobs = num2cell(cp{iHyp}); % lock on change point
        elseif lock == 2, lockobs = num2cell(dp{iHyp}); % lock on detection point
        elseif lock == 3, lockobs = num2cell(ep{iHyp}); % lock on end point
        end
        
        % Get trajectory in that window of interest
        fing = cellfun(@(p,c) Emergence_LockOnPoint(p.BarycCoord, c, xXp{lock}), ...
            D(cidx{iHyp},:), lockobs, 'uni', 0);
        
        % Convert the cells into a big 3D matrix
        fingerposwrtp{iHyp,lock} = cell2mat(reshape(fing, [1,1,size(fing)]));
    end
end

% For deterministic regularities, also express the lag in terms of the
% number of repetitions of the rule
lag2 = lag{2} ./ cellfun(@numel, dr);

%% GRAND AVERAGE OF THE POSTERIOR PROBABILITY IN THE RELEVANT HYPOTHESIS
%  =====================================================================

% Average over sequences for each type of regularity
avgfingerwrtp = cellfun(@(x) squeeze(mean(x, 3, 'OmitNaN')), fingerposwrtp, 'uni', 0);
avglag = cellfun(@(x) mean(x, 'OmitNaN'), lag, 'uni', 0);

% Average finger trajectory over subjects
avgsubtraj = cellfun(@(x) mean(x, 3, 'OmitNaN'), avgfingerwrtp, 'uni', 0);
semsubtraj = cellfun(@(x) sem( x ,3), avgfingerwrtp, 'uni', 0);

% Display trajectories in the triangular space
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Prepare a new window
figure('Position', [1 855 700 250]);

% For change-, detection- and end- points
for lock = 1:3
    subplot(1,3,lock);
    
    % Display the triangle
    Emergence_PlotTriInfo(tricc, tricol);
    
    % Display a useful custom grid on the triangle
    for k = 1:3
        lgd = Emergence_PlotGridOnTri(2, k, tricol(k,:), tricc);
        set(lgd(1,k), 'Color', [get(lgd(1,k), 'Color'), 3/4]);
        lgd = Emergence_PlotGridOnTri(3, k, tricol(k,:), tricc);
        for kk = 1:3
            set(lgd(kk,k), 'LineStyle', '--', 'Color', [get(lgd(kk,k), 'Color'), 3/4]);
        end
    end
    Emergence_PlotGridOnTri(2, [1,2], 'k', tricc);
    
    % Display the trajectories
    for iHyp = 1:2
        cc = avgsubtraj{iHyp,lock} * tricc; % cartesian coordinates
        plot(cc(:,1), cc(:,2), '.-', 'Color', tricol(iHyp,:), ...
            'LineWidth', 1, 'MarkerSize', 10);
        
        % PLot point on which trajectories are locked
        i = xXp{lock} == 0;
        plot(cc(i,1), cc(i,2), 'ko');
        
        % Make sure the triangle is equilateral and we don't see the axes
        axis('equal'); axis('off');
    end
end

% Display corresponding Barycentric coordinates
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Define subplots' width
spw = ceil(6.* (cellfun(@numel, xXp) ./ sum(cellfun(@numel, xXp))));
cspw = cat(1, [1,cumsum(spw(1:end-1))+1], cumsum(spw));

% Prepare a new window
figure('Position', [1 381 700 400]);

% For change- and detection- points
pt = {'change', 'detection', 'end'};
for lock = 1:3
    xval = xXp{lock};
    
    % For each type of regularity
    for iHyp = 1:2
        subplot(2, sum(spw), (cspw(1,lock):cspw(2,lock))+sum(spw)*(iHyp-1)); hold('on');
        
        % Draw some grid lines
        plot(xval([1,end]),    ones(1,2)./2, '-',  'Color', g, 'LineWidth', 1/2);
        plot(xval([1,end]),    ones(1,2)./3, '--', 'Color', g, 'LineWidth', 1/2);
        plot(xval([1,end]), 2.*ones(1,2)./3, '--', 'Color', g, 'LineWidth', 1/2);
        
        % Display distribution of change/detection point positions
        if     lock == 1, dist = lag{iHyp}(:);        klim = [-2,Inf];  idx = xval >= 0;
        elseif lock == 2, dist = -lag{iHyp}(:);       klim = [-Inf,2];  idx = xval <  0;
        elseif lock == 3, dist = -(N - dp{iHyp}(:));  klim = [-Inf,2];  idx = xval <= N;
        end
        [fout0,xout] = ksdensity(dist, ...          % which distribution to plot
            interp(xval(idx), 10), ...              % grid of positions
            'BandWidth',            8, ...          % bandwidth of the kernel smoothing window 
            'Support',              klim, ...       % restrict the kernel to a certain range of values
            'BoundaryCorrection',   'Reflection');	% type of correction for the boundaries
        fill([xout(1),xout,xout(end)], [0,fout0.*10,0], 'k', 'FaceColor', g);
        
        % Draw barycentric coordinates
        lt = repmat({'--'}, 1, 3);
        lt{iHyp} = '-';
        for iDim = 1:3
            plotMSEM(xval, avgsubtraj{iHyp,lock}(:,iDim), ...
                           semsubtraj{iHyp,lock}(:,iDim), ...
                0.15, tricol(iDim,:), tricol(iDim,:), 1+1*(iDim == iHyp), 1, lt{iDim}, 'none');
        end
        
        % Customize the axes
        set(gca, 'Box', 'Off', 'XLim', xval([1,end]), 'YLim', [0,1]);
        set(gca, 'YAxisLocation', 'Origin');
        if lock == 2, set(gca, 'YTickLabel', {}); end
        
        % Add some text labels
        if iHyp == 2, xlabel(sprintf('Obs.# w.r.t. %s pt', pt{lock})); end
        if lock == 1
            ylabel('Posterior beliefs p(M_i|y)');
            title(sprintf('%s regularities', proclab{iHyp}));
        end
        
        % Display the distribution of average detection points
        if lock == 1 % locked on change point
            ax = get(gca, 'Position');
            axes('Position', [ax(1) ax(2)+0.8*ax(4) ax(3) 0.2*ax(4)]);
            Emergence_PlotSubGp(avglag{iHyp}, tricol(iHyp,:));
            axis([0, 2, xval([1,end])]); axis('off');
            view([90,90]);
            Emergence_DispStatTest(avglag{iHyp});
        end
        
    end
end

%% POSTERIOR PROBABILITY IN THE RELEVANT HYPOTHESIS FOR EACH INDIVIDUAL TRIAL
%  ==========================================================================

% Get posterior probability for each individual sequence
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Prepare output variables
BEL = cell(2,2,2);
CPP = cell(2,2);

% For each type of regularity
for iHyp = 1:2
    
    % Get sequences that were correctly labelled
    detecmask = (filter{iHyp} == 1 | filter{iHyp} == 3);
    
    % Order according to change point's position and detected/undetected
    cp = cellfun(@(x) x.Jump, G(cidx{iHyp},:), 'uni', 1);
    [sortedcp, idxcp] = sortrows([cp(:), detecmask(:)], [2,1]);
    
    % Get the beliefs in the corresponding (correct) hypothesis ordered
    % according to the position the change point
    belincorhyp = cellfun(@(x) x.BarycCoord(:,iHyp)', ...
        D(cidx{iHyp},:), 'uni', 0);
    
    % Separately for sequences that were detected/undetected
    for du = 1:2
        idx = idxcp(sortedcp(:,2) == (du-1));
        bel = belincorhyp(idx)';
        cpp = num2cell(cp(idx))';
        CPP{du,iHyp} = cell2mat(cpp)';
        BEL{du,iHyp,1} = cellfun(@(p,c) [p(1:c-1), NaN(1,N-c+1)], bel, cpp, 'uni', 0);
        BEL{du,iHyp,2} = cellfun(@(p,c) [NaN(1,c-1),   p(c:end)], bel, cpp, 'uni', 0);
    end
end

% Display dynamics of beliefs in the correct hypothesis
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Combine all sequences together
BEL = cellfun(@(x) cell2mat(x'), BEL, 'uni', 0);

% Prepare a new window
figure('Position', [702 381 460 400]);

% For each type of regularity
for iHyp = 1:2 % probabilistic/deterministic
    for ba = 1:2 % before/after
        for du = 1:2 % detected/undetected
            sp = subplot(1,2,iHyp); del = 10;
            
            % Display beliefs in the correct hypothesis for all sequences
            c = BEL{du,iHyp,ba};
            x1 = (1:N)+del*(ba-1);
            ns = size(c,1);
            y1 = (1:ns) - ns*(du==1) - del*(du==1);
            imagesc(x1, y1, c, 'AlphaData', ~isnan(c)); hold('on'); 
            
            % Display change point position
            p = CPP{du,iHyp};
            for t = 0:1
                x2 = [p(1);   p;  p(end)] - 1/2 + del * t;
                y2 = [y1(1)-1; y1' ;y1(end)] + 1/2;
                stairs(x2, y2, '-', 'Color', g, 'LineWidth', 1);
            end
            
            % Customize the axes
            axis('xy'); axis('tight'); axis('on'); ScaleAxis('y');
            colormap(sp, Emergence_Colormap('Inferno'));
            
            % Add some text labels
            xlabel('Observation #');
            ylabel({'Sequence # (sorted by change point''s position)'});
            title(sprintf('%s sequences', proclab{iHyp}));
        end
    end
end

%% DISPLAY BELIEF UPDATE FOR THE TWO TYPES OF REGULARITIES
%  =======================================================

% Average update across sequences for each type of regularity
avgslope = cell2mat(cellfun(@(x) mean(x,'OmitNaN'), update, 'uni', false)');

% Perform a t-test
[h,pval,tci,stats] = ttest(diff(avgslope)');
Emergence_PrintTstats(pval,tci,stats);

% Prepare a new window
figure('Position', [702 905 120 200]);

% Display the paired difference
Emergence_PlotSubGp(avgslope', tricol(1:2,:));

% Customize the axes
set(gca, 'XLim', [0,3], 'YLim', [2,12].*1e-3,'XTick', [], 'XColor', 'none');

% Add a text label
ylabel('Belief update');
