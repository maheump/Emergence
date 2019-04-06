% This script looks at the dynamics of finger trajectory aroung change and
% detection points. We show that depending if we lock on the position of
% the change point or the detection point (when the beliefs in the true
% generative process crosses some threshold), the dynamics does not look
% the same even though detection dynamics look steeper in the case of
% deterministic regularities.
% 
% Copyright (c) 2018 Maxime Maheu

%% GET TRAJECTORIES LOCKED ON PARTICULAR POINTS AND ASSESS THEIR PROPERTIES
%  ========================================================================

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
    cp{iHyp} = cellfun(@(x) x.Jump+1/2, G(cidx{iHyp},:));
    
    % Measure important quantities
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    % Lock trajectories on change point
    lockbel = cellfun(@(x,c) Emergence_LockOnPoint(x.BarycCoord(:,iHyp),c,[0,N]), ...
            D(cidx{iHyp},:), num2cell(cp{iHyp}), 'UniformOutput', false);
    
    % Get detection point position
    lag{iHyp} = cellfun(@(p) Emergence_FindDetecPoint(p), lockbel);
    dp{iHyp} = cp{iHyp} + lag{iHyp};
    
    % Measure the build-up of beliefs along the relevant dimension
    update{iHyp} = cellfun(@(p) sum(diff(p), 'OmitNaN'), lockbel);
    
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
            D(cidx{iHyp},:), lockobs, 'UniformOutput', 0);
        
        % Convert the cells into a big 3D matrix
        fingerposwrtp{iHyp,lock} = cell2mat(reshape(fing, [1,1,size(fing)]));
    end
end

% For deterministic regularities, also express the lag in terms of the
% number of repetitions of the rule
lag2 = lag{2} ./ cellfun(@numel, dr);

%% DISPLAY AVERAGE LIKELIHOODS IN THE RELEVANT HYPOTHESIS
%  ======================================================

% Average over sequences for each type of regularity
avgfingerwrtp = cellfun(@(x) squeeze(mean(x, 3, 'OmitNaN')), fingerposwrtp, 'UniformOutput', 0);
avglag = cellfun(@(x) mean(x, 'OmitNaN'), lag, 'UniformOutput', 0);

% Average finger trajectory over subjects
avgsubtraj = cellfun(@(x) mean(x, 3, 'OmitNaN'), avgfingerwrtp, 'UniformOutput', 0);
semsubtraj = cellfun(@(x) sem( x ,3), avgfingerwrtp, 'UniformOutput', 0);

% Display trajectories in the triangular space
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Properties of the triangle
tricc  = [0, sqrt(3)/2; 1, sqrt(3)/2; 1/2, 0];
tricol = [066 146 198; 239 059 033; 065 171 093] ./ 255;

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

% Save the figure
if isfield(D{1}, 'Seq'), save2pdf(fullfile(ftapath, 'figs', 'F_Dyn_TriS.pdf'));
else, save2pdf(fullfile(ftapath, 'figs', 'F_Dyn_TriIO.pdf'));
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

% Save the figure
if isfield(D{1}, 'Seq'), save2pdf(fullfile(ftapath, 'figs', 'F_Dyn_CoordS.pdf'));
else, save2pdf(fullfile(ftapath, 'figs', 'F_Dyn_CoordIO.pdf'));
end

%% DISPLAY INDIVIDUAL TRIALS LIKELIHOODS IN THE RELEVANT HYPOTHESIS
%  ================================================================

% Prepare output variable
cp          = cell(1,2);
sortedcp    = cell(1,2);
idxcp       = cell(1,2);
belincorhyp = cell(1,2);

% For each type of regularity
for iHyp = 1:2
    
    % Get sequences that were correctly labelled
    detecmask = (filter{iHyp} == 1 | filter{iHyp} == 3);
    
    % Order according to change point's position and detected/undetected
    cp{iHyp} = cellfun(@(x) x.Jump, G(cidx{iHyp},:), 'UniformOutput', 1);
    [sortedcp{iHyp}, idxcp{iHyp}] = sortrows([cp{iHyp}(:), ...
        detecmask(:)], [2,1]);
    
    % Get the beliefs in the corresponding (correct) hypothesis ordered
    % according to the position the change point
    belincorhyp{iHyp} = cellfun(@(x) x.BarycCoord(:,iHyp)', ...
        D(cidx{iHyp},:), 'UniformOutput', 0);
end

% Prepare a new window
figure('Position', [702 381 230 400]);
cmapcol = {'Blues', 'Reds'};

% For sequences with a probabilistic/deterministic regularity
for iHyp = 1:2
    
    % Display the change in beliefs as a heatmap 
    sp = subplot(2,1,iHyp);
    bel = belincorhyp{iHyp}(idxcp{iHyp});
    imagesc(cell2mat(bel)); hold('on');
    
    % Customize the colormap
    colorbar('Location', 'EastOutside'); caxis([0,1]);
    colormap(sp, Emergence_Colormap(cmapcol{iHyp}));
    
    % Display the position of the change points
    for d = [0,1]
        x = cp{iHyp}(idxcp{iHyp}(sortedcp{iHyp}(:,2) == d));
        y = find(sortedcp{iHyp}(:,2) == d);
        stairs([x(1);x], [y(1)-1;y], 'k-');
    end
	
    % Display limits between detected and undetected sequences
    lim = find(abs(diff(sortedcp{iHyp}(:,2))) == 1) + 1/2;
    plot([0,N+1], repmat(lim, 1, 2), 'k-');
    
    % Add some text labels
    axis('xy'); set(gca, 'XTick', [1, get(gca, 'XTick')], 'YTick', [1,50]);
    xlabel('Observation #');
    ylabel({'Sequence # (sorted by', 'change point''s position)'});
    title(sprintf('%s sequences', proclab{iHyp}));
end

% Save the figure
if isfield(D{1}, 'Seq'), save2pdf(fullfile(ftapath, 'figs', 'F_Dyn_MapS.pdf'));
else, save2pdf(fullfile(ftapath, 'figs', 'F_Dyn_MapIO.pdf'));
end

%% DISPLAY BELIEF UPDATE FOR THE TWO TYPES OF REGULARITIES
%  =======================================================

% Average update across sequences for each type of regularity
avgslope = cell2mat(cellfun(@(x) mean(x,'OmitNaN'), update, 'UniformOutput', false)');

% Perform a t-test
[h,pval,tci,stats] = ttest(diff(avgslope)');
Emergence_PrintTstats(pval,tci,stats);

% Prepare a new window
figure('Position', [702 905 120 200]);

% Display the paired difference
Emergence_PlotSubGp(avgslope', tricol(1:2,:));

% Customize the axes
set(gca, 'XLim', [0,3], 'XTick', [], 'XColor', 'none');

% Add a text label
ylabel('Belief update');

% Save the figure
if isfield(D{1}, 'Seq'), save2pdf(fullfile(ftapath, 'figs', 'F_Dyn_SlopeS.pdf'));
else, save2pdf(fullfile(ftapath, 'figs', 'F_Dyn_SlopeIO.pdf'));
end
