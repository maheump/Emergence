% This script looks at the dynamics of finger trajectory aroung change and
% detection points. We show that depending if we lock on the position of
% the change point or the detection point (when the beliefs in the true
% generative process crosses some threshold), the dynamics does not look
% the same even though detection dynamics look steeper in the case of
% deterministic regularities.
% 
% Copyright (c) 2018 Maxime Maheu

%% BARYCENTRIC COORDINATES LOCKED ON DETECTION/CHANGE POINT
%  ========================================================

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
cp = cell(1,2); dp = cell(1,2); lag = cell(1,2); % change and detection points
update = cell(1,2); % measures of the trajectories
fingerposwrtp = cell(2,3); % trajectories locked on different point

% For each type of regularity
for iHyp = 1:2
    
    % Get change point positions
    cp{iHyp} = cellfun(@(x) x.Jump+1/2, G(cidx{iHyp},:));
    
    % Measure important quantities
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    % Get detection point position
    lag{iHyp} = cellfun(@(x,c) Emergence_FindDetecPoint(x.BarycCoord(c:end,iHyp)), ...
        D(cidx{iHyp},:), num2cell(cp{iHyp}));
    dp{iHyp} = cp{iHyp} + lag{iHyp};
    
    % Measure the build-up of beliefs along the relevant dimension
    update{iHyp} = cellfun(@(x,c) mean(diff(x.BarycCoord(c:end,iHyp))), ...
        D(cidx{iHyp},:), num2cell(cp{iHyp}));
    
    % Remove sequences entailing undetected regularities
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % N.B. this is computed on subject such that we keep the same sequences
    % for subjects and the ideal observer
    
    % Apply those 2 criterions to select the sequences
    detecmask = (filter{iHyp} == 3);
    cp    {iHyp}(~detecmask) = NaN;
    dp    {iHyp}(~detecmask) = NaN;
    lag   {iHyp}(~detecmask) = NaN;
    update{iHyp}(~detecmask) = NaN;
    
    % Get the trajectories around the points of interest
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    % For the important points we want to lock onto
    for lock = 1:3
        
        % Beginning and ending of the window to look in (in # of samples)
        if     lock == 1 % lock on change point
            begwin = cp{iHyp} + xcp(1);
            endwin = cp{iHyp} + xcp(end);
        elseif lock == 2 % lock on detection point
            begwin = dp{iHyp} + xdp(1);
            endwin = dp{iHyp} + xdp(end);
        elseif lock == 3 % lock on end point
            begwin = repmat(N+xep(1),   size(cp{iHyp}));
            endwin = repmat(N+xep(end), size(cp{iHyp}));
        end
        
        % Make sure the window stays in between the limits of the sequence
        begwin(begwin < 1) = 1;
        endwin(endwin > N) = N;
        begwin = num2cell(begwin);
        endwin = num2cell(endwin);
        begwin(~detecmask) = {[]};
        endwin(~detecmask) = {[]};
        
        % Get trajectory (i.e. beliefs in each hypothesis) in that window
        % of interest
        fing = cellfun(@(x,b,e) x.BarycCoord(b:e,:), D(cidx{iHyp},:), ...
            begwin, endwin, 'UniformOutput', 0);
        
        % Fill with NaNs to make the trajectories the same size of each
        % other
        nSamp = numel(xXp{lock});
        fing(cellfun(@isempty, fing)) = {NaN(nSamp,3)};
        fing = cellfun(@(x) [x; NaN(nSamp-size(x,1),3)], fing, 'UniformOutput', 0);
        
        % Convert the cells into a big 3D matrix
        fingerposwrtp{iHyp,lock} = cell2mat(reshape(fing, [1, 1, size(fing)]));
    end
end

% For deterministic regularities, also express the lag in terms of the
% number of repetitions of the rule
lag2 = lag{2} ./ cellfun(@numel, dr);

% Average over sequences for each type of regularity
avgfingerwrtp = cellfun(@(x) squeeze(mean(x, 3, 'OmitNaN')), fingerposwrtp, 'UniformOutput', 0);
avglag = cellfun(@(x) mean(x, 'OmitNaN'), lag, 'UniformOutput', 0);

% Average finger trajectory over subjects
avgsubtraj = cellfun(@(x) mean(x, 3, 'OmitNaN'), avgfingerwrtp, 'UniformOutput', 0);
semsubtraj = cellfun(@(x) sem( x ,3), avgfingerwrtp, 'UniformOutput', 0);

%% DISPLAY TRAJECTORIES LOCKED TO DIFFERENT IMPORTANT POINTS
%  =========================================================

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
