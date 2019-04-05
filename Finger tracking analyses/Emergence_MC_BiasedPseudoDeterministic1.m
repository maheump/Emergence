% FALSE ALARMS
% 
% Copyright (c) 2018 Maxime Maheu

%% INITIALIZATION
%  ==============

% Load simulations
SimuType = 'PseudoDeterministic'; % 'PseudoDeterministic' or 'BiasedPseudoDeterministic'
Emergence_MC_ModelSimulations;

%% TRIANGULAR HISTOGRAMS
%  =====================

% Perform the analysis
% ~~~~~~~~~~~~~~~~~~~~

% Select fully-stochastic parts of sequences
idxtrimap = Emergence_SelectFullyStochSeq(G, filter, 1);

% Prepare the output variables
trajmap = cell(1,nMod);
CorCoef = [];

% For each model
for iMod = 1:nMod
    
    % Compute the 2D histogram
    points = cellfun(@(x,y) x(y,:), pMgY(:,:,iMod), idxtrimap, 'UniformOutput', 0);
    [trajmap{iMod}, ~, xgrid, ygrid, mask] = Emergence_TriHist(points);
end

% Display triangular histrograms
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Prepare the window
figure('Position', [201 453 1520 250]);

% For each density map
for iMod = 1:nMod
    subplot(1,nMod,iMod);
    
    % Display the density maps of finger's position
    imagesc(xgrid, ygrid, trajmap{iMod}, 'AlphaData', mask); hold('on');
    contour(xgrid, ygrid, trajmap{iMod}, 5, 'k-', 'LineWidth', 1/2);
    image(xgrid, ygrid, repmat(~mask, [1,1,3]), 'AlphaData', ~mask);
    
    % Display an empty triangle
    fill(tricc(:,1), tricc(:,2), 'k', 'FaceColor', 'None', 'LineWidth', 2);
    
    % Overlap the grid corresponding to the resoution of the marginal
    % histrograms
    Emergence_PlotGridOnTri(10);
    
    % Customize the axes
    set(gca, 'LineWidth', 1, 'FontSize', 15);
    axis('xy'); axis('off'); axis('equal');
    axis([0,1,0,sqrt(3)/2]);
    colormap(cbrewer2('Greens', 2000));
    
    % Add some text labels
    title(options{4,iMod});
end

%%








%% DISTANCE TO AVERAGE POSITION
%  ============================

% Select fully-stochastic parts of sequences
idxtrimap = Emergence_SelectFullyStochSeq(G, filter, 1);

% 
scns = cellfun(@(x,i) mean(x(i,:)*tricc, 1)', pMgY, repmat(idxtrimap, [1,1,nMod]), 'UniformOutput', 0);
subs = cellfun(@(x,i) mean(x.BarycCoord(i,:)*tricc, 1)', G, idxtrimap, 'UniformOutput', 0);

% Transform into a matrix
scns = cell2mat(reshape(scns, [1,size(scns)]));
subs = cell2mat(reshape(subs, [1,size(subs)]));

% Average over sequences
seqavgmod = squeeze(mean(scns, 2));
seqavgsub = squeeze(mean(subs, 2));

%
figure;

%
Emergence_PlotTriInfo;

for iMod = 1:nMod
    
    m = mean(seqavgmod(:,:,iMod), 2);
    s = sem( seqavgmod(:,:,iMod), 2);
    
    plot(m(1), m(2), 'ko', 'MarkerFaceColor', modc(iMod,:));
    
end

m = mean(seqavgsub, 2);
s = sem( seqavgsub, 2);

plot(m(1) + s(1).*[-1,1], repmat(m(2),1,2), 'k-');
plot(repmat(m(1),1,2), m(2) + s(2).*[-1,1], 'k-');
plot(m(1), m(2), 'ko', 'MarkerFaceColor', 'r');

axis([1/4,3/4,0,1/2]);

%%

% Measure squared error between subjects and models
MSE = cellfun(@(x,y) sum((x-y).^2), repmat(subs, [1,1,nMod]), scns);
MSE = log(MSE);

% Average over sequences
avgMSE = squeeze(mean(MSE, 1, 'OmitNaN'));

%
groups = {1:nMod/2, nMod/2+1:nMod, nMod+1};
nG = numel(groups);
col = cbrewer2('PuRd', 5);
col = [col(2:3,:); tricol(2,:)];

% Prepare a new figure
figure('Position', [1 918 279 187]);
l = NaN(1,nG);

% 
for i = 1:nG
    
    x = orders(groups{i}) -1/6 + (1/3 .* double(i == 2));
    if isnan(x), x = nMod/2 + 2; end
    y = avgMSE(:,groups{i});
    m = mean(y);
    s = sem(y);
    
    % Display the error bars
    plot(repmat(x, [2,1]), m+s.*[-1;1], 'k-', 'LineWidth', 1/2); hold('on');
    
    % Display 
    plot(x, m, 'k-', 'LineWidth', 2);
    l(i) = plot(x, m, 'ko', 'MarkerFaceColor', col(i,:));
end 

% Customize the axes
set(gca, 'Box', 'Off', 'XTick', [min(orders):max(orders), nMod/2 + 2], ...
    'XTickLabel', cat(2, num2cell(min(orders):max(orders)), 'D'), ...
    'XLim', [min(orders)-1, nMod/2 + 3]);

% Add some text labels
legend(l(1:nG-1), {'Flat prior', 'Biased prior'}, 'Location', 'NorthWest');
xlabel('(Pseudo-) Deterministic hypothesis');
ylabel('Log-difference in average position');

%%

% Scale the colormap
toto = cellfun(@(x) x ./ max(x(:)), trajmap, 'UniformOutput', 0);

% 
%modc = hsv(nMod+1);
%modc(end,:) = tricol(2,:);

figure('Position', [1 857 330 248]); hold('on');

for iMod = 1:nMod
    contour(xgrid, ygrid, toto{iMod}, [0.5,1], '-', ...
        'EdgeColor', modc(iMod,:), 'LineWidth', 1/2);
end
for iMod = 1:nMod
    contourf(xgrid, ygrid, toto{iMod}, [0.98,1], ...
        'EdgeColor', 'k', 'FaceColor', modc(iMod,:), ...
        'LineWidth', 1);
end

% Mask what is outside the triangle
image(xgrid, ygrid, repmat(~mask, [1,1,3]), 'AlphaData', ~mask);

% Draw triangle
tr = Emergence_PlotTriInfo(tricc, tricol);
set(tr, 'LineWidth', 2);

% Customize the axes
set(gca, 'LineWidth', 1, 'FontSize', 15);
axis([0.3,0.8,0,0.6]);
axis('xy'); axis('off'); axis('equal');



%%

%
iHyp = 2;
avgPbelS = cellfun(@(p,i) mean(p.BarycCoord(i,:)), G, idxtrimap, 'UniformOutput', 0);
avgPbelIO = cellfun(@(p,i) mean(p(i,iHyp)), pMgY, repmat(idxtrimap, [1,1,nMod+1]));

figure;

m = mean(avgPbelS, 1);
bar(1, mean(m), 'FaceColor', g); hold('on');

for iOrd = 1:nMod+1
    m = mean(avgPbelIO(:,:,iOrd), 1);
    bar(2+iOrd, mean(m), 'FaceColor', modc(iOrd,:));
end
