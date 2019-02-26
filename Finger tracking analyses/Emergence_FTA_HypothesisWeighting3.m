% This script compares the relative weighting (in particular during
% fully-stochastic parts of sequences) between probabilistic and
% deterministic hypotheses from the ideal observer and the human subjects.
% In particular, it shows that the weighting between probabilistic and
% deterministic hypotheses from the ideal observer is predictive of the
% weighting from the subjects, even after taking into account several 
% confounded variables in the regression.
% 
% Copyright (c) 2018 Maxime Maheu

%% DISPLAY TRIANGULAR COLORMAPS OF PREDICTORS
%  ==========================================

% Triangle specifications
tricc = [0, sqrt(3)/2; 1, sqrt(3)/2; 1/2, 0];

% Define horizontal and vertical cartesian grids
j = 0.001;
xgrid = tricc(1,1):j:tricc(2,1);
ygrid = tricc(3,2):j:tricc(1,2);
nx = numel(xgrid);
ny = numel(ygrid);

% Deduce grid over barycentric coordinates
r = [sort(repmat(xgrid, [1, ny])); repmat(ygrid, [1, nx])];
lambda = cartes2baryc(r, tricc);
lambda(:,any(lambda < 0, 1)) = NaN;
pHpgY = lambda(1,:);
pHdgY = lambda(2,:);
pHsgY = lambda(3,:);

% Get the positions within the matrix that belong to the triangle
sub = [sort(repmat(1:nx, [1, ny])); repmat(1:ny, [1, nx])];
ind = sub2ind([nx, ny], sub(1,:), sub(2,:));

% Create maps of predictors
nMaps = 4;
Maps = repmat({NaN(nx, ny)}, 1, nMaps);
Maps{1}(ind) = ceil(pHpgY);              % b0: offset
Maps{2}(ind) = pHpgY ./ (pHpgY + pHdgY); % b1: proba./deter. ratio
Maps{3}(ind) = pHsgY;                    % b2: fully-stochastic hypothesis
Maps{4} = Maps{1} .* Maps{2};            % b3: interaction

% Define the colormaps to use for each map
prec = 100;
cmap = cell(1,nMaps);
cmap{1} = g;
cmap{2} = [flipud(cbrewer2('Reds', prec)); cbrewer2('Blues', prec)];
cmap{3} = cbrewer2('Greens', 2*prec);
cmap{4} = cbrewer2('Purples', 2*prec);

% For each predictor
figure('Name', 'Maps of predictors', 'Position', [1 974 341 130]);
for iMap = 1:nMaps
    sp = subplot(1,nMaps,iMap);
    
    % Display maps of predictors
    imagesc(xgrid, ygrid, Maps{iMap}', 'AlphaData', ~isnan(Maps{iMap})'); hold('on');
    
    % Display information about the triangle
    tr = fill(tricc(:,1), tricc(:,2), 'k', 'FaceColor', 'None', 'LineWidth', 2);
    
	% Use corresponding colormap
    colormap(sp, cmap{iMap});
    caxis([0,1]);
    
    % Customize the axes
    axis([0, 1, 0, sqrt(3)/2]);
    axis('xy'); axis('off'); axis('equal');
end

% Save the figure
save2pdf(fullfile(ftapath, 'figs', 'F_HW_PredMaps.pdf'));

%% REGRESS P/D RATIO FROM SUBJECTS AGAINST P/D RATIO FROM THE IDEAL OBSERVER
%  =========================================================================

% Define options
% ~~~~~~~~~~~~~~

% Whether to restrict to fully-stochastic parts of sequences
onlystoch = true;

% Define the number of bins to use when plotting the results
nBin = 10;

% Perform subject-specific regressions
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Select data (all or only fully-stochastic parts)
if      onlystoch, obsidx = cellfun(@(x) x.Gen == 1, G, 'UniformOutput', 0);
elseif ~onlystoch, obsidx = cellfun(@(x) true(1,N),  G, 'UniformOutput', 0);
end

% Get subjects' P/D ratio
pHpgYpHdgY_sub = cellfun(@(x,i) x.BarycCoord(i,2) ./ sum(x.BarycCoord(i,1:2), 2), ...
    G, obsidx, 'UniformOutput', 0);

% Get ideal observer's P/D ratio
pHpgYpHdgY_IO = cellfun(@(x,i) x.BarycCoord(i,2) ./ sum(x.BarycCoord(i,1:2), 2), ...
    IO, obsidx, 'UniformOutput', 0);

% Get ideal observer's belief in the fully-stochastic hypothesis
pHsgY_IO = cellfun(@(x,i) x.BarycCoord(i,3), IO, obsidx, 'UniformOutput', 0);

% Define bins' edges based on percentiles computed from the ideal
% observer's estimated probabilities
dist = cell2mat(pHpgYpHdgY_IO(:));
dist(isnan(dist)) = 1/2;
dist = dist - mean(dist, 1, 'OmitNaN');
edges = prctile(dist, linspace(0, 100, nBin+1));
binlist = num2cell(1:nBin);

% Prepare output variable
binsubn = NaN(nBin,nSub  );
subavg  = NaN(nBin,nSub,3);
ioavg   = NaN(nBin,nSub,3);
regcoef = NaN(3,   nSub,4);

% For each subject
for iSub = 1:nSub
    
    % Create the design matrix based on ideal observer's beliefs
    ioratioPD = cell2mat(pHpgYpHdgY_IO(:,iSub));
    ioratioPD(isnan(ioratioPD)) = 1/2;
    ioratioPD = ioratioPD - mean(ioratioPD, 1, 'OmitNaN');
    iobelinS  = cell2mat(pHsgY_IO(:,iSub));
    iointerac = ioratioPD .* iobelinS;
    offset    = ones(numel(ioratioPD), 1);    
    desmat    = [offset, ioratioPD, iobelinS, iointerac];
    
    % Center predictors
    desmat = desmat - [0, mean(desmat(:,2:end), 1, 'OmitNaN')];
    
    % Select variable to explain
    subratioPD = cell2mat(pHpgYpHdgY_sub(:,iSub));
    subratioPD(isnan(subratioPD)) = 1/2;
    % N.B. when the division is not possible (i.e. P+D = 0) set the ratio
    % to 1/2
    
    % For each regression model
    % (each new model entails one more predictor than the previous one)
    for iMod = 1:3
        
        % Regress subject's P/D ratio against several IO predictors
        selecpred = 1:(1+iMod);
        beta = regress(subratioPD, desmat(:,selecpred));
        regcoef(iMod,iSub,selecpred) = beta;
        
        % Remove variance from all the predictiors except the P/D ratio
        regofint = 2;
        pidx = setdiff(selecpred, regofint);
        pred = desmat(:,pidx) * beta(pidx); % predictions based on beta values
        subratioPD2 = subratioPD - pred; % pseudo-residuals
        
        % Find which observations 
        [binsubn(:,iSub),~,bins] = histcounts(desmat(:,regofint), edges);
        
        % Average IO P/D ratio and subject's residual P/D ratio
        subavg(:,iSub,iMod) = cellfun(@(i) mean(subratioPD2(bins == i)), binlist);
        ioavg(:,iSub,iMod)  = cellfun(@(i) mean(ioratioPD(  bins == i)), binlist);
    end
end

% Display distribution of regression coefficient
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Average beta coefficients over subjects
m = mean(regcoef(:,:,2), 2)';
s = sem(regcoef(:,:,2),  2)';

% Display averaged regression coefficients
figure('Position', [1 700 140 200]);
bar(1:3, m, 'FaceColor', g, 'EdgeColor', 'k'); hold('on');
plot(repmat(1:3, 2, 1), m+[-s;s], 'k-');

% Customize the axes
if onlystoch, axis([0,4,0,0.4]); end
set(gca, 'Box', 'Off', 'XTick', 1:3);
xlabel('Regression models');
ylabel('Regression coefficient for P/D beliefs');

% Save the figure
save2pdf(fullfile(ftapath, 'figs', 'F_HW_BetaCoef.pdf'));

% Display correlation between subjects' residual ratio and ideal observer's
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Regression model to display in the figure
iMod = 3; % #3 is the full model

% Test the distribution of beta values against 0
[~,pval,tci,stats] = ttest(regcoef(:,:,2)');
Emergence_PrintTstats(pval, tci, stats);

% Average over subjects
xm = mean( ioavg(:,:,iMod), 2)';
xs = sem(  ioavg(:,:,iMod), 2)';
ym = mean(subavg(:,:,iMod), 2)';
ys = sem( subavg(:,:,iMod), 2)';

% Regress average P/D ratio from subjects and the ideal observer
B = Emergence_Regress(ym, xm, 'TLS', {'beta0', 'beta1'});
confint = Emergence_Regress(ym, xm, 'TLS', 'confint');
xval = Emergence_Regress(ym, xm, 'TLS', 'confintx');

% Create colormap used to color bins according to IO belief
prec = 100;
gege = [flipud(cbrewer2('Blues', prec)); cbrewer2('Reds', prec)];
grid = linspace(-1/2, 1/2, prec*2);
[~,colidx] = min(abs(xm'-grid), [], 2);

% Create a new window
figure('Position', [142 700 200 200]);

% Display origin (because variables are centered)
plot(zeros(1,2), [-1,1], '-', 'Color', g); hold('on');
plot([-1,1], zeros(1,2), '-', 'Color', g);

% Display identity line
plot(ylim, ylim, '-', 'Color', g);

% Display the regression line
fill([xval, fliplr(xval)], [confint(1,:), fliplr(confint(2,:))], 'k', ...
       'EdgeColor', 'none', 'FaceColor', g, 'FaceAlpha', 1/2);
plot(xval, xval.*B(2) + B(1), 'k-', 'LineWidth', 3);

% Get average size of the bin
avgbinsize = mean(binsubn, 2);
dotsize = 1 + flipud(log(avgbinsize+1).*10);

% Show correlation between subjects' P/D residual ratio and IO P/D ratio
plot(xm+[-xs;xs], repmat(ym,2,1), 'k-'); % horizontal error bars
plot(repmat(xm,2,1), ym+[-ys;ys], 'k-'); % vertical   error bars
plot(xm, ym, 'k-', 'LineWidth', 1);
scatter(xm, ym, dotsize, gege(colidx,:), 'filled', 'MarkerEdgeColor', 'k');

% Customize the axes
if onlystoch, axis([-0.4,0.4,-0.1,0.1]); end
set(gca, 'Box', 'Off');

% Add some text labels
xlabel('Centered beliefs from IO');
ylabel('Residual centered beliefs from subjects');

% Save the figure
save2pdf(fullfile(ftapath, 'figs', 'F_HW_CorrRatio.pdf'));
