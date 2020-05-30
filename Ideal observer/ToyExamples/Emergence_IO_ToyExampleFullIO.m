% Toy-example script for the full ideal observer (IO) of the task. It generates 
% an example sequence that can be either (1) fully-stochastic, (2) with a
% first fully-stochastic part followed by a second part charcaterized by
% a probability bias, or (3) with a first fully-stochastic part followed by
% a second part characterized by the repetition of a given pattern. The
% full IO is presented with that sequence and assigns credence to each of
% these 3 possible generative processes (marginalizing over all possible
% change point's positions in the two latter cases). Those posterior
% beliefs are updated as new observations are received. In
% addition, the full IO also returns posterior beliefs over change point's
% positions, models' parameters, and identity of the future observation.
% 
% Copyright (c) 2020 Maxime Maheu

%% INITIALIZATION
%  ==============

% Clear the place
clear;
close('all');

% Add ideal observer functions to the MATLAB path
scriptpath = mfilename('fullpath');
ind = strfind(scriptpath,'Emergence');
folderpath = scriptpath(1:ind(end-1)+8);
addpath(genpath(folderpath));

%% GENERATE A SEQUENCE
%  ===================

% If there is no sequence already loaded in the workspace
if ~exist('sequence', 'var')

% Define the generative process of the second part of the sequence
%   >> rule = NaN; % fully-stochastic sequence
%   >> rule = [1/4 3/4]; % stochastic-to-probabilistic sequence
%   >> rule = [1 2 2 1 1 2]; % stochastic-to-deterministic sequence
rule = [3/4 3/4];

% Define the length of the sequence
N = 200;

% Define the position of the change point
Jmu = 100;
Jsigma = 15;
J = 0;
while J-Jmu < -3*Jsigma || J-Jmu > 3*Jsigma
    J = round(normrnd(Jmu, Jsigma));
end

% Generate a random sequence
if isnan(rule), sequence = GenRandSeq(N, 1/2); J = N;
elseif any(rule < 1), sequence = [GenRandSeq(J, 1/2), GenRandSeq(N-J, rule)];
else, sequence = [GenRandSeq(J, 1/2), repmat(rule, 1, N-J)];
end
end
sequence = sequence(1:N);

%% RUN THE BAYESIAN IDEAL OBSERVER
%  ===============================

% Define options for the ideal observer
pEd   = 0;                  % probability of making a memory error at each observation (deterministic rule)
pEp   = 0;                  % probability of making a memory error at each observation (probability bias)
treed = 20;                 % depth of the rules' tree to explore
stat  = 'Transition';       % statistic to be learned by the probabilistic model
p_pR  = 'Size-principle';   % the prior probability of each rule depends on its length
p_pT  = 'Bayes-Laplace';    % the prior over statistics to be learned
p_pJ  = 'Uniform';          % prior over change point's position
comp  = 'all';              % update beliefs after each observation
scale = 'log';              % scale of the model evidence
pgrid = 2e-2;               % precision of the posterior over theta
verb  = 1;                  % output some messages in the command window

% Run the observer with these options
io = Emergence_IO_FullIO(sequence, ... % binary input sequence
    pEd, pEp, treed, stat, p_pR, p_pT, p_pJ, comp, scale, pgrid, verb); % options

%% PREPARE FIGURES
%  ===============

% Define the screen on which to display the figures
res = get(0, 'MonitorPositions'); % get pixel resolution of each screen
screen = size(res,1); % display the figures on one of the available screens
sc0 = [0 + screen - 1, res(screen,2) / res(1,4)]; % offset
magfac = [res(screen,3) / res(1,3), res(screen,4) / res(1,4)]; % ratio size compared to screen 1

% Define relative size of figures
mrg = [1/20, 1/10] .* magfac; % left/right and upper/lower margins (in percent)
fw  = [0.55, 0.40] .* magfac(1); % [upper, lower] leftward* figures' width
fh  = [0.50, 0.50] .* magfac(2); % [left, right] lower* figures' height
% *: The other dimension has no degree of freedom : 1-2*margin*x

% Compute positions of the figures
fpos = [sc0(1)+mrg(1)       sc0(2)+mrg(2)+fh(1) fw(1)                    magfac(2)-2*mrg(2)-fh(1); ... % upper left
        sc0(1)+mrg(1)       sc0(2)+mrg(2)       fw(2)                    fh(1)                   ; ... % lower left
        sc0(1)+mrg(1)+fw(1) sc0(2)+mrg(2)+fh(2) magfac(1)-2*mrg(1)-fw(1) magfac(2)-2*mrg(2)-fh(2); ... % upper right
        sc0(1)+mrg(1)+fw(2) sc0(2)+mrg(2)       magfac(1)-2*mrg(1)-fw(2) fh(2)                   ];    % lower right

% Define models' name and variables to be displayed
lmlab = {'Fully-stochastic', 'Probabilistic', 'Deterministic'};
smlab = {'s', 'p', 'd'};
fname = {'Posterior model probability', ...
         'Posterior distribution over change point''s position', ...
         'Posterior distribution over model''s parameters', ...
         'Expectation-related quantities'};

% Calling function for figures
ffun = @(i) figure('Color', ones(1,3), 'Units', 'Normalized', ...
    'Name', fname{i}, 'MenuBar', 'None', 'ToolBar', 'None', ...
	'Position', fpos(i,:), 'DefaultLegendAutoUpdate', 'Off');

% Index font-size on screen resolution
if    ~ismac, scfac = 1/180; % !!! on Retina displays the resolution looks
elseif ismac, scfac = 1/120; % actually smaller than what it truly is
end
fs = (res(screen,4) - res(screen,2) + 1) * scfac * max(magfac);
fs = max(8, fs); % make sure the font size is not smaller than 8 points

% Define tickness of lines
wid1 = 1/2; % for boxes, axes, help lines
wid2 = 1;   % for line plots and relevant information

% Define various default properties
Emergence_DefaultFigureProperties;
defaxprop = @(ax) set(ax, 'FontSize', fs, 'LineWidth', wid1, ...
    'Color', 'None', 'XColor', 'k', 'YColor', 'k', ...
    'Box', 'On', 'Layer', 'Top', 'TickLabelInterpreter', 'LaTeX', ...
    'XLim', [1,N], 'XTick', [1,50:50:N], 'YDir', 'Normal');
txtopt = {'FontSize', fs, 'Interpreter', 'LaTeX', 'Rotation', 0, ...
    'HorizontalAlignment', 'Right', 'VerticalAlignment', 'Middle'};
cbropt = {'FontSize', fs, 'LineWidth', wid1, 'TickLabelInterpreter', 'LaTeX'};
cbrlabopt = {'FontSize', fs, 'Interpreter', 'LaTeX'};

% Define colours to use
tricol = [049, 130, 189; 222, 045, 038; 049, 163, 084] ./ 255;
regcol = [117, 107, 177] ./ 255;
cmap1  = hsv(N); % for the observations in the sequence
cmap2  = Emergence_Colormap('Parula'); % for the heatmaps

%% POSTERIOR DISTRIBUTIONS OVER MODELS
%  ===================================

% Prepare a new figure
ffun(1);

% Display the sequence and how it is described
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Display the observations
subplot(9,3,1:2); hold('on');
imagesc(1:N, 1.5, 1:N); alpha(1/2); colormap(cmap1);
plot(1:N, abs(sequence-3), 'k.-', 'LineWidth', wid1);
text(-5, 1.5, 'Sequence $y$', txtopt{:});

% Display the best model
for iMod = 1:3
    stricol = tricol(circshift(1:3,1),:);
    plot(find(io.Mhat == iMod), zeros(1, sum(io.Mhat == iMod)), 's', ...
        'MarkerEdgeColor', 'None', 'MarkerFaceColor', stricol(iMod,:));
    tp = io.(['dMs',smlab{iMod}]);
    plot(find(tp), zeros(1,sum(tp))-1, 's', ...
        'MarkerEdgeColor', 'None', 'MarkerFaceColor', stricol(iMod,:));
end
text(-5, 0,  'Best model $\mathcal{M}$',         txtopt{:});
text(-5, -1, 'Signif. best model $\mathcal{M}$', txtopt{:});

% Customize the axes
defaxprop(gca);
set(gca, 'XTickLabel', {}, 'YLim', [-1.5,2.5], 'Visible', 'Off');

% Display the Bayes Factors
% ~~~~~~~~~~~~~~~~~~~~~~~~~

% Create a subplot
subplot(9,3,4+[0,1,3,4]); hold('on');
limy = ceil(max(reshape(cumsum(abs([io.BFsdss; io.BFspss])), [1,N*2])));

% Set threshold (0 for logarithmic scale, 1 for linear scale)
thr = 0 + double(strcmpi(io.in.scaleme, 'lin'));

% First part of the barplot (same signs)
idx = io.BFspss < thr & io.BFsdss < thr | io.BFspss > thr & io.BFsdss > thr;
t = bar(find(idx), [io.BFspss(idx); io.BFsdss(idx)]', 1, 'Stacked');
for iMod = 1:2
    set(t(iMod), 'FaceColor', tricol(iMod,:), 'EdgeColor', 'k', 'LineWidth', wid1);
end

% Second part of the barplot (different signs)
idx = io.BFspss < thr & io.BFsdss > thr | io.BFspss > thr & io.BFsdss < thr;
for iMod = 1:2 
    bar(find(idx), io.(['BFs',smlab{iMod+1},'ss'])(idx), 1, 'EdgeColor', 'k', ...
        'FaceColor', tricol(iMod,:), 'LineWidth', wid1);
end

% Display zero level
plot([0,N+1], zeros(1,2), 'k-', 'LineWidth', wid2);

% Customize the axes
defaxprop(gca);
set(gca, 'XTickLabel', {});

% Add some text labels
if     strcmpi(scale, 'lin'), txt = '';
elseif strcmpi(scale, 'log'), txt = 'log ';
end
ylabel({[txt, ' Bayes factor'], [txt, ' $\frac{p(y_{1:K}|\mathcal{M}_{', ...
    '\mathrm{S} \rightarrow i})}{p(y_{1:K}|\mathcal{M}_{\mathrm{S}', ...
    '\rightarrow \mathrm{S}})}$']}, txtopt{:});

% Add a log/linear axis
if     strcmpi(scale, 'lin'), fun = @log;
elseif strcmpi(scale, 'log'), fun = @exp;
end
ax = {get(gca, 'YLim'), get(gca, 'YTick')};
yyaxis('right');
set(gca, 'YLim', ax{1}, 'YTick', ax{2}, 'YTickLabel', ...
    cellfun(@(x) sprintf('%1.2g', x), num2cell(fun(ax{2})), 'uni', 0));
defaxprop(gca);

% Display models' posterior probabilities the triangular arena
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Display the triangular arena
subplot(9,3,3:3:(3*9));
Emergence_PlotTrajOnTri([], J, tricol, 10, [], false, fs*1.5); alpha(3/4);
Emergence_PlotGridOnTri(3);

% Convert barycentric coordinates to cartesian ones
pMgY = cell2mat(arrayfun(@(x) io.(['pMs',smlab{x},'gY'])', [2:3,1], 'uni', 0));
tricc = [0, sqrt(3)/2; 1, sqrt(3)/2; 1/2, 0];
cartcoor = pMgY*tricc;

% Display the posterior beliefs as a trajectory
plot(cartcoor(:,1), cartcoor(:,2), 'k-', 'LineWidth', wid1);
scatter(cartcoor(:,1), cartcoor(:,2), [], cmap1, '.', 'SizeData', 100);

% Customize the axes
a = axis; defaxprop(gca); axis(a);

% Add some text labels
title({'Trajectory based on $p(\mathcal{M}_{\mathrm{S} \rightarrow i}|y_{1:K})$', ...
    ''}, 'Interpreter', 'LaTeX');

% Display the cumulative models' posterior probabilities
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Plot the cumulative models' posterior probabilities
subplot(9,3,10+[0,1,3,4]);
Emergence_PlotBarycTraj(pMgY, tricol);

% Customize the axes
defaxprop(gca); pos = get(gca, 'Position');
set(gca, 'XTickLabel', {}, 'YLim', [0,1], 'YTick', 0:1/5:1);

% Add some text labels
ylabel({'Model', 'posterior', 'probability', '$p(\mathcal{M}|y_{1:K})$'}, txtopt{:});

% Plot the amount of beliefs update
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Display the Jensen-Shannon distance
subplot(9,3,16+[0,1,3,4]); hold('on');
%plot(1:N, sqrt(io.JSpMgY), 'k.-', 'MarkerSize', 10, 'LineWidth', wid2);
plot(1:N, sqrt(io.eqAlpha), 'k.-', 'MarkerSize', 10, 'LineWidth', wid2);

% Customize the axes
defaxprop(gca);
set(gca, 'XTickLabel', {});

% Add some text labels
ylabel({'Update of', 'the posterior', 'distribution', ...
    '$\sqrt{D_\mathrm{JS}(p(\mathcal{M}|y_{1:K})}$', ...
    '$\overline{||p(\mathcal{M}|y_{1:K-1}))}$'}, txtopt{:});

% Plot the entropy of the posterior distribution over models
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% It quantifies the distance to the center of the triangle
subplot(9,3,22+[0,1,3,4]); hold('on'); limy = [0, log2(3)];
plot(1:N, io.HpMgY, 'k.-', 'MarkerSize', 10, 'LineWidth', wid2);

% Customize the axes
defaxprop(gca);
set(gca, 'YLim', limy);

% Add some text labels
xlabel('Observation ($K$)', 'Interpreter', 'LaTeX');
ylabel({'Entropy of', 'the posterior', 'distribution', ...
    '$H(p(\mathcal{M}|y_{1:K}))$'}, txtopt{:});

% Display true position of the change point in all the subplots
subplot(9,3,1:2);
plot(repmat(J,1,2), ylim, 'k-', 'LineWidth', wid1);
for sp = 4:6:22
    subplot(9,3,sp+[0,1,3,4]);
    plot(repmat(J,1,2), ylim, 'k-', 'LineWidth', wid1);
end

% Display the posterior model's probabilities at the end of the sequence
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Create an inset plot
insetw = 0.12; % width of the plot (in percent)
axes('Position', [pos(1)+pos(3)-insetw/6, pos(2)+pos(4)/2-insetw/2, insetw, insetw]);

% Display the pie chart
h = pie(pMgY(end,:));
x = find(pMgY(end,:) ~= 0);

% Customize the pie chart
for iMod = 1:numel(x)
    set(h(iMod+(iMod-1)), 'FaceColor', tricol(x(iMod),:), 'LineWidth', wid2);
    str = get(h(iMod*2), 'String');
    str = ['$' str(1:end-1) '\' str(end) '$'];
    set(h(iMod*2), 'String', str, cbrlabopt{:});
end

%% POSTERIOR DISTRIBUTIONS OVER MODELS' PARAMETERS
%  ===============================================

% Prepare a new figure
ffun(3);

% For each dimension of the posterior distibution over the probabilistic
% model's parameters (i.e. each parameter, i.e. 1 when we estimate item or
% alternation frequency, and 2^n for n-order Markov chains)
nTdim = size(io.pTgYMsp, 3);
for iTdim = 1:nTdim
    subplot(nTdim+1, 1, iTdim);
    
    % Display posterior distribution over model's parameters for the
    % probabilistic observer (i.e. transitions) marginalized over possible
    % (already observed) change point's positions
    imagesc(1:N, 0:io.in.postgrid:1, io.pTgYMsp(:,:,iTdim));
    caxis([0, max(io.pTgYMsp(:))]);
    
    % Customize the axes
    set(gca, 'YTick', 0:1/4:1);
    
    % Add some text labels
    ylabel({'Inferred', 'value of', ['$\theta_{', ...
        num2str(iTdim), '}$']}, txtopt{:});
end

% Display posterior distribution over model's parameters for the
% deterministic observer (i.e. patterns) marginalized over possible
% (already observed) change point's positions
subplot(nTdim+1, 1, nTdim+1);
imagesc(1:N, 1:io.in.nu, io.pTgYMsd);
caxis([0, max(io.pTgYMsd(:))]);
    
% Add some text labels
xlabel('Observation ($K$)', 'Interpreter', 'LaTeX');
ylabel({'Pattern', 'length', '$R_{i}$'}, txtopt{:});

% Customize each subplot
for iTdim = 1:nTdim+1
    subplot(nTdim+1, 1, iTdim); hold('on');
    
    % Display the true position of the change point
    plot(repmat(J,1,2), ylim, 'k-', 'LineWidth', wid1);
    
    % Add a colorbar
    colormap(cmap2);
    cbr = colorbar('Location', 'EastOutside', cbropt{:});
    if iTdim <= nTdim, toto = ['\theta_{', num2str(iTdim), '}'];
    else, toto = 'R';
    end
    set(cbr.Label, 'String', {'Marginal', 'posterior', 'distribution', ...
        ['$p(', toto, '|y_{1:K})$']}, 'Rotation', 0, 'HorizontalAlignment', ...
        'Left', 'VerticalAlignment', 'Middle', cbrlabopt{:});
    
    % Customize the axes
    defaxprop(gca);
end

%% POSTERIOR DISTRIBUTIONS OVER CHANGE POINT'S POSITIONS
%  =====================================================

% Choose which models to display
mod = {'gYMsp','gYMsd','gY'};
modlab = {{'Stochastic-to-probabilistic', ['model $\mathcal{M}_', ...
            '{\mathrm{S}\rightarrow \mathrm{P}}$']}, ...
          {'Stochastic-to-deterministic', ['model $\mathcal{M}_', ...
            '{\mathrm{S} \rightarrow \mathrm{D}}$']}, ...
          {'Bayesian model', 'averaging'}};
modcol = [tricol([1,2],:); regcol];

% Variables characterizing the posterior distribution over change point's position
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Choose which variables to display
var = {'', 'H', 'JS'};
varlab = {{'Posterior', 'distribution', 'over already', 'observed change', ...
            'point''s positions', '$p(j_{k}<K|y_{1:K},\mathcal{M})$'}
          {'Entropy of', 'the posterior', 'distribution', ...
            '$H(p(j_{k}|y_{1:K},\mathcal{M})$'}
          {'Update of', 'the posterior', 'distribution', ...
            '$\sqrt{D_\mathrm{JS}(p(j_{k}|y_{1:K})}$', ...
            '$\overline{||p(j_{k}|y_{1:K-1},\mathcal{M}))}$'}};

% Prepare a new figure
ffun(2);

% For each pair of model and variable
for iMod = 1:3
    for ivar = 1:3
        subplot(6,3,iMod+3*(ivar-1)); hold('on');       
        
        % Display the dynamic of that variable
        tp = io.([var{ivar}, 'pJk', mod{iMod}]);
        if ivar == 1
            tp = mat2cell(tp, N, ones(N,1));
            tp = cellfun(@(x,y) sum(x(1:y)), tp, num2cell(1:N));
        elseif ivar == 3
            tp = sqrt(tp); 
        end
        plot(1:N, tp, '-', 'Color', modcol(iMod,:), 'LineWidth', wid2);
        
        % Display other help lines
        if ivar == 1
            plot([1,N], ones(1,2)./2, 'k:', 'LineWidth', wid1); 
            plot([1,N], [0,1], 'k:', 'LineWidth', wid1);
        end
        
        % Customize the axes
        defaxprop(gca);
        if     any(ivar == [1 3]), set(gca, 'YLim', [0,1]);
        elseif ivar == 2, set(gca, 'YLim', [0,log2(N)]);
        end
        
        % Display the real position of the change point
        plot(repmat(J,1,2), ylim, 'k-', 'LineWidth', wid1);
        
        % Add some text labels
        if iMod == 1, ylabel(varlab{ivar}, txtopt{:}); end
        if ivar == 1, title(modlab{iMod}, 'Interpreter', 'LaTeX'); end
        if ivar == 4, xlabel('Observation ($K$)', 'Interpreter', 'LaTeX'); end
    end
    
    % Posterior distribution over change point's position
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    % Plot the posterior distribution of change point's position
    subplot(6,3,iMod+(9:3:3*6-iMod)); hold('on');
    tp = io.(['pJk',mod{iMod}]);
    imagesc(1:N, 1:N, tp); 
    caxis([0, max(tp(:))]);
    
    % Display limit between past and future events
    plot([1,N], [1,N], 'k-', 'LineWidth', wid1);
    
    % Display true change point's position
    plot(repmat(J,1,2), [1,N], 'k-',  'LineWidth', wid1);
    plot([1,J], repmat(J,1,2), 'k--', 'LineWidth', wid1);
    
    % Display the colorbar
    colormap(cmap2);
    cbr = colorbar('Location', 'SouthOutside', cbropt{:});
    if iMod <= 2, mathlab = ['$p(j_{k}|y_{1:K},\mathcal{M}_', ...
            '{\mathrm{S} \rightarrow \mathrm{', lmlab{iMod}(1), '}})$'];
    elseif iMod == 3, mathlab = '$p(j_{k}|y_{1:K})$';
    end
    set(cbr.Label, cbrlabopt{:}, 'String', ...
        {'Posterior distribution over', 'change point''s position', mathlab});
    
    % Customize the axes
    defaxprop(gca);
    set(gca, 'XLim', [1/2,N+1/2], 'YLim', [1/2,N+1/2], 'YTick', [1,50:50:N]);
    
    % Add some text labels
    xlabel('Observation ($K$)', 'Interpreter', 'LaTeX');
    if iMod == 1, ylabel({'Change point''s', 'position ($j_{k}$)'}, txtopt{:}); end 
end

%% MODELS' PREDICTIONS
%  ===================

% Choose which variables to display
var = {'', 'I', 'H', 'JS'};
varlab = {{'Expectation', 'regarding the', 'current observation', ...
            '$p(y_{K}=\mathrm{A}|y_{1:K-1},\mathcal{M})$'}, ...
          {'Surprise', 'evoked by the', 'current observation', ...
            '$-\log_{2}(p(y_{K}|y_{1:K-1},\mathcal{M}))$'}, ...
          {'Entropy of', 'the posterior', 'distribution', ...
            '$H(p(y_{K}|y_{1:K-1},\mathcal{M}))$'}, ...
          {'Update of', 'the posterior', 'distribution', ...
            '$\sqrt{D_\mathrm{JS}(p(\mathrm{A}|y_{1:K},\mathcal{M})}$', ...
            '$\overline{||p(\mathrm{A}|y_{1:K-1},\mathcal{M}))}$'}};
slev = [1/2, 1, 1, 0]; % level of those variables for the fully-stochastic model

% Prepare a new figure
ffun(4);

% For the 2 non-fully-stochastic models and the BMA
for iMod = 1:3
    for ivar = 1:4
        subplot(4,3,iMod+3*(ivar-1)); hold('on');
        
        % Display the same variable in the case of the fully-stochastic model
        if iMod < 3
            plot([1,N], repmat(slev(ivar),1,2), '-', ...
                'Color', tricol(end,:), 'LineWidth', wid2);
        end
        
        % Display expectation related quantities
        tp = io.([var{ivar}, 'pA', mod{iMod}]);
        if ivar == 4, tp = sqrt(tp); end
        plot(1:N, tp, '-', 'Color', modcol(iMod,:), 'LineWidth', wid2);
        
        % Customize the axes
        defaxprop(gca);
        if ivar == 2, set(gca, 'YLim', [0,4]);
        else, set(gca, 'YLim', [0,1]);
        end
        
        % Display the real position of the change point
        plot(repmat(J,1,2), ylim, 'k-', 'LineWidth', wid1);
        
        % Add some text labels
        if iMod == 1, ylabel(varlab{ivar}, txtopt{:}); end
        if ivar == 1, title(modlab{iMod}, 'Interpreter', 'LaTeX'); end
        if ivar == 4, xlabel('Observation ($K$)', 'Interpreter', 'LaTeX'); end
    end
end

%% WRAP THINGS UP
%  ==============

% Define the folder where to save the figures
try
figpath = fullfile(scriptpath(1:max(strfind(scriptpath,'/'))), 'figs/');
if exist(figpath, 'dir') ~= 7, mkdir(figpath); end

% Save the figures as separate files
figimg = arrayfun(@(x) frame2im(getframe(figure(x))), 1:4, 'uni', 0);
cellfun(@(x,y) imwrite(x, [figpath, sprintf('Emergence_IO_ToyExampleFullIO_fig%i.jpeg', ...
    y)]), figimg, num2cell(1:4), 'uni', 0);

% Save all the figures as a single files
figimg = cat(1, cat(2, figimg{1}, figimg{2}), cat(2, figimg{3}, figimg{4}));
imwrite(figimg, [figpath, 'Emergence_IO_ToyExampleFullIO.jpeg']);
catch
end

% Reoder figures on the screen
arrayfun(@figure, 4:-1:1);
