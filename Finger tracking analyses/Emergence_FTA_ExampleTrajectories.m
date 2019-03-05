% This script displays individual example and representative trajectories
% from the subjects and the ideal observer.
% 
% Copyright (c) 2018 Maxime Maheu

% N.B. To be run without subject exclusion. To do so, turn the boolean
% variable "rmbadsub" in "Emergence_FTA_LoadData" to false.

%% INITIALIZATION
%  ==============

% Make sure data from bad subjects is available
if ~exist('rmbadsub', 'var') || rmbadsub
    rmbadsub = false;
    Emergence_FTA_LoadData;
end

%% 3 EXAMPLE SEQUENCES
%  ===================

% Select example sequences to look at
sublist  = [28 27 28];
condlist = [04 16 22];
nExpSeq  = numel(sublist);

% Prepare the window
figure('Position', [1 805 800 300]);

% Display example sequences
% ~~~~~~~~~~~~~~~~~~~~~~~~~

% For each sequence
for iSeq = 1:nExpSeq
    cond = condlist(iSeq);
    sub = sublist(iSeq);
    subplot(3,1,iSeq); hold('on');
    
    % Display the generative process(es) of the sequence
    for iProc = unique(G{cond,sub}.Gen)
        if     iProc == 2 && strcmpi(G{cond,sub}.Cond(1), 'P')
            col = tricol(1,:);
            rule = sprintf('p(A|B) = %1.2f & p(B|A) = %1.2f', G{cond,sub}.Rule);
        elseif iProc == 2 && strcmpi(G{cond,sub}.Cond(1), 'D')
            col = tricol(2,:);
            rule = sprintf('[%s]^{n}', pat2str(G{cond,sub}.Rule));
        elseif iProc == 1
            col = tricol(3,:);
            rule = '';
        end
        obsidx = find(G{cond,sub}.Gen == iProc);
        fill(obsidx([1,end,end,1])+ones(1,4)/2.*[-1,1,1,-1], ...
            [0,0,1,1]+ones(1,4)/4.*[-1,-1,1,1], ...
            'k', 'FaceColor', col, 'EdgeColor', 'none');
    end
    text(mean([1, min([G{cond,sub}.Jump, N])]), -1, ...
        'Stochastic part', 'Color', tricol(3,:));
    text(mean([N, G{cond,sub}.Jump]), -1, ...
        {sprintf('%s part', G{cond,sub}.Cond), rule}, 'Color', col);
    
    % Display the position of the change point
    plot(repmat(G{cond,sub}.Jump,1,2), [-0.5,1.5], 'k-', 'LineWidth', 1);
    
    % Display the sequence
    plot(1:N, G{cond,sub}.Seq-1, 'w.-', 'LineWidth', 1/2);
    
    % Customize the axes
    axis([1,N,-1,2]); axis('off');
end

% Save the figure
save2pdf(fullfile(ftapath, 'figs', 'F_ET_ExpSeqs.pdf'));

%% SUBJECT AND IDEAL OBSERVER BELIEFS IN 3 EXAMPLE TRAJECTORIES
%  ============================================================

% Prepare a new window
figure('Position', [1 431 800 300]);

% For each sequence
for iSeq = 1:nExpSeq
    
    % For both the subject and the ideal observer
    for iObs = 1:2        
        if     iObs == 1, X = G{ condlist(iSeq),sublist(iSeq)};
        elseif iObs == 2, X = IO{condlist(iSeq),sublist(iSeq)};
        end
        
        % Plot the trajectory within the triangle
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ipos = iSeq + 3*(iObs-1) + 5*(iSeq-1);
        subplot(nExpSeq, 2*(1+2), ipos);
        Emergence_PlotTrajOnTri(X.BarycCoord, X.Jump+1/2, tricol, 6);
        Emergence_PlotGridOnTri(3);
        
        % Plot the Barycentric coordinates
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        subplot(nExpSeq, 2*(1+2), (1:2) + ipos);
        Emergence_PlotBarycTraj(X.BarycCoord, tricol);
        
        % Append help lines
        plot([1,N], repmat(1/3,1,2), '-', 'Color', g);
        plot([1,N], repmat(2/3,1,2), '-', 'Color', g);
        plot(repmat(X.Jump,1,2), [0,1], 'k-', 'LineWidth', 1);
        
        % Customize the axes
        axis([1,N,0,1]);
        set(gca, 'YTick', 0:1/3:1, 'YTickLabel', {'0','1/3','2/3','1'});
        
        % Add some text labels
        title(sprintf('Subject %1.0f, Sequence %1.0f', sublist(iSeq), condlist(iSeq)));
    end
end

% Save the figure
save2pdf(fullfile(ftapath, 'figs', 'F_ET_SubVsIO.pdf'));

%% DIVERSITY OF TRAJECTORIES
%  =========================

% Define which trajectories to display
sublist  = [08 02 15 02 20 20 06 06];
condlist = [18 02 15 22 22 18 12 18];
trajnames = {'Jolting', 'Sweeping' ...
             'Hesitation', 'Changes of mind', ...
             'Ballistic', 'Gradual', ...
             'Low confidence', 'High confidence'};
trajgp = {{'Discrete vs', 'continuous'}, ...
          {'Beliefs''', 'updating'}, ...
          {'Confidence', 'build-up'}, ...
          {'Levels of', 'confidence'}};
nExpSeq  = numel(sublist);

% Prepare a new window
figure('Position', [802 658 800 447]);

% For each sequence
for iSeq = 1:nExpSeq
    
    % Plot the trajectory within the triangle
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    % Plot the trajectory and the triangle
    ipos = (iSeq+1:iSeq+1)+2*(iSeq-1)-1;
    subplot(4, 2*(1+2), ipos);
    X = D{condlist(iSeq),sublist(iSeq)};
    Emergence_PlotTrajOnTri(X.BarycCoord, X.Jump+1/2, tricol);
    Emergence_PlotGridOnTri(3);
    if mod(iSeq,2) == 1
        text(0, 1/2, trajgp{(iSeq+1)/2}, 'FontWeight', 'Bold', 'FontSize', ...
            12, 'HorizontalAlignment', 'Right', 'Rotation', 0);
    end
    
    % Plot the Barycentric coordinates
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    % Display cumulative Barycentric coordinates
    subplot(4, 2*(1+2), (1:2) + ipos);
    Emergence_PlotBarycTraj(X.BarycCoord, tricol);
    
    % Append help lines
    plot([1,N], repmat(1/3,1,2), '-', 'Color', g);
    plot([1,N], repmat(2/3,1,2), '-', 'Color', g);
    plot(repmat(X.Jump,1,2), [0,1], 'k-', 'LineWidth', 1);
    
    % Customize the axes
    axis([1,N,0,1]);
    set(gca, 'YTick', 0:1/3:1, 'YTickLabel', {'0','1/3','2/3','1'});
    
    % Add some text labels
    title(trajnames{iSeq});
    text(7, 0.9, sprintf('Subject %i', sublist(iSeq)), ...
        'HorizontalAlignment', 'Left', 'FontSize', 8);
    text(7, 0.75, sprintf('%s', c{condlist(iSeq)}), ...
            'HorizontalAlignment', 'Left', 'FontSize', 8);
end

% Save the figure
save2pdf(fullfile(ftapath, 'figs', 'F_ET_TrajDiv.pdf'));
