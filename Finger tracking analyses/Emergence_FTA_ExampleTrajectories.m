% This script displays individual example and representative trajectories
% from the subjects and the ideal observer.
% 
% Copyright (c) 2018 Maxime Maheu

% N.B. To be run without subject exclusion. To do so, turn the boolean
% variable "rmbadsub" in "Emergence_FTA_LoadData" to false.

%% SUBJECT AND IDEAL OBSERVER BELIEFS IN 3 EXAMPLE TRAJECTORIES
%  ============================================================

% Compare with the ideal observer
subs  = [22 25 28];
conds = [10 16 22];
nseq  = numel(subs);

% Prepare a new window
figure('Position', [1 805 800 300]);

% For each sequence
for iSeq = 1:nseq
    
    % For both the subject and the ideal observer
    for iObs = 1:2        
        if     iObs == 1, X = G{conds(iSeq),subs(iSeq)};
        elseif iObs == 2, X = IO{conds(iSeq),subs(iSeq)};
        end
        
        % Plot the trajectory within the triangle
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ipos = iSeq + 3*(iObs-1) + 5*(iSeq-1);
        subplot(nseq, 2*(1+2), ipos);
        Emergence_PlotTrajOnTri(X.BarycCoord, X.Jump+1/2, tricol, 6);
        Emergence_PlotGridOnTri(3);
        
        % Plot the Barycentric coordinates
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        subplot(nseq, 2*(1+2), (1:2) + ipos);
        Emergence_PlotBarycTraj(X.BarycCoord, tricol);
        
        % Append help lines
        plot([1,N], repmat(1/3,1,2), '-', 'Color', g);
        plot([1,N], repmat(2/3,1,2), '-', 'Color', g);
        plot(repmat(X.Jump,1,2), [0,1], 'k-', 'LineWidth', 1);
        
        % Customize the axes
        axis([1,N,0,1]);
        set(gca, 'YTick', 0:1/3:1, 'YTickLabel', {'0','1/3','2/3','1'});
        
        % Add some text labels
        title(sprintf('Subject %1.0f, Sequence %1.0f', subs(iSeq), conds(iSeq)));
    end
end

% Save the figure
save2pdf('figs/F_ET_SubVsIO.pdf');

%% DIVERSITY OF TRAJECTORIES
%  =========================

% Define which trajectories to display
subs  = [08 02 15 02 20 20 06 06];
conds = [18 02 15 22 22 18 12 18];
trajnames = {'Jolting', 'Sweeping' ...
             'Hesitation', 'Changes of mind', ...
             'Ballistic', 'Gradual', ...
             'Low confidence', 'High confidence'};
trajgp = {{'Discrete vs', 'continuous'}, ...
          {'Beliefs''', 'updating'}, ...
          {'Confidence', 'build-up'}, ...
          {'Levels of', 'confidence'}};

% Prepare a new window
figure('Position', [1 284 800 447]);

% For each sequence
for iSeq = 1:numel(trajnames)
    
    % Plot the trajectory within the triangle
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    % Plot the trajectory and the triangle
    ipos = (iSeq+1:iSeq+1)+2*(iSeq-1)-1;
    subplot(4, 2*(1+2), ipos);
    X = D{conds(iSeq),subs(iSeq)};
    Emergence_PlotTrajOnTri(X.BarycCoord, X.Jump+1/2, tricol);
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
    text(7, 0.9, sprintf('Subject %i', subs(iSeq)), ...
        'HorizontalAlignment', 'Left', 'FontSize', 8);
    text(7, 0.75, sprintf('%s', c{conds(iSeq)}), ...
            'HorizontalAlignment', 'Left', 'FontSize', 8);
end

% Save the figure
save2pdf('figs/F_ET_TrajDiv.pdf');
