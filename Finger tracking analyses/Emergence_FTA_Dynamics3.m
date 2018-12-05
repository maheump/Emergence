% This script performs a correlation betweeen detection speed of
% probabilistic regularities from human subjects and the ideal observer on
% one hand and detection lag of deterministic regularities form human
% subjects and the ideal observer on the other hand.
% 
% Copyright (c) 2018 Maxime Maheu

%% RUN THE PREVIOUS SCRIPT SEPARATELY FOR SUBJECTS AND THE IO
%  ==========================================================

% Prepare output variable
Var = cell(2);

% Get IO's trajctories around change and detection points as well as
% meaningful measures such as detection delays and speed
D = IO;
Emergence_FTA_Dynamics2;
Var{1,1} = speed{1};
Var{2,1} = lag{2};

% Get subjects' trajctories around change and detection points as well as
% meaningful measures such as detection delays and speed
D = G;
Emergence_FTA_Dynamics2;
Var{1,2} = speed{1};
Var{2,2} = lag{2};

%% CORRELATION BETWEEN SUBJECTS AND IDEAL OBSERVER
%  ===============================================

% Display correlation between human subjects and the ideal observer
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Prepare two new windows
figure('Position', [1304 780 400 325]);

% For each regular hypothesis
for iHyp = 1:2
    
    % Average over subjects
    mx = mean(Var{iHyp,1}, 2, 'OmitNaN');
    sx = sem( Var{iHyp,1}, 2);
    my = mean(Var{iHyp,2}, 2, 'OmitNaN');
    sy = sem( Var{iHyp,2}, 2);
    
    % Display the correlation between the ideal observer and the
    % subjects
    confint  = Emergence_Regress(my, mx, 'TLS', 'confint');
    confintx = Emergence_Regress(my, mx, 'TLS', 'confintx');
    subplot(1,2,iHyp); hold('on');
    fill([confintx, fliplr(confintx)], [confint(1,:), fliplr(confint(2,:))], ...
        'k', 'EdgeColor', 'none', 'FaceColor', tricol(iHyp,:), 'FaceAlpha', 0.15);
    b = Emergence_Regress(my, mx, 'TLS', {'beta0', 'beta1'});
    plot(confintx([1,end]), confintx([1,end]).*b(2)+b(1), '-', ...
        'Color', tricol(iHyp,:), 'LineWidth', 3);
    
    % Display error bar (over subjects) in each condition
    plot(repmat(mx, [1,2])', (my+sy*[-1,1])', 'k-');
    plot((mx+sx*[-1,1])', repmat(my, [1,2])', 'k-');
    
    % Display average 
    for k = 1:numel(mx)
        plot(mx(k), my(k), condmkr{iHyp}(k), 'MarkerFaceColor', ...
            condmap{iHyp}(k,:), 'MarkerEdgeColor', 'k', 'MarkerSize', 10);
    end
    
    % Display identity line
    plot(xlim, xlim, '-', 'Color', g); hold('on');
    
    % Customize the axes
    set(gca, 'Box', 'Off');
    if     iHyp == 1, ylim([0,0.012]);
    elseif iHyp == 2, ylim([0,80]);
    end
	
    % Add some text labels
    xlabel('Ideal observer');
    ylabel('Subjects'); 
end

% Save the figure
save2pdf(fullfile(ftapath, 'figs', 'F_Dyn_Corr1.pdf'));

% Display distribution of subject/IO correlation across subjects
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Prepare a new window
figure('Position', [1705 780 100 325]);

% For each regular hypothesis
for iHyp = 1:2
    
    % Compute the correlation coefficient between subjects and the ideal
    % observer
    x = mat2cell(Var{iHyp,1}, size(Var{iHyp,1}, 1), ones(1,nSub));
    y = mat2cell(Var{iHyp,2}, size(Var{iHyp,2}, 1), ones(1,nSub));
    cc = cellfun(@(x,y) Emergence_Regress(y, x, 'CC', 'r'), x, y);
    
    % Test if the correlation coefficient is significantly different
    % from zero across subjects
    [~,pval,tci,stats] = ttest(cc');
    disptstats(pval,tci,stats);
    
    % Display the dispersion of correlation coefficients across
    % subjects
    subplot(2,1,iHyp);
    Emergence_PlotSubGp(cc, tricol(iHyp,:));
    plot([0,2], zeros(1,2), 'k--');
    
    % Customize the axes
    axis([0,2,-1,1]);
    set(gca, 'XTick', [], 'XColor', 'None');
end

% Save the figure
save2pdf(fullfile(ftapath, 'figs', 'F_Dyn_Corr2.pdf'));
