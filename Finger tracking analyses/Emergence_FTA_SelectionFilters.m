% This scripts shows which sequences are used in the different types of
% analyses. In particular, it shows which sequences were accurately
% classified by subjects and/or which for which sequences we were able to
% find a detection point for the subjects and the ideal observer.
% 
% Copyright (c) 2020 Maxime Maheu

% Define colors to use
cmap = Emergence_Colormap('Greys');

% Create a new window
figure('Position', [522 290 879 576]);
for iHyp = 1:3
    subplot(3,1,iHyp);
    
    % Display selection filters
    map = filter{iHyp};
    map = [map; NaN(1,size(map,2))];
    map = [map, NaN(size(map,1),1)];
    pcolor(map);
    colormap(cmap); caxis([0,3]);
    
    % Customize the axes
    axis([0,nSub+1,0,numel(cidx{iHyp}+1)])
    set(gca, 'XTick', [], 'YDir', 'Reverse');
    axis('equal'); axis('tight');
    
    % Add some text labels
    if     iHyp == 1, set(gca, 'YTick', 1:numel(pr), 'YTickLabel', pr);
    elseif iHyp == 2, set(gca, 'YTick', 1:numel(dr), 'YTickLabel', dr);
    elseif iHyp == 3, set(gca, 'YTick', []);
    end
    xlabel('Subjects');
    ylabel('Sequences');
    title(proclab{iHyp});
end

% Save the figure
save2pdf(fullfile(ftapath, 'figs', 'F_Sel_Conds.pdf'));