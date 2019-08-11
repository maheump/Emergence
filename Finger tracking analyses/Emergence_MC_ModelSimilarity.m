% This script show similarity between the inference of pairs of models.
% Different measures of similarity can be used (correlation, distance, ...
% on cartesian or barycentric coordinates).
% 
% Copyright (c) 2018 Maxime Maheu

% Load model simulations
SimuType = 'TreeDepth';
Emergence_MC_ModelSimulations;

% Compute similarity between models
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Define which metric to compute
metric = {'PC', 'EC'}; % Pearson correlation and mean squared error
nMet = numel(metric);

% Prepare output variable
dim = size(pMgY);
corr = NaN([dim, dim(end), nMet]);

% For each pair of models
for iMod1 = 1:nMod
    for iMod2 = 1:nMod
        
        % For each metric
        for iMet = 1:nMet
            
            % Measure similarity between posterior probabilities of
            % different models
            corr(:,:,iMod1,iMod2,iMet) = cellfun(@(m1, m2) ...
                Emergence_Similarity(m1, m2, metric{iMet}), ...
                pMgY(:,:,iMod1), pMgY(:,:,iMod2));
        end
    end
end

% Prepare output variable
nReg = numel(cidx)+1;
R = NaN(nSub,nMod,nMod,nReg,nMet);

% For each type of sequences
for iReg = 1:nReg
    
    % Choose which sequences to average over
    if iReg == nReg, idx = cat(2, cidx{:}); % all sequences
    else,            idx = cidx{iReg};      % sequences of one type
    end
    
    % Average over sequences
    R(:,:,:,iReg,:) = squeeze(mean(corr(idx,:,:,:,:), 1));
end

% Average over (pseudo-)subjects
AvgMat = squeeze(mean(R,1));
SemMat = squeeze(sem(R,1));

% Reorder by tree depth
[~,sidx] = sort(cell2mat(options(3,:)));
AvgMat = AvgMat(sidx,sidx,:,:);
SemMat = SemMat(sidx,sidx,:,:);

% Display results as a double-entry matrix
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Define colormaps to use
cmap = cat(3, flipud(autumn(200)), autumn(200));

% Prepare a new window
figure('Position', [661 128 600 900]);
ttl = cat(2, proclab, 'All');

% For each type of sequence
for iReg = 1:nReg
    
    % For each matrix
    for iMet = 1:nMet
        sp = subplot(nReg,nMet,iMet+nMet*(iReg-1));
        
        % Display the correlation matrix
        imagesc(1:nMod, 1:nMod, AvgMat(:,:,iReg,iMet)); hold('on');
        
        % Overlap coefficients as text labels
        text(repmat(1:nMod, 1, nMod), sort(repmat(1:nMod, 1, nMod)), ...
            cellfun(@(x) sprintf('%1.2f', x), num2cell(AvgMat(:,:,iReg,iMet)), ...
            'uni', 0), 'HorizontalAlignment', 'Center', ...
            'VerticalAlignment', 'Bottom', 'FontSize', 8);
        text(repmat(1:nMod, 1, nMod), sort(repmat(1:nMod, 1, nMod)), ...
            cellfun(@(x) sprintf('(%1.3f)', x), num2cell(SemMat(:,:,iReg,iMet)), ...
            'uni', 0), 'HorizontalAlignment', 'Center', ...
            'VerticalAlignment', 'Top', 'FontSize', 6);
        
        % Customize the colormap
        colormap(sp, cmap(:,:,iMet));
        cbr = colorbar;
        
        % Customize the axes
        axis('square'); axis('xy');
        set(gca, 'XTick', 1:nMod, 'YTick', 1:nMod);
        
        % Add some text labels
        cbr.Label.String = metric{iMet};
        xlabel('Model #1'); ylabel('Model #2');
        title(ttl{iReg});
    end
end
