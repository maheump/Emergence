% This script shows the categorisation profiles of subjects and models,
% averaged over sequences. In particular, it shows that different models
% with different representations have different categorisation profiles.
% 
% Copyright (c) 2020 Maxime Maheu

% Whether to average over sequences or not
avgseq = false;

% Get alternative model simulations
SimuType = 'DifferentPriors';
Emergence_MC_ModelSimulations;

% Add subjects to the list of models
modlab = cat(2, modlab, 'Subjects');

% For each type of generative process
for iReg = 1:3
    
    % Prepare a new window
    if (avgseq && iReg == 1) || ~avgseq
        figure('Position', [33 84 912 989]);
    end
    
    % For each model (and the subjects)
    for iMod = 1:nMod+1
        
        % Get subjects' categorisation profiles
        if iMod == nMod+1
            estgenproc = cellfun(@(x) x.Questions(2), G(cidx{iReg},:));
            estgenproc(isnan(estgenproc)) = 3;
            
            % Compute categorisation profiles
            dist = histc(estgenproc', 1:3);
            dist = dist ./ nSub;
        
        % Get models' categorisation profiles
        else
            simu = pMgY(cidx{iReg},:,iMod);
            simu = cellfun(@(p) p(end,:), simu, 'uni', 0);
            simu = mat2cell(simu, ones(numel(cidx{iReg}),1), nSub);
            simu = cellfun(@(x) cell2mat(x'), simu, 'uni', 0);
            
            % Compute categorisation profiles
            [~,dist] = cellfun(@Emergence_GetModelCategorisation, simu, 'uni', 0);
            dist = cell2mat(dist');
        end
        
        % Average over sequences
        if avgseq, dist = mean(dist, 2); end
        
        % For each sequence
        nCase = size(dist, 2);
        for iSeq = 1:nCase
            if      avgseq, nRow = 3;     iRow = iReg; % just 1 (average) sequence
            elseif ~avgseq, nRow = nCase; iRow = iSeq; % the current sequence
            end
            subplot(nRow, nMod+1, iMod+(nMod+1)*(iRow-1));
            
            % Display the pie chart
            p = pie(dist([3,1,2],iSeq));

            % Remove text labels
            delete(findobj(p, 'Type', 'text'));
        
            % Define colors
            colormap(cat(1, tricol(3,:), tricol(1:2,:)));
            
            % Add some text labels
            if iMod == 1
                if      avgseq
                    text(-4, 0, proclab{iReg}, 'Rotation', 90, 'FontWeight', 'Bold');
                elseif ~avgseq
                    if     iReg == 1, seqlab = pr{iSeq};
                    elseif iReg == 2, seqlab = dr{iSeq};
                    elseif iReg == 3, seqlab = sprintf('Random sequence #%1.0f', iSeq);
                    end
                    text(-4, 0, seqlab, 'FontWeight', 'Bold');
                end
            end
            if iRow == 1
                text(0, 2, modlab{iMod}, 'FontWeight', 'Bold');
            end
            
            % If the stat function is available
            if exist('myBinomTest', 'file') == 2
                
                % Get the result of a binomial test for correct
                % categorisation against chance level
                nSuccess = dist(iReg,iSeq) * nSub;
                chance = 1/3;
                if nSuccess > (chance * nSub)
                    pval = myBinomTest(nSuccess, nSub, chance, 'one');
                else, pval = 1;
                end
                
                % Display the significance level using stars
                lima  = [0.001 0.01 0.05, 1.1];
                stars = {'***' '**' '*', ''};
                pidx = find(pval < lima, 1, 'first');
                text(0, -1/2, stars{pidx});
            end
        end
    end
end
