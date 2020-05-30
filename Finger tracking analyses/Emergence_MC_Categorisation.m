% This script displays categorisation of each sequence for each model. It
% accounts for possible indecision between the two non-random hypotheses as
% it can happen in the models that those two hyoptheses are judged as
% equally likely at the end of the sequence.
% 
% Copyright (c) 2020 Maxime Maheu

% Run simulations
% ~~~~~~~~~~~~~~~

% Get model simulations
SimuType = 'PseudoDeterministic';
Emergence_MC_ModelSimulations;

% Get categorisation profiles
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Get subject categorisations
subc = cellfun(@(x) x.Questions(2), G);
subc(isnan(subc)) = 0; % random = 0 / stat = 1 / deter = 2

% For each model and each sequence, get best hypothesis and corresponding
% posterior probability
[bestP,bestM] = cellfun(@(x) max(x(end,:)), pMgY, 'uni', 1);
bestM(bestM == 3) = 0; % random = 0 / stat = 1 / deter = 2

% Define indecision between non-random hypotheses
prec = 1/100;
bestM(bestP >= (1/2 - prec) & bestP <= (1/2 + prec)) = 1.5; % indecision = 1.5

% Concatenate subjects and models
bestM = cat(3, subc, bestM);

% Get histogramms of categorisations
distM = mat2cell(bestM, ones(nSeq,1,1), nSub, ones(1,1,nMod+1));
distM = cell2mat(cellfun(@(x) histc(x, [0,1,1.5,2]), distM, 'uni', 0));
distM = distM ./ nSub;

% Display categorisation profiles
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Define significiance thresholds
lima  = [0.001 0.01 0.05, 1.1];
stars = {'***' '**' '*', ''};

% For each type of sequence
for iReg = 1:3
    
    % Get corresponding indices
    idx = cidx{iReg};
    if iReg == 3, idx = nSeq; end
    
    % Create a new window
    n = numel(idx);
    figure('Position', [119+560*(iReg-1) 442 560 60*n]);
    
    % For each sequence
    for iSeq = 1:n
        
        % For each model
        for iMod = 1:nMod+1
            subplot(n, nMod+1, iMod + (nMod+1) * (iSeq-1));
            
            % Get categorisation profiles
            if iReg <= 2 % display all individual sequences
                i = cidx{iReg}(iSeq);
                b = squeeze(distM(i,:,iMod))';
            elseif iReg == 3 % average over sequences
                b = mean(distM(cidx{iReg},:,iMod), 1);
            end
            
            % Display the pie chart
            p = pie(b); hold('on');
            
            % Remove text labels
            delete(findobj(p, 'Type', 'text'));
            
            % Define colors
            colormap(cat(1, tricol(3,:), tricol(1,:), [147,112,219]./255, tricol(2,:)));
            
            % If the stat function is available
            if exist('myBinomTest', 'file') == 2
                
                % Get the result of a binomial test for correct categorisation
                % against chance level
                if     iReg == 1, nSuccess = b(2) * nSub; % stats
                elseif iReg == 2, nSuccess = b(4) * nSub; % rules
                elseif iReg == 3, nSuccess = b(1) * nSub; % random
                end
                chance = 1/3;
                if nSuccess > (chance * nSub)
                    pval = myBinomTest(nSuccess, nSub, chance, 'one');
                else, pval = 1;
                end
                
                % Display the significance level using stars
                pidx = find(pval < lima, 1, 'first');
                text(0, -1/2, stars{pidx}, ...
                    'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle');
            end
            
            % Display the generative process of the sequence
            if     iReg == 1, lab = sprintf('[%1.2f,%1.2f]', r{i});
            elseif iReg == 2, lab = sprintf('%d', r{i});
            elseif iReg == 3, lab = 'Rand';
            end
            if iMod == 1
                text(-3, 0, lab, 'Rotation', 90, ...
                    'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle');
            end
            
            % Display the models' name
            if iSeq == 1
                if contains(SimuType, 'PseudoDeterministic') ...
                || contains(SimuType, 'DifferentPriors'), idx = 4;
                elseif contains(SimuType, 'Independent'), idx = 12;
                elseif contains(SimuType, 'TreeDepth'),   idx = 3;
                elseif contains(SimuType, 'Leak'),        idx = 2;
                end
                if     iMod == 1, lgd = 'Subjects';
                elseif iMod  > 1
                    x = options{idx,iMod-1};
                    if ~ischar(x), lgd = sprintf('%1.2f', x);
                    else,          lgd = x;
                    end
                end
                title(lgd, 'FontWeight', 'Normal');
            end
        end
    end
end
