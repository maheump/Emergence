% This script shows detection dynamics of probabilistic and deterministic
% regularities locked onto change point position for different models
% previously simulated.
% 
% Copyright (c) 2018 Maxime Maheu

% Run simulations
% ~~~~~~~~~~~~~~~

% Get model simulations
SimuType = 'IndependentDifferenceDiscrete';
Emergence_MC_ModelSimulations;

% Define options
% ~~~~~~~~~~~~~~

% Define the type of plot
plttype = 'dyn'; % 'tri' or 'dyn'

% Whether to restrict to sequences that were accurately classified by
% subjects
onlydetreg = true;

% Whether to display the average sequence
dispseq = false;

% Display the other (incorrect) hypotheses
dispalt = false;

% Window around change point to look into (in number of observations)
ObsWin = 0:60;
nObs = numel(ObsWin);

% Display detection dynamics
% ~~~~~~~~~~~~~~~~~~~~~~~~~~

% Define positions of figures' subplots
SpPos = {[6,7,8,10,11,12,14,15,16,2,1], [1,6,10,14,7,11,3,8,12,16]};

% Prepare output variables
AvgSeq = cell(1,2);
AvgPP  = cell(1,2);
ErrPP  = cell(1,2);

% For each type of regularity
for iReg = 1:2
    
    % Get change point position
    cp = cellfun(@(x) x.Jump + 1/2, G(cidx{iReg},:), 'uni', 0);
    
    % Select sequence likelihood around change point's position
    PP = cellfun(@(p,c) Emergence_LockOnPoint(p,c,ObsWin), ...
        pMgY(cidx{iReg},:,:), repmat(cp, [1,1,nMod]), 'uni', 0);
    PP = cell2mat(cellfun(@(x) reshape(x,[1,1,1,nObs,3]), PP, 'uni', 0));
    
    % Get average sequence
    Seq = cellfun(@(x,i) Emergence_LockOnPoint(x.Seq',i,ObsWin), ...
        G(cidx{iReg},:), cp, 'uni', 0);
    Seq = cell2mat(cellfun(@(x) reshape(x,[1,1,nObs]), Seq, 'uni', 0));
    
    % Remove regularities that were not identified by subjects
    if onlydetreg
        detecmask = (filter{2} == 3);
        nanidx = repmat(~detecmask, [1,1,nMod,nObs,3]);
        PP(nanidx) = NaN;
        nanidx = repmat(~detecmask, [1,1,nObs]);
        Seq(nanidx) = NaN;
    end
    
    % Average over (pseudo-)subjects
    AvgSeq{iReg} = squeeze(round(mean(Seq, 2, 'OmitNaN')));
    AvgPP{iReg}  = squeeze(mean(PP, 2, 'OmitNaN'));
    ErrPP{iReg}  = squeeze(sem(PP, 2));
    
    % Prepare the window
    figure('Position', [200+760*(iReg-1) 565 760 540]);
    
    % For each regularity
    for iSeq = 1:numel(cidx{iReg})
    	subplot(4,4,SpPos{iReg}(iSeq));
        
        % Display information about the triangle
        if strcmpi(plttype, 'tri')
            Emergence_PlotTriInfo;
            axis('equal'); axis('off');
        end
        
        % For each model
        for iMod = [fliplr(1:nMod-1), nMod]
            
            % Get posterior likelihood in each hypothesis
            y = squeeze(AvgPP{iReg}(iSeq,iMod,:,:));
            e = squeeze(ErrPP{iReg}(iSeq,iMod,:,:));
            
            % Define the color to use
            c = modc(iMod,:);
            
            % As a trajectory in the triangle...
            if strcmpi(plttype, 'tri')
                
                % Display trajectory
                y = y * tricc;
                plot(y(:,1), y(:,2), '-', 'Color', c);
                
                % Overlap average sequence
                if dispseq
                    obscol = [c; ones(1,3)];
                    s = scatter(y(:,1), y(:,2), 15, ...
                        obscol(AvgSeq{iReg}(iSeq,:),:), 'Filled');
                    set(s, 'MarkerEdgeColor', c);
                end
                
            % As a dynamic time-course...
            elseif strcmpi(plttype, 'dyn')
                
                % Display temporal dynamics with shaded error area
                plotMSEM(ObsWin, y(:,iReg), e(:,iReg), 0.15, c, c);
                
                % Show beliefs in the two other hypotheses
                if dispalt
                    a = setdiff(1:2, iReg);
                    plotMSEM(ObsWin, y(:,a),   e(:,a),   0.15, c, c, 1, 1, '--');
                    plotMSEM(ObsWin, y(:,end), e(:,end), 0.15, c, c, 1, 1, ':');
                end
            end
        end
        
        % Overlap the sequence
        if strcmpi(plttype, 'dyn') && dispseq
            obscol = {'k'; 'w'};
            for i = 1:2
                idx = AvgSeq{iReg}(iSeq,:) == i;
                plot(ObsWin(idx), repmat(i-1, [1,sum(idx)]), 'o', ...
                    'MarkerEdgeColor', 'k', 'MarkerFaceColor', ...
                    obscol{i}, 'MarkerSize', 5);
            end
        end
        
        % Add some text labels
        if     iReg == 1, title(pr{iSeq});
        elseif iReg == 2, title(dr{iSeq});
        end
        set(gca, 'XTick', [], 'YTick', [], 'box', 'off'); title('');
        % Customize the axes
        if     strcmpi(plttype, 'dyn'), axis([ObsWin([1,end]),0,1]); 
        elseif strcmpi(plttype, 'tri'), axis([-1,1,0,sqrt(3)/2]);
        end
    end
end
