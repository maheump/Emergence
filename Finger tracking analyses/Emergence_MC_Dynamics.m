% This script shows detection dynamics of probabilistic and deterministic
% regularities locked onto different points of intereset in the sequence
% (i.e. change/detection/end points) for the subjects vs. different models
% previously simulated. Several options are available such as the type of
% data to plot (i.e. posterior probabilities or belief difference), the
% type of plot to use (i.e. in the triangle or not), etc.
%
% Copyright (c) 2020 Maxime Maheu

%% INITIALIZATION
%  ==============

% Run simulations
% ~~~~~~~~~~~~~~~

% Get model simulations
SimuType = 'PseudoDeterministic';
Emergence_MC_ModelSimulations;

% Define options
% ~~~~~~~~~~~~~~

% On which element to 
locktype = 'detection'; % 'change' or 'detection'

% Define the type of plot
plttype = 'post'; % 'tri', 'post', or 'beldif'

% How to organise the subplots in the figure
plotar = 'matrix'; % 'col' or 'matrix'

% Whether to restrict to sequences that were accurately classified by
% subjects
onlydetreg = true;

% Whether to display the models' names
displgd = false;

%% GET DYNAMICS LOCKED ONTO POINTS OF INTERESTS
%  ============================================

% Get model names
modlab = cat(2, 'Subjects', modlab);
if ~strcmpi(plttype, 'tri')
    modlab = repmat(modlab, [2,1]);
    modlab = reshape(modlab, [(1+nMod)*2,1]);
end

% Threshold for detection
thr = 1/2;

% Define x-axes for trajectories...
if strcmpi(locktype, 'change')
    xcp = 0:65;    % ... locked on change    point
    xdp = 0;       % ... locked on detection point
    xep = -5:0;    % ... locked on end       point
elseif strcmpi(locktype, 'detection')
    xcp = 0:5;     % ... locked on change    point
    xdp = -25:25;  % ... locked on detection point
    xep = -5:0;    % ... locked on end       point
end
xXp = {xcp,xdp,xep};

% Prepare the output variables
subposwrtp = cell(2,3); % subjects
modposwrtp = cell(2,3); % models

% For each type of regularity
for iReg = 1:2
    
    % Get fixed sequence points
    % ~~~~~~~~~~~~~~~~~~~~~~~~~
    
    % Get end point positions
    ep = repmat(N, numel(cidx{iReg}), nSub);
    
    % Get change point positions
    cp = cellfun(@(x) x.Jump, G(cidx{iReg},:));
    
    % Get subjects' detection points
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    % Lock trajectories on change point
    sublockbel = cellfun(@(x,c) Emergence_LockOnPoint(x.BarycCoord(:,iReg), c), ...
            G(cidx{iReg},:), num2cell(cp), 'uni', false);
    
    % Get detection point position
    sublag = cellfun(@(p) Emergence_FindDetecPoint(p, thr), sublockbel);
    subdp = cp + sublag;
    
    % Get models' detection points
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    % Lock trajectories on change point
    modlockbel = cellfun(@(x,c) Emergence_LockOnPoint(x(:,iReg), c), ...
            pMgY(cidx{iReg},:,:), repmat(num2cell(cp), [1,1,nMod]), 'uni', false);
    
    % Get detection point position
    modlag = cellfun(@(p) Emergence_FindDetecPoint(p, thr), modlockbel);
    moddp = cp + modlag;
    
    % Remove sequences entailing undetected regularities
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if onlydetreg
        detecmask = (filter{iReg} == 3);
        cp(~detecmask) = NaN;
        subdp(~detecmask) = NaN;
        moddp(repmat(~detecmask,[1,1,nMod])) = NaN;
        ep(~detecmask) = NaN;
    end
    
    % For the important points we want to lock onto
    for iLock = 1:3
        
        % Get subjects' trajectories around the points of interest
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        % Get observations to lock onto
        if     iLock == 1, lockobs = num2cell(cp);       % lock on change point
        elseif iLock == 2, lockobs = num2cell(subdp);    % lock on detection point
        elseif iLock == 3, lockobs = num2cell(ep);       % lock on end point
        end
        
        % Get trajectory in that window of interest
        fing = cellfun(@(p,c) Emergence_LockOnPoint(p.BarycCoord, c, xXp{iLock}), ...
            G(cidx{iReg},:), lockobs, 'uni', 0);
        
        % Convert the cells into a big 3D matrix
        subposwrtp{iReg,iLock} = cell2mat(reshape(fing, [1,1,size(fing)]));
        
        % Get models' trajectories around the points of interest
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        % Models' detection points
        if     iLock == 1, lockobs = num2cell(repmat(cp, [1,1,nMod]));  % lock on change point
        elseif iLock == 2, lockobs = num2cell(moddp);                   % lock on detection point
        elseif iLock == 3, lockobs = num2cell(repmat(ep, [1,1,nMod]));  % lock on end point
        end
        
        % Get trajectory in that window of interest
        fing = cellfun(@(p,c) Emergence_LockOnPoint(p, c, xXp{iLock}), ...
            pMgY(cidx{iReg},:,:), lockobs, 'uni', 0);
        
        % Convert the cells into a big 3D matrix
        modposwrtp{iReg,iLock} = cell2mat(reshape(fing, [1,1,size(fing)]));
    end
end

%% DISPLAY DETECTION DYNAMICS
%  ==========================

% Initialization
% ~~~~~~~~~~~~~~

% Define positions of figures' subplots
if strcmpi(plotar, 'matrix')
    nrow = 4;
    ncol = 4;
    SpPos = {[6,7,8,10,11,12,14,15,16,2,1], [1,6,10,14,7,11,3,8,12,16]};
elseif strcmpi(plotar, 'col')
    nrow = 1;
    ncol = max(cellfun(@numel, cidx));
    SpPos = {1:numel(cidx{1}), 1:numel(cidx{2})};
end

% Useful plot variables
av = 0.15;
sc = zeros(1,3);

% For each type of regularity
for iReg = 1:2
	
    % Prepare the window
    if     strcmpi(plotar, 'matrix'), figure('Position', [0+960*(iReg-1) 494 960 611]);
    elseif strcmpi(plotar, 'col'),    figure('Position', [33 936-200*(iReg-1) 1856 125]);
    end
    l = cell(1,nMod+1);
    
    % For each regularity
    for iSeq = 1:numel(cidx{iReg})
    	subplot(nrow,ncol,SpPos{iReg}(iSeq));
        
        % Display information about the triangle
        if strcmpi(plttype, 'tri')
            Emergence_PlotTriInfo;
            axis('equal'); axis('off');
        end
        
        % Display dynamics locked on different points (or not)
        if     strcmpi(locktype, 'change'),    list = [1,3];
        elseif strcmpi(locktype, 'detection'), list = 1:3;
        end
        for i = 1:numel(list)
            iLock = list(i);
            
            % Get position of the trajectory on the x-axis
            if strcmpi(locktype, 'detection')
                if     iLock == 1, x = xXp{2}(1) - fliplr(1:numel(xXp{1}));
                elseif iLock == 2, x = xXp{2};
                elseif iLock == 3, x = xXp{2}(end) + (1:numel(xXp{3}));
                end
            elseif strcmpi(locktype, 'change')
                if     iLock == 1, x = xXp{1};
                elseif iLock == 3, x = xXp{1}(end) + (1:numel(xXp{3}));
                end
            end
            
            % Subjects
            % ~~~~~~~~
            
            % Get subjects belief difference
            if strcmpi(plttype, 'beldif')
                pHpgY = subposwrtp{iReg,iLock}(:,1,iSeq,:);
                pHdgY = subposwrtp{iReg,iLock}(:,2,iSeq,:);
                beldif = (pHdgY - pHpgY) ./ (pHdgY + pHpgY);
                Msub = mean(beldif, 4, 'OmitNan');
                Ssub = sem(beldif, 4);
                
            % Get subjects posterior probability
            else
                Msub = mean(subposwrtp{iReg,iLock}(:,:,iSeq,:), 4, 'OmitNan');
                Ssub = sem( subposwrtp{iReg,iLock}(:,:,iSeq,:), 4);
            end
            
            % As trajectory in the triangle
            if strcmpi(plttype, 'tri')
                m = Msub * tricc;
                l{1} = plot(m(:,1), m(:,2), '-', 'LineWidth', 2, 'Color', sc);
                
            % As temporal dynamics with shaded error area
            else
                if strcmpi(plttype, 'beldif')
                    m = Msub;
                    s = Ssub;
                elseif strcmpi(plttype, 'post')
                    m = Msub(:,iReg);
                    s = Ssub(:,iReg);
                end
                l{1} = plotMSEM(x, m, s, av/2, sc, sc, 1, 1, '--');
            end
            
            % Models
            % ~~~~~~
            
            % For each model
            for iMod = fliplr(1:nMod)
                
                % Get models' data
                if strcmpi(plttype, 'beldif')
                    pHpgY = modposwrtp{iReg,iLock}(:,1,iSeq,:,iMod);
                    pHdgY = modposwrtp{iReg,iLock}(:,2,iSeq,:,iMod);
                    beldif = (pHdgY - pHpgY) ./ (pHdgY + pHpgY);
                    Mmod = mean(beldif, 4, 'OmitNan');
                    Smod = sem(beldif, 4);
                else
                    Mmod = mean(modposwrtp{iReg,iLock}(:,:,iSeq,:,iMod), 4, 'OmitNaN');
                    Smod = sem( modposwrtp{iReg,iLock}(:,:,iSeq,:,iMod), 4);
                end
                
                % Define the color to use
                mc = modc(iMod,:);
                
                % As trajectory in the triangle
                if strcmpi(plttype, 'tri')
                    m = Mmod * tricc;
                    l{1+iMod} = plot(m(:,1), m(:,2), '-', 'Color', mc);
                
                % As temporal dynamics with shaded error area
                else
                    if strcmpi(plttype, 'beldif')
                        m = Mmod;
                        s = Smod;
                    elseif strcmpi(plttype, 'post')
                        m = Mmod(:,iReg);
                        s = Smod(:,iReg);
                    end
                    l{1+iMod} = plotMSEM(x, m, s, av, mc, mc, 2, 1, '-');
                end
            end
        end
        
        % Figure properties
        % ~~~~~~~~~~~~~~~~~
        
        % Scale the axes
        if strcmpi(plttype, 'tri')
            axis([0,1,0,sqrt(3)/2]);
        else
            if     strcmpi(plttype, 'post'),    limy = [0,1];
            elseif strcmpi(plttype, 'beldif'),  limy = [-1,1];
            end
            ylim(limy);
            if strcmpi(locktype, 'change')
                a = xXp{1}(end) + numel(xXp{3});
                xlim([0,a]);
                plot(repmat(xXp{1}(end)+1/2, 1, 2), limy, 'k-', 'LineWidth', 1);
                set(gca, 'XTick', [0,a], 'XTickLabel', {'change', 'end'}); 
            elseif strcmpi(locktype, 'detection')
                a = xXp{2}(1)   - numel(xXp{1});
                b = xXp{2}(end) + numel(xXp{3});
                xlim([a,b]);
                plot(repmat(xXp{2}(1)  -1/2, 1, 2), limy, 'k-', 'LineWidth', 1);
                plot(repmat(xXp{2}(end)+1/2, 1, 2), limy, 'k-', 'LineWidth', 1);
                plot(repmat(a,1,2), limy, 'k--', 'LineWidth', 1);
                plot(zeros(1,2),    limy, 'k--', 'LineWidth', 1);
                plot(repmat(b,1,2), limy, 'k--', 'LineWidth', 1);
                plot([a,b], repmat(thr,1,2), 'k:', 'LineWidth', 1);
                set(gca, 'XTick', [a,0,b], 'XTickLabel', {'change', 'detection', 'end'});
            end
        end
        
        % Display the number of subjects who have detected the regularity
        nDetec = sum(filter{iReg}(iSeq,:) == 3);
        text(0, limy(1), sprintf('N_{sub} = %2.0f/%2.0f', nDetec, nSub), ...
            'HorizontalAlignment', 'Left', 'VerticalAlignment', 'Bottom');
        
        % Add some text labels
        if     iReg == 1, title(pr{iSeq});
        elseif iReg == 2, title(dr{iSeq});
        end
        if strcmpi(plttype, 'tri')
            set(gca, 'XTick', [], 'YTick', []);
        else
            xlabel('Obs. # w.r.t. ... point');
            if strcmpi(plttype, 'beldif')
                ylabel('Belief difference');
            else
                ylabel('Posterior proba.');
            end
        end
        set(gca, 'Box', 'Off');
        
        % Display model names
        if SpPos{iReg}(iSeq) == 1 && displgd
            legend(modlab, 'Position', [0.1443 0.2386 0.1125 0.4324]);
        end
    end
end
