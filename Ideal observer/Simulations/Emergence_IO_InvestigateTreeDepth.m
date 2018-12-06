% This simulation script compares the posterior beliefs of the
% deterministic Bayesian observer learning repeating rules using different
% tree depth (i.e. the longest possible rule that is considered).
% 
% Copyright (c) 2018 Maxime Maheu

%% INITIALIZATION
%  ==============

% Clear the place
clear; close('all');

% Add ideal observer functions to the MATLAB path
scriptpath = mfilename('fullpath');
ind = strfind(scriptpath,'Emergence');
folderpath = scriptpath(1:ind(end-1)+8);
addpath(genpath(folderpath));

% Set default figure properties
Emergence_DefaultFigureProperties;

% Get sequences presented to subjects
Emergence_FTA_LoadData;
nSeq = size(G,1);

%% DEFINE PROPERTIES OF THE IDEAL OBSERVER
%  =======================================

% Define the depth of the trees to investigate
Nu = [unique(cellfun(@numel, dr)') 20 50 N];
nNu = numel(Nu);

% Define options of the ideal observer
pEd   = 0;                  % probability of making a memory error at each observation
pEp   = 0;                  % probability of making a memory error at each observation
treed = 10;                 % depth of the rules' tree to explore
stat  = 'Transitions';      % statistic to be learned by the probabilistic model
p_pR  = 'Size-principle';   % the prior probability of each rule depends on its length
p_pT  = 'Bayes-Laplace';    % the prior over statistics to be learnt
p_pJ  = 'Uniform';          % prior over change point's position
comp  = 'all';              % compute after each observation
scale = 'log';              % scale of the model evidence
verb  = 0;                  % do not output messages in the command window
pgrid = [];                 % do not ask for the posterior distributions

%% SIMULATION LOOP
%  ===============

% Prepare the output variable
pMdgY = NaN(N,nSub,nSeq,nNu);

% For each sequence
for iSeq = 1:nSeq
    
    % For each subject
    for iSub = 1:nSub
        
        % For each depth of the tree
        for iNu = 1:nNu
            
            % Display the status of the loop in the command window
            fprintf(['- Rule %2.0f/%2.0f, Sequence %2.0f/%2.0f, ', ...
                'Depth %1.0f/%1.0f (nu = %3.0f)... '], ...
                iSeq, nSeq, iSub, nSub, iNu, nNu, Nu(iNu));
            
            % Run the ideal observer
            io = Emergence_IO_FullIO(G{iSeq,iSub}.Seq, pEd, pEp, ...
                Nu(iNu), stat, p_pR, p_pT, p_pJ, comp, scale, pgrid, verb);
            pMdgY(:,iSub,iSeq,iNu) = io.pMsdgY;
        end
        
        % Save temporary file
        save('Emergence_IO_InvestigateTreeDepth_tmp.mat', 'pMdgY');
    end
end

% Save file containing all the simulations and delete temporary file
save('Emergence_IO_InvestigateTreeDepth.mat', 'pMdgY');
delete('Emergence_IO_InvestigateTreeDepth_tmp.mat');

%% COMPUTE SIMILARITY MATRICES
%  ===========================

% Whether to restrict the analysis to sequences that were correctly
% identified by subjects (allows comparison with other analyses in which
% data from only accurately detected regularities are shown)
RestToIdent = true;

% Get correctly labeled sequences
TrueSeqType = cellfun(@(x) find(strcmpi(x.Cond(1), {'P', 'D', 'S'})), G);
ReportedSeqType = cellfun(@(x) x.Questions(2), G);
ReportedSeqType(isnan(ReportedSeqType)) = 3;
CorrectIdent = TrueSeqType == ReportedSeqType;

% Prepare the output variables
CorMat = NaN(nNu,nNu,nSub,nSeq);
MseMat = CorMat;

% For each sequence
for iSeq = 1:nSeq
    
    % For each subject
    for iSub = 1:nSub
        
        % Remove data from sequences that were not correctly labeled by
        % human subjects
        if ~RestToIdent || RestToIdent && CorrectIdent(iSeq,iSub)
            
            % For each pair of tree depth
            for iNu1 = 1:nNu
                for iNu2 = 1:nNu
                    
                    % Get posterior beliefs of observers using different tree depth
                    x = pMdgY(2:end,iSub,iSeq,iNu1);
                    y = pMdgY(2:end,iSub,iSeq,iNu2);
                    
                    % Measure the correlation between the two
                    CorMat(iNu1,iNu2,iSub,iSeq) = corr(x, y);
                    
                    % Measure the mean squared difference between the two
                    MseMat(iNu1,iNu2,iSub,iSeq) = log(mean((x - y) .^ 2));
                end
            end
        end
    end
end

%% DISPLAY SIMULATION RESULTS
%  ==========================

% Define colormaps to use
cmap = cat(3, cbrewer2('Spectral', 101), flipud(cbrewer2('PuBuGn', 101)));

% Define text labels
labfun = @(x) sprintf('p(\\mathcal{M}_{\\mathrm{D}}|y,\\nu_{%1.0f})', x);
cbrlab = {{'Average correlation coefficient', sprintf('$\\rho(%s,%s)$', labfun(1), labfun(2))}, ...
          {'Average log-mean squared difference', sprintf('$-\\log((%s-%s)^{2})$', labfun(1), labfun(2))}};

% Prepare a new window
figure('Position', [661 128 600 900])

% For each type of sequence
for iReg = 1:3
    
    % Average over sequences
    AvgMat = cat(4, mean(CorMat(:,:,:,cidx{iReg}), 4, 'OmitNaN'), ...
                    mean(MseMat(:,:,:,cidx{iReg}), 4, 'OmitNaN'));
    
    % Average over (pseudo-)subjects
    SemMat = squeeze(sem( AvgMat, 3));
    AvgMat = squeeze(mean(AvgMat, 3));
    
    % For each matrix
    for iSp = 1:2
        sp = subplot(3,2,iSp+2*(iReg-1));
        
        % Display the correlation matrix
        imagesc(1:nNu, 1:nNu, AvgMat(:,:,iSp)); hold('on');
        
        % Overlap coefficients as text labels
        text(repmat(1:nNu, 1, nNu), sort(repmat(1:nNu, 1, nNu)), cellfun(@(x) ...
            sprintf('%1.2f', x), num2cell(AvgMat(:,:,iSp)), ...
            'UniformOutput', 0), 'HorizontalAlignment', 'Center', ...
            'VerticalAlignment', 'Bottom', 'FontSize', 8);
        text(repmat(1:nNu, 1, nNu), sort(repmat(1:nNu, 1, nNu)), cellfun(@(x) ...
            sprintf('(%1.3f)', x), num2cell(SemMat(:,:,iSp)), ...
            'UniformOutput', 0), 'HorizontalAlignment', 'Center', ...
            'VerticalAlignment', 'Top', 'FontSize', 6);
        
        % Customize the colormap
        colormap(sp, flipud(cmap(:,:,iSp)));
        if iReg == 3
            cbr = colorbar('Orientation', 'Horizontal');
            pos = get(gca, 'Position');
            set(cbr, 'Position', [pos(1), pos(2)-0.05, pos(3), 0.01]);
        end
        if iSp == 1, caxis([-1,1]); end
        
        % Customize the axes
        axis('square'); axis('xy'); set(gca, 'TickLabelInterpreter', 'LaTeX');
        set(gca, 'XTick', 1:nNu, 'XTickLabel', cellfun(@(x) sprintf('%1.0f', x), ...
            num2cell(Nu), 'UniformOutput', 0));
        set(gca, 'YTick', 1:nNu, 'YTickLabel', get(gca, 'XTickLabel'));
        
        % Add some text labels
        if iReg == 3
            cbr.Label.Interpreter = 'LaTeX';
            cbr.Label.String = cbrlab{iSp};
        end
        xlabel('$\nu_{1}$', 'Interpreter', 'LaTeX');
        ylabel('$\nu_{2}$', 'Interpreter', 'LaTeX');
        t = title(sprintf('%s sequence', proclab{iReg}), 'Interpreter', 'LaTeX');
    end
end
