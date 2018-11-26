% This simulation script compares the posterior beliefs of the
% deterministic Bayesian observer learning repeating rules using different
% tree depth (i.e. the longest possible rule that is considered).
%
% Copyright (c) 2018 Maxime Maheu

%% INITIALIZATION
%  ==============

% Clear the place
clear;
close('all');

% Add ideal observer functions to the MATLAB path
scriptpath = mfilename('fullpath');
ind = strfind(scriptpath,'Emergence');
folderpath = scriptpath(1:ind(end-1)+8);
addpath(genpath(folderpath));

% Set default figure properties
Emergence_DefaultFigureProperties;

%% DEFINE SEQUENCES' PROPERTIES
%  ============================

% Whether to run simulations on simulated sequences (opt = 1) or on
% sequences that were presented to subjects (opt = 2)
opt = 2;

% Newly simulated sequences
% ~~~~~~~~~~~~~~~~~~~~~~~~~
if opt == 1
    
    % Define the rules to use in the simulations
    dr = {'AABB', ...
          'AAABBB', 'AABABB', 'AAABAB'...
          'AAAABBBB', 'AABABABB', 'AAABBABB', ...
          'AAAAABBBBB', 'AABABABABB', 'AAABAABBAB'};
    
    % Define the number of sequences to simulate for each rule
    % N.B. The results can be interpreted as soon as with nSeq = 2 but the
    % asymptotical results are reached for nSeq = 50
    nSeq = 50;
    
    % Define the position of the change point (it is always the same)
    cp = 100;
    
    % Define the length of the sequences
    N = 200;

% Sequences presented to subjects
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
elseif opt == 2
    
    % Load data
    Emergence_FTA_LoadData;
    
    % The number of simulations equal the number of subjects
    nSeq = nSub;
end

% Get the number of deterministic regularities
nReg = numel(dr);

%% DEFINE PROPERTIES OF THE IDEAL OBSERVER
%  =======================================

% Define the depth of the trees to investigate
Nu = [unique(cellfun(@numel, dr)') 20 50 N];
nNu = numel(Nu);

% Define options of the ideal observer
pEd   = 0; % probability of making a memory error at each observation
pEp   = 0; % probability of making a memory error at each observation
treed = 10; % depth of the rules' tree to explore
stat  = 'Transitions'; % statistic to be learned by the probabilistic model
p_pR  = 'Size-principle'; % the prior probability of each rule depends on its length
p_pT  = 'Bayes-Laplace'; % the prior over statistics to be learnt
p_pJ  = 'Uniform'; % prior over change point's position
comp  = 'all'; % compute after each observation
scale = 'log'; % scale of the model evidence
verb  = 0; % do not output messages in the command window
pgrid = []; % do not ask for the posterior distributions

%% SIMULATION LOOP
%  ===============

% Prepare the output variable
pMdgY = NaN(N,nSeq,nReg,nNu);

% For each repeating rule
for iReg = 1:nReg

    % For each simulated sequence
    for iSeq = 1:nSeq

        % Get the sequence
        if opt == 1 % generate
            nrepet = ceil((N-cp)/numel(dr{iReg}));
            seq = [GenRandSeq(cp, 1/2), ...                 % first part
                   repmat(str2pat(dr{iReg}), [1,nrepet])];  % second part
        elseif opt == 2 % retrieve
            seq = G{cidx{2}(iReg),iSeq}.Seq;
        end

        % For each depth of the tree
        for iNu = 1:nNu

            % Display the status of the loop in the command window
            fprintf(['- Rule %2.0f/%2.0f (%s), Sequence %2.0f/%2.0f, ', ...
                'Depth %1.0f/%1.0f (nu = %3.0f)... '], ...
                iReg, nReg, dr{iReg}, iSeq, nSeq, iNu, nNu, Nu(iNu));

            % Run the ideal observer
            io = Emergence_IO_FullIO(seq(1:N), pEd, pEp, Nu(iNu), ...
                stat, p_pR, p_pT, p_pJ, comp, scale, pgrid, verb);
            pMdgY(:,iSeq,iReg,iNu) = io.pMsdgY;
        end
    end
end

%% COMPUTE CONFUSION MATRIX
%  ========================

% Merge rules and sequences dimensions
rs_pMdgY = reshape(pMdgY, [N, nSeq*nReg, nNu]);

% Prepare the output variables
CorMat = NaN(nNu, nNu, nSeq*nReg);
MseMat = CorMat;

% For each simulated sequence
for iSimu = 1:nSeq*nReg

    % For each pair of tree depth
    for iNu1 = 1:nNu
        for iNu2 = 1:nNu

            % Get posterior beliefs of observers using different tree depth
            x = rs_pMdgY(2:end,iSimu,iNu1);
            y = rs_pMdgY(2:end,iSimu,iNu2);

            % Measure the correlation between the two
            CorMat(iNu1,iNu2,iSimu) = corr(x, y);

            % Measure the mean squared difference between the two
            MseMat(iNu1,iNu2,iSimu) = -log(mean((x - y) .^ 2));
        end
    end
end

% Average over simulations
AvgMat = cat(3, mean(CorMat, 3), mean(MseMat, 3));
SemMat = cat(3, sem(CorMat, 3), sem(MseMat, 3));

%% DISPLAY SIMULATION RESULTS
%  ==========================

% Prepare a new window
figure('Units', 'Normalized', 'Position', [0.35 0.3 0.3 0.35]);

% Define useful variables
cmap = cat(3, cbrewer2('Spectral', 101), (cbrewer2('PuBuGn', 101)));
labfun = @(x) sprintf('p(\\mathcal{M}_{\\mathrm{D}}|y,\\nu_{%1.0f})', x);
cbrlab = {sprintf('$\\rho(%s,%s)$', labfun(1), labfun(2)), ...
          sprintf('$-\\log((%s-%s)^{2})$', labfun(1), labfun(2))};

% For each matrix
for iSp = 1:2
    sp = subplot(1,2,iSp);

    % Display the correlation matrix
    imagesc(1:nNu, 1:nNu, AvgMat(:,:,iSp)); hold('on');

    % Overlap coefficients
    text(repmat(1:nNu, 1, nNu), sort(repmat(1:nNu, 1, nNu)), cellfun(@(x) ...
        sprintf('%1.2f', x), num2cell(AvgMat(:,:,iSp)), ...
        'UniformOutput', 0), 'HorizontalAlignment', 'Center', ...
        'VerticalAlignment', 'Bottom', 'FontSize', 8);
    text(repmat(1:nNu, 1, nNu), sort(repmat(1:nNu, 1, nNu)), cellfun(@(x) ...
        sprintf('%1.2f', x), num2cell(SemMat(:,:,iSp)), ...
        'UniformOutput', 0), 'HorizontalAlignment', 'Center', ...
        'VerticalAlignment', 'Top', 'FontSize', 6);

    % Customize the colormap
    cbr = colorbar('Location', 'SouthOutside');
    colormap(sp, flipud(cmap(:,:,iSp)));
    if iSp == 1, caxis([-1,1]); end

    % Customize the axes
    axis('square'); axis('xy'); set(gca, 'TickLabelInterpreter', 'LaTeX');
    set(gca, 'XTick', 1:nNu, 'XTickLabel', cellfun(@(x) sprintf('%1.0f', x), ...
        num2cell(Nu), 'UniformOutput', 0));
    set(gca, 'YTick', 1:nNu, 'YTickLabel', get(gca, 'XTickLabel'));

    % Add some text labels
    cbr.Label.String = cbrlab{iSp};
    cbr.Label.Interpreter = 'LaTeX';
    xlabel('$\nu_{1}$', 'Interpreter', 'LaTeX');
    ylabel('$\nu_{2}$', 'Interpreter', 'LaTeX');
    if     iSp == 1, t = title('Average correlation coefficients');
    elseif iSp == 2, t = title('Average log-mean squared difference');
    end
    set(t, 'Interpreter', 'LaTeX');
end
