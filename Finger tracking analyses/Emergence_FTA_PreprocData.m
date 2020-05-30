% This script preprocesses data of each subject (e.g. organize conditions in
% the same order for each subject, project the finger's position inside the
% triangle when it is slightly outside, ...), run the ideal observer on the
% sequences that were presented to the subjects and save all the required
% information (to run the analysis scripts) in a dedicated group-level
% MATLAB file.
%
% Copyright (c) 2018 Maxime Maheu

%% INITIALIZATION
%  ==============

% Add functions to the MATLAB path
scriptpath = mfilename('fullpath');
ind = strfind(scriptpath,'Emergence');
folderpath = scriptpath(1:ind(end-1)+8);
addpath(genpath(folderpath));

% Define the location of the data
datadir = fullfile(folderpath, 'Stimulation', 'data');

% Get subjects' list
subjects = dir(fullfile(datadir, 'Subject*'));
subjects = {subjects.name};
nSub = numel(subjects);

% Order the deterministic rules
det = {[1 1 2 2], ...             % AABB       => [An,Bn] (n = 2)
       [1 1 1 2 2 2], ...         % AAABBB     => [An,Bn] (n = 3)
       [1 1 2 1 2 2], ...         % AABABB     => [A2,BAn,B2] (n = 1)
       [1 1 1 2 1 2], ...         % AAABAB     => not compressible (n = 6)
       [1 1 1 1 2 2 2 2], ...     % AAAABBBB   => [An,Bn] (n = 4)
       [1 1 2 1 2 1 2 2], ...     % AABABABB   => [A2,BAn,B2] (n = 2)
       [1 1 1 2 2 1 2 2], ...     % AAABBABB   => not compressible (n = 8)
       [1 1 1 1 1 2 2 2 2 2], ... % AAAAABBBBB => [An,Bn] (n = 5)
       [1 1 2 1 2 1 2 1 2 2], ... % AABABABABB => [A2,BAn,B2] (n = 3)
       [1 1 1 2 1 1 2 2 1 2]};    % AAABAABBAB => not compressible (n = 10)

% Order probabilistic rules
prob = {[1/3 1/3], ... % p(A|B) = 1/3 & p(B|A) = 1/3 => Rep. freq. (low)
        [1/4 1/4], ... % p(A|B) = 1/4 & p(B|A) = 1/4 => Rep. freq. (med)
        [1/5 1/5], ... % p(A|B) = 1/5 & p(B|A) = 1/5 => Rep. freq. (high)
        [2/3 2/3], ... % p(A|B) = 2/3 & p(B|A) = 2/3 => Alt. freq. (low)
        [3/4 3/4], ... % p(A|B) = 3/4 & p(B|A) = 3/4 => Alt. freq. (med)
        [4/5 4/5], ... % p(A|B) = 4/5 & p(B|A) = 4/5 => Alt. freq. (high)
        [1/3 2/3], ... % p(A|B) = 1/3 & p(B|A) = 2/3 => Item freq. (low)
        [1/4 3/4], ... % p(A|B) = 1/4 & p(B|A) = 3/4 => Item freq. (medium)
        [1/5 4/5], ... % p(A|B) = 1/5 & p(B|A) = 4/5 => Item freq. (high)
        [1/2 1/4], ... % p(A|B) = 1/2 & p(B|A) = 1/4 => Square (bottom)
        [3/4 1/2]};    % p(A|B) = 3/4 & p(B|A) = 1/2 => Square (right)
%         [1/2 3/4], ... % p(A|B) = 1/2 & p(B|A) = 3/4 => Square (top)
%         [1/4 1/2], ... % p(A|B) = 1/4 & p(B|A) = 1/2 => Square (left)
% The last two were used only to assess within-subjects reliability (they
% are discarded for the analyses).

%% PREPARE A TRIANGLE ON WHICH TO PROJECT THE TRAJECTORIES
%  =======================================================

% Get data from the first subject
subloc  = fullfile(datadir, subjects{1});
datfile = dir(fullfile(subloc, '*_Experiment_All.mat'));
subfile = load(fullfile(subloc, datfile.name));

% Screen matrix
% (0,0) is bottom left
h = subfile.h_px; % pixels
w = subfile.w_px; % pixels
bg = ones(h,w,3);
tricoord = subfile.S.Tri.CoordRev;

% Create an image with the triangle overlapping the screen matrix
fig = figure('Color', zeros(1,3), 'Visible', 'Off', 'Units', 'Normalized', 'Position', [0,0,1,1]);
imagesc(bg); hold('on');
fill(tricoord(:,1), tricoord(:,2), 'k', 'FaceColor', 'k', 'LineWidth', 1/10);
axis('equal'); axis('off');

% Capture the axis as an image
f = getframe(gca);
f = frame2im(f);
close(fig);

% Threshold and resize the image to make it match the screen's resolution
f = f(:,:,1);
f = f(f(:,1) > 0, :);
f = f == 0;
f = imresize(f, [h,w]);

% Get the x/y cartesian coordinates of each pixel composing the triangle
tripxl = find(f == 1);
[I,J] = ind2sub([h,w], tripxl);
tripxl = [J,I];

%% LOAD DATA FROM EACH SUBJECT
%  ===========================

% For each subject
for iSub = 1:nSub
    fprintf('- Preprocessing data from subject #%2.0f/%2.0f.\n', iSub, nSub);
    
    % Get data
    % ~~~~~~~~
    
    % Get subject's data
    subloc  = fullfile(datadir, subjects{iSub});
    datfile = dir(fullfile(subloc, '*_Experiment_All.mat'));
    subfile = load(fullfile(subloc, datfile.name));
    
    % Order conditions
    % ~~~~~~~~~~~~~~~~
    
    % Order conditions
    rules = cellfun(@(x) x.Cond, subfile.D, 'uni', 0)';
    
    % Stochastic condition first
    stochidx = find(strcmpi(rules, 'Stochastic'))';
    
    % Probabilistic conditions sorted per entropy level and bias type
    probar = find(strcmpi(rules, 'Probabilistic'));
    transproba = cellfun(@(x) x.Rule, subfile.D(probar), 'uni', 0);
    probaidx = NaN(1,numel(prob));
    for k = 1:numel(prob)
        dif = NaN(1,numel(transproba));
        for j = 1:numel(transproba)
            dif(j) = sum(abs(transproba{j} - prob{k}));
        end
        [~,l] = min(dif);
        probaidx(k) = probar(l);
    end
    
    % Deterministic conditions sorted by patterns' length
    deterr = find(strcmpi(rules, 'Deterministic'));
    pattern = cellfun(@(x) x.Rule, subfile.D(deterr), 'uni', 0);
    deteridx = NaN(1,numel(det));
    for k = 1:numel(det)
        dif = NaN(1,numel(pattern));
        for j = 1:numel(pattern)
            if numel(pattern{j}) == numel(det{k})
                dif(j) = sum(abs(pattern{j} - det{k}));
            end
        end
        [~,l] = min(dif);
        deteridx(k) = deterr(l);
    end
    
    % Get ordered trials' list
    ortrlist = [stochidx, probaidx, deteridx];
    nSeq = numel(ortrlist);
    
    % Prepare the output variable
    if iSub == 1
        G = cell(nSeq, nSub);
        S = cell(1, nSub);
    end
    
    % Save data in a group variable
    G(:,iSub) = subfile.D(ortrlist)';
    
    % Do not forgot to save the "S" variable
    % (0,0) is bottom left
    S{iSub} = subfile.S;
    S{iSub}.Tri.w = S{iSub}.Tri.CoordRev(2,1) - S{iSub}.Tri.CoordRev(1,1);
    S{iSub}.Tri.h = S{iSub}.Tri.CoordRev(2,2) - S{iSub}.Tri.CoordRev(3,2);
    S{iSub}.Tri.y = h - S{iSub}.Tri.y;
    S{iSub}.Tri = rmfield(S{iSub}.Tri, 'Coord');
    
    % Loop over conditions
    % ~~~~~~~~~~~~~~~~~~~~
    
    % For each condition
    for iSeq = 1:nSeq
        
        % Get the number of observations in the sequence
        N = numel(G{iSeq,iSub}.Seq);
        
        % Correct trajectories to be sure that the finger is within the
        % triangle
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        % Get the finger's movements
        traj = G{iSeq,iSub}.DiscreteMouseCoord';
        
        % Correct trajectories
        [traj, corr] = Emergence_ProjOnTri(traj, tripxl);
        G{iSeq,iSub}.CartesCoord = traj;
        G{iSeq,iSub}.Corr = corr;
        
        % Convert these x/y coordinates into barycentric ones
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        lambda = cartes2baryc(traj', tricoord)';
        
        % Some sanity steps such that we do not have to check in the
        % analyses scripts
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        % Make sure that all barycentric coordinates...
        % 1) are larger than 0 and smaller than 1
        lambda(lambda < 0) = 0;
        lambda(lambda > 1) = 1;
        % 2) sum to 1
        lambda = lambda ./ sum(lambda, 2);
        
        % Order the 3 components as follows:
        % 1: probabilistic component
        % 2: deterministic component
        % 3: stochastic component
        G{iSeq,iSub}.BarycCoord = lambda;
    end
end

% Display the average distance
corr = mat2cell(G, nSeq, ones(1,nSub));
corr = cellfun(@(y) cell2mat(cellfun(@(x) x.Corr, y, 'uni', 0)), ...
    corr, 'uni', 0);
dist = cellfun(@(x) mean(x(:,2)), corr);
fprintf('=> Average projection distance: %1.2f pixels [%1.2f, %1.2f] (%1.2f).\n', ...
    mean(dist), min(dist), max(dist), sem(dist));

%% RUN THE IDEAL OBSERVER ON SEQUENCES THAT WERE PRESENTED TO THE SUBJECTS
%  =======================================================================

% Prepare the output variable
IO = cell(size(G)); % conditions x subjects cell matrix with IO's inference

% Define options for the observer
pEd    = 0;                 % probability of making a memory error at each observation
pEp    = 0;                 % probability of making a memory error at each observation
patlen = 10;                % depth of the rules' tree to explore
stat   = 'Transitions';     % statistic to be learned by the probabilistic model
pR     = 'Size-principle';  % the prior probability of each rule depends on its length
pT     = 'Bayes-Laplace';   % the prior over statistics to be learnt
pJ     = 'Uniform';         % prior over change point's position
comp   = 'all';             % compute after each observation
scale  = 'log';             % scale of the model evidence
pgrid  = [];                % precision of the posterior over theta
verb   = 0;                 % do not output messages in the command window

% For each subject
for iSub = 1:nSub
    
    % For each condition
    for iSeq = 1:nSeq
        fprintf('- Running the IO on sequence #%2.0f/%2.0f from subject #%2.0f/%2.0f... ', ...
            iSeq, nSeq, iSub, nSub);
        
        % Run the observer with these options
        io = Emergence_IO_FullIO(G{iSeq,iSub}.Seq, ... % sequence
            pEd, pEp, patlen, stat, pR, pT, pJ, comp, scale, pgrid, verb);  % options
        
        % Append the "Jump" subfield (the IO does not know its true position)
        IO{iSeq,iSub}.Jump = G{iSeq,iSub}.Jump;
        
        % Create a variable with all models' posterior probability
        IO{iSeq,iSub}.BarycCoord = [io.pMspgY', ... % 1: probabilistic component
                                    io.pMsdgY', ... % 2: deterministic component
                                    io.pMssgY'];    % 3: stochastic    component
        
        % Create a variable with sequence likelihood under each models
        IO{iSeq,iSub}.SeqLLH = [io.pYgMsp', ... % 1: probabilistic component
                                io.pYgMsd', ... % 2: deterministic component
                                io.pYgMss'];    % 3: stochastic    component
        
        % Create a variable with posterior distribution of chance point position
        IO{iSeq,iSub}.CPbelief = io.pJkgY;
        
        % Create a variable with equivalent learning rate
        IO{iSeq,iSub}.EqLearnRate = io.eqAlpha';
        
        % Create a variable with (corrected) cartesian coordinates
        traj = round(IO{iSeq,iSub}.BarycCoord * tricoord);
        IO{iSeq,iSub}.DiscreteMouseCoord = traj;
        IO{iSeq,iSub}.CartesCoord = Emergence_ProjOnTri(traj, tripxl);
        
        % Ask the same questions to the ideal observer as those that were
        % asked to the subjects
        IO{iSeq,iSub}.Questions = NaN(1,6);
        
        % 1) Regularity?
        if io.Mhat(end) == 1
            IO{iSeq,iSub}.Questions(1) = 2; % no regularity found
            
        % 2) Which type?
        else
            IO{iSeq,iSub}.Questions(1) = 1; % a regularity was found
            IO{iSeq,iSub}.Questions(2) = io.Mhat(end) - 1;
            % 1 for probabilistic, 2 for deterministic
        end
        
        % 3) p(change)?
        pChange = io.pJkgY(:,end);
        [val,idx] = max(pChange);
        IO{iSeq,iSub}.Questions(3) = idx; % change point's position
        IO{iSeq,iSub}.Questions(4) = val; % confidence in change point's position
    end
end

% Scale the confidence such that it evolves between 0 and 100 (as subjects)
minc = min(cellfun(@(x) x.Questions(4), IO(:)));
maxc = max(cellfun(@(x) x.Questions(4), IO(:)));
for i = 1:numel(IO)
    IO{i}.Questions(4) = round(100 * ((IO{i}.Questions(4) - minc) / (maxc - minc)));
end

%% CREATE LABELS FOR THE DIFFERENT CONDITIONS
%  ==========================================

% Get the different conditions
c = cellfun(@(x) x.Cond, G(:,1), 'uni', 0);
c = cellfun(@(x) sprintf('%s', x), c, 'uni', 0);

% Get the different regularities
r = cellfun(@(x) x.Rule, G(:,1), 'uni', 0);
sr = repmat({''}, 1, sum(strcmpi(c, 'Stochastic')))';
pr = cellfun(@(x) sprintf('p(A|B) = %1.2f & p(B|A) = %1.2f', x(1), ...
    x(2)), r(strcmpi(c, 'Probabilistic')), 'uni', 0);
dr = cellfun(@pat2str, r(strcmpi(c, 'Deterministic')), 'uni', 0);

% Combine conditions and regularities
condlab = cellfun(@(x,y) sprintf('%s %s', x, y), c, [sr;pr;dr], 'uni', 0);

% Create indices variables
proclab = {'Probabilistic', 'Deterministic', 'Stochastic'};
cidx = {11:21, 22:31, 1:10};

%% WRAP THINGS UP
%  ==============

% Add useful variables to the MATLAB file in which preprocessed data will
% be saved
tricol  = [049, 130, 189; 222, 045, 038; 049, 163, 084] ./ 255;
tricc   = [0, sqrt(3)/2; 1, sqrt(3)/2; 1/2, 0];
letters = {'A','B'};
g = repmat(0.7,1,3); % grey

% Save the group data file
savefolder = fullfile(folderpath, 'Finger tracking analyses');
mkdir(savefolder);
filename = fullfile(savefolder, 'Emergence_Behaviour_GroupData.mat');
save(filename, 'G', 'IO', 'N', 'S', 'f', 'cidx', 'condlab', 'folderpath', ...
    'c', 'r', 'sr', 'pr', 'dr', 'proclab', 'subjects', 'nSub', 'nSeq', ...
    'tricol', 'tricc', 'letters', 'prob', 'det', 'g');
