function [ dist ] = Emergence_GetModelCategorisation( pMgY, prec )
% EMERGENCE_GETMODELCATEGORISATION returns the most likely hypothesis
% given the set of posterior probabilities over hypotheses. The function
% accounts for possible indecision cases between pairs of hypotheses by
% splitting the counts and attributing it in half/half proportion to each
% of the two hypotheses.
%   - "pMgY: a Nx3 matrix (where N is oftten the number of subjects)
%       specifying the posterior over hypotheses for a set of sequences.
% 
% Copyright (c) 2020 Maxime Maheu

% Get the number of subjects
nSub = size(pMgY, 1);

% Get the first and second best estimated generative processes
[val,idx] = sort(pMgY, 2);
EGP1 = idx(:,end);   % first 
EGP2 = idx(:,end-1); % second

% Get best estimated generative processes and corresponding p(H|seq)
% 1: statistical bias
% 2: deterministic rule
% 3: fully random
EstGenProc = EGP1;
MAP = val(:,end);

% Get indecision cases
if nargin < 2, prec = 0.01; end
Indecision = find(MAP >= (1/2 - prec) & MAP <= (1/2 + prec));

% Label the indecision cases
Nindec = numel(Indecision);
for i = 1:Nindec
    i = Indecision(i);
    if     (EGP1(i) == 1 & EGP2(i) == 2), EstGenProc(i) = 4; % stat/rule
    elseif (EGP1(i) == 2 & EGP2(i) == 1), EstGenProc(i) = 4; % stat/rule
    elseif (EGP1(i) == 1 & EGP2(i) == 3), EstGenProc(i) = 5; % stat/rand
    elseif (EGP1(i) == 3 & EGP2(i) == 1), EstGenProc(i) = 5; % stat/rand
    elseif (EGP1(i) == 2 & EGP2(i) == 3), EstGenProc(i) = 6; % rule/rand
    elseif (EGP1(i) == 3 & EGP2(i) == 2), EstGenProc(i) = 6; % rule/rand
    end
end

% Get the frequency of each 
dist = histc(EstGenProc, 1:6);

% Indecision between statistical bias and deterministic rule
for i = 4:6
    Nindec = dist(i) ./ 2;
    if i == 4
        dist(1) = dist(1) + Nindec; % statistical bias
        dist(2) = dist(2) + Nindec; % deterministic rule
    elseif i == 5
        dist(1) = dist(1) + Nindec; % statistical bias
        dist(3) = dist(3) + Nindec; % fully random
    elseif i == 6
        dist(2) = dist(2) + Nindec; % deterministic rule
        dist(3) = dist(3) + Nindec; % fully random
    end
end

% Remove indecision cases
dist = dist(1:3);

% Normalize such that the distribution sums to 1
dist = dist ./ nSub;

end