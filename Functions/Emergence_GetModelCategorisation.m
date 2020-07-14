function [ EstGenProc, dist, odist ] = Emergence_GetModelCategorisation( pMgY, prec )
% EMERGENCE_GETMODELCATEGORISATION returns the most likely hypothesis
% given the set of posterior probabilities over hypotheses. The function
% accounts for possible indecision cases between pairs of hypotheses by
% splitting the counts and attributing it in half/half proportion to each
% of the two hypotheses.
%   - "pMgY: a Nx3 matrix (where N is often the number of subjects)
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

% Deal with indecision between 2 hypotheses
if nargin < 2, prec = 0.001; end
Indecision2 = find(MAP >= (1/2 - prec) & MAP <= (1/2 + prec));
Nindec = numel(Indecision2);
for i = 1:Nindec
    i = Indecision2(i);
    if     (EGP1(i) == 2 & EGP2(i) == 3), EstGenProc(i) = 4; % rule/rand
    elseif (EGP1(i) == 3 & EGP2(i) == 2), EstGenProc(i) = 4; % rule/rand
    elseif (EGP1(i) == 1 & EGP2(i) == 3), EstGenProc(i) = 5; % stat/rand
    elseif (EGP1(i) == 3 & EGP2(i) == 1), EstGenProc(i) = 5; % stat/rand
    elseif (EGP1(i) == 1 & EGP2(i) == 2), EstGenProc(i) = 6; % stat/rule
    elseif (EGP1(i) == 2 & EGP2(i) == 1), EstGenProc(i) = 6; % stat/rule
    end
end

% Deal with indecision between 3 hypotheses
Indecision3 = all(MAP >= (1/3 - prec) & MAP <= (1/3 + prec), 2);
EstGenProc(Indecision3) = 7;

% If we have to export the distribution
if nargout > 1
    
    % Get the frequency of each cell in the reponse/sequence matrix
    odist = histc(EstGenProc, 1:7);
    
    % Indecision between statistical bias and deterministic rule
    for i = 4:7
        if     i  < 7, Nindec = odist(i) ./ 2; % between 2 hypotheses
        elseif i == 7, Nindec = odist(i) ./ 3; % between 3 hypotheses
        end
        if i == 4
            odist(2) = odist(2) + Nindec; % deterministic rule
            odist(3) = odist(3) + Nindec; % fully random
        elseif i == 6
            odist(1) = odist(1) + Nindec; % statistical bias
            odist(2) = odist(2) + Nindec; % deterministic rule
        elseif i == 5
            odist(1) = odist(1) + Nindec; % statistical bias
            odist(3) = odist(3) + Nindec; % fully random
        elseif i == 7
            odist = odist + Nindec; % all 3 hypotheses
        end
    end
    
    % Remove indecision cases
    dist = odist(1:3);
    
    % Normalize such that the distribution sums to 1
    dist = dist ./ nSub;
end

end