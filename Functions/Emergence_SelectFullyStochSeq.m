function randidx = Emergence_SelectFullyStochSeq(data, filter, restopt)
% EMERGENCE_SELECTFULLYSTOCHSEQ selects fully-stochastic parts of sequences
%   - "data": the cell matrix of group data structures.
%   - "filter": a matrix specifying which fully-stochastic sequences were
%       accurately classified by the subjects.
%   - "restopt": a scalar specifying 
%       * 1: all fully-stochastic parts
%       * 2: only fully-stochastic sequences
%       * 3: only fully-stochastic sequences that were correctly labeled
%       * 4: only first-part of stochastic-to-regular sequences
% 
% Copyright (c) 2018 Maxime Maheu

% Select moment in which sequences do not entail any regularities
randidx = cellfun(@(x) (x.Gen == 1)', data, 'uni', 0);

% Get fully-stochastic sequences
seqfs = all(cellfun(@(x) strcmpi(x.Cond(1), 'S'), data), 2);

% Get fully-stochastic sequences accurately classified by subjects
subacclas = (filter{3} == 1);

% For sequences to remove, fill with falses
N = numel(data{1}.Gen);
fill = false(N,1);

% Restric to fully-stochastic sequences
if restopt == 2 || restopt == 3
    randidx(~seqfs,:) = {fill};
    
    % Restric to those that were correctly labeled
    if restopt == 3
        detecmask = zeros(size(data));
        detecmask(seqfs,:) = subacclas;
        randidx(~detecmask) = {fill};
    end
    
% Restric to first part of stochastic-to-regular sequences
elseif restopt == 4
    randidx(seqfs,:) = {fill};
end

end