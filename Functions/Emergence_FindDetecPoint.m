function dp = Emergence_FindDetecPoint(bel)
% EMERGENCE_FINDDETECPOINT returns position of the detection point.
%   - "bel": the 1xN array of posterior probabilities as a function of
%       time.
% 
% Copyright (c) 2020 Maxime Maheu

% The belief threshold to exceed in order to "detect" a regularity
threshold = 1/2;

% CONSTRAINT 1: larger than threshold
idx1 = (bel >= threshold);

% CONSTRAINT 2: it should cross the threshold in a positive manner
% i.e. the one before should be smaller than the threshold
% CONSTRAINT 3: it should not be already above theshold at the time of
% change point
idx2 = [false; (bel(1:end-1) < threshold)];
logicdet = (idx1 & idx2);

% The detection point is the first threshold crossing
dp = find(logicdet, 1, 'First') - 1;

% If no detection points has been found, return a NaN
if isempty(dp), dp = NaN; end

end