function H = Emergence_MarkovEntropy( pAgB, pBgA )
% EMERGENCE_MARKOVENTROPY computes Shannon entropy of a set of transition
% probabilities.
%   - "pAgB": a scalar specifying the value of p(A|B).
%   - "pBgA": a scalar specifying the value of p(B|A).
%
% Copyright (c) 2020 Maxime Maheu

% Compute frequency of all transitions
pAgA = 1 - pBgA;
pBgB = 1 - pAgB;

% Compute resulting p(A) and p(B)
pA = pAgB ./ (pAgB + pBgA);
pB = 1 - pA;

% Compute the 4 different parts of the entropy
H = [pA .* pAgA .* (log2(pA) + log2(pAgA)), ...
     pA .* pBgA .* (log2(pA) + log2(pBgA)), ...
     pB .* pAgB .* (log2(pB) + log2(pAgB)), ...
     pB .* pBgB .* (log2(pB) + log2(pBgB))];

% Get transitions for which the probability is null
na = [pAgA, pBgA, pAgB, pBgB] == 0;

% Replace the corresponding (uncomputablte) part in the entropy by the
% limit 0*log2(0) = 0
H(na) = 0;
 
% Compute Shannon entropy by summing over the 4 different parts
H = -sum(H);

end