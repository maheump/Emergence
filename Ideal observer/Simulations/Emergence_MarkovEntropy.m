function H = Emergence_MarkovEntropy( pAgB, pBgA )
%EMERGENCE_MARKOVENTROPY computes Shannon entropy of a set of transition
%probabilities.
%   - "pAgB": a scalar specifying the value of p(A|B).
%   - "pBgA": a scalar specifying the value of p(B|A).
%
% Copyright (c) 2018 Maxime Maheu

% Compute resulting p(A) and p(B)
pA = pAgB ./ (pAgB + pBgA);
pB = 1 - pA;

% Compute entropy
H = -(pA.*(1-pBgA).*(log(pA) + log(1-pBgA)) ...
    + pA.*pBgA.*(log(pA) + log(pBgA)) ...
    + pB.*pAgB.*(log(pB) + log(pAgB)) ...
    + pB.*(1-pAgB).*(log(pB) + log(1-pAgB)));

end