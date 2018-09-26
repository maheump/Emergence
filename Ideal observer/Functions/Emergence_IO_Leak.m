function leak = Emergence_IO_Leak( pErr, N )
% EMERGENCE_IO_LEAK returns a 1xN vector of exponentialy decaying weights.
% This can be use to implement a leaky memory. Technicaly speaking, it uses
% probability of substitution instead of genuine exponential leak (as what
% we used in Meyniel, Maheu & Dehaene, PCB, 2016).
%   - "pErr": the probability of making a mistake at each observation put
%       in memory (a scalar E [0,1]). 
%   - "N": the number of observations in the sequence (an integer > 0).
% 
% Example:
%   >> leak = ComputeLeak(1/5, 10);
%   => [.51 .51 .51 .52 .54 .56 .61 .68 .80  1]
%        #1  #2  #3  #4  #5  #6  #7  #8  #9 #10
%        past observations ...<==============| current observation
%        as we remember the past, the weights of each observation (in the
%        current estimates) decrease such that the most recent observations
%        have more influence in the current estimate.
% 
% Copyright (c) 2018 Maxime Maheu

% Prepare the output variable
leak = NaN(1,N);

% For each observation
for k = 1:N
    
    % Probability that there was no error up to element k
    pnochange = (1-pErr) ^ (N-k);
    
    % Probability that there was an error up to element k
    % (or that this error was corrected)
    j = 1:((N-k)/2);
    Nj = numel(j); C = NaN(1,Nj); % loop over all past positions
    for i = 1:Nj, C(i) = BinCoef(N-k, 2*j(i)); end % faster than "arrayfun"
    pchange = sum(C .* pErr.^(2*j) .* (1-pErr).^(N-2*j-k));
    
    % Combine both probabilities
    leak(k) = pnochange + pchange;
end

% Custom nchoosek function (see the MATLAB function)
function c = BinCoef( n, k )

nums = (n-k+1):n;
dens = 1:k;
nums = nums./dens;
c = round(prod(nums));

end

end