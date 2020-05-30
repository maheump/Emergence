function Emergence_PrintFstats( RMtbl )
% EMERGENCE_PRINTFSTATS displays in the MATLAB command window the results
% of a repeated measures ANOVA.
%   - "RMtbl": a table from ranova.m.
% 
% Copyright (c) 2020 Maxime Maheu

% Get factors of interest
facidx = find(contains(RMtbl.Properties.RowNames, '(Intercept):'));
I = numel(facidx);

% Define significiance thresholds
lima  = [0.001 0.01 0.05 1];
stars = {'***' '**' '*' 'ns'};

% Compute estimates of effect size for each factor of interest
nt              = size(RMtbl, 1) / 2 - 1;
eta2            = NaN((nt+1)*2,1);
omega2          = NaN((nt+1)*2,1);
absent          = repmat([false; true], [nt+1,1]);
SStreatment     = RMtbl.SumSq(1:2:end);
SStotal         = sum(RMtbl.SumSq);
eta2(1:2:end)   = SStreatment ./ SStotal;
RMtbl.eta2      = internal.stats.DoubleTableColumn(eta2, absent);
DFtreatment     = RMtbl.DF(1:2:end);
MSerror         = RMtbl.MeanSq(2:2:end);
omega2(1:2:end) = (SStreatment - DFtreatment .* MSerror) ./ (SStotal + MSerror);
RMtbl.omega2    = internal.stats.DoubleTableColumn(omega2, absent);

% For each factor of interest
for i = 1:I
    
    % Print stats, dof and effect size
    name = RMtbl(facidx(i),:).Row{1}(13:end);
    efsz = RMtbl(facidx(i),:).omega2;
    df1  = RMtbl(facidx(i),:).DF;
    df2  = RMtbl(2,:).DF;
    stat = RMtbl(facidx(i),:).F;
    fprintf('%s: w2 = %1.2f, F(%1.0f,%1.0f) = %1.2f, ', name, efsz, df1, df2, stat);
    
    % Print p-value
    pval = RMtbl(facidx(i),:).pValue;
    idx = find(pval < lima, 1, 'first');
    if     pval >= 0.001, fprintf('p = %1.3f %s\n', pval, stars{idx});
    elseif pval <  0.001, fprintf('p = %1.2d %s\n', pval, stars{idx});
    end
end

end
