% Display the triangle as well as remarkable boundaries and locations of
% the triangle.
%
% Copyright (c) 2020 Maxime Maheu

%% INITIALIZATION
%  ==============

% Clear the place
clear;
close('all');

% Add functions to the MATLAB path
scriptpath = mfilename('fullpath');
ind = strfind(scriptpath,'Emergence');
folderpath = scriptpath(1:ind(end-1)+8);
addpath(genpath(folderpath));

% Set default figure properties
Emergence_DefaultFigureProperties;

% Coordinates of the triangles limits (cartesian coordinates)
tricc = [0, sqrt(3)/2; 1, sqrt(3)/2; 1/2, 0];

% Colors associated to each vertex of the triangle
tricol = [049, 130, 189; 222, 045, 038; 049, 163, 084] ./ 255;

% Label associated to each vertex of the triangle
m = {'P', 'D', 'S'};

% Boundary regarding posterior probabilities to highlight
thr = {'1/3', '1/2', '2/3'};

%% DRAW AN EXAMPLE TRIANGLE
%  ========================

% Whether or not to display landmarks
plotlandm = true;

% Whether or not to display an example sequence
plotexpseq = true;

% Prepare the window
figure;

% Display the triangle and its background
Emergence_PlotTrajOnTri(ones(1,3), NaN, tricol, 10, {'v', 'p', '^'}, false);

% Display landmarks
% ~~~~~~~~~~~~~~~~~
if plotlandm
    
    % For each side of the triangle
    for iSide = 1:3
        
        % For each boundary to display
        for iBnd = 1:numel(thr)
            
            % Specify that the first dimension should be equal to the currently
            % considered threshold
            x1 = eval(thr{iBnd});
            
            % Deduce what degrees of freedom we have left for the other
            % dimensions
            rem = 1-x1;
            
            % Precise how many scenarios we should look at. Since the boundary
            % is a straight line, we just need two points to draw it.
            dt = rem;
            % N.B. Use 0.01 if you want to convince you that it is a straight line!
            
            % Deduce the barycentric coordinates the other two dimensions can
            % take
            x2 = 0:dt:rem;
            x3 = rem:-dt:0;
            
            % Combine the three barycentric coordinates
            n = numel(x2);
            if     iSide == 1, bc = [repmat(x1, [n,1]), x2', x3'];
            elseif iSide == 2, bc = [x2', repmat(x1, [n,1]), x3'];
            elseif iSide == 3, bc = [x2', x3', repmat(x1, [n,1])];
            end
            
            % Convert barycentric to cartesian coordinates
            cc = bc*tricc;
            
            % Display the boundary
            alpha(2/3);
            plot(cc(:,1), cc(:,2), '-', 'Color', tricol(iSide,:), 'LineWidth', 2);
            
            % Specify where the put the text labels
            x = cc(1,1);
            y = cc(1,2);
            if     iSide == 1, loc = 'Right';
            elseif iSide == 2, loc = 'Left';
            elseif iSide == 3
                x = mean(cc(:,1));
                y = mean(cc(:,2));
                loc = 'Center';
            end
            
            % Display remarkable boundaries of the triangle
            text(x, y, ['$p(\mathcal{M}_{\mathrm{', m{iSide}, '}}) = ', ...
                thr{iBnd}, '$'], 'Interpreter', 'LaTeX', 'Color', tricol(iSide,:), ...
                'HorizontalAlignment', loc, 'VerticalAlignment', 'Top');
        end
    
        % Display remarkable dots of the triangle
        for i = 1:3
            for j = 1:3
                mx = mean(tricc([i,j],1));
                my = mean(tricc([i,j],2));
                col = mean(tricol([i,j],:),1);
                plot(mx, my, 'ko', 'MarkerFaceColor', col, 'MarkerSize', 15);
            end
        end
        mid = ones(1,3)/3*tricc;
        plot(mid(1), mid(2), 'ko', 'MarkerFaceColor', mean(tricol,1), 'MarkerSize', 15);
    end
end

% Display example trajectory
% ~~~~~~~~~~~~~~~~~~~~~~~~~~
if plotexpseq
    
    % Make sure data from bad subjects is available
    rmbadsub = false;
    Emergence_FTA_LoadData;
    
    % Which trajectory to use
    sublist  = [28 27 28];
    condlist = [04 16 22];
    expidx = 3;
    sub = sublist(expidx); % subject number
    seq = condlist(expidx); % sequence number    
    
    % Display example dot
    expc = G{seq,sub}.BarycCoord * tricc;
    alpha(0);
    plot(expc(:,1), expc(:,2), 'k-');
    scatter(expc(:,1), expc(:,2), 50+ones(200,1), hsv(200), 'filled', 'MarkerEdgeColor', 'k');
    
    % Display finger position at the time of change point
    cp = G{seq,sub}.Jump + 1/2;
    plot(expc(cp,1), expc(cp,2), 'k.-', 'MarkerSize', 10);
end
