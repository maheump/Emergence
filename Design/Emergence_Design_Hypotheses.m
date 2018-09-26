% This script displays the theoretical differences between the different
% conditions (fully-stochastic, probabilistic and deterministic ones) in
% terms of Shannon's entropy.
%
% Copyright (c) 2018 Maxime Maheu

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

%% DISPLAY THE DIFFERENT CONDITIONS ON THE ENTROPY FUNCTION
%  ========================================================

% Define the precision of the grid onto which the entropy function will be
% evaluated
p = 0:0.001:1;
entarr = arrayfun(@Emergence_IO_Entropy, p);

% Prepare a new window
figure('Position', [1 945 280 160]); lgd = NaN(1,3);

% Display a line at the maximum entropy level
plot(ones(1,2)./2, [0,1], '-', 'Color', ones(1,3)./3); hold('on');

% Display the entropy function
plot(p, entarr, '-', 'Color', tricol(1,:), 'LineWidth', 5);

% Overlap the fully-stochastic case (1/2)
[~,i] = min(abs(p - median(p)));
plot(p(i), entarr(i), 'ko', 'MarkerFaceColor', tricol(3,:), 'MarkerSize', 10);

% Overlap deterministic cases (0 and 1)
plot(0, 0, 'ko', 'MarkerFaceColor', tricol(2,:), 'MarkerSize', 10);
plot(1, 0, 'ko', 'MarkerFaceColor', tricol(2,:), 'MarkerSize', 10);

% Customize the axes
set(gca, 'XTick', 0:0.2:1, 'YTick', 0:0.2:1, 'Box', 'Off', 'Layer', 'Bottom');

% Add some text labels
xlabel('p(A)');
ylabel({'Shannon', 'entropy', '(bits)'}, 'Rotation', 0, 'HorizontalAlignment', 'Right');

%% DISPLAY THE ENTROPY CONTINUUM ON ONE DIMENSION
%  ==============================================

% Look at only half of the probability continuum this time
p = 1/2:0.001:1;
entarr = arrayfun(@Emergence_IO_Entropy, p);

% Prepare a new window
figure('Position', [282 1055 500 50]);

% Display a colormap indexed on the entropy levels
imagesc(p, 1, entarr);
colormap(cbrewer2('Blues')); caxis([0,1]);
axis('off');
