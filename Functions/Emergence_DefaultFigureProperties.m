% This script defines default options for the figures
% 
% Copyright (c) 2020 Maxime Maheu

set(groot, 'DefaultAxesLineWidth',              1);                     % axes
set(groot, 'DefaultBarLineWidth',               1);                     % bars
set(groot, 'DefaultColorbarLineWidth',          1);                     % colorbar
set(groot, 'DefaultLineLineWidth',              1);                     % line plots
set(groot, 'DefaultFigureColor',                ones(1,3));             % remove the uggly gray surrounding
set(groot, 'DefaultAxesColor',                  'None');                % transparent background
set(groot, 'DefaultTextColor',                  zeros(1,3));            % black text
set(groot, 'DefaultAxesLineWidth',              1);                     % axes
set(groot, 'DefaultBarLineWidth',               1);                     % bars
set(groot, 'DefaultColorbarLineWidth',          1);                     % colorbar
set(groot, 'DefaultLineLineWidth',              1);                     % line plots
set(groot, 'DefaultAxesFontName',               'Helvetica');           % Helvetica font
set(groot, 'DefaultTextFontName',               'Helvetica');           % Helvetica font
set(groot, 'DefaultAxesFontUnits',              'Points');              % unit
set(groot, 'DefaultAxesFontSize',               10);                    % default font-size
set(groot, 'DefaultTextVerticalAlignment',      'Middle');              % center the text
set(groot, 'DefaultTextHorizontalAlignment',    'Center');              % center the text
set(groot, 'DefaultAxesBox',                    'off');                 % no upper and right axes
set(groot, 'DefaultAxesLayer',                  'Top');                 % axes on top of the plot
set(groot, 'DefaultAxesTickDir',                'In');                  % ticks oriented outside the plot
set(groot, 'DefaultAxesTickLength',             repmat(0.01, [1 2]));   % length of ticks
set(groot, 'DefaultAxesTickDirMode',            'Auto');                % otherwise the previous commands do not work