function Emergence_Stimulation_MouseCoordCallback( timerObj, event, trial, txtfile ) 
% EMERGENCE_STIMULATION_MOUSECOORDCALLBACK tracks finger position (mouse
% position actually) in the background in order to reach high sampling
% rate.
% 
% Copyright (c) 2020 Maxime Maheu

% Read the coordinates
coords = get(0, 'PointerLocation');

% Get relevant information
t = GetSecs;
x = coords(1);

% In PTB, the origin (0,0) is the top left
% In MATLAB, the origin (0,0) is the bottom left
y = coords(2); % abs(h_px - coords(2));

% Write in a text file
fid = fopen(txtfile, 'a');
fprintf(fid, '%i,%i,%i,%i\n', trial, t, x, y);
fclose(fid);

end