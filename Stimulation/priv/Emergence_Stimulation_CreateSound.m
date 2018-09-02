function stim = Emergence_Stimulation_CreateSound( freq, rise, dur, audfreq )
% EMERGENCE_STIMULATION_CREATESOUND creates auditory stimuli.
% 
% Same technical details as El Karoui et al., Cereb. Cortex, 2014.
% Each sound was 50-ms long and composed of 3 sinusoidal tones 
% (350, 700, and 1400 Hz, sound A; or 500, 1000, and 2000 Hz, sound B). 
% All tones were prepared with 7-ms rise and 7-ms fall times.
% 
% Copyright (c) 2018 Maxime Maheu

% Fill the inputs
if nargin < 4
    dv = PsychPortAudio('GetDevices');
    audfreq = dv(1).DefaultSampleRate;
    if nargin < 3
        dur = 0.050; % second
        if nargin < 2
            rise = 0.007; % second
            if nargin < 1
                freq = [350, 700, 1400]; % Hertz
            end
        end
    end
end

% Make the stimulus
stim = zeros(1, round(dur*audfreq));
for k = 1:3
    stim = stim + sin(2*pi*freq(k)/audfreq*(1:round(dur*audfreq)));
end

% Add waning and waxing 
stim(1:round(rise*audfreq))         = linspace(0, 1, round(rise*audfreq)) .* stim(1:round(rise*audfreq));
stim(end-round(rise*audfreq)+1:end) = linspace(1, 0, round(rise*audfreq)) .* stim(end-round(rise*audfreq)+1:end);

% Clip range to [-1 1]
stim = stim / max(abs(stim));

end