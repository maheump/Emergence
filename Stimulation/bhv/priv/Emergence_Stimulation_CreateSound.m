function S1 = Emergence_Stimulation_CreateSound( stim_freq, stim_rise, stim_dur, aud_freq )
% From El Karoui et al Cereb Cortex 2014.
% Each sound was 50-ms long and composed of 3 sinusoidal tones 
% (350, 700, and 1400 Hz, sound A; or 500, 1000, and 2000 Hz, sound B). 
% All tones were prepared with 7-ms rise and 7-ms fall times.

if nargin < 4
    dv = PsychPortAudio('GetDevices');
    aud_freq = dv(1).DefaultSampleRate;
    if nargin < 3
        stim_dur = 0.05; % ms
        if nargin < 2
            stim_rise = 0.007; % ms
            if nargin < 1
                stim_freq = [350, 700, 1400]; % Hz
            end
        end
    end
end

% make S1, add waning and waxing and clip range to [-1 1];
S1 = zeros(1, round(stim_dur*aud_freq));
for k = 1:3
    S1 = S1 + sin(2*pi*stim_freq(k)/aud_freq*(1:round(stim_dur*aud_freq)));
end
S1(1:round(stim_rise*aud_freq)) = linspace(0, 1, round(stim_rise*aud_freq)) .* S1(1:round(stim_rise*aud_freq));
S1(end-round(stim_rise*aud_freq)+1:end) = linspace(1, 0, round(stim_rise*aud_freq)) .* S1(end-round(stim_rise*aud_freq)+1:end);
S1 = S1 / max(abs(S1));

end