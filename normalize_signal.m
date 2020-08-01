
function [normed_signal] = normalize_signal(data, sampling_frequency)
line_frequency = 60;
if(line_frequency == 60)
    if(sampling_frequency > 27000)
        bandpass_poles = 3;

    elseif(sampling_frequency > 10000)
        bandpass_poles = 4;
    else
        bandpass_poles = 5;
    end
elseif(line_frequency == 50)
    bandpass_poles = 3; %fix
else
    bandpass_poles = 3;
end
num_harmonics = min(floor(sampling_frequency/(5*line_frequency)), 5);
low_bound = line_frequency/2;
high_bound = line_frequency * (num_harmonics + 1);
% [b, a] = butter(3, 2/sampling_frequency, 'high');
 [b, a] = butter(bandpass_poles, 2 * [low_bound high_bound]/sampling_frequency, 'bandpass');
normed_signal = filtfilt(b, a, data);
amp_signal = median(abs(normed_signal));
normed_signal = data ./ amp_signal;
end
