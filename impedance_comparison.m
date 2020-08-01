function [impedance_models, impedance_scores, impedance_models_bandpass, impedance_scores_bandpass, real_amplitude] = impedance_comparison(raw_data, sampling_frequency, impedance_frequency, impedance_amplitude, cycles_per_subsample, number_of_subsamples, cycles_per_template)
    %seconds_in_sample = numel(raw_data)/sampling_frequency;
    if(cycles_per_subsample == 0)
        cycles_per_subsample = 100;
    end
    length_per_subsample = round(cycles_per_subsample * sampling_frequency/impedance_frequency);
    time_window_starts = floor(0:number_of_subsamples * length_per_subsample);
    impedance_scores_bandpass = zeros(number_of_subsamples, 1);
    impedance_scores = zeros(number_of_subsamples, 1);
    impedance_models_bandpass = zeros(length_per_subsample, number_of_subsamples);
    impedance_models = zeros(length_per_subsample, number_of_subsamples);
    impedance_data = zeros(length_per_subsample, number_of_subsamples);
    

    temp = 1:numel(raw_data);
    num_poles = 4;
    if(impedance_frequency == 20)
        num_poles = 3;
    end
    [b, a] = butter(num_poles, 2* [min(impedance_frequency - 5, impedance_frequency * 0.95), max(impedance_frequency * 1.05, impedance_frequency + 5)]/sampling_frequency, 'bandpass');
    unnoised_filtered_data = filtfilt(b, a, raw_data);
    if(impedance_amplitude == 0)
        impedance_data(:, 3) = unnoised_filtered_data(time_window_starts(3) + 1:time_window_starts(3) + length_per_subsample);
        median_of_data = median(impedance_data(:, 3));
        base_amplitude = median(abs(impedance_data(:, 3) - median_of_data));
        real_amplitude = base_amplitude * 10;
    else
        real_amplitude = impedance_amplitude;
    end
    noised_data = raw_data + real_amplitude * sin(2*pi*temp* impedance_frequency/sampling_frequency)';
    filtered_data = filtfilt(b, a, noised_data);
    filtered_data = filtered_data - mean(filtered_data);
    for i = 1:number_of_subsamples
        impedance_data(:, i) = noised_data(time_window_starts(i) + 1: time_window_starts(i) + length_per_subsample);
        impedance_models_bandpass(:, i) = filtered_data(time_window_starts(i) + 1: time_window_starts(i) + length_per_subsample);
        impedance_bandpass_median = median(impedance_models_bandpass(:, i));
        impedance_scores_bandpass(i) = median(abs(impedance_models_bandpass(:, i) - impedance_bandpass_median));
        [impedance_scores(i), impedance_models(:, i)] = impedance_modeler(impedance_data(:, i), sampling_frequency, impedance_frequency, cycles_per_template, 1, 1);
    end
end