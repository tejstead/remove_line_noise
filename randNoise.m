

function [noised_signal, noise, r] = randNoise(signal, sampling_frequency, line_frequency, noise_ratio)

%TODO: remove the noise first?
    samps_per_cycle = sampling_frequency/line_frequency;
    MIN_SAMPS_PER_CYCLE = 60;
    if(samps_per_cycle > MIN_SAMPS_PER_CYCLE)
        interp_mode = "linear";
    else
        interp_mode = "spline";
    end
    sample_length = numel(signal);
    cycles_in_data = sample_length/samps_per_cycle;
    int_samps_per_cycle =  MIN_SAMPS_PER_CYCLE; 
    upsampled_sample_length = round(int_samps_per_cycle * cycles_in_data);
    upsampled_sampling_frequency = int_samps_per_cycle * line_frequency;
    upsample_spline_x = linspace(1, sample_length, upsampled_sample_length);    
    if(strcmp(interp_mode, "spline"))
        upsampled_raw_data = spline(1:sample_length, signal, upsample_spline_x);
    else %linear - faster, but inaccurate if original sampling frequency is low
        upsampled_raw_data = lin_interp(1:sample_length, signal, upsample_spline_x);
    end
%     [b,a] = butter(3, 2 * line_frequency * 8 / sampling_frequency, 'low');
%     upsampled_raw_data = filtfilt(b,a,upsampled_raw_data);
    upsampled_normalized_data = normalize_signal(upsampled_raw_data, upsampled_sampling_frequency);
    denoised_data = upsampled_normalized_data';
    %denoised_data = score_2020(upsampled_normalized_data', sampling_frequency, line_frequency, 0,0,0,1,0);
    temp = (1:upsampled_sample_length);
    noise = zeros(upsampled_sample_length, 1);
    r = rand(upsampled_sample_length * 2, 1);
    [b, a] = butter(2, 2 * 0.02 / sampling_frequency, 'low');
    r = filtfilt(b, a, r);
    len_r = numel(r);
    r = r(floor(len_r / 4):floor(len_r / 4) + upsampled_sample_length - 1);
    r = r - min(r);
    r = (r / (range(r))) * 1.2;
    r = r - 0.1;
    r(r > 1) = 1;
    r(r < 0) = 0;

    NOISE_CONSTANTS = [1 0.056 0.1322 0.01193 0.2005]; 
    for i = 1:numel(NOISE_CONSTANTS)
                noise = (noise) + (NOISE_CONSTANTS(i) .* sin(temp * i * (2 * pi * line_frequency / upsampled_sampling_frequency)))';
    end

    %noise = noise .* r;
    noise = noise * noise_ratio / median(abs(noise));
    upsampled_noised_signal = (denoised_data + noise);
    %figure(111); pwelch(upsampled_noised_signal, upsampled_sampling_frequency, [], upsampled_sampling_frequency, upsampled_sampling_frequency);
    if(strcmp(interp_mode, "spline"))
        noised_signal = spline(upsample_spline_x, upsampled_noised_signal, 1:sample_length)';
        noise = spline(upsample_spline_x, noise, 1:sample_length)';
    else %linear - faster, but inaccurate if original sampling frequency is low
        noised_signal = lin_interp(upsample_spline_x, upsampled_noised_signal, 1:sample_length)';
        noise = lin_interp(upsample_spline_x, noise, 1:sample_length)';
    end
end


