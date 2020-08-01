function [denoised_data, noise_score, predicted_score] = score_2020(raw_data, sampling_frequency, line_frequency, cycles_per_template, calc_score_flag, remove_noise_flag, threshold_score)
    SCORE_LOOKUP_TABLE = [ 00 01 02 02 02 02 02 02 02 02 02 03 03 03 03 03 ... 
                        03 03 03 03 03 03 03 03 03 03 03 03 03 03 03 03 ...
                        03 03 03 03 03 03 03 03 03 03 03 04 04 04 04 04 ...
                        04 04 04 04 04 04 04 04 04 04 04 04 04 04 04 04 ...
                        04 04 05 05 05 05 05 05 05 05 05 05 05 05 05 05 ...
                        05 05 05 06 06 06 06 06 06 06 06 06 06 06 06 06 ...
                        06 07 07 07 07 07 07 07 07 07 07 07 07 08 08 08 ...
                        08 08 08 08 08 08 08 08 08 09 09 09 09 09 09 09 ...
                        09 09 09 09 10 10 10 10 10 10 10 10 10 10 10 10 ...
                        10 10 11 11 11 11 11 11 11 11 11 11 11 11 12 12 ...
                        12 12 12 12 12 12 12 12 12 12 12 13 13 13 13 13 ...
                        13 13 13 13 13 13 13 13 14 14 14 14 14 14 14 14 ...
                        14 14 14 15 15 15 15 15 15 15 15 15 15 15 15 15 ...
                        16 16 16 16 16 16 16 16 16 16 16 16 16 17 17 17 ...
                        17 17 17 17 17 17 17 17 17 17 17 17 17 18 28 18 ...
                        18 18 18 18 18 18 18 18 18 18 18 19 19 19 19 20]; 
                    %based on empirical data, used for predicted_score
                    
    BANDPASS_CUTOFFS_60 = [];
    MIN_SAMPS_PER_CYCLE = 60; %60 = LCM(1,2,3,4,5,6) - this has nothing to do with 60Hz
                                     
SCORE_SCALING_CONSTANT = 1;
if(sampling_frequency <= 0)
    error("Sampling frequency set to %f, invalid value.", sampling_frequency);
end
if(line_frequency == 0)
    fprintf("Line frequency == 0, returning data unchanged.");
    denoised_data = raw_data; noise_score = 0; predicted_score = 0;
    return;
end

if(cycles_per_template == 0)
    cycles_per_template = 60;
end
sample_length = numel(raw_data);
samps_per_cycle = sampling_frequency/line_frequency;
if(samps_per_cycle > MIN_SAMPS_PER_CYCLE)
    interp_mode = "spline";
else
    interp_mode = "spline";
end

int_samps_per_cycle =  MIN_SAMPS_PER_CYCLE; %resample to 60
cycles_in_data = sample_length/samps_per_cycle;
cycles_in_grid = floor(cycles_in_data);
upsampled_sampling_frequency = sampling_frequency * int_samps_per_cycle / samps_per_cycle;
upsampled_sample_length = round(int_samps_per_cycle * cycles_in_data);
nyquist = sampling_frequency / 2; %Nyquist limit is half the sampling frequency

num_harmonics = min(floor(sampling_frequency/(2*line_frequency)) - 1, 5);
low_bound = line_frequency - 10;
high_bound = line_frequency * num_harmonics + 30;
if(high_bound > nyquist)
    high_bound = nyquist - 1;
end
reasonable = sampling_frequency/5;
if high_bound > reasonable
   fprintf("Warning: Sampling frequency is %f and filter upper bound is %f, consider a lower number of harmonics.\n", sampling_frequency, high_bound);
end
    highpass_poles = 3;


upsample_spline_x = linspace(1, sample_length, upsampled_sample_length);
if(strcmp(interp_mode, "spline"))
    upsampled_raw_data = spline(1:sample_length, raw_data, upsample_spline_x);
else %linear - faster, but inaccurate if original sampling frequency is low
    upsampled_raw_data = lin_interp(1:sample_length, raw_data, upsample_spline_x);
end
[b, a] = butter(highpass_poles, 2/upsampled_sampling_frequency, 'high');
upsampled_highpassed_data = filtfilt(b, a, upsampled_raw_data);
if(line_frequency == 60)
    if(upsampled_sampling_frequency > 30000)
        bandpass_poles = 3;
        SCORE_SCALING_CONSTANT = SCORE_SCALING_CONSTANT * 13/12;
    elseif(upsampled_sampling_frequency > 14000)
        bandpass_poles = 4;
    else
        bandpass_poles = 5;
        SCORE_SCALING_CONSTANT = SCORE_SCALING_CONSTANT * 11/12;
    end
elseif(line_frequency == 50)
    bandpass_poles = 3; %fix
else
    bandpass_poles = 3;
end

if(high_bound < upsampled_sampling_frequency / 2)
    [b, a] = butter(bandpass_poles, 2 * [low_bound high_bound]/upsampled_sampling_frequency, 'bandpass');
else
    [b, a] = butter(highpass_poles, 2 * low_bound/upsampled_sampling_frequency, 'high');
end

upsampled_bandpassed_data = filtfilt(b, a, upsampled_highpassed_data);


median_grid = zeros(cycles_in_grid, int_samps_per_cycle);

%put bandpassed data into median grid
filtered_data_index = 1;
for i = 1:cycles_in_grid
    median_grid(i, :) = upsampled_bandpassed_data(filtered_data_index:(filtered_data_index + int_samps_per_cycle - 1));
    filtered_data_index = filtered_data_index + int_samps_per_cycle;
end

%median filter median grid
for i = 1:int_samps_per_cycle
    median_grid(:, i) = quantfilt1(median_grid(:, i), 0.5, cycles_per_template, 'extrapolate');
end

if(0) %change to 1 for debugging - this will show the model of noise generated
    median_grid_temporal = zeros(upsampled_sample_length, 1);
    filtered_data_index = 1;
    for i = 1:cycles_in_grid
        median_grid_temporal(filtered_data_index:(filtered_data_index + int_samps_per_cycle - 1)) = median_grid(i, :);
        filtered_data_index = filtered_data_index + int_samps_per_cycle;
    end
    figure(999); plot(upsample_spline_x, median_grid_temporal);
end

noise_score = 0;
if(calc_score_flag)
    scores = zeros(cycles_in_grid, 1);
    filtered_data_index = 1;
    SCALING_CONSTANT = 600 * SCORE_SCALING_CONSTANT;
    for i = 1:cycles_in_grid
        median_of_cycle = median(upsampled_highpassed_data(filtered_data_index:(filtered_data_index + int_samps_per_cycle - 1)));
        baseline_waveform = abs(upsampled_highpassed_data(filtered_data_index:(filtered_data_index + int_samps_per_cycle - 1)) - median_of_cycle);
        waveform_size = median(baseline_waveform);
        template_size = median(abs(median_grid(i, :)));
        size_ratio = template_size^2/(waveform_size + log(template_size + 1));
        scores(i) = round(size_ratio * SCALING_CONSTANT);
        filtered_data_index = filtered_data_index + int_samps_per_cycle;
    end
        noise_score = min(round(median(scores)), 255);
        predicted_score = SCORE_LOOKUP_TABLE(noise_score + 1); %add 1 because MATLAB is 1-indexed, remove in C
else
    if(threshold_score ~= 0)
        fprintf("Threshold score is not used, since calc_score_flag = 0.");
    end
end
passes_threshold = (calc_score_flag == 0) || (noise_score > threshold_score);

denoised_data = raw_data;

if(remove_noise_flag && passes_threshold)
    %downsample the model of noise
    median_grid_temporal = zeros(upsampled_sample_length, 1);
    num_remaining_points = upsampled_sample_length - (cycles_in_grid * int_samps_per_cycle);
    filtered_data_index = 1;
    for i = 1:cycles_in_grid
        median_grid_temporal(filtered_data_index:(filtered_data_index + int_samps_per_cycle - 1)) = median_grid(i, :);
        filtered_data_index = filtered_data_index + int_samps_per_cycle;
    end
    median_grid_temporal(filtered_data_index:upsampled_sample_length) = median_grid_temporal((filtered_data_index - int_samps_per_cycle): (filtered_data_index - int_samps_per_cycle + num_remaining_points - 1));
    %downsample_spline_x = linspace(1, upsampled_sample_length, sample_length);
    if(strcmp(interp_mode, "spline"))
        downsampled_template = spline(upsample_spline_x, median_grid_temporal, 1:sample_length);
    else
        downsampled_template = lin_interp(upsample_spline_x, median_grid_temporal, 1:sample_length);
    end
    figure(994); plot(1:sample_length, downsampled_template, 'k');
    %plot(raw_data);
    denoised_data = raw_data - downsampled_template';
end
