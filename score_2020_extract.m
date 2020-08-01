function [desampled_template] = score_2020_template(raw_data, sampling_frequency, line_frequency, cycles_per_template, calc_score_flag, remove_noise_flag, threshold_score)
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
                    %based on empirical data, used for predicted_score -
                    %FIX
                    
    NUM_HARMONICS_SCALING_CONSTANTS = [0.95 0.9 0.9 0.95 1];
    MIN_SAMPS_PER_CYCLE = 60; %60 = LCM(1,2,3,4,5,6) - this has nothing to do with 60Hz
                                     
SCORE_SCALING_CONSTANT = 1;
if(sampling_frequency <= 0)
    error("Sampling frequency set to %f, invalid value.", sampling_frequency);
end
if(line_frequency == 0 || isempty(line_frequency))
    fprintf("Line frequency == 0, returning data unchanged.");
    denoised_data = raw_data; noise_score = 0; predicted_score = 0;
    return;
end

if(cycles_per_template == 0 || isempty(cycles_per_template))
    cycles_per_template = 120;
end
sample_length = numel(raw_data);
samps_per_cycle = sampling_frequency/line_frequency;
if(samps_per_cycle > MIN_SAMPS_PER_CYCLE)
    interp_mode = "downsample";
else
    interp_mode = "upsample";
end

if(line_frequency == 60)
    if(sampling_frequency > 27000)
        bandpass_poles = 3;
        SCORE_SCALING_CONSTANT = SCORE_SCALING_CONSTANT * 0.99;
    elseif(sampling_frequency > 10000)
        bandpass_poles = 4;
    else
        bandpass_poles = 5;
        SCORE_SCALING_CONSTANT = SCORE_SCALING_CONSTANT * 1.01;
    end
elseif(line_frequency == 50)
    bandpass_poles = 3; %fix
else
    bandpass_poles = 3;
end

int_samps_per_cycle = MIN_SAMPS_PER_CYCLE; %force 60
cycles_in_data = sample_length/samps_per_cycle;
cycles_in_grid = floor(cycles_in_data);
% resampled_sampling_frequency = sampling_frequency * int_samps_per_cycle / samps_per_cycle;
resampled_sample_length = round(int_samps_per_cycle * cycles_in_data);
num_harmonics = min(floor(sampling_frequency/(5*line_frequency)), 5);
SCORE_SCALING_CONSTANT = SCORE_SCALING_CONSTANT * NUM_HARMONICS_SCALING_CONSTANTS(num_harmonics);
low_bound = line_frequency/2;
high_bound = line_frequency * (num_harmonics + 1);
[b, a] = butter(bandpass_poles, 2 * [low_bound high_bound]/sampling_frequency, 'bandpass');
bandpassed_data = filtfilt(b, a, raw_data);
bandpassed_data = bandpassed_data - mean(bandpassed_data);

resample_spline_x = linspace(1, sample_length, resampled_sample_length);
% figure(800); plot(bandpassed_data); hold on;
%fprintf("Reached Line 75");
if(strcmp(interp_mode, "upsample"))
    resampled_bandpassed_data = spline(1:sample_length, bandpassed_data, resample_spline_x);
    q1 = lin_interp(resample_spline_x, resampled_bandpassed_data, 1:sample_length);
else %linear - faster, but inaccurate if original sampling frequency is low
    resampled_bandpassed_data = lin_interp(1:sample_length, bandpassed_data, resample_spline_x);
    q2 = spline(resample_spline_x, resampled_bandpassed_data, 1:sample_length);
end

% all good to here (q1 & q2 are garbage)

median_grid = zeros(cycles_in_grid, int_samps_per_cycle);
%put bandpassed data into median grid
filtered_data_index = 1;
for i = 1:cycles_in_grid
    median_grid(i, :) = resampled_bandpassed_data(filtered_data_index:(filtered_data_index + int_samps_per_cycle - 1));
    filtered_data_index = filtered_data_index + int_samps_per_cycle;
end
%fprintf("Reached Line 93");
%median filter median grid
for i = 1:int_samps_per_cycle
    median_grid(:, i) = quantfilt1(median_grid(:, i), 0.5, cycles_per_template, 'extrapolate');
end

if(0) %change to 1 for debugging - this will show the model of noise generated
    median_grid_temporal = zeros(resampled_sample_length, 1);
    filtered_data_index = 1;
    for i = 1:cycles_in_grid
        median_grid_temporal(filtered_data_index:(filtered_data_index + int_samps_per_cycle - 1)) = median_grid(i, :);
        filtered_data_index = filtered_data_index + int_samps_per_cycle;
    end
    figure(996); plot(resample_spline_x, median_grid_temporal);
end
%fprintf("Reached Line 108");
noise_score = 0;
median_grid_temporal = zeros(resampled_sample_length, 1);
num_remaining_points = resampled_sample_length - (cycles_in_grid * int_samps_per_cycle);
filtered_data_index = 1;
for i = 1:cycles_in_grid % reshape()
    median_grid_temporal(filtered_data_index:(filtered_data_index + int_samps_per_cycle - 1)) = median_grid(i, :);
    filtered_data_index = filtered_data_index + int_samps_per_cycle;
end
%fprintf("Reached Line 117");
 median_grid_temporal(filtered_data_index:resampled_sample_length) = median_grid_temporal((filtered_data_index - int_samps_per_cycle): (filtered_data_index - int_samps_per_cycle + num_remaining_points - 1));


% figure(2);
% plot(resample_spline_x, median_grid_temporal, 'r'); hold on; plot(resample_spline_x, resampled_bandpassed_data, 'k'); hold off;
% pause;
%downsample_spline_x = linspace(1, upsampled_sample_length, sample_length);

if(strcmp(interp_mode, "downsample"))
    desampled_template = spline(resample_spline_x, median_grid_temporal, 1:sample_length);
else
    desampled_template = lin_interp(resample_spline_x, median_grid_temporal, 1:sample_length);
end
%  figure(800); hold on; plot(desampled_template, 'r'); plot(bandpassed_data, 'b'); plot(noise, 'k'); hold off;
% median(abs(bandpassed_data))
% median(abs(bandpassed_data - desampled_template'))
% median(abs(desampled_template))
 figure(800); plot(1:sample_length, desampled_template);
%fprintf("Reached Line 125");
if(calc_score_flag)
    scores = zeros(cycles_in_grid, 1);
    waveform_sizes = zeros(cycles_in_grid, 1);
    template_sizes = zeros(cycles_in_grid, 1);
    filtered_data_index = 1;
    SCALING_CONSTANT = 200 * SCORE_SCALING_CONSTANT;
    for i = 1:cycles_in_grid 
        start_idx = floor(filtered_data_index);
        end_idx = floor(filtered_data_index + samps_per_cycle - 1);
        baseline_waveform = bandpassed_data(start_idx:end_idx);
        baseline_waveform = baseline_waveform - desampled_template(start_idx:end_idx)';
        waveform_size = median(abs(baseline_waveform));
        template_size = median(abs(desampled_template(start_idx:end_idx)));
        waveform_sizes(i) = waveform_size;
        template_sizes(i) = template_size;
%         size_ratio = template_size^2/(waveform_size + log(template_size + 1));
        size_ratio = (template_size)/(waveform_size);
        scores(i) = round(size_ratio * SCALING_CONSTANT);
        filtered_data_index = filtered_data_index + samps_per_cycle;
    end
%     figure(999); plot(waveform_sizes , 'k'); hold on; plot(template_sizes, 'r'); hold off;
        noise_score = max(min(round(median(scores)), 254),0);
        predicted_score = SCORE_LOOKUP_TABLE(noise_score + 1); %add 1 because MATLAB is 1-indexed, remove in C %FIX
else
    if(threshold_score ~= 0)
        fprintf("Threshold score is not used, since calc_score_flag = 0.");
    end
end
%fprintf("Reached Line 150");
passes_threshold = (calc_score_flag == 0) || (noise_score > threshold_score);

denoised_data = raw_data;

if(remove_noise_flag && passes_threshold)
    %downsample the model of noise
    denoised_data = raw_data - desampled_template';
end
