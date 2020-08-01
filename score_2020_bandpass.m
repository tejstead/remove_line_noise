function [desampled_template, denoised_data, noise_score] = score_2020_bandpass(raw_data, sampling_frequency, line_frequency, cycles_per_template, calc_score_flag, remove_noise_flag, threshold_score)

                    
NUM_HARMONICS_SCALING_CONSTANTS_60 = [0.95 0.9 0.9 0.95 1]; % constants that make score consistent given varying # of harmonics
NUM_HARMONICS_SCALING_CONSTANTS_50 = [0.97 0.92 0.92 0.95 1];
MIN_SAMPS_PER_CYCLE = 60; % 60 = LCM(1,2,3,4,5,6) - this has nothing to do with 60Hz
                                     
SCORE_SCALING_CONSTANT = 1; % constant that will be adjusted for # harmonics and # bandpass poles

% set default values
if(sampling_frequency <= 0)
    error("Sampling frequency set to %f, invalid value.", sampling_frequency);
end
if(line_frequency == 0 || isempty(line_frequency))
    fprintf("Line frequency == 0, returning data unchanged.");
    denoised_data = raw_data; noise_score = 0; 
    return;
end
if(cycles_per_template == 0 || isempty(cycles_per_template))
    cycles_per_template = 2 * line_frequency; % 2 seconds 
end

sample_length = numel(raw_data);
samps_per_cycle = sampling_frequency/line_frequency;

% spline is used for upsampling, linear for downsampling
if(samps_per_cycle > MIN_SAMPS_PER_CYCLE)
    interp_mode = "downsample";
else
    interp_mode = "upsample";
end

% Restrict number of harmonics to [1,5], but subtract 2 before doing so because of bandpass filter's high bound
num_harmonics = max(1, min(floor(sampling_frequency/(5*line_frequency)) - 2, 5));

% adjust filter poles
if(line_frequency == 60) 
    if(sampling_frequency > 27000)
        bandpass_poles = 3;
        SCORE_SCALING_CONSTANT = 0.96;
    elseif(sampling_frequency > 10000)
        bandpass_poles = 4;
        SCORE_SCALING_CONSTANT = 0.97;
    else
        bandpass_poles = 5;
        SCORE_SCALING_CONSTANT = 0.99 * NUM_HARMONICS_SCALING_CONSTANTS_60(num_harmonics);
    end
elseif(line_frequency == 50)
    if(sampling_frequency > 22000)
        bandpass_poles = 3;
        SCORE_SCALING_CONSTANT = 1.05;
    elseif(sampling_frequency > 8000)
        bandpass_poles = 4;
        SCORE_SCALING_CONSTANT = 1.05;
    else
        bandpass_poles = 5;
        SCORE_SCALING_CONSTANT = 1.07 * NUM_HARMONICS_SCALING_CONSTANTS_50(num_harmonics);
    end
else % non-standard line frequency - this shouldn't be used
    bandpass_poles = 3;
end

int_samps_per_cycle = MIN_SAMPS_PER_CYCLE; % force 60 samples per cycle
cycles_in_data = sample_length/samps_per_cycle;
cycles_in_grid = floor(cycles_in_data);
resampled_sample_length = round(int_samps_per_cycle * cycles_in_data);


% bandpass filter data
low_bound = line_frequency/2;
high_bound = line_frequency * (num_harmonics + 2);
[b, a] = butter(bandpass_poles, 2 * [low_bound high_bound]/sampling_frequency, 'bandpass');
bandpassed_data = filtfilt(b, a, raw_data);
bandpassed_data = bandpassed_data - mean(bandpassed_data);

% resample data
resample_spline_x = linspace(1, sample_length, resampled_sample_length);
if(strcmp(interp_mode, "upsample"))
    resampled_bandpassed_data = spline(1:sample_length, bandpassed_data, resample_spline_x);
else %linear - faster, but inaccurate if original sampling frequency is low
    resampled_bandpassed_data = lin_interp(1:sample_length, bandpassed_data, resample_spline_x);
end

% put bandpassed data into median grid
median_grid = zeros(cycles_in_grid, int_samps_per_cycle);
filtered_data_index = 1;
for i = 1:cycles_in_grid
    median_grid(i, :) = resampled_bandpassed_data(filtered_data_index:(filtered_data_index + int_samps_per_cycle - 1));
    filtered_data_index = filtered_data_index + int_samps_per_cycle;
end

% median filter median grid
for i = 1:int_samps_per_cycle
    median_grid(:, i) = quantfilt1(median_grid(:, i), 0.5, cycles_per_template, 'extrapolate');
end


noise_score = 0;
median_grid_temporal = zeros(resampled_sample_length, 1);
num_remaining_points = resampled_sample_length - (cycles_in_grid * int_samps_per_cycle);
filtered_data_index = 1;
for i = 1:cycles_in_grid % reshape median grid to be in time order
    median_grid_temporal(filtered_data_index:(filtered_data_index + int_samps_per_cycle - 1)) = median_grid(i, :);
    filtered_data_index = filtered_data_index + int_samps_per_cycle;
end

% handle leftover points by copying a portion of the last cycle
median_grid_temporal(filtered_data_index:resampled_sample_length) = median_grid_temporal((filtered_data_index - int_samps_per_cycle): (filtered_data_index - int_samps_per_cycle + num_remaining_points - 1));



if(strcmp(interp_mode, "downsample")) %interp_modes are reversed here
    desampled_template = spline(resample_spline_x, median_grid_temporal, 1:sample_length)/1.014;
else
    desampled_template = lin_interp(resample_spline_x, median_grid_temporal, 1:sample_length)/1.014;
end
% 
% for i = 1:(sampling_frequency):sample_length %1-second blocks
%     start_idx = i;
%     end_idx = i + (sampling_frequency) - 1;
%     
%     % amplitude correction using least absolute deviations regression line
%     [b, m] = lad_reg2(desampled_template(start_idx:end_idx), bandpassed_data(start_idx:end_idx));
%     desampled_template(start_idx:end_idx) = (desampled_template(start_idx:end_idx) * m) + b;
% end

%calculate score
if(calc_score_flag)
    scores = zeros(cycles_in_grid, 1);
    filtered_data_index = 1;
    SCALING_CONSTANT = 280 * SCORE_SCALING_CONSTANT;
    for i = 1:cycles_in_grid 
        start_idx = floor(filtered_data_index);
        end_idx = floor(filtered_data_index + samps_per_cycle - 1);
        baseline_waveform = bandpassed_data(start_idx:end_idx);
        baseline_waveform = baseline_waveform - desampled_template(start_idx:end_idx)';
        waveform_size = median(abs(baseline_waveform));
        template_size = median(abs(desampled_template(start_idx:end_idx)));
        size_ratio = (template_size)/(waveform_size);
        scores(i) = round(size_ratio * SCALING_CONSTANT);
        filtered_data_index = filtered_data_index + samps_per_cycle;
    end
        noise_score = max(min(round(median(scores)) - 30, 294),0); %restricts score to range [0,254]
else
    if(threshold_score ~= 0)
        fprintf("Threshold score is not used, since calc_score_flag = 0.");
    end
end

figure(5); plot(bandpassed_data); hold on; 
passes_threshold = (calc_score_flag == 0) || (noise_score >= threshold_score); 

if(remove_noise_flag && passes_threshold)
    denoised_data = raw_data - desampled_template';
else 
    denoised_data = raw_data;
end
