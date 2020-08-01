function [sampling_frequencies, linear_scores, spline_scores] = determine_cutoff(raw_data,orig_sampling_frequency,line_frequency,cycles_per_template,num_harmonics,calc_score_flag,remove_noise_flag,threshold_score)
    sampling_frequencies = 400:100:20000;
    linear_scores = zeros(numel(sampling_frequencies), 1);
    spline_scores = zeros(numel(sampling_frequencies), 1);
    for i = 1:numel(sampling_frequencies)
        base_x = resample(raw_data, sampling_frequencies(i), orig_sampling_frequency);
        [~, linear_scores(i), ~] = score_2020(base_x, sampling_frequencies(i), line_frequency,cycles_per_template,num_harmonics,calc_score_flag,0,threshold_score);
        %[~, spline_scores(i), ~] = score_2020(base_x, sampling_frequencies(i), line_frequency,cycles_per_template,num_harmonics,calc_score_flag,0,threshold_score);
    end
end