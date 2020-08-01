function [sampling_frequencies, score_at_max_noise] = max_score_by_sampling_frequency(raw_data,orig_sampling_frequency,line_frequency,cycles_per_template,num_harmonics,calc_score_flag,remove_noise_flag,threshold_score)
    sampling_frequencies = 400:100:32600;
    score_at_max_noise = zeros(numel(sampling_frequencies), 1);
    for i = 1:numel(sampling_frequencies)
        base_x = 
        if(sampling_frequencies(i) > 15000)
        [~, score_at_max_noise(i), ~] = score_2020(base_x, sampling_frequencies(i), line_frequency,cycles_per_template,num_harmonics,calc_score_flag,remove_noise_flag,threshold_score)
    end
end