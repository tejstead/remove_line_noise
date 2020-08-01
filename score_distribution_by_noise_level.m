function [noise_levels, noise_scores] = score_distribution_by_noise_level(raw_data,sampling_frequency,line_frequency,cycles_per_template)
    noise_levels = 0:254;
    noise_scores = zeros(255,1);
%     denoised_noise_scores = zeros(255, 1);
    for i = 0:254
        i
        [noised_data] = randNoise(raw_data, sampling_frequency, line_frequency, i/254);
        [~, noise_scores(i + 1)] = score_2020_bandpass(noised_data,sampling_frequency,line_frequency,cycles_per_template,1,0,0);
%         [~, denoised_noise_scores(i + 1)] = score_2020_bandpass(dnx,sampling_frequency,line_frequency,cycles_per_template,1,0,0);
    end 
end
