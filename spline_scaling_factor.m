function [sampling_frequencies, scaling_factors] = spline_scaling_factor()

    sampling_frequencies = 4000:1000:50000;
    scaling_factors = zeros(numel(sampling_frequencies), 1);
    for j = 1:numel(sampling_frequencies)
        sampling_frequency = sampling_frequencies(j)
        sca_fac = zeros(10, 1);
        for k = 1:10
            white_noise = rand(sampling_frequency * 200, 1);
            [b,a] = butter(3, 2/sampling_frequency, 'high');
            white_noise = fft(white_noise);
            white_noise = white_noise .* [1./(1:(sampling_frequency*100)) 1./((sampling_frequency*100):-1:1)]';
            white_noise = ifft(white_noise);
            white_noise = abs(white_noise);
            white_noise = filtfilt(b,a,white_noise);
            [raw_data, normed_noise] = randNoise(white_noise, sampling_frequency, 60, 1);
            desampled_template = score_2020_bandpass(raw_data,sampling_frequency,60,0,1,1,0);
            sca_fac(k) = median(abs(desampled_template))/median(abs(normed_noise));
        end
        scaling_factors(j) = median(sca_fac);
    end
end