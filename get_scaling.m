function [ fit1 fit2 ] = get_scaling(nx, n, sampling_freq)

    for i = 1:10
        disp(i);
        [~, denoised_x, ~, ~] = new_score(nx(1:i:length(nx)), sampling_freq / i, 60, 1, 1, 60, 5);
        f = robustfit(n(1:i:length(nx)), nx(1:i:length(nx)) - denoised_x');
        fit1(i) = f(1);
        fit2(i) = f(2);
    end

end

