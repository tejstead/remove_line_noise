function [dists_400, dists_5000, dists_32556] = score_algo_tester(data_400, data_5000, data_32556)
    dists_400 = zeros(255, 1);
    dists_5000 = zeros(255, 10);
    dists_32556 = zeros(255, 4);
    for i = 1:1 %handle data_400
        i
        [~, dists_400(:, i)] = score_distribution_by_noise_level(data_400(:, i), 400,60,0);
    end
    
    for i = 1:10 %handle data_5000
        i
        [~, dists_5000(:, i)] = score_distribution_by_noise_level(data_5000(:, i), 5000,60,0);
    end
    
    for i = 1:4 %handle data_32556
        i
        [~, dists_32556(:, i)] = score_distribution_by_noise_level(data_32556(:, i), 32556,60,0);
    end
end