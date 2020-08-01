function [mean_r,median_r,stdev_r,min_r,max_r,r] = custom_min(data_1, data_2, data_3)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%Every time you process a dataset, store scores[](raw, unrounded) and fen
    %Process Data_1 - 600 seconds at 5000Hz sampling frequency
    data_1_scores = zeros(3e6,10);
    data_2_scores = zeros(1e6,10);
    data_3_scores = zeros(6511200,3);
    data_1_envelopes = data_1_scores;
    data_2_envelopes = data_2_scores;
    data_3_envelopes = data_3_scores;
    for j = 1:10     % take each column in  Data_1 ( each dataset )
        x = data_1(:,j);
        nx = randNoise(x,5000,60); % remove noise, and add random noise
        [scores, dnx, ~] = new_score(nx,5000,60,1,1,60,4); % take noised_x and use new_score to generate scores and dnx
        q = nx - dnx;
        [b, a] = butter(4, 2 * 100 / 5000, 'low');
        q = filtfilt(b,a,q);
        [crits, signs] = find_crits(q);
        crits = crits(signs == 1);             
        en = spline(crits(2:end-1), q(crits(2:end-1)), [1:length(dnx)]);
        data_1_scores(:,j) = scores;
        data_1_envelopes(:,j) = en;
%         plot(q, 'b'); hold on; plot(en, 'r'); hold off;
        % connect peaks only
        % lowpass filter below 65hz on noised_x - dnx
%         [b, a] = butter(4, 2 * 110 / 5000, 'low');
%         size(nx)
%         size(dnx)
%         fdnx = filtfilt(b, a, q);
%         plot(q, 'b'); hold on; plot(fdnx, 'r'); hold off;
    end
    % Same process on Data_2 - 200 seconds at 5000Hz sampling frequency
    for j = 1:10     % take each column in  Data_1 ( each dataset )
        x = data_2(:,j);
        nx = randNoise(x,5000,60); % remove noise, and add random noise
        [scores, dnx, ~] = new_score(nx,5000,60,1,1,60,4); % take noised_x and use new_score to generate scores and dnx
        q = nx - dnx;
        [b, a] = butter(4, 2 * 100 / 5000, 'low');
        q = filtfilt(b,a,q);
        [crits, signs] = find_crits(q);
        crits = crits(signs == 1);             
        en = spline(crits(2:end-1), q(crits(2:end-1)), [1:length(dnx)]);
        data_2_scores(:,j) = scores;
        data_2_envelopes(:,j) = en;
    end
    %Same (?) process on Data_3 - 200 seconds at 32556Hz sampling frequency
    for j = 1:3     % take each column in  Data_1 ( each dataset )
        x = data_3(:,j);
        nx = randNoise(x,32556,60); % remove noise, and add random noise
        [scores, dnx, ~] = new_score(nx,32556,60,1,1,60,4); % take noised_x and use new_score to generate scores and dnx
        q = nx - dnx;
        [b, a] = butter(4, 2 * 100 / 32556, 'low');
        q = filtfilt(b,a,q);
        [crits, signs] = find_crits(q);
        crits = crits(signs == 1);             
        en = spline(crits(2:end-1), q(crits(2:end-1)), [1:length(dnx)]);
        data_3_scores(:,j) = scores;
        data_3_envelopes(:,j) = en;
    end
    r = zeros(23,1);
    k = 1;
    for j = 1:10
        temp = corrcoef((data_1_scores(:,j)), data_1_envelopes(:,j));
        r(k) = temp(2,1);
        k = k + 1;
    end
    for j = 1:10
        temp = corrcoef(data_2_scores(:,j),data_2_envelopes(:,j));
        r(k) = temp(2,1);
        k = k + 1;
    end
    for j = 1:3
        temp = corrcoef(data_3_scores(:,j),data_3_envelopes(:,j));
        r(k) = temp(2,1);
        k = k + 1;
    end
    mean_r = mean(r);
    stdev_r = std(r);
    median_r = median(r);
    min_r = min(r);
    max_r = max(r);
end
