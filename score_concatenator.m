function [ old_scores, new_scores ] = score_concatenator( data_1, data_2, data_3 )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%I don't know if I spelled that right

old_scores = zeros(3e6*60*10/5000 + 1e6*60*10/5000 + 6511200*60*3/32556,1);
new_scores = zeros(3e6*60*10/5000 + 1e6*60*10/5000 + 6511200*60*3/32556,1);
score_counter = 1;
for j = 1:10     % take each column in  Data_1 ( each dataset )
        x = data_1(:,j);
        nx = randNoise(x,5000,60); % remove noise, and add random noise
        [scores, dnx, ~] = new_score(nx,5000,60,1,1,60,4); % take noised_x and use new_score to generate scores and dnx
        old_scores(score_counter:score_counter+3e6*60/5000 - 1) = scores;
        [scores,~,~] = new_score(dnx,5000,60,1,1,60,4);
        new_scores(score_counter:score_counter+3e6*60/5000 - 1) = scores;
        score_counter = score_counter + 3e6*60/5000;
end
for j = 1:10     % take each column in  Data_2 ( each dataset )
        x = data_2(:,j);
        nx = randNoise(x,5000,60); % remove noise, and add  noised_x and use new_score to generate scores and dnx
        [scores, dnx, ~] = new_score(nx,5000,60,1,1,60,4); % take noised_x and use new_score to generate scores and dnx
        old_scores(score_counter:score_counter+1e6*60/5000 - 1) = scores;
        [scores,~,~] = new_score(dnx,5000,60,1,1,60,4);
        new_scores(score_counter:score_counter+1e6*60/5000 - 1) = scores;
        score_counter = score_counter + 1e6*60/5000;
end
for j = 1:3
        x = data_3(:,j);
        nx = randNoise(x,32556,60); % remove noise, and add random noise
        [scores, dnx, ~] = new_score(nx,32556,60,1,1,60,4); % take noised_x and use new_score to generate scores and dnx
        old_scores(score_counter:score_counter+6511200*60/32556 - 1) = scores;
        [scores,~,~] = new_score(dnx,32556,60,1,1,60,4);
        new_scores(score_counter:score_counter+6511200*60/32556 - 1) = scores;
        score_counter = score_counter + 6511200*60/32556;
end

end

