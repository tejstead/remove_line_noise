    function [] = empiric_constant_finder( data_1, data_2, data_3 )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

data_1 = round(data_1);
data_2 = round(data_2);
data_3 = round(data_3);'
data_1_templates = zeros(size(data_1));
data_2_templates = zeros(size(data_2));
data_3_templates = zeros(size(data_3));

for i = 1:10
    current_data = data_1(:,i);
    [~,~,~,~,mg] = Copy_of_new_score(current_data,5000,60,0,);
    
end
end

