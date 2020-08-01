function [ new_y ] = lin_interp(old_x, old_y, new_x)
    %Input is old Y-values, which are assumed to correspond to
    %1:numel(old_y) if no old x-values are specified. The function them
    %interpolates to have new_len number of elements. X values must be
    %evenly spaced, start at 1, and end at the same value..
    new_len = numel(new_x);
    old_len = numel(old_y);
    if(new_len < 2 || old_len < 2)
        new_ys = old_ys; return;
    end
    if(numel(old_x) ~= old_len)
        error("Y and X are mismatched.");
    end
    if isempty(old_x)
        old_x = 1:numel(old_y);
    end
    new_y = zeros(numel(new_x), 1);
    old_x_spacing = old_x(2) - old_x(1);
    new_x_spacing = new_x(2) - new_x(1);
    x = 1;
    for i = 1:new_len
        bot_x_index = floor((x - 1)/old_x_spacing) + 1;
        top_x_index = min(ceil((x - 1)/old_x_spacing) + 1, old_len);
        y_diff = old_y(top_x_index) - old_y(bot_x_index);
        x_diff = x - old_x(bot_x_index);
        new_y(i) = old_y(bot_x_index) + x_diff*y_diff/old_x_spacing;
        x = x + new_x_spacing;
    end
    new_y = new_y';
end

