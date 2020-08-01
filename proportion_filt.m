function px = proportion_filt(x, proportion, span)

    % set up
    x = x(:);
    len = length(x);
    px = zeros(size(x));
    wind_mat = zeros(span + 3, 3);
    VALUE_COL = 1;
    PREV_IDX_COL = 2;
    NEXT_IDX_COL = 3;
    IN_IDX_COL = 4;
    
    % build initial window
    temp_wind = zeros(span, IN_IDX_COL);
    temp_wind(1:span, VALUE_COL) = x(1:span);
    temp_wind(:, IN_IDX_COL) = [1:span];
    temp_wind = sortrows(temp_wind, VALUE_COL);
    temp_wind(:, PREV_IDX_COL) = [0 ; temp_wind(1:span - 1, IN_IDX_COL)];
    temp_wind(:, NEXT_IDX_COL) = [temp_wind(2:span, IN_IDX_COL) ; 0];
    min_idx = temp_wind(1, IN_IDX_COL);
    prop_idx = temp_wind(round(span * proportion), IN_IDX_COL);
    max_idx = temp_wind(span, IN_IDX_COL);
    temp_wind = sortrows(temp_wind, IN_IDX_COL);
    wind_mat = temp_wind(1:span, VALUE_COL:NEXT_IDX_COL);
    clear temp_wind;
    wind_mat(span + 2, VALUE_COL) = -inf;
    wind_mat(span + 2, NEXT_IDX_COL) = min_idx;
    wind_mat(min_idx, PREV_IDX_COL) = span + 2;
    wind_mat(span + 3, VALUE_COL) = inf;
    wind_mat(span + 3, PREV_IDX_COL) = max_idx;
    wind_mat(max_idx, NEXT_IDX_COL) = span + 3;
    
    % slide window
    prop_val = wind_mat(prop_idx, VALUE_COL);
    px(1:round(span / 2)) = prop_val;
    out_idx = round(span / 2) + 1;
    oldest_idx = 1;
    newest_idx = span;
    span_plus_1 = span + 1;
    blank_idx = span_plus_1;
    in_idx = span_plus_1;
    while in_idx <= len
        
        % insert new value into blank
        new_val = x(in_idx);
        wind_mat(blank_idx, VALUE_COL) = new_val;
        newest_val = wind_mat(newest_idx, VALUE_COL);
        curr_idx = newest_idx;
        if new_val > prop_val
            prop_shift = 0.5;
        elseif new_val < prop_val
            prop_shift = -0.5;
        else
            prop_shift = 0;
        end
        if new_val > newest_val
            % search forward
            while 1
                next_idx = wind_mat(curr_idx, NEXT_IDX_COL);
                if new_val <= wind_mat(next_idx, VALUE_COL) 
                    break;
                end
                curr_idx = next_idx;
            end
            if prop_shift == 0
                prop_shift = val_equals_prop(curr_idx);
                if prop_shift == 0
                    prop_shift = 0.5;  % because we insert after curr_idx
                end
            end
            % insert after curr_idx
            wind_mat(blank_idx, PREV_IDX_COL) = curr_idx;
            wind_mat(blank_idx, NEXT_IDX_COL) = next_idx;
            wind_mat(curr_idx, NEXT_IDX_COL) = blank_idx;
            wind_mat(next_idx, PREV_IDX_COL) = blank_idx;
        else
            % search backward
            while 1
                prev_idx = wind_mat(curr_idx, PREV_IDX_COL);
                if new_val >= wind_mat(prev_idx, VALUE_COL) 
                    break;
                end
                curr_idx = prev_idx;
            end
            if prop_shift == 0
                prop_shift = val_equals_prop(curr_idx);
                if prop_shift == 0
                    prop_shift = -0.5;  % because we insert before curr_idx
                end
            end
            % insert before curr_idx
            wind_mat(blank_idx, NEXT_IDX_COL) = curr_idx;
            wind_mat(blank_idx, PREV_IDX_COL) = prev_idx;
            wind_mat(curr_idx, PREV_IDX_COL) = blank_idx;
            wind_mat(prev_idx, NEXT_IDX_COL) = blank_idx;
        end
        
        % update proportion value
        old_val = wind_mat(oldest_idx, VALUE_COL);
        if old_val > prop_val
            prop_shift = prop_shift - 0.5;
        elseif old_val < prop_val
            prop_shift = prop_shift + 0.5;
        else
            t_shift = val_equals_prop(oldest_idx);
            if t_shift == 0
                if prop_shift > 0
                    prop_shift = 1;
                else
                    prop_shift = -1;
                end
            else
                prop_shift = prop_shift - t_shift;
            end
        end
        if prop_shift == 1
            prop_idx = wind_mat(prop_idx, NEXT_IDX_COL);
            prop_val = wind_mat(prop_idx, VALUE_COL);
        elseif prop_shift == -1
            prop_idx = wind_mat(prop_idx, PREV_IDX_COL);
            prop_val = wind_mat(prop_idx, VALUE_COL);
        end
        px(out_idx) = prop_val;
        
        % remove oldest value
        prev_idx = wind_mat(oldest_idx, PREV_IDX_COL);
        next_idx = wind_mat(oldest_idx, NEXT_IDX_COL);
        wind_mat(prev_idx, NEXT_IDX_COL) = next_idx;
        wind_mat(next_idx, PREV_IDX_COL) = prev_idx;
                               
        % update indices
        oldest_idx = oldest_idx + 1;
        if oldest_idx > span_plus_1
            oldest_idx = 1;
        end
        newest_idx = newest_idx + 1;
        if newest_idx > span_plus_1
            newest_idx = 1;
        end
        blank_idx = blank_idx + 1;
        if blank_idx > span_plus_1
            blank_idx = 1;
        end
        in_idx = in_idx + 1;
        out_idx = out_idx + 1;
        
    end
        
    % fill in tail
    px(out_idx:len) = prop_val;
    
    % start nested functions
    
    function shift = val_equals_prop(idx)
        % find whether passed index is before or after proportion value index
        if idx == prop_idx
            shift = 0;
        else
            % search forward
            n_idx = idx;
            while 1
                n_idx = wind_mat(n_idx, NEXT_IDX_COL);
                if n_idx == prop_idx 
                    shift = -0.5;
                    break;
                end
                if wind_mat(n_idx, VALUE_COL) ~= prop_val
                    shift =  0.5;
                    break;
                end
            end
        end
    end

    % end nested functions    

end

