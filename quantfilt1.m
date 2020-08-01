function qx = quantfilt1(x, quantile, span, varargin)

    %
    %    qx = quantfilt1(x, quantile, span, [tail_option])
    %
    %    quantfilt1() is a span-order quantile filter for one-dimensional arrays
    %    quantile is the proportion of the window data whose value is returned at the center of each window
    %    span is the filter order (window size)
    %
    %    quantile range: [0 to 1]       Examples:    0.0 == local minimum filter
    %                                                0.5 == local median filter
    %                                                1.0 == local maximum filter
    %
    %    For odd span, qx(i) is the quantile value of x((i - ((span - 1) / 2)) : (i + ((span - 1) / 2)))
    %    For even span, qx(i) is the quantile value of x((i - ((span - 2) / 2)) : (i + (span / 2)))
    %
    %    Tail half-windows are treated with the following options:
    %        'truncate': fill with the quantile values of progressively shrinking windows (default)
    %        'extrapolate': fill with the quantile values of closest full window
    %        'zeropad': fill with zeros
    %
    
    % check arguments
    DEFAULT_TAIL_OPTION = 'truncate';
    if nargin == 3
        tail_option = DEFAULT_TAIL_OPTION;
    elseif nargin == 4
        tail_option = varargin{1};
        if (strcmp(tail_option, 'truncate') || strcmp(tail_option, 'extrapolate') || strcmp(tail_option, 'zeropad')) == 0
            disp('quantfilt1(): unrecognized tail option');
            return;
        end
    else
        disp('quantfilt1(): improper number of input arguments');
        return;
    end        
    if nargout ~= 1
        disp('quantfilt1(): improper number of output arguments');
        return;
    end        
    if span < 0
        disp('quantfilt1(): span must be a positive integer');
        return;
    end        
    if (quantile < 0) || (quantile > 1)
        disp('quantfilt1(): quantile out of range');
        return;
    end    
    x = x(:);
    if size(x, 2) ~= 1
        disp('quantfilt1(): input is not a one-dimensional array');
        return;
    end
    len = length(x);
    if len < 1
        qx = [];
        return;
    end
    if len == 1
        qx = x;
        return;
    end
    if span > len 
        span = len;
    end
        
    % set up
    nodes = zeros(span + 3, 3);
    qx = zeros(len, 1);    
    VAL = 1;
    PREV_IDX = 2;
    NEXT_IDX = 3;
    span_plus_1 = span + 1;
    INIT_BLANK_IDX = span_plus_1;
    HEAD_IDX = span + 2;
    TAIL_IDX = span + 3;
        
    % build initial window (tail_option == 'truncate')
    nodes(HEAD_IDX, VAL) = -realmax;
    nodes(HEAD_IDX, NEXT_IDX) = 1;
    nodes(TAIL_IDX, VAL) = realmax;
    nodes(TAIL_IDX, PREV_IDX) = 1;
    nodes(1, VAL) = x(1);
    nodes(1, NEXT_IDX) = TAIL_IDX;
    nodes(1, PREV_IDX) = HEAD_IDX;
    out_idx = 2;
    odd_span = mod(span, 2);
    if odd_span
    	qx(1) = x(1);
        prev_new_idx = 1;
        new_idx = 2;
        in_idx = 2;
    else % even span
        nodes(2, VAL) = x(2);
        if x(2) < x(1)
            nodes(2, NEXT_IDX) = 1;
            nodes(2, PREV_IDX) = HEAD_IDX;
            nodes(1, PREV_IDX) = 2;
            nodes(HEAD_IDX, NEXT_IDX) = 2;
        else
            nodes(2, NEXT_IDX) = TAIL_IDX;
            nodes(2, PREV_IDX) = 1;
            nodes(1, NEXT_IDX) = 2;
            nodes(TAIL_IDX, PREV_IDX) = 2;
        end
        qx(1) = ((1 - quantile) * x(1)) + (quantile * x(2));
        prev_new_idx = 2;
        new_idx = 3;
        in_idx = 3;
    end

    while in_idx < span
        
        % insert a new node
        prev_new_val = nodes(prev_new_idx, VAL);
        new_val = x(in_idx);
        nodes(new_idx, VAL) = new_val;
        curr_idx = prev_new_idx;
        if new_val >= prev_new_val
                % search forward
                while 1
                    next_idx = nodes(curr_idx, NEXT_IDX);
                    if new_val < nodes(next_idx, VAL)
                        break;
                    end
                    curr_idx = next_idx;
                end
                % insert new_idx in front of curr_idx
                nodes(new_idx, PREV_IDX) = curr_idx;
                nodes(new_idx, NEXT_IDX) = next_idx;
                nodes(curr_idx, NEXT_IDX) = new_idx;
                nodes(next_idx, PREV_IDX) = new_idx;
        else % new_val < prev_new_val
                % search backward
                while 1
                    prev_idx = nodes(curr_idx, PREV_IDX);
                    if new_val >= nodes(prev_idx, VAL)
                        break;
                    end
                    curr_idx = prev_idx;
                end
                % insert new_idx behind curr_idx
                nodes(new_idx, NEXT_IDX) = curr_idx;
                nodes(new_idx, PREV_IDX) = prev_idx;
                nodes(curr_idx, PREV_IDX) = new_idx;
                nodes(prev_idx, NEXT_IDX) = new_idx;
        end
        
        % insert another new node
        prev_new_idx = new_idx;
        new_idx = new_idx + 1;
        in_idx = in_idx + 1;
        prev_new_val = nodes(prev_new_idx, VAL);
        new_val = x(in_idx);
        nodes(new_idx, VAL) = new_val;
        curr_idx = prev_new_idx;
        if new_val >= prev_new_val       
                % search forward
                while 1
                    next_idx = nodes(curr_idx, NEXT_IDX);
                    if new_val < nodes(next_idx, VAL)
                        break;
                    end
                    curr_idx = next_idx;
                end
                % insert before curr_idx
                nodes(new_idx, PREV_IDX) = curr_idx;
                nodes(new_idx, NEXT_IDX) = next_idx;
                nodes(curr_idx, NEXT_IDX) = new_idx;
                nodes(next_idx, PREV_IDX) = new_idx;
        else % new_val < prev_new_val
                % search backward
                while 1
                    prev_idx = nodes(curr_idx, PREV_IDX);
                    if new_val >= nodes(prev_idx, VAL)
                        break;
                    end
                    curr_idx = prev_idx;
                end
                % insert behind curr_idx
                nodes(new_idx, NEXT_IDX) = curr_idx;
                nodes(new_idx, PREV_IDX) = prev_idx;
                nodes(curr_idx, PREV_IDX) = new_idx;
                nodes(prev_idx, NEXT_IDX) = new_idx;
        end

        % calculate output
        if quantile ~= 1
            temp_idx = (quantile * (in_idx - 1)) + 1;
            low_q_idx = floor(temp_idx);
            high_val_q = temp_idx - low_q_idx;
            low_val_q = 1.0 - high_val_q;
            curr_idx = nodes(HEAD_IDX, NEXT_IDX);
            for i = 2:low_q_idx
                curr_idx = nodes(curr_idx, NEXT_IDX);
            end
            low_q_idx = curr_idx;
            low_q_val = nodes(low_q_idx, VAL);
            high_q_idx = nodes(low_q_idx, NEXT_IDX);
            high_q_val = nodes(high_q_idx, VAL);
            qx(out_idx) = (low_q_val * low_val_q) + (high_q_val * high_val_q);
        else % quantile == 1
            low_q_idx = nodes(TAIL_IDX, PREV_IDX);
            low_q_val = nodes(low_q_idx, VAL);
            qx(out_idx) = low_q_val;
        end
        
        % update for next loop
        prev_new_idx = new_idx;
        new_idx = new_idx + 1;
        in_idx = in_idx + 1;
        out_idx = out_idx + 1;
    end
    
    % handle other tail options (initial window)
    if strcmp(tail_option, 'extrapolate')
        qx(1:(out_idx - 2)) = qx(out_idx - 1);
    else
        if strcmp(tail_option, 'zeropad')
            qx(1:(out_idx - 2)) = 0;
        end
    end
    
    % slide window
    oldest_idx = 1;
    new_idx = INIT_BLANK_IDX;
    low_q_val = nodes(low_q_idx, VAL);
    while in_idx <= len

        % insert new value into blank node
        new_val = x(in_idx);
        nodes(new_idx, VAL) = new_val;
        prev_new_val = nodes(prev_new_idx, VAL);
        curr_idx = prev_new_idx;
        if new_val >= prev_new_val
                % search forward
                while 1
                    next_idx = nodes(curr_idx, NEXT_IDX);
                    if new_val < nodes(next_idx, VAL)
                        break;
                    end
                    curr_idx = next_idx;
                end
                % insert after curr_idx
                nodes(new_idx, PREV_IDX) = curr_idx;
                nodes(new_idx, NEXT_IDX) = next_idx;
                nodes(curr_idx, NEXT_IDX) = new_idx;
                nodes(next_idx, PREV_IDX) = new_idx;
        else % new_val < prev_new_val
                % search backward
                while 1
                    prev_idx = nodes(curr_idx, PREV_IDX);
                    if new_val >= nodes(prev_idx, VAL)
                        break;
                    end
                    curr_idx = prev_idx;
                end
                % insert before curr_idx
                nodes(new_idx, NEXT_IDX) = curr_idx;
                nodes(new_idx, PREV_IDX) = prev_idx;
                nodes(curr_idx, PREV_IDX) = new_idx;
                nodes(prev_idx, NEXT_IDX) = new_idx;
        end
        if new_val >= low_q_val
            q_shift = 0.5;
        else % new_val < low_q_val
            q_shift = -0.5;
        end

        % update q node
        oldest_val = nodes(oldest_idx, VAL);
        if oldest_val > low_q_val
            q_shift = q_shift - 0.5;
        elseif oldest_val < low_q_val
            q_shift = q_shift + 0.5;
        else % oldest_val == low_q_val
            if oldest_idx == low_q_idx
                q_shift = q_shift * 2;
            else % oldest_idx ~= low_q_idx
                q_shift = q_shift + 0.5;
            end
        end
                
        % remove oldest node
        prev_idx = nodes(oldest_idx, PREV_IDX);
        next_idx = nodes(oldest_idx, NEXT_IDX);
        nodes(prev_idx, NEXT_IDX) = next_idx;
        nodes(next_idx, PREV_IDX) = prev_idx;  

        if q_shift == 1
            low_q_idx = nodes(low_q_idx, NEXT_IDX);
        else
            if q_shift == -1
                low_q_idx = nodes(low_q_idx, PREV_IDX);
            end
        end

        % output new q value
        if quantile ~= 1
            low_q_val = nodes(low_q_idx, VAL);
            high_q_idx = nodes(low_q_idx, NEXT_IDX);
            high_q_val = nodes(high_q_idx, VAL);
            qx(out_idx) = (low_q_val * low_val_q) + (high_q_val * high_val_q);
        else
            low_q_idx = nodes(TAIL_IDX, PREV_IDX);
            low_q_val = nodes(low_q_idx, VAL);
            qx(out_idx) = low_q_val;
        end
               
        % update rotating indices
        oldest_idx = oldest_idx + 1;
        if oldest_idx > span_plus_1
            oldest_idx = 1;
        end
        prev_new_idx = new_idx;
        new_idx = new_idx + 1;
        if new_idx > span_plus_1
            new_idx = 1;
        end
        in_idx = in_idx + 1;
        out_idx = out_idx + 1;

    end
    
    % build terminal window (tail_option == 'truncate')
    last_sliding_out_idx = out_idx;
    while out_idx < len
        
        % remove oldest node
        prev_idx = nodes(oldest_idx, PREV_IDX);
        next_idx = nodes(oldest_idx, NEXT_IDX);
        nodes(prev_idx, NEXT_IDX) = next_idx;
        nodes(next_idx, PREV_IDX) = prev_idx;
        oldest_idx = oldest_idx + 1;
        if oldest_idx > span_plus_1
            oldest_idx = 1;
        end
        
        % remove next oldest node
        prev_idx = nodes(oldest_idx, PREV_IDX);
        next_idx = nodes(oldest_idx, NEXT_IDX);
        nodes(prev_idx, NEXT_IDX) = next_idx;
        nodes(next_idx, PREV_IDX) = prev_idx;
        oldest_idx = oldest_idx + 1;
        if oldest_idx > span_plus_1
            oldest_idx = 1;
        end

        % calculate output
        span = span - 2;
        if quantile ~= 1
            temp_idx = (quantile * (span - 1)) + 1;
            low_q_idx = floor(temp_idx);
            high_val_q = temp_idx - low_q_idx;
            low_val_q = 1.0 - high_val_q;
            curr_idx = nodes(HEAD_IDX, NEXT_IDX);
            for i = 2:low_q_idx
                curr_idx = nodes(curr_idx, NEXT_IDX);
            end
            low_q_idx = curr_idx;
            low_q_val = nodes(low_q_idx, VAL);
            high_q_idx = nodes(low_q_idx, NEXT_IDX);
            high_q_val = nodes(high_q_idx, VAL);
            qx(out_idx) = (low_q_val * low_val_q) + (high_q_val * high_val_q);
        else
            low_q_idx = nodes(TAIL_IDX, PREV_IDX);
            low_q_val = nodes(low_q_idx, VAL);
            qx(out_idx) = low_q_val;
        end
        
        % update for next loop
        out_idx = out_idx + 1;
    end
    qx(len) = x(len);

    % handle other tail options (terminal window)
    if strcmp(tail_option, 'extrapolate')
        qx(last_sliding_out_idx:len) = qx(last_sliding_out_idx - 1);
    else
        if strcmp(tail_option, 'zeropad')
            qx(last_sliding_out_idx:len) = 0;
        end
    end
    
end
    


