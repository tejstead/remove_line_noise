function [b, m] = lad_reg2(x, y)

    %  least absolute differences linear regression
    %  fit: y = mx + b

    y = y(:);
    x = x(:);
    n = length(y);
    
    min_y = min(y);
    max_y = max(y);
    min_x = min(x);
    max_x = max(x);
    max_m = (max_y - min_y) / (max_x - min_x);
    min_m = -max_m;
        
    upper_m = max_m;
    lower_m = min_m;
    safe_eps = eps * 1000;
    thresh = safe_eps * 10;
    d = max_m - min_m;
    while d > thresh  
        m = (upper_m + lower_m) / 2;
        t = y - (m * x);
        b = median(t);
        lad = sum(abs(t - b));
        m_eps = m + safe_eps;
        t = y - (m_eps * x);
        b_eps = median(t);
        lad_eps = sum(abs(t - b_eps));
        test_m = lad_eps - lad;
        if test_m > 0
            upper_m = m;
        else
            if test_m < 0
                lower_m = m;
            else
                break;
            end
        end
        d = upper_m - lower_m;
    end
    
end

