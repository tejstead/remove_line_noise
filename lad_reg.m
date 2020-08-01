function [b, m] = lad_reg(y)

    %  least absolute differences linear regression
    %  assumes x to be 1:length(x)
    %  fit: y = mx + b

    y = y(:);
    n = length(y);
    x = [1:n]';
    
    min_y = min(y);
    max_y = max(y);
    min_m = (min_y - max_y) / (n - 1);
    max_m = -min_m;
    
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

