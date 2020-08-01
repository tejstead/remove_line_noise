function [crits, signs] = find_crits(x)

    len = length(x); 
    crits = zeros(len, 1);
    signs = zeros(len, 1);
    crits(1) = 1;
    crit_count = 1;
    
    i = 1;
    j = 2;
    while j <= len
        if x(1) == x(j)
            j = j + 1;
        else
            break;
        end
    end
    
    PEAK = 1;
    TROUGH = -1;
    if x(1) < x(j)
        mode = PEAK;
    else
        mode = TROUGH;
    end
    signs(1) = -mode;
    
    i = j - 1;
    while j <= len
        if mode == PEAK
            while j <= len
                if x(j) > x(i)
                    j = j + 1;
                    i = j - 1;
                elseif x(i) == x(j)
                    j = j + 1;
                else
                    break;
                end
            end
            mode = TROUGH;
        else % mode == TROUGH
            while j <= len
                if x(j) < x(i)
                    j = j + 1;
                    i = j - 1;
                elseif x(i) == x(j)
                    j = j + 1;
                else
                    break;
                end
            end
            mode = PEAK;
        end
        crit_count = crit_count + 1;        
        if i == (j - 1)
            crits(crit_count) = i;
        else
            crits(crit_count) = round((i + j) / 2);
            i = j - 1;    
        end
        signs(crit_count) = -mode;
    end

    if crits(crit_count) ~= len
        crit_count = crit_count + 1;
        crits(crit_count) = len;
        signs(crit_count) = mode;
    end

    crits = crits(1:crit_count);
    signs = signs(1:crit_count);
    
end
