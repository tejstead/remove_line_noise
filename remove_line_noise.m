function [ fx ] = remove_line_noise(x, template , sampling_frequency, line_frequency)

    template_len = numel(template);
    template = [template ; template(1)];
    f_template_len = sampling_frequency / line_frequency;
    template_step = f_template_len / template_len;
    tx = 1:template_step:(f_template_len + 1);
    lx = numel(x);
    fx = zeros(size(x));
        
    i_start = 1;
    i_end = floor(tx(end));
    while i_end < lx
        q = spline(tx, template, (i_start:i_end));
        fx(i_start:i_end) = round(x(i_start:i_end) - q');
        tx = tx + f_template_len;
        i_start = i_end + 1;
        i_end = floor(tx(end));
    end
    i_end = lx;
    q = spline(tx, template, (i_start:i_end));
    fx(i_start:i_end) = round(x(i_start:i_end) - q');

end

