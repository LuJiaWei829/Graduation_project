function ncc = compute_ncc(x, y)
    % ��չƽ
    x = x(:);
    y = y(:);
    
    % �����ֵ
    n = length(x);
    mean_x = mean(x);
    mean_y = mean(y);

    sigma_xy = sum((x - mean_x) .* (y - mean_y)) / n;
    d_x = sum(power(x - mean_x, 2)) / n;
    d_y = sum(power(y - mean_y, 2)) / n;
    
    ncc = sigma_xy / sqrt(d_x * d_y);
end