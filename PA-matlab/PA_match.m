% 
% This function is the algorithm implementation of the paper "A novel dense 
% descriptor based on structure tensor voting for multi-modal image
% matching".
% 
% Parameter
% ---------
%   theta_l : 2D array, uint8 type
%           tensor orientation of the template image defined in the paper
%   theta_r : 2D array, uint8 type
%           tensor orientation of the base image defined in the paper
%   points : 2D array
%            points to be matched, each row contains x coordinate, y coordinate in the img_l.
%   offset : double
%            search region radius(relative to the points parameter)
%   win_size : double
%              template window size
% 
% Return
% ------
%   matched_points : 2D array
%           the matched_points, each row contains x coordinate, y coordinate, corresponding points index in the template image.
%
function matched_points = PA_match(theta_l, theta_r, points, offset, win_size)    
    pNum = size(points, 1);
    [rows_l, cols_l] = size(theta_l);
    [rows_r, cols_r] = size(theta_r);
    matched_points = [];
    % for each point
    for i = 1:pNum    
        pt = points(i, :);
        x = pt(1); y = pt(2);
        
        x_l = x - fix(win_size / 2);
        x_h = x + fix(win_size / 2);
        y_l = y - fix(win_size / 2);
        y_h = y + fix(win_size / 2);
        % when out of boundary
        if x_l < 1 || x_h > rows_l || y_l < 1 || y_h > cols_l
            continue;
        end
        
        % It is the O(t) in the paper
        descriptor_l = theta_l(x_l:x_h, y_l:y_h);
        
        % search for the maximum similarity and the corresponding position
        % in the search region
        max_similar = 0;
        for j = (x - offset):(x + offset)
            for k = (y - offset):(y + offset)
                j_l = j - fix(win_size / 2);
                j_h = j + fix(win_size / 2);
                k_l = k - fix(win_size / 2);
                k_h = k + fix(win_size / 2);
                
                % when out of boundary ? 
                if j_l < 1 || j_h > rows_r || k_l < 1 || k_h > cols_r
                    continue;
                end
                
                % It is the O(w) in the paper
                descriptor_r = theta_r(j_l:j_h, k_l:k_h);
                
                similarity = compute_similariry(descriptor_l, descriptor_r);
                % choose the larger similariry
                if similarity > max_similar
                   max_similar = similarity;
                   max_x = j;
                   max_y = k;
                end
            end
        end
        
        matched_points = [matched_points; max_x, max_y, i];
    end
end