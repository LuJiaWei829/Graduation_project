% 
% the implementation of the equation (10) in the paper.
% 
function theta = calculate_tensor_orientation(tensor)
    [rows, cols, ~] = size(tensor);
    theta = zeros(rows, cols);
    for i = 1:rows
        for j = 1:cols
            t11 = tensor(i, j, 1);
            t22 = tensor(i, j, 4);
            t12 = tensor(i, j, 3);
            theta(i, j) = 1 / 2 * atan(2 * t12 / (t11 - t22));
        end
    end
end