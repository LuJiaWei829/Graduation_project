% 
% the implementation of the equation (5) in the paper.
% 
function T = compute_structure_tensor(dx, dy)
    [rows, cols] = size(dx);
    T = zeros(rows, cols, 4);
    T(:, :, 1) = dx .* dx;
    T(:, :, 2) = dx .* dy;
    T(:, :, 3) = dx .* dy;
    T(:, :, 4) = dy .* dy;
end