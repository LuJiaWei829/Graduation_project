% 
% the implementation of tensor voting described in Sectino 3.2.2 in the
% paper. We define the effective neighborhood as circular region rather than
% squared one.
% 
function T = do_tensor_voting(tensor, theta, sigma, radius)
    radius = round(radius);
    [rows, cols, ~] = size(tensor);
    T = zeros(size(tensor));
    for i = 1:rows
       for j = 1:cols
           alpha = theta(i, j);
           for w_i = (i - radius):(i + radius)
               for w_j = (j - radius):(j + radius)
                   % if (w_i, w_j) is out of image boundary
                   if w_i < 1 || w_i > rows || w_j < 1 || w_j > cols
                       continue;
                   end              
                   r = sqrt(power(w_i - i, 2) + power(w_j - j, 2));
                   % if (w_i, w_j) is out of effective neighbor 
                   if r > radius
                       continue;
                   end
                   
                   beta = theta(w_i, w_j);
                   m = exp(-r * r / (2 * sigma * sigma)) * exp(1 - (pi / 2) / abs(pi / 2 - abs(alpha - beta) + eps));
                   T(i, j, :) = T(i, j, :) + m * tensor(w_i, w_j, :);
               end
           end
       end
    end
end