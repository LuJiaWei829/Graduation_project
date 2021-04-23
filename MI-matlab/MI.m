% 
% This is function compute mutual information between two gray images of the same size.
% 
% MI(A, B) = H(A) + H(B) - H(AB)
% 
% Parameter
% ---------
%   img_l: left image
%   ima_r: right image
% 
% Return
% ------
%   val : computed mutual information
% 

function val = MI(img_l, img_r)
    ent_l = entropy(img_l);    % H(A)
    ent_r = entropy(img_r);    % H(B)
    temp = zeros(256, 256);
    [rows, cols] = size(img_l);
    % compute joint probability distribution
    for i = 1:rows
       for j = 1:cols
           pixel_l = img_l(i, j);
           pixel_r = img_r(i, j);
           temp(pixel_l + 1, pixel_r + 1) = temp(pixel_l + 1, pixel_r + 1) + 1;
       end
    end
    temp = temp ./ sum(sum(temp));   % calculate the probability
    epsilon = 1e-6;    % prevent temp from being zero, the logarithm function requires the variable to be greater than 0
    ent_cross = sum(sum(-temp .* log2(temp + epsilon)));    % calculate joint entropy
    val = ent_l + ent_r - ent_cross;
end