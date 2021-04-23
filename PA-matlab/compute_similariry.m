%
% the implementation of the equation (11) in the paper.
% 
function s = compute_similariry(descriptor_l, descriptor_r)
    temp = abs(pi / 2 - abs(descriptor_l - descriptor_r));
    s = sum(sum(temp));
end