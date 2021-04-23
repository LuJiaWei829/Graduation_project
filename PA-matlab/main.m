% read SAR and optical images
tic
img_l = imread('E:\MATLAB\R2016b\bin\bishe\HOPC\data\optical_ref.png');
% img_l_temp = img_l;
if ndims(img_l) == 3
    img_l = rgb2gray(img_l);
end
img_r = imread('E:\MATLAB\R2016b\bin\bishe\HOPC\data\SAR_sen.png');
% img_r_temp = img_r;
if ndims(img_r) == 3
    img_r = rgb2gray(img_r);
end
CP_Check_file = 'E:\MATLAB\R2016b\bin\bishe\HOPC\data\OpticaltoSAR_CP.txt';
% arguments of the algorithm
sigma = 0.6;
sigma_ = 1.2;
C = 0;    % the number of correct match
CM = 0;    % the number of total match
C_e = 0;    % the number of mismatch
tranFlag = 3;
disthre = 15;

%% 


% compute gaussian derivatives.
% note that we let x-axis points to the down side and y-axis points to
% the left side, so that x corresponds to i, y corresponds to j.
x = round(-2 * sigma) : (round(2 * sigma));
dx = -x .* exp(-x .* x / (2 * sigma * sigma)) / sqrt(2 * pi * power(sigma, 4));%加不加sqrt怎么差不多？
dy = dx';
img_l_dx = conv2(img_l, dx', 'same');
img_l_dy = conv2(img_l, dy', 'same');
img_r_dx = conv2(img_r, dx', 'same');
img_r_dy = conv2(img_r, dy', 'same');

% compute orientation angle and structure tensor
theta_l = atan2(img_l_dy, img_l_dx);
tensor_l = compute_structure_tensor(img_l_dx, img_l_dy);
theta_r = atan2(img_r_dy, img_r_dx);
tensor_r = compute_structure_tensor(img_r_dx, img_r_dy);

% do tensor voting and calculate tensor orientation
neighbor = 3 * sigma_;
tensor_l = do_tensor_voting(tensor_l, theta_l, sigma_, neighbor);
tensor_r = do_tensor_voting(tensor_r, theta_r, sigma_, neighbor);
theta_l = calculate_tensor_orientation(tensor_l);
theta_r = calculate_tensor_orientation(tensor_r);

%% start matching
% points_to_match = [200, 200; 100, 100];
points_to_match = harris(img_l);
matched_points = PA_match(theta_l, theta_r, points_to_match, 20, 101);

%% draw the matching result
% imshow(img_l);
% n = size(matched_points, 1);
% for i = 1:n
%     index = matched_points(i, 3);
%     x = points_to_match(index, 1);
%     y = points_to_match(index, 2);
% %     rectangle('Position', [y - 50, x - 50, 100, 100], 'EdgeColor', 'r');
%     hold on;
%     scatter(y, x, 6, 'filled', 'r');
% end
% 
% figure;
% imshow(img_r);
% for i = 1:n
%     x = matched_points(i, 1);
%     y = matched_points(i, 2);
% %     rectangle('Position', [y - 50, x - 50, 100, 100], 'EdgeColor', 'r');
%     hold on;
%     scatter(y, x, 6, 'filled', 'r');
% end

figure;
[rows_l, cols_l] = size(img_l);
[rows_r, cols_r] = size(img_r);
max_rows = max(rows_l, rows_r);
img = zeros(max_rows, cols_l + cols_r);
img(1:rows_l, 1:cols_l) = img_l;
img(1:rows_r, (cols_l + 1):end) = img_r;
img = uint8(img);
imshow(img);
hold on;
n = size(matched_points, 1);
match_l = zeros(n,2);
match_r = zeros(n,2);
for i = 1:n
    index = matched_points(i, 3);
    x_l = points_to_match(index, 1);
    y_l = points_to_match(index, 2);
    match_l(i,1) = x_l;
    match_l(i,2) = y_l;
    scatter(y_l, x_l, 40, 'filled', 'y');
    
    x_r = matched_points(i, 1);
    y_r = matched_points(i, 2);
    match_r(i,1) = x_r;
    match_r(i,2) = y_r;
    scatter(y_r + cols_l, x_r, 40, 'filled', 'y');
    
%     plot([y_l, y_r + cols_l], 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'y', 'MarkerSize', 5);
%     plot([x_l, x_r], 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'y', 'MarkerSize', 5);
    plot([y_l, y_r + cols_l], [x_l, x_r]);
end
[CMR,mat_diff] = Caculate_CMR(CP_Check_file,n,match_l,match_r,tranFlag,disthre);
x= sprintf('the correct match ratio is %4.3f', CMR);
disp(x)
toc