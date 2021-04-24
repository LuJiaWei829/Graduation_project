clear; 
clc;

%% 读取影像
img_ref_origin = imread('..\data\optical_ref.png');    % reference image（Optical）
img_sen_origin = imread('..\data\SAR_sen.png');    % sensed image（SAR）
% img_ref1 = imread('.\lena.bmp');
% img_sen1 = imread('.\lena.bmp');

if (ndims(img_ref_origin) == 3)    % guarantee image is gray-image.
    img_ref = rgb2gray(img_ref_origin);
else
    img_ref = img_ref_origin;
end
img_ref = double(img_ref);
if (ndims(img_sen_origin) == 3)
    img_sen = rgb2gray(img_sen_origin);
else
    img_sen = img_sen_origin;
end
img_sen = double(img_sen);

%% 参数设定
tranFlag = 3;    % 变换模型
disthre = 5;    % 坐标差阈值
targ_rad = 15;    % 模板窗口半径
search_radius = 20;    % 搜索区半径
[ref_r, ref_c] = size(img_ref);
[sen_r, sen_c] = size(img_sen);

%% 提取角点
harris_points = Harris(img_ref);          %Harris detection
pnum = size(harris_points, 1);

% exclude the points out of the edge of template
exclude_points_num = 0;
exclude_points_no = zeros(pnum, 1);
for i = 1:pnum
    x = harris_points(i, 1);
    y = harris_points(i, 2);
    % 如果超限
    if(x <= targ_rad || y <= targ_rad || x >= ref_r - targ_rad || y >= ref_c - targ_rad)
        exclude_points_num = exclude_points_num + 1;
        exclude_points_no(i) = 1;
    end
end
points_to_match = zeros(pnum - exclude_points_num, 2);
[ir, ~] = find(exclude_points_no == 0);

% 下边的代码等价于: points_to_match = harris_points(ir, :);
pnum = pnum - exclude_points_num;
for i = 1:pnum
    points_to_match(i, 1) = harris_points(ir(i), 1);
    points_to_match(i, 2) = harris_points(ir(i), 2);
end

%% 画图
figure(1);
imshow(img_ref_origin);
hold on;
title('Reference Image');
for k = 1:pnum
    scatter(points_to_match(k, 2), points_to_match(k, 1), 40, 'filled', 'y');
end

figure(2);
imshow(img_sen_origin);
title('Sensed Image');
hold on;

%% 匹配
tic;
for i = 1:pnum
    % 当前点
    x = points_to_match(i, 1);
    y = points_to_match(i, 2);
    % reference image上邻域窗口
    x_l = x - targ_rad;
    x_r = x + targ_rad;
    y_l = y - targ_rad;
    y_r = y + targ_rad;
    ref_template = img_ref(x_l:x_r, y_l:y_r);
    
    ncc = zeros(2 * search_radius + 1, 2 * search_radius + 1);
    % 开始搜索
    for j = -search_radius:search_radius
        for k = -search_radius:search_radius
            j_l = j + x - targ_rad;
            j_h = j + x + targ_rad;
            k_l = k + y - targ_rad;
            k_h = k + y + targ_rad;

            % 如果超限
            if j_l < 1 || j_h > sen_r || k_l < 1 || k_h > sen_c
                continue;
            end
            
            % sensed image的模板窗口
            sen_template = img_sen(j_l:j_h, k_l:k_h); 
            
            ncc_ = compute_ncc(sen_template, ref_template);
            ncc(j + search_radius + 1, k + search_radius + 1) = ncc_;
        end
        
    end
    max_ncc = max(max(ncc));
    [max_x, max_y] = find(ncc == max(max(ncc))); %找到最大匹配值的位置
    if size(max_x, 1) ~= 1
        max_x = max_x(1);
        max_y = max_y(1);
    end
    max_xy(i, 1) = max_x - search_radius + 1 + x;
    max_xy(i, 2) = max_y - search_radius + 1 + y;
    scatter(max_xy(i, 2), max_xy(i, 1), 40, 'filled', 'y')
end

CP_Check_file = '..\data\OpticaltoSAR_CP.txt';
[CMR, mat_diff] = Caculate_CMR(CP_Check_file, pnum, points_to_match, max_xy, tranFlag, disthre);
x= sprintf('the correct match ratio is %4.3f', CMR);
disp(x)
toc;
