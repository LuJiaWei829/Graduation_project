%% pre-process
clear;tic;
img_ref1 = imread('E:\MATLAB\R2016b\bin\bishe\HOPC\data\infrared_sen.tif');%read images and check control points file path
img_sen1 = imread('E:\MATLAB\R2016b\bin\bishe\HOPC\data\visible_ref.tif');
% img_ref1 = imread('.\lena.bmp');%read images and check control points file path
% img_sen1 = imread('.\lena.bmp');
CP_Check_file = 'E:\MATLAB\R2016b\bin\bishe\HOPC\data\VisibletoInfrared_CP.txt';
if (size(img_ref1,3)==3)    %guarantee image is gray-image.
    img_ref = double(rgb2gray(img_ref1));
else
    img_ref = double(img_ref1);
end
if (size(img_sen1,3)==3)
    img_sen = double(rgb2gray(img_sen1));
else
    img_sen = double(img_sen1);
end
C = 0;    % the number of correct match
CM = 0;    % the number of total match
C_e = 0;    % the number of mismatch
tranFlag = 3;
disthre = 10;
targ_rad = 15; 
[ref_r,ref_c] = size(img_ref);
[sen_r,sen_c] = size(img_sen);
%% 
pos = Harris(img_ref);          %Harris detection
pnum = size(pos,1);
%% exclude the points out of the edge of template
exclude_points_num = 0;
exclude_points_no = zeros(pnum,1);
for i=1:pnum
    if(pos(i,1)<=targ_rad || pos(i,2)<=targ_rad || pos(i,1)>=ref_r-targ_rad || pos(i,2)>=ref_c-targ_rad)
        exclude_points_num = exclude_points_num + 1;
        exclude_points_no(i) = 1;
    end
end
match_points = zeros(pnum-exclude_points_num,2);
[ir,~] = find(exclude_points_no==0);
pnum1 = pnum-exclude_points_num;
for i=1:pnum1
    match_points(i,1) = pos(ir(i),1);
    match_points(i,2) = pos(ir(i),2);
end
%% 

figure(1);
imshow(img_ref1);
hold on;
title('Reference Image');
for k=1:pnum
    scatter(pos(k,2),pos(k,1),40, 'filled', 'y');
end

[image_x,image_y] = size(img_sen);
avgray_tar =double(0);
figure(2)
imshow(img_sen1);
title('Sensed Image');
hold on;
for i=1:pnum1
   tar_gray = img_ref(match_points(i,1)-targ_rad:match_points(i,1)+targ_rad,match_points(i,2)-targ_rad:match_points(i,2)+targ_rad);
   [target_x,target_y] = size(tar_gray);
   xy = target_x *target_y;
   for x=1:target_x
       for y=1:target_y
           avgray_tar = avgray_tar + tar_gray(x,y);%计算模板的总灰度值
       end
   end
   avgray_tar = avgray_tar/xy;%模板平均灰度值
   for w =1 : (image_x-target_x-1)
       disp(w);
       for h = 1:(image_y-target_y-1) 
           avgray_im=double(0);
           tempt_image = img_sen(w:w+target_x-1,h:h+target_y-1);%裁剪中模板大小的区域
           
           for x  =1:target_x
               for y = 1:target_y
                   avgray_im = avgray_im + tempt_image(x,y);%计算裁剪区域的总灰度值
               end
           end
           avgray_im = avgray_im/xy;%裁剪部分的平均灰度值
           err_tar =double(0);%模板的灰度误差
           err_img =double(0);%裁剪区域的灰度误差
           err_all = double(0);%二者相乘
           var_img =double(0);%裁剪区域标准差
           var_tar = double(0);%模板标准差
           for x=1:target_x
               for y=1:target_y
                   err_tar = (tar_gray(x,y) - avgray_tar);
                   err_img = (tempt_image(x,y)-avgray_im);
                   err_all = err_all+err_tar*err_img;
                   var_tar = power((tar_gray(x,y)-avgray_tar),2)+var_tar;
                   var_img =power((tempt_image(x,y)-avgray_im),2)+var_img;
               end
           end
           ncc(w,h) = ((err_all)/(sqrt(var_tar)*sqrt(var_img)))/(xy-1);%计算NCC值
       end
   end
   max_ncc = max(max(ncc));
   [max_x,max_y] = find(ncc == max(max(ncc))); %找到最大匹配值的位置
   if size(max_x,1)~=1
       max_x=max_x(1);
       max_y=max_y(1);
   end
   max_xy(i,1) = max_x;
   max_xy(i,2) = max_y;
   scatter(max_y+target_y/2,max_x+target_x/2,40, 'filled', 'y')
end
[CMR,mat_diff] = Caculate_CMR(CP_Check_file,pnum1,match_points,max_xy,tranFlag,disthre);
x= sprintf('the correct match ratio is %4.3f', CMR);
disp(x)
toc
