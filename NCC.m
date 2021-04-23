%% pre-process
clc;tic;
img_ref1 = imread('E:\MATLAB\R2016b\bin\bishe\HOPC\data\optical_ref.png');%read images and check control points file path
img_sen1 = imread('E:\MATLAB\R2016b\bin\bishe\HOPC\data\SAR_sen.png');
cp_check_file = 'E:\MATLAB\R2016b\bin\bishe\HOPC\data\OpticaltoSAR_CP.txt';
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
%% 
pos = Harris(img_ref);          %Harris detection
[ref_r,ref_c] = size(img_ref);
[sen_r,sen_c] = size(img_sen);

targ_rad = 1;               %define template window size
pnum = size(pos,1);
total_gray = zeros(pnum,1);
av_gray = zeros(pnum,1);
traget_window = zeros(pnum,1);
for i=1:pnum         %caculate average target gray values in template window
    target_window = img_ref(pos(i,1)-targ_rad:pos(i,1)+targ_rad,pos(i,2)-targ_rad:pos(i,2)+targ_rad);
    total_gray(i) = sum(sum(target_window));
    av_gray(i) = total_gray(i) / (2*targ_rad+1)^2;
end

[targ_r,targ_c] = size(target_window);
err_tar =double(0);%ģ��ĻҶ����
err_img =double(0);%�ü�����ĻҶ����     
err_all = zeros(pnum,1);%�������       
var_img = zeros(pnum,1);%�ü������׼��        
var_tar = zeros(pnum,1);%ģ���׼��       
ncc = double(0);
for h=1:(sen_r-targ_r-1)
    disp (h);
    for w=1:(sen_c-targ_c-1)
        temp_image = img_sen(h:h+targ_r-1,w:w+targ_c-1);
        total_temp_gray = sum(sum(temp_image));
        av_temp_gray = total_temp_gray / (targ_r * targ_c);      
        for i=1:pnum
           for x=1:targ_r          
               for y=1:targ_c               
                   err_img(x,y) = temp_image(x,y) - av_temp_gray;              
                   target_window = img_ref(pos(i,1)-targ_rad:pos(i,1)+targ_rad,pos(i,2)-targ_rad:pos(i,2)+targ_rad);
                   err_tar(x,y,i) = target_window(x,y) - av_gray(i);        
                   err_all(i) = err_all(i) + err_tar(x,y,i) * err_img(x,y);
                   var_tar(i) = power((target_window(x,y)-av_gray(i)),2) + var_tar(i);
                   var_img(i) = power((temp_image(x,y)-av_gray(i)),2) + var_img(i);
               end               
           end
           ncc(h,w,i) = (err_all(i)/(sqrt(var_tar(i)) * sqrt(var_img(i)))) / (targ_r * targ_c - 1); 
        end
    end
end
each_ncc = zeros(sen_r-targ_r-1,sen_c-targ_c-1);
max_ncc = zeros(pnum,1);
max_x = zeros(pnum,1);
max_y = zeros(pnum,1);
for j=1:pnum
   each_ncc = ncc(:,:,j);
   max_ncc(j) = max(max(each_ncc));
   [max_x(j),max_y(j)] = find(each_ncc == max_ncc(j));
end
figure(1)
imshow(img_ref1);
hold on;
title('Reference Image');
for k=1:pnum
    scatter(pos(k,1),pos(k,2),40, 'filled', 'y');
end
figure(2)
imshow(img_sen1);
hold on;
title('Sensed Image');
for k=1:pnum
    scatter(max_x(k),max_y(k),40, 'filled', 'y');
end
toc
% tic;
% im_gray = double(imread('.\lena.bmp'));
% tar_gray = double(imread('.\lena1.jpg'));
% image = imread('.\lena.bmp');
% target = imread('.\lena1.jpg');
% im_gray = double(rgb2gray(image));
% tar_gray = double(rgb2gray(target));
% [target_x,target_y] = size(tar_gray);
% [image_x,image_y] = size(im_gray);
% xy = target_x *target_y;
% avgray_tar =double(0);
% 
% for x =1:target_x
%     for y = 1: target_y
%         avgray_tar = avgray_tar + tar_gray(x,y);%����ģ����ܻҶ�ֵ
%     end
% end
% avgray_tar = avgray_tar/xy;%ģ��ƽ���Ҷ�ֵ
% for w =1 : (image_x-target_x-1)
%     disp(w);
%     for h = 1:(image_y-target_y-1) 
%         avgray_im=double(0);
%         tempt_image = im_gray(w:w+target_x-1,h:h+target_y-1);%�ü���ģ���С������
%        
%         for x  =1:target_x
%             for y = 1:target_y
%                 avgray_im = avgray_im + tempt_image(x,y);%����ü�������ܻҶ�ֵ
%             end
%         end
%         avgray_im = avgray_im/xy;%�ü����ֵ�ƽ���Ҷ�ֵ
%         err_tar =double(0);%ģ��ĻҶ����
%         err_img =double(0);%�ü�����ĻҶ����
%         err_all = double(0);%�������
%         var_img =double(0);%�ü������׼��
%         var_tar = double(0);%ģ���׼��
%         for x=1:target_x
%             for y=1:target_y
%                 err_tar = (tar_gray(x,y) - avgray_tar);
%                 err_img = (tempt_image(x,y)-avgray_im);
%                 err_all = err_all+err_tar*err_img;
%                 var_tar = power((tar_gray(x,y)-avgray_tar),2)+var_tar;
%                 var_img =power((tempt_image(x,y)-avgray_im),2)+var_img;
%             end
%         end
%         ncc(w,h) = ((err_all)/(sqrt(var_tar)*sqrt(var_img)))/(xy-1);%����NCCֵ
%     end
% end
% max_ncc = max(max(ncc));
% [max_x,max_y] = find(ncc == max(max(ncc))); %�ҵ����ƥ��ֵ��λ��
% figure(1)
% imshow(uint8(target));
% ncc=mapminmax(ncc,-1,1);%��һ����[-1,1]
% figure(2)
% imshow(uint8(image));%��ԭͼ
% hold on;
% show = imagesc(ncc);    %���ݾ����Ϊ����ͼ
% set(show,'AlphaData',0.6); %����͸����
% rectx =[max_y,max_x,target_y,target_x]; %���þ�������
% rectangle('Position',rectx,'Edgecolor','r','LineWidth',2);%�����Σ���shai����2
% toc
% disp(toc)