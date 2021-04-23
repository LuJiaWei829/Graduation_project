clear;tic
img_l = imread('E:\MATLAB\R2016b\bin\bishe\HOPC\data\infrared_sen.tif');
% img_l = imread('.\lena.bmp');
img_l_copy = img_l;
if ndims(img_l) == 3
    img_l = rgb2gray(img_l);
end
[l_r,l_c] = size(img_l);
img_r = imread('E:\MATLAB\R2016b\bin\bishe\HOPC\data\visible_ref.tif');
% img_r = imread('.\lena.bmp');
img_r_copy = img_r;
if ndims(img_r) == 3
    img_r = rgb2gray(img_r);
end
[r_r,r_c] = size(img_r);
CP_Check_file = 'E:\MATLAB\R2016b\bin\bishe\HOPC\data\VisibletoInfrared_CP.txt';

points_to_match = harris(img_r);
temp_rad = 10;
tranFlag = 3;
disthre = 10;
C = 0;    % the number of correct match
CM = 0;    % the number of total match
C_e = 0;    % the number of mismatch
e = 0.0000001;    % avoid divided by zero
pnum = size(points_to_match,1);

%% exclude the points out of the edge of template
exclude_points_num = 0;
exclude_points_no = zeros(pnum,1);
for i=1:pnum
    if(points_to_match(i,1)<temp_rad || points_to_match(i,2)<temp_rad || points_to_match(i,1)>r_r-temp_rad || points_to_match(i,2)>r_c-temp_rad)
        exclude_points_num = exclude_points_num + 1;
        exclude_points_no(i) = 1;
    end
end
match_points = zeros(pnum-exclude_points_num,2);
[ir,~] = find(exclude_points_no==0);
pnum1 = pnum-exclude_points_num;
for i=1:pnum1
    match_points(i,1) = points_to_match(ir(i),1);
    match_points(i,2) = points_to_match(ir(i),2);
end

%% Caculate MI
margin = 20;
figure(1);
imshow(img_r_copy);
title('Reference Image');
hold on;
for i=1:pnum1
   scatter(match_points(i,2),match_points(i,1),30,'filled','y') ;   
end
figure(2);
imshow(img_l_copy);
title('Sensed Image');
hold on;
rc = zeros(pnum1,2);
disp(pnum1)
for i=1:pnum1
    template = img_r(match_points(i,1)-temp_rad/2+1:match_points(i,1)+temp_rad/2,match_points(i,2)-temp_rad/2+1:match_points(i,2)+temp_rad/2);
    ET = entropy(template);   %%Define template window and caculate the entropy of template window
    [m,n] = size(template);
    Val_MI = zeros(l_r-temp_rad,l_c-temp_rad);
    disp(i)
    for j=max(1,match_points(i,1)-temp_rad/2+1-margin):min(l_r-temp_rad,match_points(i,1)+temp_rad/2+margin)%j=1:l_r-temp_rad%
        for k=max(1,match_points(i,2)-temp_rad/2+1-margin):min(l_r-temp_rad,match_points(i,2)+temp_rad/2+margin)% k=1:l_c-temp_rad%
            img_s = img_l(j:j+temp_rad-1,k:k+temp_rad-1);  %%define search window
            Val_MI(j,k) = MI(template,img_s);
%             ES = entropy(img_s);  %%caculate the entropy of search window
%             histq = zeros(256,256);
%             for s=1:m
%                 for t=1:n
%                     x = template(s,t) + 1;
%                     y = img_s(s,t) + 1;
%                     histq(x,y) = histq(x,y) + 1;
%                 end
%             end
%             p = histq./sum(sum(histq));  %%caculate Joint probability density
%             EST = -sum(sum(p.*log2(p+eps)));
%             MI(j,k) = ES + ET - EST;
        end
    end
    mi_max=max(max(Val_MI));
    [r,c]=find(mi_max==Val_MI);
    if size(r,1)~=1 || size(c,1) ~=1
        r = r(1);
        c = c(1);
    end
    rc(i,1) = r;rc(i,2) = c;
    scatter(c+temp_rad/2,r+temp_rad/2,30,'filled','y');
end

% checkPt = textread(CP_Check_file);
% refpt = [checkPt(:, 1), checkPt(:, 2)]; %the check points in the referencing image
% senpt = [checkPt(:, 3), checkPt(:, 4)]; %the check points in the sensed image
% % solve the geometric tranformation parameter
% % tran 0:affine, 1: projective, 2: Quadratic polynomial,3: cubic polynomial,the default is 3
% tform = [];
% if tranFlag == 0
%     tform = cp2tform(refpt, senpt, 'affine');
%     T = tform.tdata.T;
% elseif tranFlag == 1
%     tform = cp2tform(refpt, senpt, 'projective');
%     T = tform.tdata.T;
% else
%     T = solvePoly(refpt, senpt, tranFlag);
% end
% H = T';   % the geometric transformation parameters from im_Ref to im_Sen
% for n = 1:pnum1
%     %transform the (x,y) of reference image to sensed image by the geometric relationship of check points
%     %to determine the search region
%     tempCo1 = transferTo(tform, match_points(n,:), H, tranFlag);
%     % tranformed coordinate (X_Sen_c, Y_Sen_c)
%     X_Sen_c = tempCo1(1);
%     Y_Sen_c = tempCo1(2);
%     Y_match = rc(n,2);
%     X_match = rc(n,1);
%     % calculate the match errors
%     diffY = abs(Y_match - Y_Sen_c);
%     diffX = abs(X_match - X_Sen_c);
%     diff = sqrt(diffX .^ 2 + diffY .^ 2);
%     % calculate the numbers of correct match, mismatch and total match
%     if diff <= disthre
%         C = C + 1; % the number of correct matches
%     else
%         C_e = C_e + 1;        
%     end
%     CM = CM + 1;
% end
% %the correct ratio
% CMR = C / CM;
[CMR,mat_diff] = Caculate_CMR(CP_Check_file,pnum1,match_points,rc,tranFlag,disthre);
x= sprintf('the correct match ratio is %4.3f', CMR);
disp(x)
toc