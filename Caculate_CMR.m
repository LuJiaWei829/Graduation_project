function [CMR,mat_diff]=Caculate_CMR(CP_Check_file,points_number,points_ref,points_sen,tranFlag,disthre)
match_points = points_ref;
rc = points_sen;
pnum1 = points_number;
C = 0;    % the number of correct match
CM = 0;    % the number of total match
C_e = 0;    % the number of mismatch
checkPt = textread(CP_Check_file);
refpt = [checkPt(:, 1), checkPt(:, 2)]; %the check points in the referencing image
senpt = [checkPt(:, 3), checkPt(:, 4)]; %the check points in the sensed image
% solve the geometric tranformation parameter
% tran 0:affine, 1: projective, 2: Quadratic polynomial,3: cubic polynomial,the default is 3
tform = [];
if tranFlag == 0
    tform = cp2tform(refpt, senpt, 'affine');
    T = tform.tdata.T;
elseif tranFlag == 1
    tform = cp2tform(refpt, senpt, 'projective');
    T = tform.tdata.T;
else
    T = solvePoly(refpt, senpt, tranFlag);
end
H = T';   % the geometric transformation parameters from im_Ref to im_Sen
mat_diff = zeros(pnum1,1);
for n = 1:pnum1
    %transform the (x,y) of reference image to sensed image by the geometric relationship of check points
    %to determine the search region
    tempCo1 = transferTo(tform, match_points(n,:), H, tranFlag);
    % tranformed coordinate (X_Sen_c, Y_Sen_c)
    X_Sen_c = tempCo1(1);
    Y_Sen_c = tempCo1(2);
    Y_match = rc(n,2);
    X_match = rc(n,1);
    % calculate the match errors
    diffY = abs(Y_match - Y_Sen_c);
    diffX = abs(X_match - X_Sen_c);
    diff = sqrt(diffX .^ 2 + diffY .^ 2);
    mat_diff(n) = diff;
    % calculate the numbers of correct match, mismatch and total match
    if diff <= disthre
        C = C + 1; % the number of correct matches
    else
        C_e = C_e + 1;        
    end
    CM = CM + 1;
end
%the correct ratio
CMR = C / CM;
end