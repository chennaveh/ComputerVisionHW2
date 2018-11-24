function [P] = stereo_list(ps1,ps2, ML,MR)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%ps1 and ps2 are m× 2 where m is the number of points

%turn ps1 and ps2 to homogenous by adding a column of 1
p_l = [ps1 ones(size(ps1,1),1)];
p_r = [ps2 ones(size(ps2,1),1)];

%calc COP_R
Null_vector= null(MR);
COP_R = Null_vector(1:3, :) / Null_vector(4);

%calc COP_L
Null_vector= null(ML);
COP_L = Null_vector(1:3, :) / Null_vector(4);

S_L = p_l*pinv(ML)';
S_R = p_r*pinv(MR)';

%change U to be just (X Y Z) cuz the point is in inf
S_L(:,4)= double(S_L(:,4)==zeros(size(S_L(2),1))) + S_L(:,4); %in case of a point in infinity
S_R(:,4)= double(S_R(:,4)==zeros(size(S_R(2),1))) + S_R(:,4); %in case of a point in infinity

S_L_nh = S_L(:,1:3)./repmat(S_L(:,4),1,3);
S_R_nh = S_R(:,1:3)./repmat(S_R(:,4),1,3);

U_L = S_L_nh-COP_L';
U_R = S_R_nh-COP_R';

%b = COP_L-COP_R

dominator = (COP_L-COP_R);
lambda=[(-1)*U_L;U_R]'\dominator ;% TODO - make sure this is true!

%based on lambda recover the point P
x=1;
for i=1:2:size(lambda)
    PL(x,:) = lambda(i)*COP_L' + (1-lambda(i))*(S_L_nh(x,:));
    PR(x,:) = lambda(i+1)*COP_R' + (1-lambda(i+1))*(S_R_nh(x,:));
    x=x+1;
end

%P=PL;
P = ((PR + PL)/2);%find the avarage between reconstruction of left and right
%P = (P(:,1:3)/P(:,4))';

end

