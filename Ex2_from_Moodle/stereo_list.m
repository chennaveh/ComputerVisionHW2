function [P] = stereo_list(ps1,ps2, ML,MR)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%?? and ?? are ?× 2 where m is the number of points

%turn ps1 and ps2 to homogenous
p_l = [ps1 ones(size(ps1,1),1)];
p_r = [ps2 ones(size(ps2,1),1)];

%calc COP_R
COP_R= null(MR);
%COP_R = Null_vector(1:3, :) / Null_vector(4);

%calc COP_L
COP_L= null(ML);
%COP_L = Null_vector(1:3, :) / Null_vector(4);

%todo need to fix dimenstino of COP
U_L = (p_l*pinv(ML)')'-COP_L;
U_R = (p_r*pinv(MR)')'-COP_R;

%b = COP_L-COP_R

dominator = (COP_L-COP_R);
lambda=[U_L U_R]\dominator ;

%based on lambda recover the point P
PR = lambda(1)*COP_R' + (1-lambda(1))*(p_r*pinv(MR)');
PL = lambda(2)*COP_L' + (1-lambda(2))*(p_l*pinv(ML)');

P = (PR + PL)/2;%find the avarage between reconstruction of left and right
P = (P(:,1:3)/P(:,4))';


end

