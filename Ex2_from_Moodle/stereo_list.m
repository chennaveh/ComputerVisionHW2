function [P] = stereo_list(ps1,ps2, ML,MR)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%?? and ?? are ?× 2 where m is the number of points

%turn ps1 and ps2 to homogenous
p_l = [ps1 ones(size(ps1,1),1)];
p_r = [ps2 ones(size(ps2,1),1)];

%calc COP_R
Null_vector = null(M_R);
COP_R = Null_vector(1:3, :) / Null_vector(4);

%calc COP_L
Null_vector = null(M_L);
COP_L = Null_vector(1:3, :) / Null_vector(4);

U_L = pinv(ML).*p_l-COP_L;
U_R = pinv(MR).*p_r-COP_R;

%b = COP_L-COP_R

lambda=[U_L U_R]\(COP_L-COP_R);

%based on lambda recover the point P
PR = lambda*COP_R + (1-lambda)*pinv(MR).*p_r;
PL = lambda*COP_L + (1-lambda)*pinv(ML).*p_l;

P = (PR + PL)\2;%find the avarage between reconstruction of left and right
end

