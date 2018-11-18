function [P] = stereo_list(ps1,ps2, ML,MR)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%?? and ?? are ?× 2 where m is the number of points

%turn ps1 and ps2 to homogenous
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

%U_L = (S_L/S_L(:,4))-COP_L;%change U to be just (X Y Z) cuz the point is in inf

S_L(:,4)= double(S_L(:,4)==0) + S_L(:,4); %in case of a point in infinity
S_R(:,4)= double(S_R(:,4)==0) + S_R(:,4); %in case of a point in infinity

S_L_nh = S_L(:,1:3)/S_L(:,4);
S_R_nh = S_R(:,1:3)/S_R(:,4);

U_L = S_L_nh-COP_L';
U_R = S_R_nh-COP_R';

%b = COP_L-COP_R

dominator = (COP_L-COP_R);
lambda=[(-1)*U_L;U_R]'\dominator ;

%based on lambda recover the point P
PR = lambda(2)*COP_R' + (1-lambda(2))*(S_R_nh);
PL = lambda(1)*COP_L' + (1-lambda(1))*(S_L_nh);

P = ((PR + PL)/2)';%find the avarage between reconstruction of left and right
%P = (P(:,1:3)/P(:,4))';

p_L=proj(ML,P);
p_R=proj(MR,P);

end

