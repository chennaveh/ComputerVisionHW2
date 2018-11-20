function [D] = disparityCalc(im1, im2,Sx,Sy,d_min, d_max)
%DISPARITYCALC Summary of this function goes here
%   Detailed explanation goes here

width = size(im1,2);
height = size(im1,1);

new_height = height+Sy-mod(Sy,2);
new_width = width+Sx-mod(Sx,2);

im1_padded = zeros(new_height,new_width); % pad zero frame
im1_padded (Sy/2:end-Sy/2, Sx/2:end-Sx/2) = im1;

im2_padded = zeros(new_height,new_width); % pad zero frame
im2_padded(Sy/2:end-Sy/2, Sx/2:end-Sx/2) = im2;
D = zeros(size(im1));

%loop over epipolar lines and for each pair calc the disparity
Sy_half_flr = floor(Sy/2);
Sx_half_flr = floor(Sx/2);

for i=ceil(Sx/2)+d_min:new_width-ceil(Sx/2)+1-d_max
    for j=ceil(Sy/2):new_height-ceil(Sy/2)+1
        %need to add minimum on the range dmin:dmax
        v1 = reshape (im1_padded(j-Sy_half_flr:j+Sy_half_flr,i-Sx_half_flr :i+Sx_half_flr),[1,Sx*Sy]);
        v2 = reshape (im2_padded(j-Sy_half_flr:j+Sy_half_flr,i-Sx_half_flr+d_max :i+Sx_half_flr+d_max ),[1,Sx*Sy]);
 
        distance = (v1*v2')/(sqrt(sum(v1.^2))*sqrt(sum(v2.^2)));        
        D(j-Sy_half_flr,i-Sx_half_flr) = distance;
%         if strcmpi(distanceAlgo,'Cosine')
%             distance = (v1*v2')/(sqrt(sum(v1.^2))*sqrt(sum(v2.^2)));
%         elseif strcmpi(distanceAlgo,'Gradient')
%             [G_X,G_Y] = gradient(v1,v2);
%             distance = (G_X.^2 + G_Y.^2).^0.5;
%         else
%             distance = 0;
%         end 
%         D(j-Sy_half_flr,i-Sx_half_flr) = distance;
    end
end
end

