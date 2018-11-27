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
%mem = zeros(size(im2_padded))-1;
D = zeros(size(im1));

%loop over epipolar lines and for each pair calc the disparity
Sy_half_flr = floor(Sy/2);
Sx_half_flr = floor(Sx/2);

for i=ceil(Sx/2):new_width-ceil(Sx/2)-d_min
    for j=ceil(Sy/2):new_height-ceil(Sy/2)+1
        %need to add minimum on the range dmin:dmax
        v1 = reshape (im1_padded(j-Sy_half_flr:j+Sy_half_flr,i-Sx_half_flr :i+Sx_half_flr),[1,Sx*Sy]);
        x=1;
        for d=i+d_min:i+d_max
            if d < new_width-ceil(Sx/2)+1
                 %   if mem(i,d)==-1 %memoization for each pair option i (first pic) and d (second pic)
                v2 = reshape (im2_padded(j-Sy_half_flr:j+Sy_half_flr,d-Sx_half_flr :d+Sx_half_flr ),[1,Sx*Sy]);
                distance(x) = (v1*v2')/(sqrt(sum(v1.^2))*sqrt(sum(v2.^2)));
                x=x+1;
            end
        end
        [~, index] = min(distance);
        D(j-Sy_half_flr,i-Sx_half_flr) = index;
    end
end
end

