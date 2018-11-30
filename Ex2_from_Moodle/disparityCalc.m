function [D] = disparityCalc(im1, im2,Sx,Sy,d_min, d_max)
%DISPARITYCALC Summary of this function goes here
%   Detailed explanation goes here

width = size(im1,2);
height = size(im1,1);

new_height = height+Sy-mod(Sy,2);
new_width = width+Sx-mod(Sx,2);

%loop over epipolar lines and for each pair calc the disparity
Sy_half_flr = floor(Sy/2);
Sx_half_flr = floor(Sx/2);

im1_padded = zeros(new_height,new_width); % pad zero frame
im1_padded (Sy_half_flr+1:end-Sy_half_flr, Sx_half_flr+1:end-Sx_half_flr) = im1;

im2_padded = zeros(new_height,new_width); % pad zero frame
im2_padded(Sy_half_flr+1:end-Sy_half_flr, Sx_half_flr+1:end-Sx_half_flr) = im2;

D = zeros(size(im1));


for i = 1:height 
    for j = 1:width
        patch1 = im1_padded(i:i+Sy-1, j:j+Sx-1);
        start_idx = min(j+d_min,width);
        end_idx   = min(j+d_max,width);
        results_vec = zeros(1, end_idx - start_idx + 1);
        for k = start_idx:end_idx
            patch2 = im2_padded(i:i+Sy-1, k:k+Sx-1);
            results_vec(k - start_idx + 1) = sum(patch1.*patch2,'all') / sqrt((sum(patch1.^2, 'all')*sum(patch2.^2, 'all')));
        end
        [~, index] = min(results_vec);
        D(i,j) = index;
    end
end

end