function [] = draw_epipolar_lines(im_L, im_R, F, p, f1, f2)
% This function takes point p from im_L and draw the corresponding epipolar
% line im im_R using F

p_h = [p ; 1]; % homogeneous coordinates
l = F * p_h;   % epipolar lines
l = l / l(3);

point = lineToBorderPoints(l', size(im_R)); % adjust lines to image dimensions

figure(f2);

line(point([1,3])',point([2,4])');

end

