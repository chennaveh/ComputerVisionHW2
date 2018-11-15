function [output] = proj(M, P)
% P - the 3D point
% M - projection  matrix
% output - the  2D coordinates of the projected point.
P(4) = 1; % homogeneous coordinates
output = M*P;
output = output / output(3);
output = output(1:2);
end

