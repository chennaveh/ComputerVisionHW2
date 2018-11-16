function [Pairs_clean] = remove_incorrect_matches(Pairs, Points_L, Points_R, F, th)
% This function removes incorrect matches using epipolar geometry. 

% points_L_h = [Points_L.Location ones(size(Points_L.Location, 1), 1)]; % homogeneous coordinates
% lines = (F * points_L_h')';

x1 = Points_L.Location(Pairs(:,1), :)';
x2 = Points_R.Location(Pairs(:,2), :)';

dist_vec = sampsonDistance(F, x1, x2);

Pairs_clean = Pairs(dist_vec < th, :);

end

