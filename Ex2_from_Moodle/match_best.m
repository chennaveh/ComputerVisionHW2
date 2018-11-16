function [index_R] = match_best(Features_L, Features_R, ind_L)
% finds the most similar feature to Features_L(ind_L)

feature_ref = Features_L(ind_L, :);

[~, i] = min(vecnorm(Features_R - feature_ref, 2, 2)); % calc 2nd norm on rows
index_R = i;

end

