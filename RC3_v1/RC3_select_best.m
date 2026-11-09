function [idx, score] = RC3_select_best(Obj)
%RC3_SELECT_BEST  Best compromise solution of a Pareto set by the fuzzy
%satisfying method: mu = (f_max - f)/(f_max - f_min), score = sum of mu,
%normalised over all solutions; the largest score is selected.
fmin = min(Obj, [], 1);  fmax = max(Obj, [], 1);
rng_ = fmax - fmin;  rng_(rng_ == 0) = 1;
mu   = (fmax - Obj) ./ rng_;
score = sum(mu, 2) / sum(mu(:));
[~, idx] = max(score);
end
