function [Obj, Con, X] = RC3_evaluate(X, P)
%RC3_EVALUATE  The single evaluation path used by every solver:
%   1  repair every row (RC3_repair, deterministic and idempotent)
%   2  objectives and constraint violations through the MP-Opt-Model
%      object (RC3_fitness -> RC3_om_eval -> RC3_om_fcn -> RC3_objectives)
%
%   [Obj, Con, X] = RC3_evaluate(X, P)
%   X   : n x P.D decision vectors (returned repaired)
%   Obj : n x P.nObj objectives (m = 1: weighted single objective)
%   Con : n x P.nCon constraint violations (feasible if all <= 0)
%RC3_PlatEMO performs the same two steps in PROBLEM.Evaluation (CalDec =
%RC3_repair, CalObj/CalCon = RC3_fitness), so both solvers return identical
%values for identical decision vectors (see RC3_check_solvers).
X = RC3_repair(X, P);
n = size(X, 1);
Obj = zeros(n, P.nObj);  Con = zeros(n, P.nCon);
for i = 1:n
    [Obj(i, :), Con(i, :)] = RC3_fitness(X(i, :), P);
end
end
