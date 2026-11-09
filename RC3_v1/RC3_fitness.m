function [F, C] = RC3_fitness(x, P)
%RC3_FITNESS  Objective vector and constraint violations of one (repaired)
%decision vector, evaluated through the MP-Opt-Model object.
%
%   [F, C] = RC3_fitness(x, P)
%   F : 1 x P.nObj   m = 1: total model cost  sum_k w_k F_k / Fref_k
%                    m > 1: the named costs f_<name> of the user objectives
%                    (P.objNames, e.g. f1..f5) in the order of opts.objectives
%   C : 1 x P.nCon   [nonlinear set 'net', one aggregated violation per
%                    linear set]; feasible if all C <= 0
%
%   RC3_fitness('reset', maxFE)   clear the evaluation log
%   H = RC3_fitness('get')        log: H.F (FE x nObj), H.CV (FE x 1)
persistent LOG
if ischar(x)
    switch x
        case 'reset'
            n = 1e4;  if nargin > 1, n = P; end
            LOG = struct('F', nan(n, 0), 'CV', nan(n, 1), 'n', 0);
        case 'get'
            if isempty(LOG), F = struct('F', [], 'CV', [], 'n', 0); return; end
            F = LOG;  F.F = F.F(1:F.n, :);  F.CV = F.CV(1:F.n);
    end
    return
end
x = reshape(x, 1, []);
if P.m == 1
    F = RC3_om_eval(P, 'costall', x);
else
    F = zeros(1, P.M);
    for k = 1:P.M
        F(k) = RC3_om_eval(P, 'cost', x, ['f_' P.objNames{k}]);
    end
end
if nargout < 2, return; end
g  = RC3_om_eval(P, 'nln', x);
A  = P.lin.A;
if isempty(A)
    lv = zeros(1, 0);
else
    Ax = A * x';
    vr = max(Ax - P.lin.u, 0) + max(P.lin.l - Ax, 0);           % row violations
    vr(vr < 1e-7) = 0;
    lv = cellfun(@(r) sum(vr(r)), P.lin.rows);
end
C = [g(:)', lv];
if ~isempty(LOG)                                                % evaluation log
    if LOG.n >= size(LOG.CV, 1)
        LOG.CV = [LOG.CV; nan(size(LOG.CV))];  LOG.F = [LOG.F; nan(size(LOG.F))];
    end
    if size(LOG.F, 2) ~= numel(F), LOG.F = nan(size(LOG.CV, 1), numel(F)); LOG.n = 0; end
    LOG.n = LOG.n + 1;
    LOG.F(LOG.n, :) = F;  LOG.CV(LOG.n) = sum(max(C, 0));
end
end
