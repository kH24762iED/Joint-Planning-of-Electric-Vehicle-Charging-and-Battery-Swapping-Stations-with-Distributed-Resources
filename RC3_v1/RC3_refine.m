function out = RC3_refine(P, x0, opt)
%RC3_REFINE  Optional local polish of the continuous set-points of a solution
%(single objective, m = 1) with the NLP solver of MP-Opt-Model, all
%integer/binary decisions (sites, counts, types, switches, taps, steps)
%fixed at their values in x0.
%
%   out = RC3_refine(P, x0)
%   out = RC3_refine(P, x0, opt)    opt.alg ('MIPS'), opt.max_it (40), opt.rho (1e3)
%
%   Sub-model (mp.opt_model / opt_model):
%     variables  z in [0, 1]^nc, x_c = lo + z .* (up - lo)  (continuous blocks)
%     cost       f(z) = F_1(x) + rho * sum(max(G(x), 0))    (exact penalty)
%     gradient   forward differences;  Hessian  (|grad|_inf / 0.05) * I
%   The refined point is kept only if it is feasible and better than x0.
if nargin < 3, opt = struct(); end
alg  = getf(opt, 'alg', 'MIPS');  mxit = getf(opt, 'max_it', 40);  rho = getf(opt, 'rho', 1e3);
out  = struct('x', x0, 'improved', false, 'f0', NaN, 'f', NaN, 'eflag', NaN, 'msg', '');
ic   = find(P.enc == 1);
if isempty(ic) || ~any(strcmp(P.om_api, {'mp', 'legacy'}))
    out.msg = 'no continuous variables or MP-Opt-Model not available';  return
end
lo = P.lower(ic);  up = P.upper(ic);  rg = max(up - lo, 1e-12);
z0 = (x0(ic) - lo) ./ rg;
[f0, C0] = RC3_fitness(x0, P);
out.f0 = f0 + rho * sum(max(C0, 0));
fz = @(z) pen_cost(z, x0, ic, lo, rg, P, rho);
try
    if strcmp(P.om_api, 'mp')
        mm = mp.opt_model();
        mm.var.add('z', numel(ic), z0(:), zeros(numel(ic), 1), ones(numel(ic), 1));
        mm.nlc.add(mm.var, 'pen', 1, fz);
    else
        mm = opt_model();
        mm.add_var('z', numel(ic), z0(:), zeros(numel(ic), 1), ones(numel(ic), 1));
        mm.add_nln_cost('pen', 1, fz);
    end
    sopt = struct('alg', alg, 'verbose', 0);
    sopt.mips = struct('max_it', mxit, 'step_control', 1);
    [z, f, eflag] = mm.solve(sopt);
    x = x0;  x(ic) = lo + z(:)' .* rg;
    x = RC3_repair(x, P);
    [F1, C1] = RC3_fitness(x, P);
    out.f = F1 + rho * sum(max(C1, 0));  out.eflag = eflag;
    if all(C1 <= 1e-9) && F1 < f0 - 1e-9
        out.x = x;  out.improved = true;
    end
    out.msg = sprintf('%s: %.6g -> %.6g (exit flag %d)', alg, out.f0, f, eflag);
catch err
    out.msg = ['refinement skipped: ' err.message];
end
end

function [f, df, d2f] = pen_cost(z, x0, ic, lo, rg, P, rho)
z = z(:)';
x = x0;  x(ic) = lo + z .* rg;
f = pen(x, P, rho);
if nargout > 1
    n  = numel(z);  df = zeros(n, 1);  h = 1e-5;
    for j = 1:n
        zp = z;  zp(j) = min(z(j) + h, 1);  if zp(j) == z(j), zp(j) = z(j) - h; end
        xp = x0;  xp(ic) = lo + zp .* rg;
        df(j) = (pen(xp, P, rho) - f) / (zp(j) - z(j));
    end
    if nargout > 2
        d2f = speye(n) * max(norm(df, inf) / 0.05, 1e-6);
    end
end
end

function v = pen(x, P, rho)
[F, C] = RC3_fitness(x, P);
v = F + rho * sum(max(C, 0));
end

function v = getf(s, name, default)
if isfield(s, name) && ~isempty(s.(name)), v = s.(name); else, v = default; end
end
