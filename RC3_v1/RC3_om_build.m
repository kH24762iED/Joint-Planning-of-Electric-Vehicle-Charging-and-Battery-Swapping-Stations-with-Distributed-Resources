function [om, api] = RC3_om_build(P)
%RC3_OM_BUILD  MP-Opt-Model formulation of the problem defined by RC3_setup.
%
%   [om, api] = RC3_om_build(P)
%
%   Variables     one named set per decision block (P.blocks): initial value
%                 = base case, bounds, type 'C' | 'I' | 'B'
%   Linear cons.  l <= A x <= u : device-count limits, tap ramp limits, DR
%                 energy neutrality and curtailment budget, BSS thresholds
%   Nonlinear     g(x) <= 0, set 'net': radiality, voltage, thermal,
%   constraints   operations, DG penetration, siting, BSS energy, EV
%                 service, islanding, frequency, reverse power (P.conNl)
%   Costs         one general nonlinear cost f_<name> per user objective
%                 (opts.objectives: basic, weighted sum such as f2, or handle);
%                 with m = 1 the costs are w_k * F_k / Fref_k, so the total
%                 cost of the model is the weighted single objective.
%
%   api = 'mp'     MP-Opt-Model >= 5.0 (mp.opt_model, MATPOWER >= 8.1)
%         'legacy' MP-Opt-Model 2.x-4.x (opt_model)
%         'none'   MP-Opt-Model not on the path: same data in a struct
%
%   The load flow inside the cost/constraint functions is solved by
%   RC3_objectives, so the functions are black boxes; they are evaluated by
%   the metaheuristic (PlatEMO or RC3_nsga2) through this model, and their
%   finite-difference derivatives are available for RC3_refine.

api = 'none';
try
    om = mp.opt_model();  api = 'mp';
catch
    try
        om = opt_model();  api = 'legacy';
    catch
        warning('RC3_om_build:none', ['MP-Opt-Model was not found (add MATPOWER to the path); ' ...
                'the model is kept in an equivalent internal struct.']);
        om = struct('var', struct('name', {}, 'N', {}, 'v0', {}, 'vl', {}, 'vu', {}, 'vt', {}), ...
                    'lin', struct('name', {}, 'A', {}, 'l', {}, 'u', {}), ...
                    'nli', struct('name', {}, 'N', {}, 'fcn', {}), ...
                    'nlc', struct('name', {}, 'fcn', {}));
    end
end
vtc = 'CI-B';                                        % enc 1 -> C, 2 -> I, 4 -> B
D = P.D;

%% ---- variables ----------------------------------------------------------------------
for b = P.blocks
    N  = numel(b.idx);
    v0 = P.xbase(b.idx)';  vl = b.lo';  vu = b.up';
    vt = repmat(vtc(b.enc(1)), 1, N);
    switch api
        case 'mp',     om.var.add(b.name, N, v0, vl, vu, vt);
        case 'legacy', om.add_var(b.name, N, v0, vl, vu, vt);
        otherwise,     om.var(end+1) = struct('name', b.name, 'N', N, 'v0', v0, 'vl', vl, 'vu', vu, 'vt', vt);
    end
end

%% ---- linear constraints ---------------------------------------------------------------
for c = P.linsets
    switch api
        case 'mp',     om.lin.add(om.var, c.name, c.A, c.l, c.u);
        case 'legacy', om.add_lin_constraint(c.name, c.A, c.l, c.u);
        otherwise,     om.lin(end+1) = struct('name', c.name, 'A', c.A, 'l', c.l, 'u', c.u);
    end
end

%% ---- nonlinear inequality constraints g(x) <= 0 --------------------------------------
gfcn = @(x) RC3_om_fcn('g', x, P);
hess = @(x, lam) sparse(D, D);                       % not used (derivative-free)
switch api
    case 'mp',     om.nli.add(om.var, 'net', P.nConNl, gfcn, hess);
    case 'legacy', om.add_nln_constraint('net', P.nConNl, 0, gfcn, hess);
    otherwise,     om.nli(end+1) = struct('name', 'net', 'N', P.nConNl, 'fcn', gfcn);
end

%% ---- nonlinear costs: one per selected objective ---------------------------------------
for k = 1:P.M
    if P.m == 1, sc = P.weights(k) / P.Fref(k); else, sc = 1; end
    name = ['f_' P.objNames{k}];
    ffcn = @(x) RC3_om_fcn('f', x, P, k, sc);
    switch api
        case 'mp',     om.nlc.add(om.var, name, 1, ffcn);
        case 'legacy', om.add_nln_cost(name, 1, ffcn);
        otherwise,     om.nlc(end+1) = struct('name', name, 'fcn', ffcn);
    end
end
end
