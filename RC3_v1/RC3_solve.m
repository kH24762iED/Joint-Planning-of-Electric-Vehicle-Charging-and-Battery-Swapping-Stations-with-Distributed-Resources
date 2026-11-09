function R = RC3_solve(T, opts, alg, N, maxFE)
%RC3_SOLVE  Build the model (RC3_setup), solve it and post-process the result.
%
%   R = RC3_solve(T, opts)
%   R = RC3_solve(T, opts, alg, N, maxFE)          (backward compatible call)
%
%   opts (besides the RC3_setup options):
%     .solver    'auto' (PlatEMO if on the path, else built-in) | 'platemo' | 'builtin'
%     .alg       PlatEMO algorithm: handle (@LMOCSO), name ('LMOCSO') or cell with
%                parameters ({@FDV, 0.8, 0.4, 5}); default @GA (m = 1),
%                @NSGAII (m = 2, 3), @NSGAIII (m >= 4)
%     .N, .maxFE population size and evaluation budget (100, 20000)
%     .seed      seed of the initial population and of the solver (1)
%     .validate  re-solve the best solution with MATPOWER runpf (true if available)
%     .refine    m = 1: polish the continuous set-points with MP-Opt-Model (false)
%     .verbose   print progress (true)
%
%   Both solvers start from the same initial population (RC3_init_pop, seed
%   opts.seed) and evaluate through the same function (repair + MP-Opt-Model
%   costs and constraints). After the run the final population is
%   evaluated again with RC3_evaluate; R.check reports the largest difference
%   from the values returned by the solver (0 when consistent).
%
%   R.P          problem data (P.om = MP-Opt-Model object)
%   R.F_base     user objectives of the base case,  R.out_base  all its indicators
%   R.Dec, R.Obj feasible non-dominated set (m > 1) or best solutions (m = 1)
%   R.ibest      best compromise (fuzzy satisfying method) or best solution
%   R.x_best, R.F_best, R.out_best   decision vector, objectives, indicators
%   R.Fuser_pareto, R.Fall_pareto, R.Flib_pareto   user, basic and library
%                objectives of every solution of R.Dec
%   R.nfeas      feasible solutions in the final population
%   R.hist       evaluation log (objectives, constraint violation per FE)
%   R.val        MATPOWER validation table of the best solution
%   R.time       run time (s)
if nargin < 2 || isempty(opts), opts = struct(); end
if nargin >= 3 && ~isempty(alg),   opts.alg = alg; end
if nargin >= 4 && ~isempty(N),     opts.N = N; end
if nargin >= 5 && ~isempty(maxFE), opts.maxFE = maxFE; end
t0 = tic;
P = RC3_setup(T, opts);
N = getopt(opts, 'N', 100);  maxFE = getopt(opts, 'maxFE', 20000);
verb = getopt(opts, 'verbose', true);

%% base case
[~, ~, R.out_base] = RC3_objectives(P.xbase, P);
R.F_base = R.out_base.F;  R.Flib_base = R.out_base.Flib;

%% solver
solver = lower(getopt(opts, 'solver', 'auto'));
if strcmp(solver, 'auto')
    if exist('platemo', 'file') == 2, solver = 'platemo'; else, solver = 'builtin'; end
end
if verb
    fprintf('\n[%s] T = %d | devices: %s | objectives: %s | m = %d | island: %s | D = %d | solver: %s\n', ...
            datestr(now, 'HH:MM:SS'), T, strjoin(P.devices, '+'), strjoin(P.objNames, ','), ...
            P.m, P.isl.mode, P.D, solver); %#ok<TNOW1,DATST>
end
RC3_fitness('reset', maxFE + 2 * N);
here = pwd;                                   % platemo() changes the current folder
addpath(fileparts(mfilename('fullpath')));    % keep the RC3_* files on the path
switch solver
    case 'platemo'
        alg = getopt(opts, 'alg', []);
        if isempty(alg)
            if P.m == 1, alg = @GA; elseif P.M <= 3, alg = @NSGAII; else, alg = @NSGAIII; end
        end
        % accepted forms: @LMOCSO, 'LMOCSO', {@FDV, 0.8, 0.4, 5} or {'FDV', 0.8, 0.4, 5}
        if ischar(alg) || isstring(alg), alg = str2func(char(alg)); end
        if iscell(alg) && (ischar(alg{1}) || isstring(alg{1})), alg{1} = str2func(char(alg{1})); end
        try
            [Dec, Obj0, Con0] = platemo('problem', {@RC3_PlatEMO, P}, 'algorithm', alg, ...
                                        'N', N, 'maxFE', maxFE);
        catch err
            cd(here);  rethrow(err);
        end
        cd(here);
        if iscell(alg), R.alg = func2str(alg{1}); else, R.alg = func2str(alg); end
    case 'builtin'
        rng(P.seed, 'twister');               % same stream as RC3_PlatEMO.Initialization
        [Dec, Obj0, Con0] = RC3_nsga2(P, N, maxFE, struct('verbose', verb));
        R.alg = 'RC3_nsga2';
    otherwise
        error('RC3_solve: unknown solver %s', solver);
end
R.hist = RC3_fitness('get');

%% re-evaluate the final population with the common evaluator
[Obj, Con, Dec] = RC3_evaluate(Dec, P);
rel = @(a, b) max([0; abs(a(:) - b(:)) ./ max(1, abs(b(:)))]);
R.check = struct('maxdObj', rel(Obj, Obj0), 'maxdCon', rel(Con, Con0), 'solver', solver);
if R.check.maxdObj > 1e-9 || R.check.maxdCon > 1e-9
    warning('RC3_solve:consistency', ['Re-evaluated objectives/constraints differ from the solver ' ...
            'values (max rel. diff %.3g / %.3g).'], R.check.maxdObj, R.check.maxdCon);
end

%% feasible set, non-dominated filter, best solution
CV   = sum(max(Con, 0), 2);
feas = CV <= 0;
R.nfeas = nnz(feas);
if ~any(feas)
    warning('RC3_solve:infeasible', ['No feasible solution: the least-violating one is ' ...
            'reported. Increase maxFE or relax the limits.']);
    feas = CV == min(CV);
end
Dec = Dec(feas, :);  Obj = Obj(feas, :);
[Obj, iu] = unique(Obj, 'rows');  Dec = Dec(iu, :);
if P.m > 1
    nd = true(size(Obj, 1), 1);
    for i = 1:size(Obj, 1)
        nd(i) = ~any(all(Obj <= Obj(i, :), 2) & any(Obj < Obj(i, :), 2));
    end
    Dec = Dec(nd, :);  Obj = Obj(nd, :);
    [~, o] = sort(Obj(:, 1));  R.Dec = Dec(o, :);  R.Obj = Obj(o, :);
    [R.ibest, R.score] = RC3_select_best(R.Obj);
else
    [~, o] = sort(Obj(:, 1));  R.Dec = Dec(o, :);  R.Obj = Obj(o, :);
    R.ibest = 1;  R.score = [];
end
R.x_best = R.Dec(R.ibest, :);

%% optional local refinement of the continuous set-points (m = 1)
if P.m == 1 && getopt(opts, 'refine', false)
    R.refine = RC3_refine(P, R.x_best);
    if R.refine.improved, R.x_best = R.refine.x; end
end
[~, ~, R.out_best] = RC3_objectives(R.x_best, P);
R.F_best = R.out_best.F;  R.Flib_best = R.out_best.Flib;
np = size(R.Dec, 1);
R.Fuser_pareto = zeros(np, P.M);  R.Fall_pareto = zeros(np, P.nAll);
R.Flib_pareto  = zeros(np, numel(P.objLib));
for i = 1:np
    [Fu, ~, ~, Fa] = RC3_objectives(R.Dec(i, :), P);
    R.Fuser_pareto(i, :) = Fu;  R.Fall_pareto(i, :) = Fa;
    R.Flib_pareto(i, :)  = RC3_objfun(Fa, P, P.objLib);
end

%% MATPOWER validation of the best solution
R.val = [];
if getopt(opts, 'validate', true) && exist('runpf', 'file') == 2
    try
        R.val = RC3_validate_solution(P, R.out_best);
    catch err
        warning('RC3_solve:validate', 'MATPOWER validation skipped: %s', err.message);
    end
end
R.P = P;  R.opts = opts;  R.solver = solver;  R.N = N;  R.maxFE = maxFE;
R.time = toc(t0);
if verb
    fprintf('  done in %.1f s: %d feasible in the final population, %d non-dominated / best\n', ...
            R.time, R.nfeas, size(R.Dec, 1));
end
end

function v = getopt(s, name, default)
if isfield(s, name) && ~isempty(s.(name)), v = s.(name); else, v = default; end
end
