function C = RC3_check_solvers(T, opts, nTest)
%RC3_CHECK_SOLVERS  Check that PlatEMO and the built-in solver evaluate any
%decision vector identically, and compare one run of each.
%
%   C = RC3_check_solvers()                    default case (T = 24)
%   C = RC3_check_solvers(T, opts, nTest)
%
%   1  Evaluation test: nTest decision vectors (the common initial
%      population plus raw random vectors, i.e. not repaired) are evaluated
%        a) by the built-in path        RC3_evaluate(X, P)
%        b) by PlatEMO's own path       RC3_PlatEMO(...).Evaluation(X)
%      The repaired vectors, objectives and constraint violations must be
%      identical (C.eval: maximum absolute differences, expected 0).
%   2  Idempotence: evaluating the repaired vectors again gives the same
%      vectors and values (C.idem).
%   3  Initial population: RC3_PlatEMO.Initialization and RC3_nsga2 start from
%      the same population for the same seed (C.init).
%   4  Runs: RC3_solve with solver 'builtin' and 'platemo', same seed, N and
%      maxFE (C.runs: feasible solutions, best objectives, re-evaluation
%      differences). The two algorithms are different implementations, so
%      their final populations differ, but both start from the same feasible
%      population and use the same evaluation.
%   Steps 1b, 3 and 4 need PlatEMO on the MATLAB path.
if nargin < 1 || isempty(T), T = 24; end
if nargin < 2 || isempty(opts)
    opts = struct('devices', {{'ONR','DG','SC','TAP','DSTATCOM','BSS','DRM','ISL'}}, ...
                  'objectives', {{'f1','f2','f5'}}, 'island_mode', 'outage', ...
                  'N', 40, 'maxFE', 2000, 'seed', 1);
end
if nargin < 3 || isempty(nTest), nTest = 40; end
hasPE = exist('platemo', 'file') == 2 && exist('PROBLEM', 'class') == 8;
P = RC3_setup(T, opts);
C = struct('eval', [], 'idem', [], 'init', [], 'runs', []);

%% 1) the same vectors through both evaluation paths
rng(P.seed, 'twister');
X0 = RC3_init_pop(P, ceil(nTest / 2));
Xr = P.lower + rand(nTest - size(X0, 1), P.D) .* (P.upper - P.lower);   % raw, not repaired
X  = [X0; Xr];
[O1, C1, X1] = RC3_evaluate(X, P);
fprintf('Built-in evaluation: %d of %d vectors feasible\n', nnz(all(C1 <= 0, 2)), size(X, 1));
if hasPE
    Pr  = RC3_PlatEMO('parameter', {P}, 'N', size(X, 1), 'maxFE', 1e6);
    Pop = Pr.Evaluation(X);
    C.eval = struct('dDec', max(abs(Pop.decs - X1), [], 'all'), ...
                    'dObj', max(abs(Pop.objs - O1), [], 'all'), ...
                    'dCon', max(abs(Pop.cons - C1), [], 'all'));
    fprintf('PlatEMO vs built-in evaluation: max |dDec| %.3g, |dObj| %.3g, |dCon| %.3g\n', ...
            C.eval.dDec, C.eval.dObj, C.eval.dCon);
else
    fprintf('PlatEMO not on the path: evaluation test b) skipped.\n');
end

%% 2) idempotence of repair + evaluation
[O2, C2, X2] = RC3_evaluate(X1, P);
C.idem = struct('dDec', max(abs(X2 - X1), [], 'all'), 'dObj', max(abs(O2 - O1), [], 'all'), ...
                'dCon', max(abs(C2 - C1), [], 'all'));
fprintf('Re-evaluation of repaired vectors: max |dDec| %.3g, |dObj| %.3g, |dCon| %.3g\n', ...
        C.idem.dDec, C.idem.dObj, C.idem.dCon);

%% 3) common initial population
if isfield(opts, 'N') && ~isempty(opts.N), N = opts.N; else, N = 100; end
rng(P.seed, 'twister');  Xb = RC3_init_pop(P, 2 * ceil(N / 2));          % as RC3_nsga2
if hasPE
    Pr  = RC3_PlatEMO('parameter', {P}, 'N', 2 * ceil(N / 2), 'maxFE', 1e6);
    rng('shuffle');                                                       % as platemo()
    Pop = Pr.Initialization();
    C.init = max(abs(Pop.decs - Xb), [], 'all');
    fprintf('Initial populations PlatEMO vs built-in: max |dDec| %.3g (%d feasible of %d)\n', ...
            C.init, nnz(all(Pop.cons <= 0, 2)), size(Xb, 1));
end

%% 4) one run of each solver
solvers = {'builtin'};
if hasPE, solvers{end+1} = 'platemo'; end
rows = cell(numel(solvers), 1);  vals = [];
for k = 1:numel(solvers)
    o = opts;  o.solver = solvers{k};  o.validate = false;  o.verbose = false;
    R = RC3_solve(T, o);
    rows{k} = solvers{k};
    vals(k, :) = [R.nfeas, size(R.Dec, 1), R.F_best, R.check.maxdObj, R.check.maxdCon, R.time]; %#ok<AGROW>
end
vn = [{'FeasibleFinalPop', 'nNonDominated'}, strcat('best_', P.objNames), {'ReEvalDiffObj', 'ReEvalDiffCon', 'Time_s'}];
C.runs = array2table(vals, 'RowNames', rows, 'VariableNames', vn);
disp(C.runs);
end
