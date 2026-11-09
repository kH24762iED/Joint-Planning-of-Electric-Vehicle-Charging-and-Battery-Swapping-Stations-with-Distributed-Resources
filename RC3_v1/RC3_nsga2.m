function [Dec, Obj, Con, hist] = RC3_nsga2(P, N, maxFE, opt)
%RC3_NSGA2  Built-in mixed-integer NSGA-II (Deb et al. 2002) with
%constrained domination (Deb 2000); used when PlatEMO is not on the path
%(opts.solver = 'builtin'). With one objective (m = 1) it is an elitist
%real/integer/binary-coded genetic algorithm.
%
%   [Dec, Obj, Con, hist] = RC3_nsga2(P, N, maxFE)
%   [...] = RC3_nsga2(P, N, maxFE, opt)   opt.pc (0.9) opt.etac (20) opt.etam (20)
%                                        opt.verbose (true)
%   Variables: continuous -> SBX + polynomial mutation; integer -> the same,
%   rounded; binary -> uniform crossover + bit flip. The initial population
%   is RC3_init_pop (the same as RC3_PlatEMO for the same seed) and every
%   population is evaluated with RC3_evaluate (repair + MP-Opt-Model costs
%   and constraints), the same path as PlatEMO. Returns the final population.
%   Call rng(P.seed) before (RC3_solve does) for reproducible runs.
if nargin < 4, opt = struct(); end
pc   = getf(opt, 'pc', 0.9);  etac = getf(opt, 'etac', 20);  etam = getf(opt, 'etam', 20);
verb = getf(opt, 'verbose', true);
N    = 2 * ceil(N / 2);
D    = P.D;  lo = P.lower;  up = P.upper;  enc = P.enc;
isB  = enc == 4;  isR = ~isB;  isI = enc == 2;
pm   = 1 / D;

%% common initial population (base case + local and global seeds)
[Obj, Con, Dec] = RC3_evaluate(RC3_init_pop(P, N), P);  FE = N;
[Front, Crowd] = rank_pop(Obj, Con);
hist = struct('FE', [], 'best', [], 'nfeas', []);
gen = 0;
while FE < maxFE
    gen = gen + 1;
    %% mating: binary tournament on (front, -crowding)
    M  = tournament(N, Front, Crowd);
    P1 = Dec(M(1:2:end), :);  P2 = Dec(M(2:2:end), :);
    n2 = size(P1, 1);
    %% SBX on real/integer variables
    mu   = rand(n2, D);
    beta = zeros(n2, D);
    beta(mu <= 0.5) = (2 * mu(mu <= 0.5)).^(1 / (etac + 1));
    beta(mu > 0.5)  = (2 - 2 * mu(mu > 0.5)).^(-1 / (etac + 1));
    beta = beta .* (-1).^randi([0 1], n2, D);
    beta(rand(n2, D) < 0.5) = 1;
    beta(repmat(rand(n2, 1) > pc, 1, D)) = 1;
    O1 = (P1 + P2) / 2 + beta .* (P1 - P2) / 2;
    O2 = (P1 + P2) / 2 - beta .* (P1 - P2) / 2;
    %% uniform crossover on binary variables
    sw = rand(n2, D) < 0.5 & repmat(rand(n2, 1) <= pc, 1, D);
    B1 = P1;  B2 = P2;  B1(sw) = P2(sw);  B2(sw) = P1(sw);
    O1(:, isB) = B1(:, isB);  O2(:, isB) = B2(:, isB);
    Off = [O1; O2];
    %% polynomial mutation (real/integer) and bit flip (binary)
    no   = size(Off, 1);
    Lo   = repmat(lo, no, 1);  Up = repmat(up, no, 1);
    mut  = rand(no, D) < pm & repmat(isR, no, 1) & Up > Lo;
    mu   = rand(no, D);
    Off  = min(max(Off, Lo), Up);
    rng_ = Up - Lo;  rng_(rng_ == 0) = 1;
    d1   = (Off - Lo) ./ rng_;  d2 = (Up - Off) ./ rng_;
    t1   = mut & mu <= 0.5;  t2 = mut & mu > 0.5;
    Off(t1) = Off(t1) + rng_(t1) .* ((2 .* mu(t1) + (1 - 2 .* mu(t1)) .* (1 - d1(t1)).^(etam + 1)).^(1 / (etam + 1)) - 1);
    Off(t2) = Off(t2) + rng_(t2) .* (1 - (2 .* (1 - mu(t2)) + 2 .* (mu(t2) - 0.5) .* (1 - d2(t2)).^(etam + 1)).^(1 / (etam + 1)));
    flip = rand(no, D) < pm & repmat(isB, no, 1);
    Off(flip) = 1 - Off(flip);
    iI = repmat(isI, no, 1) & rand(no, D) < pm;                   % extra integer jumps
    Off(iI) = Lo(iI) + floor(rand(nnz(iI), 1) .* (Up(iI) - Lo(iI) + 1));
    nOff = min(no, maxFE - FE);  Off = Off(1:nOff, :);
    [OObj, OCon, Off] = RC3_evaluate(Off, P);  FE = FE + nOff;
    %% environmental selection
    Dec = [Dec; Off];  Obj = [Obj; OObj];  Con = [Con; OCon];
    [Front, Crowd] = rank_pop(Obj, Con);
    [~, ord] = sortrows([Front, -Crowd]);
    keep = ord(1:N);
    Dec = Dec(keep, :);  Obj = Obj(keep, :);  Con = Con(keep, :);
    Front = Front(keep);  Crowd = Crowd(keep);
    %% history
    feas = all(Con <= 0, 2);
    if any(feas), b = min(Obj(feas, :), [], 1); else, b = nan(1, size(Obj, 2)); end
    hist.FE(end+1, 1) = FE;  hist.best(end+1, :) = b;  hist.nfeas(end+1, 1) = nnz(feas);
    if verb && (mod(gen, 10) == 0 || FE >= maxFE)
        fprintf('  gen %4d  FE %6d  feasible %3d  best %s\n', gen, FE, nnz(feas), mat2str(b, 5));
    end
end
end

% =========================================================================================
function [Front, Crowd] = rank_pop(Obj, Con)
% constrained non-dominated sorting: feasible solutions first (by Pareto
% fronts), then infeasible ones by total constraint violation
n  = size(Obj, 1);  CV = sum(max(Con, 0), 2);  feas = CV <= 0;
Front = inf(n, 1);  Crowd = zeros(n, 1);
idx = find(feas);
if ~isempty(idx)
    Front(idx) = nd_sort(Obj(idx, :));
end
mf = max([0; Front(feas)]);
idx = find(~feas);
if ~isempty(idx)
    [~, o] = sort(CV(idx));
    Front(idx(o)) = mf + (1:numel(o))';
end
for f = unique(Front(feas))'
    k = find(Front == f);
    Crowd(k) = crowding(Obj(k, :));
end
end

function F = nd_sort(Obj)
n = size(Obj, 1);  F = inf(n, 1);  rest = true(n, 1);  f = 0;
while any(rest)
    f = f + 1;  r = find(rest);
    dom = false(numel(r), 1);
    for a = 1:numel(r)
        oa = Obj(r(a), :);
        dom(a) = any(all(Obj(r, :) <= oa, 2) & any(Obj(r, :) < oa, 2));
    end
    F(r(~dom)) = f;  rest(r(~dom)) = false;
end
end

function C = crowding(Obj)
[n, m] = size(Obj);  C = zeros(n, 1);
if n <= 2, C(:) = Inf; return; end
for j = 1:m
    [v, o] = sort(Obj(:, j));
    C(o([1 end])) = Inf;
    rg = v(end) - v(1);
    if rg > 0
        C(o(2:end-1)) = C(o(2:end-1)) + (v(3:end) - v(1:end-2)) / rg;
    end
end
end

function M = tournament(N, Front, Crowd)
a = randi(numel(Front), N, 1);  b = randi(numel(Front), N, 1);
better = Front(a) < Front(b) | (Front(a) == Front(b) & Crowd(a) > Crowd(b));
M = b;  M(better) = a(better);
end

function v = getf(s, name, default)
if isfield(s, name) && ~isempty(s.(name)), v = s.(name); else, v = default; end
end
