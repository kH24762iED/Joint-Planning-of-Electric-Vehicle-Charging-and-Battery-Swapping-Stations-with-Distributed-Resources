function X = RC3_repair(X, P)
%RC3_REPAIR  Deterministic, idempotent repair of decision vectors (one per row
%of X). It is the first step of every evaluation (RC3_evaluate, and
%RC3_PlatEMO.CalDec), and the repaired vector is stored by the solver
%(Lamarckian repair), so PlatEMO and RC3_nsga2 always evaluate the same
%vector in the same way. A vector that is already valid is not changed.
%
%   1  bounds and integrality of every variable
%   2  radial ONR configuration in every slot (loop encoding: keep the chosen
%      open branches if they give a spanning tree, otherwise a depth-first
%      search over the loops takes, for each loop, the candidate nearest to the
%      chosen one that keeps the network connected and backtracks when a later
%      loop has no admissible branch; binary encoding: spanning tree that keeps
%      as many closed switches as possible)
%   3  PCC transitions <= island.max_trans (flip the shortest free run of
%      slots that reduces the number of transitions)
%   4  minimum number of units  sum_k y_dk >= nmin_d
%   5  one unit per bus / branch (a duplicate site moves to the nearest free
%      candidate)
%   6  DG set-points <= available renewable power, installed DG capacity
%      <= pen * total load (set-points scaled down)
%   7  tap changes per slot within the OLTC / transformer ramp limits
%   8  DR energy neutrality (projection) and curtailment budget (scaling)
%   9  BSS rule thresholds  alpha_dis(t) >= alpha_ch(t) + gap
lo = P.lower;  up = P.upper;
X  = min(max(X, lo), up);                    % NaN -> lower bound
ii = P.enc ~= 1;
X(:, ii) = round(X(:, ii));
for i = 1:size(X, 1)
    x = X(i, :);
    x = rep_onr(x, P);
    x = rep_pcc(x, P);
    x = rep_count(x, P);
    x = rep_sites(x, P);
    x = rep_dg(x, P);
    x = rep_ramp(x, P, 'TAP_oltc', 1, P.oltc.ramp);
    x = rep_ramp(x, P, 'TX_tap', P.tx.nmax, P.tx.ramp);
    x = rep_dr(x, P);
    x = rep_bss(x, P);
    X(i, :) = x;
end
end

% =========================================================================================
function x = rep_onr(x, P)
T = P.T;
if isfield(P.bidx, 'ONR_open')
    nt = numel(P.tie);  V = reshape(x(P.bidx.ONR_open), nt, T);
    for t = 1:T, V(:, t) = radial_loop(V(:, t), P); end
    x(P.bidx.ONR_open) = V(:)';
elseif isfield(P.bidx, 'ONR_sw')
    V = reshape(x(P.bidx.ONR_sw), P.nl, T);
    for t = 1:T, V(:, t) = radial_tree(V(:, t), P); end
    x(P.bidx.ONR_sw) = V(:)';
end
end

function pos = radial_loop(pos, P)
% loop encoding: one open branch per fundamental loop. A radial vector is
% kept; otherwise a depth-first search over the loops (candidates of each
% loop in order of distance from the chosen position: 0, +1, -1, +2, ...)
% opens one branch per loop while the network stays connected, and
% backtracks when a later loop has no admissible branch. The base
% configuration (tie switches open) is radial, so a solution always exists.
closed = true(P.nl, 1);
for k = 1:numel(P.loops), closed(P.loops{k}(pos(k))) = false; end
if connected(P, closed) && nnz(~closed) == numel(P.loops), return; end   % already radial
[ok, new] = loop_dfs(1, true(P.nl, 1), pos, pos, P);
if ok, pos = new; end
end

function [ok, new] = loop_dfs(k, closed, pos, new, P)
ok = false;
if k > numel(P.loops), ok = true;  return; end
L = P.loops{k};  n = numel(L);  p0 = pos(k);
d = [0, reshape([1:n; -(1:n)], 1, [])];                  % 0, +1, -1, +2, -2, ...
cand = p0 + d;  cand = cand(cand >= 1 & cand <= n);
for q = cand
    b = L(q);
    if ~closed(b), continue; end
    closed(b) = false;
    if connected(P, closed)
        new(k) = q;
        [ok, new] = loop_dfs(k + 1, closed, pos, new, P);
        if ok, return; end
    end
    closed(b) = true;
end
end

function s = radial_tree(s, P)
% binary encoding: Kruskal spanning tree preferring the closed switches
closed = s(:) > 0.5;
if connected(P, closed) && nnz(closed) == P.nb - 1, return; end
order = [find(closed); find(~closed)];
root  = 1:P.nb;  sel = false(P.nl, 1);
for e = order'
    a = findroot(root, P.fb(e));  b = findroot(root, P.tb(e));
    if a ~= b, root(a) = b;  sel(e) = true; end
end
s = double(sel);
end

function r = findroot(root, a)
r = a;
while root(r) ~= r, r = root(r); end
end

function tf = connected(P, closed)
[~, ~, ~, ~, nReach] = RC3_tree(P.fb, P.tb, closed, P.nb);
tf = nReach == P.nb;
end

function x = rep_pcc(x, P)
if ~isfield(P.bidx, 'PCC_u'), return; end
free = isnan(P.ugForced);
ug = P.ugForced;  u = x(P.bidx.PCC_u);  ug(free) = u(free);
n  = ntrans(ug);
while n > P.isl.max_trans
    d  = [true, diff(ug) ~= 0];  st = find(d);  en = [st(2:end) - 1, numel(ug)];
    best = [];  bestn = n;  bestlen = Inf;
    for r = 1:numel(st)
        sl = st(r):en(r);
        if ~all(free(sl)), continue; end
        ug2 = ug;  ug2(sl) = 1 - ug2(sl);  n2 = ntrans(ug2);
        if n2 < n && (n2 < bestn || (n2 == bestn && numel(sl) < bestlen))
            best = sl;  bestn = n2;  bestlen = numel(sl);
        end
    end
    if isempty(best), break; end
    ug(best) = 1 - ug(best);  n = bestn;
end
u(free) = ug(free);  x(P.bidx.PCC_u) = u;
end

function n = ntrans(ug)
n = sum(abs(diff([1, ug])));
end

function x = rep_count(x, P)
cnt = {'DG_y','dg'; 'SC_y','sc'; 'TX_y','tx'; 'ST_y','st'; 'UP_y','up'; 'EV_y','ev'; 'BSS_y','bss'};
for k = 1:size(cnt, 1)
    if ~isfield(P.bidx, cnt{k, 1}), continue; end
    ix = P.bidx.(cnt{k, 1});  y = x(ix);  nmin = P.(cnt{k, 2}).nmin;
    if sum(y) < nmin
        z = find(y < 0.5);  y(z(1:min(numel(z), nmin - sum(y)))) = 1;
        x(ix) = y;
    end
end
end

function x = rep_sites(x, P)
bus = {'DG_y','DG_bus','dg'; 'SC_y','SC_bus','sc'; 'ST_y','ST_bus','st'; ...
       'EV_y','EV_bus','ev'; 'BSS_y','BSS_bus','bss'};
for k = 1:size(bus, 1)
    if ~isfield(P.bidx, bus{k, 1}), continue; end
    x = place(x, P.bidx.(bus{k, 1}), P.bidx.(bus{k, 2}), P.(bus{k, 3}).cand, false(P.nb, 1));
end
usedBr = false(P.nl, 1);                                 % TX and UPFC share the branches
if isfield(P.bidx, 'TX_y'), [x, usedBr] = place(x, P.bidx.TX_y, P.bidx.TX_br, P.tx.cand, usedBr); end
if isfield(P.bidx, 'UP_y'), x = place(x, P.bidx.UP_y, P.bidx.UP_br, P.up.cand, usedBr); end
end

function [x, used] = place(x, iy, is, cand, used)
% installed units in index order; a unit on a used site moves to the
% nearest candidate (in the candidate list) whose site is free
n = numel(cand);
for u = 1:numel(iy)
    if x(iy(u)) < 0.5, continue; end
    c = x(is(u));
    if used(cand(c))
        for dlt = 1:n
            c2 = [c + dlt, c - dlt];  c2 = c2(c2 >= 1 & c2 <= n);
            c2 = c2(~used(cand(c2)));
            if ~isempty(c2), c = c2(1);  break; end
        end
        x(is(u)) = c;
    end
    used(cand(c)) = true;
end
end

function x = rep_dg(x, P)
if ~isfield(P.bidx, 'DG_P'), return; end
n = P.dg.nmax;  T = P.T;
y  = x(P.bidx.DG_y) > 0.5;  y = y(:);
ty = P.dg.allowed(x(P.bidx.DG_type));  ty = ty(:);
Pm  = reshape(x(P.bidx.DG_P), n, T);
prf = P.dg.avail(ty, :);
isR = ~strcmp(P.dg.profile(ty), 'const');  isR = isR(:);
Pm(isR, :) = min(Pm(isR, :), P.dg.Pmax * prf(isR, :));  % set-point <= available power
ratio = Pm ./ max(prf, eps);  ratio(prf < 1e-3) = 0;
cap = min(max(ratio, [], 2), P.dg.Pmax);  cap(~y) = 0;   % capacity, as in RC3_decode
lim = P.dg.pen * P.Ptot;
if sum(cap) > lim
    Pm(y, :) = Pm(y, :) * (lim / sum(cap) * (1 - 1e-9));
end
x(P.bidx.DG_P) = Pm(:)';
end

function x = rep_ramp(x, P, name, n, R)
if ~isfield(P.bidx, name) || P.T < 2 || ~isfinite(R), return; end
V = reshape(x(P.bidx.(name)), n, P.T);
for t = 2:P.T
    V(:, t) = min(max(V(:, t), V(:, t-1) - R), V(:, t-1) + R);
end
x(P.bidx.(name)) = V(:)';
end

function x = rep_dr(x, P)
if isfield(P.bidx, 'DR_shift')
    ix = P.bidx.DR_shift;
    x(ix) = RC3_dr_project(x(ix), P.dr_w, P.lower(ix), P.upper(ix));
end
if isfield(P.bidx, 'DR_il')
    ix = P.bidx.DR_il;  B = P.drm.il_frac * sum(P.dr_w);
    e = x(ix) * P.dr_w';
    if e > B * (1 + 1e-12), x(ix) = x(ix) * (B / e); end
end
end

function x = rep_bss(x, P)
if ~isfield(P.bidx, 'BSS_ach'), return; end
ia = P.bidx.BSS_ach;  id = P.bidx.BSS_adis;
x(ia) = min(x(ia), P.upper(id) - P.bss.gap);
x(id) = min(max(x(id), x(ia) + P.bss.gap), P.upper(id));
end
