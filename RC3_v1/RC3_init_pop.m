function X = RC3_init_pop(P, N)
%RC3_INIT_POP  Common initial population of BOTH solvers (RC3_PlatEMO and
%RC3_nsga2), so that PlatEMO and the built-in solver start from exactly the
%same N decision vectors for the same seed.
%
%   X = RC3_init_pop(P, N)      uses the global random stream: call
%                              rng(P.seed) first (RC3_solve and
%                              RC3_PlatEMO.Initialization do this)
%
%   row 1        the base case P.xbase (feasible reference solution)
%   even rows    local seeds: about 30 % of the decision blocks of the base
%                case are re-drawn at random
%   odd rows     global seeds: every block is drawn at random
%   Distributions: installation flags 0/1 with probability 0.5; sites,
%   types and chargers uniform; continuous set-points uniform in their
%   bounds (DG set-points scaled to the penetration limit, DR shares to
%   half their range); taps start in [-2, 2] and change by at most one step
%   per slot; the PCC (optimal islanding) opens for at most one contiguous
%   period of up to T/6 slots. Every row is then repaired (RC3_repair), so
%   all ONR configurations are radial and the DG penetration, tap ramp,
%   siting and DR constraints hold.
X  = repmat(P.xbase, N, 1);
B  = P.blocks;
for i = 2:N
    local = mod(i, 2) == 0;
    for j = 1:numel(B)
        if local && rand > 0.3, continue; end
        b = B(j);  n = numel(b.idx);  lo = b.lo;  up = b.up;
        switch b.name
            case {'TAP_oltc', 'TX_tap'}
                nu = b.dims(1);  T = b.dims(2);
                V = zeros(nu, T);  V(:, 1) = randi([-2 2], nu, 1);
                for t = 2:T, V(:, t) = V(:, t-1) + randi([-1 1], nu, 1); end
                v = min(max(V(:)', lo), up);
            case 'PCC_u'
                v = ones(1, n);
                if rand < 0.5
                    L = randi([1, max(1, round(n / 6))]);  s = randi([1, n - L + 1]);
                    v(s:s+L-1) = 0;
                end
            case 'ISL_cut'
                v = (rand(1, n) < 0.5) .* randi([1, P.nl], 1, n);
            case 'DG_P'
                sc = min(1, P.dg.pen * P.Ptot / max(P.dg.nmax * P.dg.Pmax, eps));
                v  = lo + rand(1, n) .* (up - lo) * sc;
            case {'DR_shift', 'DR_il'}
                v = 0.5 * (lo + rand(1, n) .* (up - lo));
            otherwise
                if all(b.enc == 4)
                    v = double(rand(1, n) < 0.5);
                elseif all(b.enc == 2)
                    v = lo + floor(rand(1, n) .* (up - lo + 1));
                else
                    v = lo + rand(1, n) .* (up - lo);
                end
        end
        X(i, b.idx) = v;
    end
end
X = RC3_repair(X, P);
end
