function r = RC3_pf_fbs(net, V0, S0, Ysh, np, nq, pv, up, tol, maxit, Vinit)
%RC3_PF_FBS  Forward-backward sweep load flow of one radial component in the
%matrix (BIBC/BCBV) form of Teng (2003), extended with
%   - ideal tap-changing transformers (OLTC / in-line regulators): all
%     quantities are referred to the root side, V = c .* Vt, I = It ./ c,
%     Zt = Z ./ zc.^2 (exact for real tap ratios, MATPOWER branch model)
%   - UPFC series voltage sources vse = s * r * exp(j*theta) * V_from (Noroozian
%     et al. 1997), the shunt converter drawing the series active power
%   - voltage-controlled (PV) nodes for D-STATCOMs in 'pv' mode with
%     reactive limits (compensation method, Shirmohammadi et al. 1988)
%
%   net : .n nodes (local; node 1 = root) .K path matrix .par local parent
%         .Zt referred branch impedance feeding each node (p.u., Zt(1) = 0)
%         .c  cumulative tap ratio of each node
%   V0  : root voltage (p.u., complex)
%   S0  : load at 1 p.u. voltage (p.u., load convention), n x 1
%   Ysh : shunt admittance (p.u.), n x 1
%   np, nq : load voltage exponents  S = P0*|V|^np + j*Q0*|V|^nq
%   pv  : .node .Vset .Qmin .Qmax  (injection, p.u.; unique nodes)
%   up  : .node (node fed by the UPFC branch) .ft (from-terminal node)
%         .k = s*r*exp(j*theta) (s = +1 if the from-bus is the parent) .Qsh
%   Vinit : optional warm start (referred voltages)
%
%   r.V actual voltages, r.Vt referred voltages, r.Ib referred current of
%   the branch feeding each node, r.Ss / r.Sr sending / receiving branch
%   power, r.loss complex losses, r.Sroot power delivered by the root
%   source (root load included), r.Q PV-node injections, r.Sse UPFC series
%   power, r.it iterations, r.success.
n  = net.n;  K = net.K;  Zt = net.Zt;  c = net.c;  c2 = c.^2;
P0 = real(S0);  Q0 = imag(S0);
if nargin > 10 && ~isempty(Vinit) && numel(Vinit) == n, V = Vinit; else, V = V0 * ones(n, 1); end
npv = numel(pv.node);  Q = zeros(npv, 1);  act = true(npv, 1);
if npv
    Kp  = K(pv.node, :);
    Xpv = full((Kp .* imag(Zt).') * Kp.');           % dVt/dQ of the PV nodes
    Xpv = c(pv.node) .* Xpv;                          % actual voltage sensitivity
end
nup = numel(up.node);  Pse = zeros(nup, 1);  vse = zeros(n, 1);
r.success = false;
for it = 1:maxit
    I   = load_current(V, c, c2, P0, Q0, np, nq, Ysh, pv.node, Q, up, Pse, n);
    Ib  = K.' * I;                                    % backward sweep (BIBC)
    if nup, vse(up.node) = up.k .* V(up.ft); end
    Vn  = V0 - K * (Zt .* Ib - vse);                  % forward sweep (BCBV)
    err = max(abs(Vn - V));  V = Vn;
    if nup
        Psn = real(vse(up.node) .* conj(Ib(up.node)));
        err = max(err, max(abs(Psn - Pse)));  Pse = Psn;
    end
    if npv                                            % PV-node compensation
        e   = pv.Vset - c(pv.node) .* abs(V(pv.node));
        act = act | (Q >= pv.Qmax & e < 0) | (Q <= pv.Qmin & e > 0);
        dq  = zeros(npv, 1);
        if any(act), dq(act) = Xpv(act, act) \ e(act); end
        Qn  = min(max(Q + dq, pv.Qmin), pv.Qmax);
        err = max([err; abs(Qn - Q); abs(e(act))]);
        act = act & Qn > pv.Qmin & Qn < pv.Qmax;
        Q   = Qn;
    end
    if err < tol, r.success = true; break; end
end
I  = load_current(V, c, c2, P0, Q0, np, nq, Ysh, pv.node, Q, up, Pse, n);
Ib = K.' * I;
if nup, vse(up.node) = up.k .* V(up.ft); end
pp = net.par;  pp(1) = 1;
r.Vt   = V;  r.V = c .* V;  r.Ib = Ib;  r.Q = Q;  r.it = it;
r.Ss   = V(pp) .* conj(Ib);  r.Ss(1) = 0;            % sending end (power invariant)
r.Sr   = V .* conj(Ib);      r.Sr(1) = 0;            % receiving end
r.loss = sum(abs(Ib).^2 .* Zt);                       % complex losses (p.u.)
r.Sse  = zeros(nup, 1);
if nup, r.Sse = vse(up.node) .* conj(Ib(up.node)); end
r.Sroot = V0 * conj(I(1) + sum(Ib(net.par == 1)));    % root source (incl. root load)
end

function I = load_current(V, c, c2, P0, Q0, np, nq, Ysh, pvn, Q, up, Pse, n)
% referred current drawn at each node: voltage-dependent load + shunt
% - PV injection + UPFC shunt converter (draws Pse, injects Qsh)
Vm = c .* abs(V);
S  = complex(P0 .* Vm.^np, Q0 .* Vm.^nq);
if ~isempty(pvn), S(pvn) = S(pvn) - 1j * Q; end
if ~isempty(up.node)
    S = S + complex(accumarray(up.ft(:), Pse(:), [n 1]), ...
                    -accumarray(up.ft(:), up.Qsh(:), [n 1]));
end
I = conj(S ./ V) + c2 .* Ysh .* V;
end
