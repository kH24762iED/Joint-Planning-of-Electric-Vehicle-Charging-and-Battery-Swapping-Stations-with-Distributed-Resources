function V = RC3_validate_solution(P, out, slots, mpopt)
%RC3_VALIDATE_SOLUTION  Re-solve time slots of an evaluated solution with
%MATPOWER's runpf (Newton-Raphson by default) and compare with the
%forward-backward sweep results of RC3_objectives.
%
%   V = RC3_validate_solution(P, out)              all slots
%   V = RC3_validate_solution(P, out, slots, mpopt)
%
%   UPFCs are added through their injection model (Norton equivalent of
%   the series source, Noroozian et al. 1997): with vse = r*exp(j*theta)*V_f,
%       S_f = -V_f * conj(vse/Z),   S_t = V_t * conj(vse/Z)
%   plus the shunt converter (draws Re(vse*conj(I_ft)), injects Qsh); the
%   injections depend on the voltages, so runpf is repeated until they
%   converge.
%
%   V is a table with, per slot: max |V_FBS - V_MP| (p.u.), losses (kW),
%   grid import (MW), reference-generator mismatch (MW), runpf iterations,
%   UPFC fixed-point iterations and the success flag.
if nargin < 3 || isempty(slots), slots = 1:P.T; end
if nargin < 4 || isempty(mpopt)
    mpopt = mpoption('verbose', 0, 'out.all', 0, 'pf.tol', 1e-10, ...
                     'pf.enforce_q_lims', double(strcmp(P.st.mode, 'pv')));
    np = P.lim.zip_np;  nq = P.lim.zip_nq;
    if np ~= 0 || nq ~= 0
        zp = [0 0 0];  zp(min(max(round(np), 0), 2) + 1) = 1;
        zq = [0 0 0];  zq(min(max(round(nq), 0), 2) + 1) = 1;
        mpopt = mpoption(mpopt, 'exp.sys_wide_zip_loads.pw', zp, ...   % [P I Z] shares
                                'exp.sys_wide_zip_loads.qw', zq);
    end
end
ns = numel(slots);
[dV, lossF, lossM, PsF, PsM, dRef, its, upIt] = deal(nan(ns, 1));  ok = false(ns, 1);
for i = 1:ns
    t = slots(i);
    [mpc0, info] = RC3_mpc_solution(P, out, t);
    U = info.upfc;  nU = numel(U.br);
    Vb = out.V(:, t);  Vb(isnan(Vb)) = 1;                % start from the FBS voltages
    for k = 1:30
        mpc = mpc0;
        for u = 1:nU                                     % UPFC injection model
            b = U.br(u);
            if mpc.branch(b, 11) == 0, continue; end
            f = P.fb(b);  tt = P.tb(b);  Z = P.Zbr(b);
            vse = U.r(u) * exp(1j * U.th(u)) * Vb(f);
            Is  = vse / Z;  Ift = (Vb(f) + vse - Vb(tt)) / Z;
            Sf  = -Vb(f) * conj(Is);  St = Vb(tt) * conj(Is);   % injections (p.u.)
            Pse = real(vse * conj(Ift));
            mpc.bus(f, 3)  = mpc.bus(f, 3)  + (-real(Sf) + Pse) * P.baseMVA;
            mpc.bus(f, 4)  = mpc.bus(f, 4)  + (-imag(Sf) - U.Qsh(u) / P.baseMVA) * P.baseMVA;
            mpc.bus(tt, 3) = mpc.bus(tt, 3) - real(St) * P.baseMVA;
            mpc.bus(tt, 4) = mpc.bus(tt, 4) - imag(St) * P.baseMVA;
        end
        [res, success] = runpf(mpc, mpopt);
        Vn = res.bus(:, 8) .* exp(1j * res.bus(:, 9) * pi / 180);
        chg = max(abs(Vn(info.energized) - Vb(info.energized)));
        Vb = Vn;
        if ~success || nU == 0 || chg < 1e-9, break; end
    end
    e = info.energized;
    ok(i)    = success;  upIt(i) = k;
    if isfield(res, 'iterations'), its(i) = res.iterations; end
    dV(i)    = max(abs(abs(out.V(e, t)) - res.bus(e, 8)));
    lossF(i) = out.loss_kW(t);
    lossM(i) = 1000 * sum(real(get_losses(res)));
    PsF(i)   = out.Psub(t);
    PsM(i)   = res.gen(1, 2) * (res.gen(1, 8) > 0);
    if out.ug(t) < 0.5                                   % reference generators of the islands
        ng = size(res.gen, 1);  nref = numel(out.roots{t});
        dRef(i) = max(abs(res.gen(ng-nref+1:ng, 2)));
    else
        dRef(i) = 0;
    end
end
V = table(slots(:), ok, dV, lossF, lossM, PsF, PsM, dRef, its, upIt, ...
    'VariableNames', {'Slot', 'Success', 'MaxDeltaV_pu', 'Loss_FBS_kW', 'Loss_MP_kW', ...
                      'Psub_FBS_MW', 'Psub_MP_MW', 'RefGenMismatch_MW', 'NR_iter', 'UPFC_iter'});
end
