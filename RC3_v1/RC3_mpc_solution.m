function [mpc, info] = RC3_mpc_solution(P, out, t)
%RC3_MPC_SOLUTION  MATPOWER case of time slot t of an evaluated solution.
%
%   [mpc, info] = RC3_mpc_solution(P, out, t)      out from RC3_objectives
%
%   The case contains the device settings of slot t:
%     branch(:, BR_STATUS)  switches incl. ONR, PCC and microgrid cuts
%     branch(:, TAP)        in-line transformers (MATPOWER ratio = 1 / a)
%     bus(:, PD/QD)         served loads (after DR and shedding) + EV charging
%     bus(:, BS)            shunt capacitors (MVAr at 1 p.u.)
%     gen                   substation (REF, VG = OLTC voltage, status = PCC),
%                           DGs, BSS pools, D-STATCOMs ('pq': fixed Q, 'pv':
%                           PV bus with Q limits) and, in islanded slots, one
%                           reference generator (PG = QG = 0 set-point, VG = droop
%                           voltage) at the root bus of every energized microgrid
%   Buses of dead islands are type NONE. UPFCs are not native MATPOWER
%   devices: info.upfc lists them; RC3_validate_solution adds their
%   injection model (Noroozian et al. 1997) iteratively.
%
%   info.energized (nb x 1 logical), info.ref (root buses), info.ngen0
S    = out.S;  nb = P.nb;
bus  = P.mpc.bus;  branch = P.mpc.branch;  gen0 = P.mpc.gen(1, :);
ug   = out.ug(t) > 0.5;
%% branches
branch(:, 11) = double(out.X(:, t));
branch(:, 9)  = 0;
for k = 1:numel(S.tx.br), branch(S.tx.br(k), 9) = 1 / S.tx.a(k, t); end
%% loads (served) and shunts
Pl = out.Pl(:, t);  Ql = out.Ql(:, t);
ks = ones(nb, 1);  ip = Pl > 0;  ks(ip) = 1 - out.shed(ip, t) ./ Pl(ip);
Pev = accumarray(S.ev.bus(:), out.ev.P(:, t) - out.evShed(:, t), [nb 1]) / 1000;
bus(:, 3) = Pl .* ks + Pev;
bus(:, 4) = Ql .* ks + Pev * P.tan_ev;
bus(:, 6) = bus(:, 6) + accumarray(S.sc.bus(:), S.sc.Q(:, t), [nb 1]);
bus(:, 2) = 1;  bus(:, 8) = 1;  bus(:, 9) = 0;
%% generators
mk  = @(b, pg, qg, vg) [b, pg, qg, qg, qg, vg, P.baseMVA, 1, 999, -999, zeros(1, 11)];
gen = gen0;  gen(1, 8) = double(ug);  gen(1, 6) = 1 + P.oltc.step * S.oltc.tap(t);
for k = 1:numel(S.dg.bus)
    if out.Pdg(k, t) ~= 0 || out.Qdg(k, t) ~= 0
        gen(end+1, :) = mk(S.dg.bus(k), out.Pdg(k, t), out.Qdg(k, t), 1); %#ok<AGROW>
    end
end
for k = 1:numel(S.bss.bus)
    Pk = -out.bss.P(k, t) / 1000;  Qk = -out.Qbss(k, t);           % injection
    if Pk ~= 0 || Qk ~= 0, gen(end+1, :) = mk(S.bss.bus(k), Pk, Qk, 1); end %#ok<AGROW>
end
stPV = strcmp(P.st.mode, 'pv');
for k = 1:numel(S.st.bus)
    if stPV
        if S.st.Q(k, t) > 0
            g = mk(S.st.bus(k), 0, 0, P.st.vset);  g(4) = S.st.Q(k, t);  g(5) = -S.st.Q(k, t);
            gen(end+1, :) = g;  bus(S.st.bus(k), 2) = 2; %#ok<AGROW>
        end
    elseif S.st.Q(k, t) ~= 0
        gen(end+1, :) = mk(S.st.bus(k), 0, S.st.Q(k, t), 1); %#ok<AGROW>
    end
end
%% reference buses and dead buses
energized = ~isnan(out.V(:, t));
info.ref = [out.roots{t}.bus];
bus(1, 2) = 3;
if ~ug
    bus(1, 2) = 4;
    for k = 1:numel(out.roots{t})
        rb = out.roots{t}(k).bus;
        bus(rb, 2) = 3;
        gen(end+1, :) = mk(rb, 0, 0, abs(out.roots{t}(k).V)); %#ok<AGROW>
        gen(end, 4) = 999;  gen(end, 5) = -999;
    end
end
dead = ~energized;  dead(1) = dead(1) || ~ug;
bus(dead, 2) = 4;  bus(dead, 3:4) = 0;
gen(ismember(gen(:, 1), find(dead)), 8) = 0;
% all generators of a regulated bus share its voltage set-point
if stPV
    for k = 1:numel(S.st.bus), gen(gen(:, 1) == S.st.bus(k), 6) = P.st.vset; end
end
if ug
    gen(gen(:, 1) == 1, 6) = gen(1, 6);
else
    for k = 1:numel(out.roots{t})
        gen(gen(:, 1) == out.roots{t}(k).bus, 6) = abs(out.roots{t}(k).V);
    end
end
mpc = struct('version', '2', 'baseMVA', P.baseMVA, 'bus', bus, 'gen', gen, 'branch', branch);
info.energized = energized;  info.ngen0 = 1;
info.upfc = struct('br', S.up.br(:), 'r', S.up.r(:, t), 'th', S.up.th(:, t), 'Qsh', S.up.Qsh(:, t));
end
