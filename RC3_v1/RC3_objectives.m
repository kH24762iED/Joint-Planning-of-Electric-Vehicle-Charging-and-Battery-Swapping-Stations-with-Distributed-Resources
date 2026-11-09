function [F, G, out, Fall] = RC3_objectives(x, P)
%RC3_OBJECTIVES  Evaluate ONE decision vector: device models, slot-by-slot
%load flow of every energized component, all basic objectives, the user
%objective functions and the nonlinear constraints.
%
%   [F, G, out, Fall] = RC3_objectives(x, P)
%
%   F    : 1 x numel(P.objdef), the user objective functions requested in
%          opts.objectives, in that order (RC3_objfun): basic objectives
%          (raw), weighted sums of normalised basic objectives such as
%          f2 = 0.3*LOSS + 0.3*VD + 0.4*VSI, or function handles
%   G    : 1 x numel(P.conNl) violations (feasible if all <= 0)
%          g_rad g_V g_S g_ops g_pen g_site g_soc g_ev g_isl g_f g_rev
%   out  : all indicators, slot results and device settings; out.F = F,
%          out.Fall = Fall, out.Flib = library objectives f1..f5,
%          out.terms / out.libTerms = breakdown of the weighted sums
%   Fall : 1 x numel(P.allObj), all basic objectives (order of P.allObj);
%          load-flow based objectives are P.bigM if a slot is not radial
%          or a load flow does not converge
%       SWCOST   switching, PCC and tap operation cost (USD/day)
%       INVCOST  capital cost of DG, SC, TX, D-STATCOM, UPFC, EVCS, BSS (USD)
%       COST     SWCOST + INVCOST (USD)
%       DEGCOST  BSS degradation cost (USD/day)
%       OPCOST   SWCOST + DEGCOST + energy + fuel + DR + lost load - revenue (USD/day)
%       TAC      CRF*INVCOST + annual operating cost (USD/yr)
%       ENS SAIDI SAIFI   reliability (kWh/yr, h/cust.yr, int/cust.yr)
%       LOSS     mean active loss (kW)        VD   mean sum((1-V)^2) (p.u.^2)
%       VSI      mean 1/min VSI               PEAK max grid import (MW)
%       LVAR     variance of grid import      PAR  peak-to-average ratio
%       EMISSION CO2 of grid import and DGs (kg/day)
%       EVUNMET  EV + swap energy not served (kWh/day)   EVWAIT mean wait (min)
%       EVCSUN   EV charging energy not served (kWh/day)
%       BSSUN    battery swap energy not served (kWh/day)
%       DRCOST   DR incentive cost (USD/day)  DISCOMF  share of enrolled energy moved
%       LSHED    load shed (kWh/day)          WLSHED   priority-weighted load shed
%       FDEV     mean |frequency deviation| of the islands (Hz)
%       SAI      1 / min supply adequacy of the islands (0 if never islanded)
%       SAD      mean supply-adequacy deficit of the islands, in [0, 1]
%
%   Grid-connected slot (PCC closed): one radial network fed by the
%   substation (bus 1, V0 set by the OLTC), forward-backward sweep (RC3_pf_fbs).
%   Islanded slot (PCC open): every component that contains a grid-forming
%   unit (MT DG or BSS) is a microgrid solved by the droop load flow
%   RC3_pf_island (with load shedding); components without one are dead.
%   The result depends only on x and P: every solver calls this function
%   through RC3_evaluate / RC3_fitness, so a decision vector gets the same
%   objectives and constraints whatever the solver.

persistent cache lastX lastV lastKey
if isempty(cache), cache = containers.Map('KeyType', 'char', 'ValueType', 'any'); end
x = reshape(x, 1, []);
if nargout < 3 && ~isempty(lastX) && strcmp(lastKey, P.key) && isequal(x, lastX)
    Fall = lastV(1:P.nAll);  G = lastV(P.nAll+1:end);  F = userobj(Fall, P);  return
end
key = hashkey(P.key, x);
if nargout < 3 && isKey(cache, key)
    v = cache(key);  lastX = x;  lastV = v;  lastKey = P.key;
    Fall = v(1:P.nAll);  G = v(P.nAll+1:end);  F = userobj(Fall, P);  return
end

S    = RC3_decode(x, P);
nb   = P.nb;  nl = P.nl;  T = P.T;  dt = P.dt;  base = P.baseMVA;  w = P.w;
lim  = P.lim;  isl = P.isl;
bsum = @(b, v) accumarray(b(:), v(:), [nb 1]);
nEV = numel(S.ev.n);  nDG = numel(S.dg.bus);  nST = numel(S.st.bus);
nUP = numel(S.up.br); nTX = numel(S.tx.br);   nBS = numel(S.bss.n);
stPV = strcmp(P.st.mode, 'pv');

%% ---- loads after demand response, EV stations, BSS state --------------------------
Pl0  = P.Pd .* P.mult .* S.dr.fac;                       % MW, nb x T
Ql0  = P.Qd .* P.mult .* S.dr.fac;
LDR  = sum(S.dr.kappa .* P.Pd) * P.mult;                 % enrolled load (MW), 1 x T
ev   = RC3_dev_ev(S, P);
Lnet = sum(Pl0, 1) + sum(ev.P, 1) / 1000 - sum(S.dg.P, 1);
bss  = RC3_dev_bss('init', [], S, P, Lnet);

%% ---- slot-by-slot load flow ------------------------------------------------------
g_rad = 0;  g_V = 0;  g_S = 0;
[ENS_t, SAIDI_t, SAIFI_t, VD_t, loss, Psub, nIt] = deal(zeros(1, T));
VSI_t = inf(1, T);  hourOK = true(1, T);
Vc = nan(nb, T);  Sbr = zeros(nl, T);
Qst = zeros(nST, T);  Sse = zeros(nUP, T);  Ssh = zeros(nUP, T);
Pdg = zeros(nDG, T);  Qdg = zeros(nDG, T);  shedL = zeros(nb, T);  evShed = zeros(nEV, T);
curt = zeros(1, T);  dfT = zeros(1, T);  SAT = nan(1, T);  nIslT = zeros(1, T);
demIsl = zeros(1, T);  deadT = zeros(1, T);  phiT = zeros(1, T);
Qbs = zeros(nBS, T);  roots = cell(1, T);  Vws = [];  lastC = [];  netG = [];  keyG = [];
for t = 1:T
    bss    = RC3_dev_bss('swap', bss, t);
    closed = S.X(:, t);  ugt = S.ug(t) > 0.5;
    if ~isequal(closed, lastC), Vws = []; lastC = closed; end   % new topology: cold start
    tau    = ones(nl, 1);
    for k = 1:nTX, tau(S.tx.br(k)) = 1 / S.tx.a(k, t); end
    [comp, nc] = components(P, closed);
    nE    = accumarray(comp(P.fb(closed)), 1, [nc 1]);
    nV    = accumarray(comp, 1, [nc 1]);
    loops = max(nE - (nV - 1), 0);
    if any(loops), g_rad = g_rad + sum(loops);  hourOK(t) = false; end
    Pev = bsum(S.ev.bus, ev.P(:, t)) / 1000;  Qev = Pev * P.tan_ev;
    Ysh = P.Ysh0 + 1j * bsum(S.sc.bus, S.sc.Q(:, t)) / base;
    Qpq = zeros(nb, 1);
    if ~stPV && nST, Qpq = bsum(S.st.bus, S.st.Q(:, t)); end   % D-STATCOM injection
    if ugt
        %% grid-connected: one network fed by the substation
        [bss, Pb] = RC3_dev_bss('grid', bss, t);                 % kW, + = charge
        Pbs = bsum(S.bss.bus, Pb) / 1000;
        Pdg(:, t) = S.dg.P(:, t);  Qdg(:, t) = Pdg(:, t) .* S.dg.tanphi;
        Pg  = bsum(S.dg.bus, Pdg(:, t));  Qg = bsum(S.dg.bus, Qdg(:, t));
        S0  = complex(Pl0(:, t) + Pev + Pbs - Pg, ...
                      Ql0(:, t) + Qev + Pbs * P.tan_bss - Qg - Qpq) / base;
        if ~isequal(keyG, [closed; tau]), netG = build_net(P, closed, 1, tau);  keyG = [closed; tau]; end
        net = netG;  o = net.order;                               % reused while topology/taps repeat
        [pvL, ipv] = pv_local(P, S, net, t);  [upL, iup] = up_local(P, S, net, t);
        V0  = 1 + P.oltc.step * S.oltc.tap(t);
        r   = RC3_pf_fbs(net, V0, S0(o), Ysh(o), lim.zip_np, lim.zip_nq, pvL, upL, ...
                        lim.fbs_tol, lim.fbs_maxit, Vws);
        bss = RC3_dev_bss('commit', bss, t, Pb);
        Qbs(:, t) = Pb / 1000 * P.tan_bss;
        roots{t} = struct('bus', 1, 'V', V0);
        if r.success, Vws = r.Vt; end
        Psub(t) = real(r.Sroot) * base;
        dead = find(comp ~= comp(1));
        if ~isempty(dead)                                       % unintended islands
            g_rad = g_rad + numel(dead);  hourOK(t) = false;
            shedL(dead, t) = Pl0(dead, t);
            evShed(ismember(S.ev.bus, dead), t) = ev.P(ismember(S.ev.bus, dead), t);
            Pdg(ismember(S.dg.bus, dead), t) = 0;  Qdg(ismember(S.dg.bus, dead), t) = 0;
        end
        cres = {struct('net', net, 'r', r, 'Pl', Pl0(o, t), 'Pev', Pev(o), 'pv', ipv, 'up', iup, 'upQ', upL.Qsh)};
    else
        %% islanded: microgrids formed around grid-forming units
        [bss, cap] = RC3_dev_bss('island', bss, t);
        Pb   = zeros(nBS, 1);  cres = {};  roots{t} = struct('bus', {}, 'V', {});
        Pres = bsum(S.dg.bus(~S.dg.gf), S.dg.P(~S.dg.gf, t));        % renewable / non-GF
        Qres = bsum(S.dg.bus(~S.dg.gf), S.dg.P(~S.dg.gf, t) .* S.dg.tanphi(~S.dg.gf));
        for c = setdiff(1:nc, comp(1))
            buses = find(comp == c);
            iDG = find(S.dg.gf & S.dg.cap > 1e-6 & comp(S.dg.bus) == c);
            iBS = zeros(0, 1);
            if P.bss.gf && nBS, iBS = find(comp(S.bss.bus) == c); end
            gP  = [S.dg.cap(iDG);  cap.Pdis(iBS) / 1000];          % MW
            gA  = [zeros(numel(iDG), 1);  cap.Pch(iBS) / 1000];
            qB  = zeros(0, 1);
            if ~isempty(iBS), qB = bss.Pdism(iBS) / 1000; end
            gQ  = isl.q_frac * [S.dg.cap(iDG);  qB];
            gB  = [S.dg.bus(iDG);  S.bss.bus(iBS)];
            dem = sum(Pl0(buses, t)) + sum(Pev(buses));
            if dem > 1e-9
                SAT(t) = min(SAT(t), (sum(gP) + sum(Pres(buses))) / dem);
                demIsl(t) = demIsl(t) + sum(Pl0(buses, t));
            end
            inC = ismember(S.ev.bus, buses);  dgC = ismember(S.dg.bus, buses);
            if sum(gP) <= 1e-6                                     % dead component
                shedL(buses, t) = Pl0(buses, t);  evShed(inC, t) = ev.P(inC, t);
                Pdg(dgC, t) = 0;  Qdg(dgC, t) = 0;
                deadT(t) = deadT(t) + nnz(Pl0(buses, t) > 0);
                continue
            end
            [~, ir] = max(gP);  og = [ir, setdiff(1:numel(gP), ir)];
            net = build_net(P, closed, gB(ir), tau);  o = net.order;
            gfL = struct('node', net.loc(gB(og)), 'Pmax', gP(og) / base, ...
                         'Pabs', gA(og) / base, 'Qmax', max(gQ(og), 1e-6) / base);
            L = struct('Pl', Pl0(o, t) / base, 'Ql', Ql0(o, t) / base, ...
                       'Pev', Pev(o) / base, 'Qev', Qev(o) / base, ...
                       'Pres', Pres(o) / base, 'Qres', Qres(o) / base, ...
                       'Pfix', zeros(numel(o), 1), 'Qfix', -Qpq(o) / base, 'w', isl.w(o));
            [pvL, ipv] = pv_local(P, S, net, t);  [upL, iup] = up_local(P, S, net, t);
            [r, il] = RC3_pf_island(net, gfL, L, Ysh(o), pvL, upL, P);
            if ~il.ok, r.success = false; end
            shedL(o, t) = il.shed * base;
            evShed(inC, t) = il.fracEV * ev.P(inC, t);
            nG = ~S.dg.gf & dgC;
            Pdg(nG, t) = (1 - il.fracCurt) * S.dg.P(nG, t);  Qdg(nG, t) = Pdg(nG, t) .* S.dg.tanphi(nG);
            curt(t) = curt(t) + il.fracCurt * sum(Pres(buses));
            Pgf = zeros(numel(gP), 1);  Qgf = Pgf;
            Pgf(og) = il.Pgf * base;  Qgf(og) = il.Qgf * base;    % back to unit order
            nd = numel(iDG);
            Pdg(iDG, t) = Pgf(1:nd);  Qdg(iDG, t) = Qgf(1:nd);
            Pb(iBS) = -1000 * Pgf(nd+1:end);  Qbs(iBS, t) = -Qgf(nd+1:end);
            dfT(t) = max(dfT(t), abs(il.df));  phiT(t) = max(phiT(t), il.phi);
            nIslT(t) = nIslT(t) + 1;
            roots{t}(end+1) = struct('bus', gB(ir), 'V', il.Vroot);
            PlS = Pl0(o, t) - shedL(o, t);
            cres{end+1} = struct('net', net, 'r', r, 'Pl', PlS, 'Pev', (1 - il.fracEV) * Pev(o), ...
                                 'pv', ipv, 'up', iup, 'upQ', upL.Qsh); %#ok<AGROW>
        end
        bss  = RC3_dev_bss('commit', bss, t, Pb);
        Psub(t) = 0;
    end
    %% indicators of the energized components
    if isempty(cres), VSI_t(t) = NaN; end
    for ci = 1:numel(cres)
        cr = cres{ci};  net = cr.net;  r = cr.r;  o = net.order;  n = net.n;
        if ~r.success, g_V = g_V + 1;  hourOK(t) = false;  continue; end
        nIt(t) = max(nIt(t), r.it);
        Vm = abs(r.V);  Vc(o, t) = r.V;
        g_V = g_V + sum(max(P.Vmin - Vm, 0)) + sum(max(Vm - P.Vmax, 0));
        loss(t) = loss(t) + real(r.loss) * base * 1000;
        VD_t(t) = VD_t(t) + sum((1 - Vm).^2);
        if n > 1
            j  = (2:n)';  br = net.br(j);
            Sb = max(abs(r.Ss(j)), abs(r.Sr(j))) * base;  Sbr(br, t) = Sb;
            g_S = g_S + sum(max(Sb - P.Smax(br), 0));
            V1 = abs(r.Vt(net.par(j)));  R = real(net.Zt(j));  X = imag(net.Zt(j));
            Pr = real(r.Sr(j));  Qr = imag(r.Sr(j));
            VSI_t(t) = min(VSI_t(t), min(V1.^4 - 4*(Pr.*X - Qr.*R).^2 - 4*(Pr.*R + Qr.*X).*V1.^2));
            Pdown = net.K.' * ((cr.Pl + cr.Pev) * 1000);           % kW, subtree sums
            Ndown = net.K.' * P.Ncust(o);
            lr    = P.lambda(br) .* P.r_rep(br);
            ENS_t(t)   = ENS_t(t) + sum(lr .* Pdown(j));
            SAIDI_t(t) = SAIDI_t(t) + sum(lr .* Ndown(j)) / sum(P.Ncust);
            SAIFI_t(t) = SAIFI_t(t) + sum(P.lambda(br) .* Ndown(j)) / sum(P.Ncust);
        end
        if ~isempty(cr.pv.unit)                                   % D-STATCOM 'pv' mode
            Qst(cr.pv.unit, t) = r.Q(cr.pv.ic) .* cr.pv.share * base;
        end
        if ~isempty(cr.up.unit)
            Sse(cr.up.unit, t) = abs(r.Sse) * base;
            Ssh(cr.up.unit, t) = abs(complex(real(r.Sse), cr.upQ)) * base;
        end
    end
    if ~stPV && nST, Qst(:, t) = S.st.Q(:, t); end
    if ~isfinite(VSI_t(t)), VSI_t(t) = NaN; end
end
bss = RC3_dev_bss('final', bss);
evShedE   = sum(evShed(:)) * dt;                                  % kWh/day
ev.unmet  = ev.unmet + evShedE;  ev.Edel = ev.Edel - evShedE;

%% ---- operation, investment and energy costs ----------------------------------------
ops    = abs(double(S.X) - double([P.x0, S.X(:, 1:end-1)]));
ops(P.pcc, :) = 0;
nTrans = sum(abs(diff([1, S.ug > 0.5])));
tapOps = abs(diff([0, S.oltc.tap]));
txOps  = abs(diff([zeros(nTX, 1), S.tx.tap], 1, 2));
SWCOST = sum(P.c_sw .* sum(ops, 2)) + P.pccd.c_op * nTrans ...
       + P.oltc.c_op * sum(tapOps) + P.tx.c_op * sum(txOps(:));
cdg  = P.dg.c_cap(S.dg.type);   cdg = cdg(:);
INV.DG   = sum(cdg .* S.dg.cap);
INV.SC   = numel(S.sc.bus) * P.sc.c_fix + P.sc.c_var * sum(S.sc.cap);
INV.TX   = nTX * P.tx.c_fix;
INV.DSTATCOM = P.st.c_var * sum(S.st.cap);
Sup  = max(Sse, [], 2) + max(Ssh, [], 2);                         % MVA rating of each UPFC
if nUP == 0, Sup = zeros(0, 1); end
cc   = P.up.c_curve;
INV.UPFC = sum(1000 * Sup .* max(cc(1) * Sup.^2 + cc(2) * Sup + cc(3), 0));
cev  = P.ev.c_ch(S.ev.type);    cev = cev(:);
INV.EVCS = sum(S.ev.n .* cev) + P.ev.c_site * nEV;
cbs  = P.bss.c_ch(S.bss.type);  cbs = cbs(:);
INV.BSS  = sum(S.bss.n .* cbs) + P.bss.c_site * nBS + P.bss.c_bat * sum(S.bss.Ecap);
INVCOST  = INV.DG + INV.SC + INV.TX + INV.DSTATCOM + INV.UPFC + INV.EVCS + INV.BSS;
EGCOST = 1000 * dt * sum(P.price .* (max(Psub, 0) - P.econ.sell_frac * max(-Psub, 0)));
cf   = P.dg.c_fuel(S.dg.type);  cf = cf(:);
FUEL = dt * sum(cf .* sum(Pdg, 2));
DRCOST = dt * (P.drm.c_shift * sum(max(-S.dr.shift, 0) .* LDR) + P.drm.c_il * sum(S.dr.il .* LDR));
if P.use.DRM && P.drm.site_select, DRCOST = DRCOST + P.drm.c_enrol * sum(S.dr.z); end
LSHED  = 1000 * dt * sum(shedL(:));
WLSHED = 1000 * dt * sum(isl.w' * shedL);
SHEDCOST = isl.voll * LSHED;
REV    = P.ev.tariff * ev.Edel + P.bss.tariff * bss.Eswap;
nSW    = 1 + 364 * (T > 1);                                       % daily schedule repeats

%% ---- indicators -------------------------------------------------------------------
ok = all(hourOK);
val.SWCOST   = SWCOST;   val.INVCOST = INVCOST;   val.COST = SWCOST + INVCOST;
val.DEGCOST  = bss.deg;
val.OPCOST   = SWCOST + bss.deg + EGCOST + FUEL + DRCOST + SHEDCOST - REV;
val.TAC      = P.crf * INVCOST + nSW * SWCOST + 365 * (bss.deg + EGCOST + FUEL + DRCOST + SHEDCOST - REV);
val.ENS      = sum(w .* ENS_t);
val.SAIDI    = sum(w .* SAIDI_t);
val.SAIFI    = sum(w .* SAIFI_t);
val.LOSS     = sum(w .* loss);
val.VD       = sum(w .* VD_t);
vs = VSI_t(isfinite(VSI_t));
if isempty(vs), val.VSI = 0; else, val.VSI = mean(1 ./ max(vs, 1e-6)); end
val.PEAK     = max(Psub);
mP = sum(w .* Psub);
val.LVAR     = sum(w .* (Psub - mP).^2);
val.PAR      = max(Psub) / max(mP, 1e-3);
efd = P.dg.ef(S.dg.type);  efd = efd(:);
val.EMISSION = dt * sum(P.econ.ef_grid * max(Psub, 0)) + dt * sum(efd .* sum(Pdg, 2));
val.EVUNMET  = ev.unmet + bss.unmet;
val.EVWAIT   = ev.wait;
val.EVCSUN   = ev.unmet;
val.BSSUN    = bss.unmet;
val.DRCOST   = DRCOST;
Een = sum(LDR) * dt;
if Een > 0
    val.DISCOMF = dt * (0.5 * sum(abs(S.dr.shift) .* LDR) + sum(S.dr.il .* LDR)) / Een;
else
    val.DISCOMF = 0;
end
val.LSHED    = LSHED;  val.WLSHED = WLSHED;
hI = S.ug < 0.5;
if any(nIslT > 0), val.FDEV = mean(dfT(nIslT > 0)); else, val.FDEV = 0; end
sa = SAT(hI & isfinite(SAT));
if isempty(sa)
    val.SAI = 0;  val.SAD = 0;
else
    val.SAI = 1 / max(min(sa), 1e-3);                      % 1 / worst adequacy
    val.SAD = mean(max(0, 1 - sa));                        % mean adequacy deficit in [0, 1]
end

Fall = zeros(1, P.nAll);
for m = 1:P.nAll
    name = P.allObj{m};
    if ok || any(strcmp(name, P.pfFree)), Fall(m) = val.(name); else, Fall(m) = P.bigM; end
end

%% ---- remaining constraints -------------------------------------------------------------
g_ops  = sum(max(sum(ops, 2) - P.onr.max_ops, 0)) + max(sum(tapOps) - P.oltc.max_ops, 0) ...
       + sum(max(sum(txOps, 2) - P.tx.max_ops, 0));
g_pen  = max(sum(S.dg.cap) - P.dg.pen * P.Ptot, 0);
g_site = dup(S.dg.bus) + dup(S.sc.bus) + dup(S.st.bus) + dup(S.ev.bus) + dup(S.bss.bus) ...
       + dup([S.tx.br; S.up.br]);
Edem   = ev.Edem + bss.Edem;
g_ev   = 0;
if Edem > 0, g_ev = max(P.ev.serve_min - (1 - val.EVUNMET / Edem), 0); end
g_isl  = max(nTrans - isl.max_trans, 0);
if any(hI) && sum(demIsl) > 0
    served = 1 - sum(sum(shedL(:, hI))) / sum(demIsl(hI));
    g_isl  = g_isl + max(isl.min_serve - served, 0);
end
g_f    = sum(max(dfT - isl.f_lim, 0));
g_rev  = sum(max(-Psub - P.pccd.Pexp_max, 0));
G = [g_rad, g_V, g_S, g_ops, g_pen, g_site, bss.g_soc, g_ev, g_isl, g_f, g_rev];

if cache.Count > 2e4, remove(cache, keys(cache)); end
v = [Fall, G];  cache(key) = v;  lastX = x;  lastV = v;  lastKey = P.key;
F = userobj(Fall, P);

if nargout > 2
    out = val;  out.Fall = Fall;  out.G = G;  out.feasible = all(G <= 0);  out.ok = ok;
    out.F = F;
    if ~isempty(F), [~, out.terms] = RC3_objfun(Fall, P); else, out.terms = []; end
    if isfield(P, 'objLib') && ~isempty(P.objLib)
        [out.Flib, out.libTerms] = RC3_objfun(Fall, P, P.objLib);
    else
        out.Flib = [];  out.libTerms = [];
    end
    out.EGCOST = EGCOST;  out.FUEL = FUEL;  out.SHEDCOST = SHEDCOST;  out.REVENUE = REV;
    out.INV = INV;  out.nTrans = nTrans;
    out.S = S;  out.X = S.X;  out.ops = ops;  out.tapOps = tapOps;  out.txOps = txOps;
    out.ENS_t = ENS_t;  out.SAIDI_t = SAIDI_t;  out.SAIFI_t = SAIFI_t;
    out.loss_kW = loss;  out.VD_t = VD_t;  out.VSI_t = VSI_t;  out.VSImin = min(VSI_t);
    out.Psub = Psub;  out.Lnet = Lnet;  out.V = Vc;  out.Vm = abs(Vc);
    out.Vmin = min(abs(Vc), [], 1);  out.Sbr = Sbr;
    out.hourOK = hourOK;  out.fbs_it = nIt;  out.ev = ev;  out.bss = bss;
    out.Pl = Pl0;  out.Ql = Ql0;  out.LDR = LDR;  out.shed = shedL;  out.evShed = evShed;
    out.Pdg = Pdg;  out.Qdg = Qdg;  out.Qst = Qst;  out.Sse = Sse;  out.Ssh = Ssh;  out.Sup = Sup;
    out.Qbss = Qbs;  out.curt = curt;  out.df = dfT;  out.phi = phiT;  out.SA = SAT;
    out.nIsl = nIslT;  out.dead = deadT;  out.roots = roots;  out.ug = S.ug;
    out.open = arrayfun(@(t) find(~S.X(:, t))', 1:T, 'UniformOutput', false);
end
end

% =========================================================================================
function F = userobj(Fall, P)
% user objective functions (empty while RC3_setup computes the base case)
if isfield(P, 'objdef') && ~isempty(P.objdef), F = RC3_objfun(Fall, P); else, F = zeros(1, 0); end
end

function k = hashkey(pkey, x)
% exact cache key of a decision vector (MD5 of its bytes, or all digits)
persistent useJava
if isempty(useJava), useJava = usejava('jvm'); end
if useJava
    md = java.security.MessageDigest.getInstance('MD5');
    h  = typecast(md.digest(typecast(x(:), 'int8')), 'uint8');
    k  = [pkey '|' sprintf('%02x', h)];
else
    k  = [pkey '|' sprintf('%.17g,', x)];
end
end

function n = dup(v)
% number of repeated entries in v
n = numel(v) - numel(unique(v));
end

function [comp, nc] = components(P, closed)
% connected components of the network formed by the closed branches
nb = P.nb;  comp = zeros(nb, 1);  nc = 0;
adj = cell(nb, 1);
for e = find(closed(:))'
    adj{P.fb(e)}(end+1) = P.tb(e);  adj{P.tb(e)}(end+1) = P.fb(e);
end
for s = 1:nb
    if comp(s), continue; end
    nc = nc + 1;  comp(s) = nc;  q = s;  h = 1;
    while h <= numel(q)
        u = q(h);  h = h + 1;
        for v = adj{u}
            if ~comp(v), comp(v) = nc;  q(end+1) = v; end %#ok<AGROW>
        end
    end
end
end

function net = build_net(P, closed, root, tau)
% radial network of the component containing 'root': BFS order, path
% matrix and the tap-ratio referral of every node (see RC3_pf_fbs)
[par, pbr, order] = RC3_tree(P.fb, P.tb, closed, P.nb, root);
n = numel(order);  loc = zeros(P.nb, 1);  loc(order) = 1:n;
parL = zeros(n, 1);  br = zeros(n, 1);  fwd = true(n, 1);
if n > 1
    parL(2:n) = loc(par(order(2:n)));
    br(2:n)   = pbr(order(2:n));
    fwd(2:n)  = P.fb(br(2:n)) == par(order(2:n));
end
K = zeros(n);
c = ones(n, 1);  zc = ones(n, 1);
for j = 2:n
    K(j, :) = K(parL(j), :);  K(j, j) = 1;
    tb = tau(br(j));
    if fwd(j), c(j) = c(parL(j)) / tb;  zc(j) = c(j);       % tap at the sending end
    else,      c(j) = c(parL(j)) * tb;  zc(j) = c(parL(j)); % tap at the receiving end
    end
end
Zt = zeros(n, 1);
if n > 1, Zt(2:n) = P.Zbr(br(2:n)) ./ zc(2:n).^2; end
net = struct('order', order, 'loc', loc, 'n', n, 'par', parL, 'br', br, 'fwd', fwd, ...
             'K', K, 'c', c, 'zc', zc, 'Zt', Zt);
end

function [pv, info] = pv_local(P, S, net, t)
% D-STATCOMs in 'pv' mode of this component as voltage-controlled nodes
pv = struct('node', zeros(0, 1), 'Vset', zeros(0, 1), 'Qmin', zeros(0, 1), 'Qmax', zeros(0, 1));
info = struct('unit', zeros(0, 1), 'ic', zeros(0, 1), 'share', zeros(0, 1));
if ~strcmp(P.st.mode, 'pv') || isempty(S.st.bus), return; end
u = find(net.loc(S.st.bus) > 1 & S.st.Q(:, t) > 0);          % not at the root node
if isempty(u), return; end
[nodes, ~, ic] = unique(net.loc(S.st.bus(u)));
qmax = accumarray(ic, S.st.Q(u, t)) / P.baseMVA;
pv = struct('node', nodes(:), 'Vset', P.st.vset * ones(numel(nodes), 1), 'Qmin', -qmax, 'Qmax', qmax);
info = struct('unit', u(:), 'ic', ic(:), 'share', S.st.Q(u, t) / P.baseMVA ./ qmax(ic));
end

function [up, info] = up_local(P, S, net, t)
% UPFCs whose branch is closed in this component (series source + shunt converter)
up = struct('node', zeros(0, 1), 'ft', zeros(0, 1), 'k', zeros(0, 1), 'Qsh', zeros(0, 1));
info = struct('unit', zeros(0, 1));
if isempty(S.up.br), return; end
for u = 1:numel(S.up.br)
    j = find(net.br == S.up.br(u), 1);
    if isempty(j), continue; end
    if net.fwd(j), ft = net.par(j); s = 1; else, ft = j; s = -1; end
    up.node(end+1, 1) = j;  up.ft(end+1, 1) = ft;
    up.k(end+1, 1)    = s * S.up.r(u, t) * exp(1j * S.up.th(u, t));
    up.Qsh(end+1, 1)  = S.up.Qsh(u, t) / P.baseMVA;
    info.unit(end+1, 1) = u;
end
end
