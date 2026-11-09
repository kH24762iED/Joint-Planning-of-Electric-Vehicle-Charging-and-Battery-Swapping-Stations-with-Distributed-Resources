function S = RC3_decode(x, P)
%RC3_DECODE  Decision vector -> settings of the installed devices.
%
%   Time-independent (site / count / technology) decisions decide which
%   units exist and where; time-varying decisions give one set-point per
%   slot (columns = slots). Only installed units (y = 1) are returned.
%
%   S.X  (nl x T) branch status incl. ONR, PCC and microgrid cuts (true = closed)
%   S.Xonr        switch status set by ONR only;  S.ug (1 x T) PCC status
%   S.cut (ncut x T) branches opened to form microgrids (0 = none)
%   S.dg   .unit .bus .type .P (MW, n x T) .cap (MW) .avail (MW, n x T) .gf .tanphi
%   S.sc   .unit .bus .Q (MVAr at 1 p.u., n x T) .steps .cap
%   S.oltc .tap (1 x T)
%   S.tx   .unit .br .tap .a (regulation ratio, n x T)
%   S.st   .unit .bus .Q (n x T; injection 'pq' / limit 'pv') .cap
%   S.up   .unit .br .r .th .Qsh (n x T)
%   S.ev   .unit .bus .n .type .Pch (kW/charger) .capf (p.u., n x T)
%   S.bss  .unit .bus .n .type .Pch .Nbat .Ecap (kWh) .u | .ach .adis (n x T)
%   S.dr   .shift .il (1 x T) .z (nb x 1) .kappa (nb x 1) .fac (nb x T)
x = reshape(x, 1, []);
T = P.T;  nb = P.nb;
e = zeros(0, 1);  eT = zeros(0, T);
v = struct();
for b = P.blocks
    u = x(b.idx);
    if b.enc(1) ~= 1, u = round(u); end
    u = min(max(u, b.lo), b.up);
    v.(b.name) = reshape(u, b.dims);
end

%% ---- switches, PCC and microgrid cuts --------------------------------------------
X = repmat(P.x0, 1, T);
if isfield(v, 'ONR_open')
    nt = numel(P.tie);  X = true(P.nl, T);
    for t = 1:T
        for k = 1:nt
            X(P.loops{k}(v.ONR_open(k, t)), t) = false;
        end
    end
elseif isfield(v, 'ONR_sw')
    X = v.ONR_sw > 0.5;
end
S.Xonr = X;
ug = P.ugForced;
if isfield(v, 'PCC_u'), d = isnan(ug); ug(d) = v.PCC_u(d); end
ug(isnan(ug)) = 1;
S.ug = ug;
X(P.pcc, :) = X(P.pcc, :) & (ug > 0.5);
S.cut = zeros(P.ncut, T);
if isfield(v, 'ISL_cut')
    for t = find(ug < 0.5)
        c = v.ISL_cut(:, t);  c = c(c > 0 & c ~= P.pcc);
        X(c, t) = false;  S.cut(1:numel(c), t) = c;
    end
end
S.X = X;

%% ---- DG -------------------------------------------------------------------------
S.dg = struct('unit', e, 'bus', e, 'type', e, 'P', eT, 'cap', e, 'avail', eT, 'gf', false(0, 1), 'tanphi', e);
if P.use.DG
    on = find(v.DG_y > 0.5);
    if ~isempty(on)
        ty  = P.dg.allowed(v.DG_type(on));  ty = ty(:);
        prf = P.dg.avail(ty, :);                                 % p.u. availability
        xs  = v.DG_P(on, :);
        isR = ~strcmp(P.dg.profile(ty), 'const');  isR = isR(:);
        Pdg = xs;
        Pdg(isR, :) = min(xs(isR, :), P.dg.Pmax * prf(isR, :));   % renewable: <= available
        ratio = Pdg ./ max(prf, eps);  ratio(prf < 1e-3) = 0;
        cap = max(ratio, [], 2);  cap = min(cap, P.dg.Pmax);
        S.dg.unit = on(:);  S.dg.bus = P.dg.cand(v.DG_bus(on))';  S.dg.type = ty;
        S.dg.P = Pdg;  S.dg.cap = cap;  S.dg.avail = cap .* prf;
        S.dg.gf = logical(P.dg.gf(ty))' & P.isl.enabled;
        S.dg.tanphi = P.dg.tanphi(ty)';
    end
end

%% ---- SC, OLTC, TX, D-STATCOM, UPFC -------------------------------------------------
S.sc = struct('unit', e, 'bus', e, 'steps', eT, 'Q', eT, 'cap', e);
if P.use.SC
    on = find(v.SC_y > 0.5);
    S.sc.unit = on(:);  S.sc.bus = P.sc.cand(v.SC_bus(on))';
    S.sc.steps = v.SC_step(on, :);  S.sc.Q = S.sc.steps * P.sc.Qstep;
    S.sc.cap = max(S.sc.Q, [], 2);
    if isempty(on), S.sc.Q = eT; S.sc.steps = eT; S.sc.cap = e; end
end
S.oltc.tap = zeros(1, T);
if P.use.TAP, S.oltc.tap = v.TAP_oltc; end
S.tx = struct('unit', e, 'br', e, 'tap', eT, 'a', eT);
if P.use.TX
    on = find(v.TX_y > 0.5);
    if ~isempty(on)
        S.tx.unit = on(:);  S.tx.br = P.tx.cand(v.TX_br(on))';
        S.tx.tap = v.TX_tap(on, :);  S.tx.a = 1 + P.tx.step * S.tx.tap;
    end
end
S.st = struct('unit', e, 'bus', e, 'Q', eT, 'cap', e);
if P.use.DSTATCOM
    on = find(v.ST_y > 0.5);
    if ~isempty(on)
        S.st.unit = on(:);  S.st.bus = P.st.cand(v.ST_bus(on))';
        S.st.Q = v.ST_Q(on, :);  S.st.cap = max(abs(S.st.Q), [], 2);
    end
end
S.up = struct('unit', e, 'br', e, 'r', eT, 'th', eT, 'Qsh', eT);
if P.use.UPFC
    on = find(v.UP_y > 0.5);
    if ~isempty(on)
        S.up.unit = on(:);  S.up.br = P.up.cand(v.UP_br(on))';
        S.up.r = v.UP_r(on, :);  S.up.th = v.UP_th(on, :);  S.up.Qsh = v.UP_Qsh(on, :);
    end
end

%% ---- EV charging stations and battery swapping stations ------------------------------
S.ev = struct('unit', e, 'bus', e, 'n', e, 'type', e, 'Pch', e, 'capf', eT);
if P.use.EVCS
    on = find(v.EV_y > 0.5);
    if ~isempty(on)
        S.ev.unit = on(:);  S.ev.bus = P.ev.cand(v.EV_bus(on))';
        S.ev.n = v.EV_nch(on);  S.ev.type = v.EV_type(on);
        S.ev.Pch = P.ev.rates(S.ev.type)';  S.ev.capf = v.EV_cap(on, :);
        S.ev.n = S.ev.n(:);  S.ev.type = S.ev.type(:);
    end
end
S.bss = struct('unit', e, 'bus', e, 'n', e, 'type', e, 'Pch', e, 'Nbat', e, 'Ecap', e, ...
               'u', eT, 'ach', eT, 'adis', eT);
if P.use.BSS
    on = find(v.BSS_y > 0.5);
    if ~isempty(on)
        S.bss.unit = on(:);  S.bss.bus = P.bss.cand(v.BSS_bus(on))';
        S.bss.n = v.BSS_nch(on);  S.bss.n = S.bss.n(:);
        S.bss.type = v.BSS_type(on);  S.bss.type = S.bss.type(:);
        S.bss.Pch  = P.bss.rates(S.bss.type)';
        S.bss.Nbat = S.bss.n * P.bss.bpc;  S.bss.Ecap = S.bss.Nbat * P.bss.Ebat;
        if strcmp(P.bss.control, 'rule')
            S.bss.ach  = v.BSS_ach(on, :);
            S.bss.adis = max(v.BSS_adis(on, :), S.bss.ach + P.bss.gap);
        else
            S.bss.u = v.BSS_u(on, :);
        end
    end
end

%% ---- demand response ---------------------------------------------------------------
S.dr.shift = zeros(1, T);  S.dr.il = zeros(1, T);
S.dr.z = double(P.drm.part(:) > 0);
if P.use.DRM
    if isfield(v, 'DR_shift')                                 % same projection as RC3_repair
        S.dr.shift = RC3_dr_project(v.DR_shift, P.dr_w, -P.drm.rho_max, P.drm.rho_max);
    end
    S.dr.il = v.DR_il;
    B = P.drm.il_frac * sum(P.dr_w);  E = S.dr.il * P.dr_w';
    if E > B * (1 + 1e-12), S.dr.il = S.dr.il * (B / E); end
    if isfield(v, 'DR_z')
        S.dr.z = zeros(nb, 1);  S.dr.z(P.drm.cand) = v.DR_z(:);
    end
else
    S.dr.z = zeros(nb, 1);
end
S.dr.kappa = P.drm.part(:) .* S.dr.z;
S.dr.fac   = 1 + S.dr.kappa * (S.dr.shift - S.dr.il);              % nb x T load factor
end
