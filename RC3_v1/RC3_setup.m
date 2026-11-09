function P = RC3_setup(T, opts)
%RC3_SETUP  Problem data, decision vector and MP-Opt-Model formulation for the
%multi-period planning/operation of a distribution network with any
%combination of the devices below (IEEE 33-bus by default, RC3_case33).
%
%   P = RC3_setup(T)          T = 1 (peak-hour snapshot) or T = 24 (hourly)
%   P = RC3_setup(T, opts)
%
%   opts.devices   any combination of
%       'ONR'      network reconfiguration (switch status, hourly)
%       'DG'       distributed generators: count, bus, type, P(t)
%       'SC'       shunt capacitors:       count, bus, steps(t)
%       'TAP'      substation OLTC:        tap(t)
%       'TX'       in-line transformers:   count, branch, tap(t)
%       'DSTATCOM' D-STATCOMs:             count, bus, Q(t)
%       'UPFC'     UPFCs:                  count, branch, r(t), theta(t), Qsh(t)
%       'EVCS'     EV charging stations:   count, bus, chargers, type, cap(t)
%       'BSS'      battery swapping:       count, bus, chargers, type, u(t)
%       'DRM'      demand response:        shift(t), curtailment(t) [, enrolment]
%       'ISL'      islanding decisions:    PCC status u(t), microgrid cuts(t)
%   opts.objectives  cell array of user objective functions, each one of
%       'TAC'                           a basic objective (P.allObj), raw value
%       'f2'                            a function of the library RC3_objlib
%       '0.3*LOSS + 0.3*VD + 0.4*VSI'   weighted sum of normalised basic objectives
%       'tech = 0.5*LOSS + 0.5*VD'      the same with a name
%       struct('name',..,'expr',..,'norm',true|false,'label',..)
%       @(V, Vn) ...  or struct('name',..,'fn',@(V,Vn) ..)   any function of the
%                       raw (V) and normalised (Vn) basic objectives
%   opts.objref      struct of reference values overriding P.objRef, e.g.
%                    struct('LOSS', 202.68)
%   opts.m           1 = single objective (weighted sum of the user objectives
%                    normalised by their base-case values, opts.weights), or
%                    m = numel(objectives) (multi-objective)
%   opts.island_mode 'off' | 'outage' | 'optimal' | 'islanded'
%   opts.mpc         MATPOWER case with device data (default RC3_case33())
%   opts.set         {'field.sub', value; ...} overrides of mpc fields
%   opts.dg_types    e.g. {'PV'} or {'PV','MT'}: DG types allowed
%   opts.load_mult   load multiplier(s) overriding the profile
%   opts.seed        seed of the common initial population and of the
%                    solver's random stream (default 1)
%
%   Site / count / technology variables are time-independent; all
%   operating set-points (powers, steps, taps, ratios) have one value per
%   time slot, so their dimension grows with T.

persistent nSetup
if isempty(nSetup), nSetup = 0; end
if nargin < 1 || isempty(T), T = 24; end
if nargin < 2 || isempty(opts), opts = struct(); end
validateattributes(T, {'numeric'}, {'scalar', 'integer', 'positive'}, mfilename, 'T', 1);
row = @(v) reshape(v, 1, []);

%% ---- MATPOWER case and overrides -----------------------------------------------------
mpc = getopt(opts, 'mpc', []);
if isempty(mpc), mpc = RC3_case33(); elseif ischar(mpc), mpc = RC3_case33(mpc); end
ov = getopt(opts, 'set', {});
for k = 1:size(ov, 1)
    mpc = setpath(mpc, ov{k, 1}, ov{k, 2});
end
P.mpc = mpc;
P.seed = getopt(opts, 'seed', 1);

%% ---- network (MATPOWER column indices are used as literals) ---------------------------
nb = size(mpc.bus, 1);  nl = size(mpc.branch, 1);
assert(isequal(mpc.bus(:, 1), (1:nb)'), 'RC3_setup: buses must be numbered 1..nb (use ext2int).');
assert(mpc.bus(1, 2) == 3, 'RC3_setup: bus 1 must be the substation (REF) bus.');
P.baseMVA = mpc.baseMVA;  P.baseKV = mpc.bus(1, 10);
P.fb  = mpc.branch(:, 1);  P.tb = mpc.branch(:, 2);
P.Zbr = complex(mpc.branch(:, 3), mpc.branch(:, 4));          % p.u.
P.Pd  = mpc.bus(:, 3);  P.Qd = mpc.bus(:, 4);                  % MW, MVAr
P.Ysh0 = complex(mpc.bus(:, 5), mpc.bus(:, 6)) / P.baseMVA;    % existing shunts (p.u.)
rate = mpc.branch(:, 6);  rate(rate <= 0) = Inf;  P.Smax = rate;  % MVA
x0 = mpc.branch(:, 11) > 0;                                    % base: closed branches
P.nb = nb;  P.nl = nl;  P.T = T;  P.dt = 24 / T;
P.x0 = x0;  P.sec = find(x0);  P.tie = find(~x0);
P.pcc = mpc.pcc.branch;  P.pccd = mpc.pcc;
P.Ptot = sum(P.Pd);
P.lim = mpc.limits;
P.Vmin = mpc.limits.Vmin;  P.Vmax = mpc.limits.Vmax;
P.bigM = 1e9;

%% ---- devices and island mode ------------------------------------------------------------
allDev = {'ONR','DG','SC','TAP','TX','DSTATCOM','UPFC','EVCS','BSS','DRM','ISL'};
P.devices = reshape(upper(cellstr(getopt(opts, 'devices', {'ONR'}))), 1, []);
assert(all(ismember(P.devices, allDev)), 'devices must be from: %s', strjoin(allDev, ', '));
P.use = cell2struct(num2cell(ismember(allDev, P.devices)), allDev, 2);
mode = lower(getopt(opts, 'island_mode', mpc.island.mode));
assert(ismember(mode, {'off','outage','optimal','islanded'}), ...
       'island_mode must be off, outage, optimal or islanded.');
if P.use.ISL && strcmp(mode, 'off'), mode = 'optimal'; end
P.isl = mpc.island;  P.isl.mode = mode;  P.isl.enabled = ~strcmp(mode, 'off');

%% ---- basic objectives (all are computed for every solution) --------------------------------
info = {'SWCOST',   'USD/day',      'switching, PCC and tap operation cost'; ...
        'INVCOST',  'USD',          'investment cost of the installed units'; ...
        'COST',     'USD',          'switching + investment cost'; ...
        'DEGCOST',  'USD/day',      'BSS battery degradation cost'; ...
        'OPCOST',   'USD/day',      'daily operating cost'; ...
        'TAC',      'USD/yr',       'total annual cost'; ...
        'ENS',      'kWh/yr',       'energy not supplied'; ...
        'SAIDI',    'h/cust.yr',    'system average interruption duration index'; ...
        'SAIFI',    'int/cust.yr',  'system average interruption frequency index'; ...
        'LOSS',     'kW',           'mean active power loss'; ...
        'VD',       'p.u.^2',       'mean voltage deviation sum((1-V)^2)'; ...
        'VSI',      '-',            'mean of 1/min voltage stability index'; ...
        'PEAK',     'MW',           'peak grid import'; ...
        'LVAR',     'MW^2',         'variance of the grid import'; ...
        'PAR',      '-',            'peak-to-average ratio of the grid import'; ...
        'EMISSION', 'kg/day',       'CO2 emission of grid import and DGs'; ...
        'EVUNMET',  'kWh/day',      'EV charging + battery swap energy not served'; ...
        'EVWAIT',   'min',          'mean EV waiting time (lost EVs included)'; ...
        'DRCOST',   'USD/day',      'demand-response incentive cost'; ...
        'DISCOMF',  '-',            'DR discomfort (share of enrolled energy moved)'; ...
        'LSHED',    'kWh/day',      'load shed'; ...
        'WLSHED',   'kWh/day',      'priority-weighted load shed'; ...
        'FDEV',     'Hz',           'mean |frequency deviation| of the islands'; ...
        'SAI',      '-',            '1 / minimum supply adequacy of the islands'; ...
        'SAD',      '-',            'mean supply-adequacy deficit of the islands'; ...
        'EVCSUN',   'kWh/day',      'EV charging energy not served'; ...
        'BSSUN',    'kWh/day',      'battery swap energy not served'};
P.allObj  = info(:, 1)';  P.objUnit = info(:, 2)';  P.objDesc = info(:, 3)';
P.nAll    = numel(P.allObj);
P.pfFree  = {'SWCOST','INVCOST','COST','DEGCOST','EVUNMET','EVWAIT','DRCOST','DISCOMF','EVCSUN','BSSUN'};
P.scaleOnly = {'FDEV','SAI','SAD','DISCOMF','EVWAIT'};   % always normalised by their scale
P.objdef = [];  P.objLib = [];  P.objRef = ones(1, P.nAll);

%% ---- time slots and profiles -------------------------------------------------------------
pr = mpc.profile;
if T == 1
    hrs = pr.peak_hour;
    pick = @(p) p(hrs);
    P.mult = 1;  avail = struct('pv', 1, 'wt', 1, 'const', 1);   % rated snapshot
    P.ev_arr  = pr.ev(hrs) / sum(pr.ev);
    P.bss_arr = pr.swap(hrs) / sum(pr.swap);
else
    hrs  = linspace(1, 24, T);
    pick = @(p) interp1(1:24, row(p), hrs);
    P.mult = pick(pr.load);
    avail = struct('pv', pick(pr.pv), 'wt', pick(pr.wt), 'const', ones(1, T));
    P.ev_arr  = pick(pr.ev / sum(pr.ev));   P.ev_arr  = P.ev_arr  / (sum(P.ev_arr)  * P.dt);
    P.bss_arr = pick(pr.swap / sum(pr.swap)); P.bss_arr = P.bss_arr / (sum(P.bss_arr) * P.dt);
end
lm = getopt(opts, 'load_mult', []);
if ~isempty(lm), P.mult = row(lm) .* ones(1, T); end
P.hours = hrs;
P.price = pick(pr.price);
P.w     = ones(1, T) / T;

%% ---- device data -----------------------------------------------------------------------
P.onr  = mpc.onr;  P.onr.encoding = lower(getopt(opts, 'encoding', mpc.onr.encoding));
P.c_sw = P.onr.c_sec * x0 + P.onr.c_tie * ~x0;
P.dg   = mpc.dg;   P.sc = mpc.sc;   P.oltc = mpc.oltc;  P.tx = mpc.tx;
P.st   = mpc.dstatcom;  P.up = mpc.upfc;  P.ev = mpc.evcs;  P.bss = mpc.bss;
P.drm  = mpc.drm;  P.rel = mpc.rel;  P.econ = mpc.econ;
dgt = getopt(opts, 'dg_types', {});
dgp = lower(getopt(opts, 'dg_profile', ''));             % backward compatibility
if isempty(dgt) && strcmp(dgp, 'pv'), dgt = {'PV'}; end
if isempty(dgt) && strcmp(dgp, 'dispatchable'), dgt = {'MT'}; end
if ~isempty(dgt), [~, P.dg.allowed] = ismember(upper(cellstr(dgt)), upper(P.dg.types)); end
assert(all(P.dg.allowed > 0), 'Unknown DG type in opts.dg_types.');
nty = numel(P.dg.types);
P.dg.avail = zeros(nty, T);                               % availability (p.u. of rating)
for k = 1:nty, P.dg.avail(k, :) = avail.(P.dg.profile{k}); end
P.dg.tanphi = tan(acos(P.dg.pf));
P.st.mode = lower(P.st.mode);
P.bss.control = lower(P.bss.control);  P.bss.signal = lower(P.bss.signal);
P.tan_ev = tan(acos(P.ev.pf));  P.tan_bss = tan(acos(P.bss.pf));
P.crf = P.econ.ir * (1 + P.econ.ir)^P.econ.life / ((1 + P.econ.ir)^P.econ.life - 1);
P.lambda = P.rel.lambda(:) .* ones(nl, 1);  P.r_rep = P.rel.r_rep(:) .* ones(nl, 1);
P.Ncust  = P.rel.Ncust(:);
P.isl.w  = P.isl.w(:);
if (P.use.EVCS || P.use.BSS || P.use.DRM) && T == 1
    warning('RC3_setup:T1', 'EVCS/BSS/DRM with T = 1 is a peak-hour snapshot; use T = 24 for scheduling.');
end
P.dr_w = P.mult * P.dt;                                    % DR energy-neutrality weights

%% ---- PCC status: forced values (NaN = decision variable) -------------------------------
out = false(1, T);
if numel(P.isl.outage) == 2
    out = hrs >= P.isl.outage(1) & hrs <= P.isl.outage(2) + 1 - 1e-9;
end
switch mode
    case 'off',      ug = ones(1, T);
    case 'outage',   ug = double(~out);
    case 'islanded', ug = zeros(1, T);
    case 'optimal',  ug = nan(1, T);  ug(out) = 0;
end
P.ugForced = ug;

%% ---- fundamental loops of the base network ------------------------------------------------
[par, pbr] = RC3_tree(P.fb, P.tb, x0, nb);
P.loops = cell(1, numel(P.tie));
for k = 1:numel(P.tie)
    b = P.tie(k);
    P.loops{k} = [setxor(climb(P.fb(b), par, pbr), climb(P.tb(b), par, pbr))'; b];
end
assert(~any(cellfun(@(L) ismember(P.pcc, L), P.loops)), 'The PCC branch must not lie in a loop.');

%% ---- decision vector: one block per variable set ----------------------------------------
% enc: 1 = continuous ('C'), 2 = integer ('I'), 4 = binary ('B')  (PlatEMO / MP-Opt-Model)
B = struct('name', {}, 'idx', {}, 'lo', {}, 'up', {}, 'enc', {}, 'dims', {}, 'desc', {});
nt = numel(P.tie);  o = @(n) ones(1, n);
if P.use.ONR
    if strcmp(P.onr.encoding, 'loop')
        B = addblock(B, 'ONR_open', o(nt*T), repmat(cellfun(@numel, P.loops), 1, T), 2, [nt T], ...
                     'open branch of each fundamental loop (position in loop), per slot');
    else
        B = addblock(B, 'ONR_sw', zeros(1, nl*T), o(nl*T), 4, [nl T], 'status of every switch, per slot');
    end
end
if strcmp(mode, 'optimal') && any(isnan(ug))
    B = addblock(B, 'PCC_u', zeros(1, T), o(T), 4, [1 T], 'PCC breaker status (1 = grid-connected), per slot');
end
P.ncut = 0;
if P.use.ISL && P.isl.enabled && P.isl.ncut > 0
    P.ncut = P.isl.ncut;
    B = addblock(B, 'ISL_cut', zeros(1, P.ncut*T), nl * o(P.ncut*T), 2, [P.ncut T], ...
                 'extra open branch forming microgrids (0 = none), per islanded slot');
end
P.dg.cand = row(P.dg.cand);  P.sc.cand = row(P.sc.cand);  P.st.cand = row(P.st.cand);
P.ev.cand = row(P.ev.cand);  P.bss.cand = row(P.bss.cand);
P.tx.cand = row(P.tx.cand);  P.up.cand = row(P.up.cand);  P.drm.cand = row(P.drm.cand);
if P.use.DG
    n = P.dg.nmax;
    B = addblock(B, 'DG_y',    zeros(1, n), o(n), 4, [n 1], 'DG unit installed');
    B = addblock(B, 'DG_bus',  o(n), numel(P.dg.cand) * o(n), 2, [n 1], 'DG bus (index in candidate list)');
    B = addblock(B, 'DG_type', o(n), numel(P.dg.allowed) * o(n), 2, [n 1], 'DG technology (index in allowed types)');
    B = addblock(B, 'DG_P',    zeros(1, n*T), P.dg.Pmax * o(n*T), 1, [n T], 'DG active power set-point (MW), per slot');
end
if P.use.SC
    n = P.sc.nmax;
    B = addblock(B, 'SC_y',    zeros(1, n), o(n), 4, [n 1], 'capacitor bank installed');
    B = addblock(B, 'SC_bus',  o(n), numel(P.sc.cand) * o(n), 2, [n 1], 'capacitor bus');
    B = addblock(B, 'SC_step', zeros(1, n*T), P.sc.nsteps * o(n*T), 2, [n T], 'capacitor steps in service, per slot');
end
if P.use.TAP
    B = addblock(B, 'TAP_oltc', P.oltc.tap_min * o(T), P.oltc.tap_max * o(T), 2, [1 T], 'OLTC tap position, per slot');
end
if P.use.TX
    n = P.tx.nmax;
    B = addblock(B, 'TX_y',   zeros(1, n), o(n), 4, [n 1], 'in-line transformer installed');
    B = addblock(B, 'TX_br',  o(n), numel(P.tx.cand) * o(n), 2, [n 1], 'transformer branch (index in candidate list)');
    B = addblock(B, 'TX_tap', P.tx.tap_min * o(n*T), P.tx.tap_max * o(n*T), 2, [n T], 'transformer tap position, per slot');
end
if P.use.DSTATCOM
    n = P.st.nmax;
    if strcmp(P.st.mode, 'pv'), qlo = 0; else, qlo = -P.st.Qmax; end
    B = addblock(B, 'ST_y',   zeros(1, n), o(n), 4, [n 1], 'D-STATCOM installed');
    B = addblock(B, 'ST_bus', o(n), numel(P.st.cand) * o(n), 2, [n 1], 'D-STATCOM bus');
    B = addblock(B, 'ST_Q',   qlo * o(n*T), P.st.Qmax * o(n*T), 1, [n T], ...
                 'D-STATCOM reactive injection (pq) or limit (pv) (MVAr), per slot');
end
if P.use.UPFC
    n = P.up.nmax;
    B = addblock(B, 'UP_y',   zeros(1, n), o(n), 4, [n 1], 'UPFC installed');
    B = addblock(B, 'UP_br',  o(n), numel(P.up.cand) * o(n), 2, [n 1], 'UPFC branch');
    B = addblock(B, 'UP_r',   zeros(1, n*T), P.up.r_max * o(n*T), 1, [n T], 'series voltage magnitude ratio, per slot');
    B = addblock(B, 'UP_th',  -pi * o(n*T), pi * o(n*T), 1, [n T], 'series voltage angle (rad), per slot');
    B = addblock(B, 'UP_Qsh', -P.up.Qsh_max * o(n*T), P.up.Qsh_max * o(n*T), 1, [n T], ...
                 'shunt converter reactive injection (MVAr), per slot');
end
if P.use.EVCS
    n = P.ev.nmax;
    B = addblock(B, 'EV_y',    zeros(1, n), o(n), 4, [n 1], 'EV station installed');
    B = addblock(B, 'EV_bus',  o(n), numel(P.ev.cand) * o(n), 2, [n 1], 'EV station bus');
    B = addblock(B, 'EV_nch',  o(n), P.ev.nch_max * o(n), 2, [n 1], 'number of chargers');
    B = addblock(B, 'EV_type', o(n), numel(P.ev.rates) * o(n), 2, [n 1], 'charger type');
    B = addblock(B, 'EV_cap',  P.ev.cap_min * o(n*T), o(n*T), 1, [n T], 'charging-power cap (p.u.), per slot');
end
if P.use.BSS
    n = P.bss.nmax;
    B = addblock(B, 'BSS_y',    zeros(1, n), o(n), 4, [n 1], 'BSS installed');
    B = addblock(B, 'BSS_bus',  o(n), numel(P.bss.cand) * o(n), 2, [n 1], 'BSS bus');
    B = addblock(B, 'BSS_nch',  o(n), P.bss.nch_max * o(n), 2, [n 1], 'number of chargers');
    B = addblock(B, 'BSS_type', o(n), numel(P.bss.rates) * o(n), 2, [n 1], 'charger type');
    if strcmp(P.bss.control, 'rule')
        B = addblock(B, 'BSS_ach',  P.bss.ach(1) * o(n*T),  P.bss.ach(2) * o(n*T),  1, [n T], 'charge threshold, per slot');
        B = addblock(B, 'BSS_adis', P.bss.adis(1) * o(n*T), P.bss.adis(2) * o(n*T), 1, [n T], 'discharge threshold, per slot');
    else
        B = addblock(B, 'BSS_u', -o(n*T), o(n*T), 1, [n T], 'pool power set-point (+ charge / - B2G, p.u.), per slot');
    end
end
if P.use.DRM
    if T > 1
        B = addblock(B, 'DR_shift', -P.drm.rho_max * o(T), P.drm.rho_max * o(T), 1, [1 T], ...
                     'shifted share of the enrolled load, per slot');
    end
    B = addblock(B, 'DR_il', zeros(1, T), P.drm.chi_max * o(T), 1, [1 T], 'curtailed share of the enrolled load, per slot');
    if P.drm.site_select
        n = numel(P.drm.cand);
        B = addblock(B, 'DR_z', zeros(1, n), o(n), 4, [n 1], 'bus enrolled in the DR programme');
    end
end
assert(~isempty(B), 'Select at least one device.');
P.blocks = B;
P.lower = [B.lo];  P.upper = [B.up];  P.enc = [B.enc];  P.D = numel(P.lower);
P.bidx = struct();
for b = B, P.bidx.(b.name) = b.idx; end

%% ---- base-case decision vector ------------------------------------------------------------
xb = zeros(1, P.D);
for b = B
    switch b.name
        case 'ONR_open', xb(b.idx) = repmat(cellfun(@numel, P.loops), 1, T);   % tie switches open
        case 'ONR_sw',   xb(b.idx) = repmat(double(x0'), 1, T);
        case 'PCC_u',    xb(b.idx) = 1;
        case {'EV_cap'}, xb(b.idx) = 1;
        case {'DG_bus','DG_type','SC_bus','TX_br','ST_bus','UP_br','EV_bus','EV_nch', ...
              'EV_type','BSS_bus','BSS_nch','BSS_type','BSS_ach','BSS_adis'}
            xb(b.idx) = b.lo;                  % irrelevant: unit not installed
        otherwise,       xb(b.idx) = 0;        % y = 0, set-points, taps, steps, DR = 0
    end
end
P.xbase = xb;

%% ---- linear constraint sets (registered in MP-Opt-Model) -----------------------------------
L = struct('name', {}, 'A', {}, 'l', {}, 'u', {}, 'desc', {}, 'repaired', {});
cnt = {'DG','DG_y','dg';  'SC','SC_y','sc';  'TX','TX_y','tx';  'DSTATCOM','ST_y','st'; ...
       'UPFC','UP_y','up'; 'EVCS','EV_y','ev'; 'BSS','BSS_y','bss'};
for k = 1:size(cnt, 1)
    if P.use.(cnt{k, 1})
        d = P.(cnt{k, 3});
        A = sparse(1, P.bidx.(cnt{k, 2}), 1, 1, P.D);
        L(end+1) = struct('name', ['cnt_' cnt{k, 3}], 'A', A, 'l', d.nmin, 'u', d.nmax, ...
                          'desc', sprintf('number of %s units', cnt{k, 1}), 'repaired', true); %#ok<AGROW>
    end
end
if isfield(P.bidx, 'TAP_oltc') && T > 1
    A = ramp_matrix(P.bidx.TAP_oltc, 1, T, P.D);
    r = P.oltc.ramp * ones(size(A, 1), 1);
    L(end+1) = struct('name', 'ramp_oltc', 'A', A, 'l', -r, 'u', r, ...
                      'desc', 'OLTC tap change per slot', 'repaired', true);
end
if isfield(P.bidx, 'TX_tap') && T > 1
    A = ramp_matrix(P.bidx.TX_tap, P.tx.nmax, T, P.D);
    r = P.tx.ramp * ones(size(A, 1), 1);
    L(end+1) = struct('name', 'ramp_tx', 'A', A, 'l', -r, 'u', r, ...
                      'desc', 'transformer tap change per slot', 'repaired', true);
end
if isfield(P.bidx, 'DR_shift')
    A = sparse(1, P.bidx.DR_shift, P.dr_w, 1, P.D);
    L(end+1) = struct('name', 'dr_neutral', 'A', A, 'l', 0, 'u', 0, ...
                      'desc', 'energy neutrality of load shifting', 'repaired', true);
end
if isfield(P.bidx, 'DR_il')
    A = sparse(1, P.bidx.DR_il, P.dr_w, 1, P.D);
    L(end+1) = struct('name', 'dr_il_budget', 'A', A, 'l', -Inf, 'u', P.drm.il_frac * sum(P.dr_w), ...
                      'desc', 'daily curtailed-energy budget', 'repaired', true);
end
if isfield(P.bidx, 'BSS_ach')
    n = P.bss.nmax * T;
    A = sparse(1:n, P.bidx.BSS_adis, 1, n, P.D) - sparse(1:n, P.bidx.BSS_ach, 1, n, P.D);
    L(end+1) = struct('name', 'bss_gap', 'A', A, 'l', P.bss.gap * ones(n, 1), 'u', Inf(n, 1), ...
                      'desc', 'alpha_dis - alpha_ch >= gap', 'repaired', true);
end
P.linsets = L;

%% ---- constraint bookkeeping and cache key ------------------------------------------------
P.conNl = {'g_rad','g_V','g_S','g_ops','g_pen','g_site','g_soc','g_ev','g_isl','g_f','g_rev'};
P.nConNl = numel(P.conNl);
P.conNames = [P.conNl, cellfun(@(s) ['lin_' s], {L.name}, 'UniformOutput', false)];
P.nCon = numel(P.conNames);
nSetup = nSetup + 1;
P.key  = sprintf('%d_%d_%.0f', nSetup, P.D, 1e6 * posixtime(datetime('now')));

%% ---- normalisation references of the basic objectives -------------------------------------
%  ref_k = |F_k(x0)| (base case) when it is a usable value, otherwise the
%  default scale S_k of default_scales(); opts.objref overrides both.
P.objScale = default_scales(P);
[~, ~, ~, Fb] = RC3_objectives(P.xbase, P);               % all basic objectives, base case
ref  = P.objScale;
good = isfinite(Fb) & abs(Fb) > 1e-6 * P.objScale & Fb < P.bigM & ~ismember(P.allObj, P.scaleOnly);
ref(good) = abs(Fb(good));
ovr = getopt(opts, 'objref', struct());
for f = fieldnames(ovr)'
    k = find(strcmpi(P.allObj, f{1}), 1);
    assert(~isempty(k), 'opts.objref: unknown objective %s', f{1});
    ref(k) = ovr.(f{1});
end
P.objRef = ref;  P.Fbase_all = Fb;

%% ---- user objective functions (opts.objectives) -------------------------------------------
P.objSpec = getopt(opts, 'objectives', {'LOSS', 'VD'});
P.objdef  = RC3_objdef(P.objSpec, P);
lib = RC3_objlib();
P.objLib  = RC3_objdef({lib.name}, P);                     % f1..f5, always reported
P.objNames  = {P.objdef.name};
P.objLabels = {P.objdef.label};
P.objectives = P.objNames;                                % names used by the solvers/tables
P.M = numel(P.objdef);
P.m = getopt(opts, 'm', P.M);
assert(P.m == 1 || P.m == P.M, ['opts.m must be 1 (weighted single objective) or ' ...
       'equal to numel(opts.objectives) = %d.'], P.M);
P.nObj = P.m;
w = row(getopt(opts, 'weights', ones(1, P.M) / P.M));
assert(numel(w) == P.M, 'opts.weights needs one weight per objective.');
P.weights = w;
Fu = RC3_objfun(Fb, P);                                    % user objectives, base case
fr = abs(Fu);  fr(~isfinite(fr) | fr < 1e-9 | fr >= P.bigM) = 1;
P.Fref = fr;                                              % normalisation for m = 1

%% ---- MP-Opt-Model formulation ---------------------------------------------------------------
[P.om, P.om_api] = RC3_om_build(P);
[A, l, u] = RC3_om_eval(P, 'lin');                       % linear constraints read back
P.lin.A = A;  P.lin.l = l;  P.lin.u = u;
P.lin.rows = cell(1, numel(L));  r0 = 0;
for k = 1:numel(L)
    P.lin.rows{k} = r0 + (1:size(L(k).A, 1));  r0 = r0 + size(L(k).A, 1);
end
end

% =========================================================================================
function S = default_scales(P)
%DEFAULT_SCALES  Scale S_k of every basic objective (same unit as F_k), used
%to normalise an objective whose base-case value is zero or not meaningful.
mu   = mean(P.mult);
Eday = 1000 * P.dt * sum(P.Ptot * P.mult);                  % daily load energy (kWh)
Cday = 1000 * P.dt * sum(P.price .* P.Ptot .* P.mult);      % daily energy cost (USD)
Eev  = 0;  Esw = 0;
if P.use.EVCS, Eev = P.ev.N * P.ev.E; end                   % EV energy demand (kWh/day)
if P.use.BSS,  Esw = P.bss.N * P.bss.Ebat * (P.bss.soc_full - P.bss.soc_arr); end
Lenr = P.mult * sum(P.drm.part(:) .* P.Pd);                 % enrolled DR load (MW)
v.TAC      = 365 * Cday;
v.INVCOST  = v.TAC / P.crf;     v.COST = v.INVCOST;
v.SWCOST   = Cday;              v.OPCOST = Cday;
v.DEGCOST  = P.bss.c_bat * P.bss.nmax * P.bss.nch_max * P.bss.bpc * P.bss.Ebat / (365 * P.bss.deg.cal_life);
v.ENS      = sum(P.lambda .* P.r_rep) * 1000 * P.Ptot * mu;
v.SAIDI    = sum(P.lambda .* P.r_rep);
v.SAIFI    = sum(P.lambda);
v.LOSS     = 0.05 * 1000 * P.Ptot * mu;
v.VD       = (P.nb - 1) * 0.05^2;
v.VSI      = 1;
v.PEAK     = P.Ptot * max(P.mult);
v.LVAR     = max(var(P.Ptot * P.mult, 1), (0.01 * P.Ptot)^2);
v.PAR      = 1;
v.EMISSION = P.econ.ef_grid * P.Ptot * sum(P.mult) * P.dt;
v.EVUNMET  = Eev + Esw;   v.EVCSUN = Eev;   v.BSSUN = Esw;
v.EVWAIT   = P.ev.w_pen;
v.DRCOST   = P.dt * sum((P.drm.c_shift * P.drm.rho_max + P.drm.c_il * P.drm.chi_max) * Lenr);
v.DISCOMF  = 1;
v.LSHED    = Eday;
v.WLSHED   = 1000 * P.dt * sum(P.mult) * sum(P.isl.w .* P.Pd);
v.FDEV     = P.isl.df_max;
v.SAI      = 1;
v.SAD      = 1;
S = cellfun(@(n) v.(n), P.allObj);
S(~(S > 0) | ~isfinite(S)) = 1;
end

function B = addblock(B, name, lo, up, enc, dims, desc)
if isempty(B), s = 0; else, s = B(end).idx(end); end
B(end+1) = struct('name', name, 'idx', s + (1:numel(lo)), 'lo', lo, 'up', up, ...
                  'enc', enc * ones(1, numel(lo)), 'dims', dims, 'desc', desc);
end

function A = ramp_matrix(idx, n, T, D)
% rows: x(k,t) - x(k,t-1), k = 1..n, t = 2..T  (block stored unit-major per slot)
r = 0;  I = [];  J = [];  V = [];
for t = 2:T
    for k = 1:n
        r = r + 1;
        I = [I, r, r];  J = [J, idx((t-1)*n + k), idx((t-2)*n + k)];  V = [V, 1, -1]; %#ok<AGROW>
    end
end
A = sparse(I, J, V, r, D);
end

function v = getopt(s, name, default)
if isfield(s, name) && ~isempty(s.(name)), v = s.(name); else, v = default; end
end

function s = setpath(s, path, val)
parts = strsplit(path, '.');
s = subsasgn(s, struct('type', repmat({'.'}, 1, numel(parts)), 'subs', parts), val);
end

function p = climb(a, par, pbr)
p = [];
while a ~= 1, p(end+1) = pbr(a); a = par(a); end %#ok<AGROW>
end
