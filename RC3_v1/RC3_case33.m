function mpc = RC3_case33(casename)
%RC3_CASE33  IEEE 33-bus feeder in MATPOWER case format, extended with the
%data of every optional device in extra fields of the mpc struct.
%
%   mpc = RC3_case33()            embedded IEEE 33-bus data (= case33bw)
%   mpc = RC3_case33('case33bw')  network read with MATPOWER's loadcase, then
%                                the device fields below are added
%
%   Standard MATPOWER fields (MATPOWER User's Manual, Appendix B):
%     version, baseMVA, bus, gen, branch, gencost
%     12.66 kV, 10 MVA base, Baran & Wu (1989) data, loads in MW / MVAr,
%     impedances in p.u.; branches 33-37 are the normally open tie lines
%     (BR_STATUS = 0), bus 1 is the substation (REF) bus.
%
%   Extra fields (ignored by runpf/runopf, preserved by loadcase/savecase):
%     mpc.profile   24-h profiles: load, PV, wind, price, EV arrivals, swaps
%     mpc.onr       optimal network reconfiguration (switches and costs)
%     mpc.pcc       point of common coupling (main-grid breaker)
%     mpc.island    islanded operation (droop, frequency, priorities, VOLL)
%     mpc.dg        distributed generators: PV, wind (WT), microturbine (MT)
%     mpc.sc        switched shunt capacitor banks
%     mpc.oltc      on-load tap changer of the substation transformer
%     mpc.tx        in-line tap-changing transformers (step voltage regulators)
%     mpc.dstatcom  distribution static compensators
%     mpc.upfc      unified power flow controllers
%     mpc.evcs      electric-vehicle charging stations
%     mpc.bss       battery swapping stations (battery pool = equivalent BESS)
%     mpc.drm       demand response management (load shifting + interruptible)
%     mpc.rel       reliability data (failure rate, repair time, customers)
%     mpc.econ      economic and environmental data
%     mpc.limits    operating limits and load-flow settings
%
%   Every value can be changed before calling RC3_setup, e.g.
%       mpc = RC3_case33();  mpc.dg.nmax = 4;  mpc.island.outage = [17 22];
%       opts.mpc = mpc;     R = RC3_solve(24, opts);
%   or through opts.set = {'dg.nmax', 4; 'island.outage', [17 22]}.

if nargin < 1 || isempty(casename)
    mpc = case33_embedded();
else
    mpc = loadcase(casename);                 % requires MATPOWER on the path
end
nb = size(mpc.bus, 1);  
nl = size(mpc.branch, 1);
Pd = mpc.bus(:, 3);                           % MW
loadBus = find(Pd > 0);

%% ---- 24-h profiles (hour 1 = 00:00-01:00) -----------------------------------------
mpc.profile.load  = [0.55 0.50 0.47 0.45 0.45 0.48 0.55 0.65 0.72 0.76 0.80 0.83 ...
                     0.85 0.84 0.82 0.80 0.82 0.88 0.96 1.00 0.97 0.88 0.76 0.64];
mpc.profile.pv    = [0 0 0 0 0 0.02 0.12 0.30 0.50 0.68 0.82 0.92 ...
                     0.95 0.90 0.80 0.64 0.44 0.22 0.05 0 0 0 0 0];
mpc.profile.wt    = [0.62 0.66 0.70 0.71 0.69 0.65 0.58 0.51 0.44 0.38 0.34 0.31 ...
                     0.30 0.32 0.35 0.40 0.47 0.53 0.58 0.62 0.65 0.67 0.66 0.64];
mpc.profile.price = [0.08*ones(1,7) 0.12*ones(1,10) 0.20*ones(1,5) 0.12*ones(1,2)];  % $/kWh
mpc.profile.ev    = [1.0 0.8 0.6 0.5 0.5 1.0 2.5 4.5 5.5 5.0 4.5 4.5 ...
                     5.0 5.0 5.0 5.5 6.5 8.0 8.5 7.5 6.0 4.5 3.0 1.6];   % relative EV arrivals
mpc.profile.swap  = [2.0 1.5 1.0 1.0 1.0 2.0 4.0 6.0 6.5 5.5 5.0 5.0 ...
                     5.0 5.0 5.0 5.5 6.0 6.5 7.0 6.5 5.5 4.5 3.5 2.5];   % relative swap requests
mpc.profile.peak_hour = 20;                   % hour used when T = 1 (peak snapshot)

%% ---- ONR: sectionalizing switches (closed) and tie switches (open) -----------------
mpc.onr.c_sec    = 20;                        % $ per sectionalizing-switch operation
mpc.onr.c_tie    = 25;                        % $ per tie-switch operation
mpc.onr.max_ops  = Inf;                       % operations per switch and day
mpc.onr.encoding = 'loop';                    % 'loop' (one open branch per loop) or 'binary'

%% ---- PCC: breaker between the main grid and the feeder ------------------------------
mpc.pcc.branch   = 1;                         % branch 1-2 is the PCC breaker
mpc.pcc.c_op     = 50;                        % $ per opening/closing of the PCC
mpc.pcc.Pexp_max = Inf;                       % max. export to the main grid (MW)

%% ---- islanded operation ------------------------------------------------------------
mpc.island.mode      = 'off';                 % 'off' | 'outage' | 'optimal' | 'islanded'
mpc.island.outage    = [18 21];               % main grid unavailable from hour 18 to 21
mpc.island.ncut      = 2;                     % extra open switches to form microgrids
mpc.island.max_trans = 2;                     % PCC openings + closings per day
mpc.island.f0        = 50;                    % nominal frequency (Hz)
mpc.island.df_max    = 0.5;                   % P-f droop: deviation at full GF loading (Hz)
mpc.island.f_lim     = 0.5;                   % allowed frequency deviation (Hz)
mpc.island.V_nom     = 1.0;                   % no-load voltage of the Q-V droop (p.u.)
mpc.island.dV_max    = 0.05;                  % Q-V droop: deviation at rated Q (p.u.)
mpc.island.q_frac    = 0.75;                  % reactive capability of GF units (x rating)
mpc.island.voll      = 5.0;                   % value of lost load ($/kWh)
mpc.island.min_serve = 0;                     % min. served share of the islanded demand
mpc.island.w = ones(nb, 1);                   % load priority weights
if nb == 33
    mpc.island.w([7 14 24 30])  = 10;         % critical loads
    mpc.island.w([8 25 31 32])  = 5;          % semi-critical loads
end
mpc.island.w(1) = 0;

%% ---- DG: PV, wind turbine (WT), microturbine (MT, dispatchable, grid-forming) -------
mpc.dg.nmax    = 10;                           % candidate units (0..nmax are installed)
mpc.dg.nmin    = 0;
mpc.dg.cand    = 2:nb;                        % candidate buses
mpc.dg.Pmax    = 1.0;                         % MW per unit
mpc.dg.types   = {'PV', 'WT', 'MT'};
mpc.dg.profile = {'pv', 'wt', 'const'};       % availability profile of each type
mpc.dg.allowed = [1 2 3];                     % types the optimizer may choose
mpc.dg.pf      = [1.0 1.0 0.9];               % power factor (lagging, injection)
mpc.dg.gf      = [0 0 1];                     % grid-forming capability in islands
mpc.dg.c_cap   = [1.00e6 1.30e6 0.90e6];      % $ per MW installed
mpc.dg.c_fuel  = [0 0 85];                    % $ per MWh generated
mpc.dg.ef      = [0 0 600];                   % kg CO2 per MWh generated
mpc.dg.pen     = 0.6;                         % max. installed DG / total load

%% ---- SC: switched capacitor banks ------------------------------------------------------
mpc.sc.nmax   = 3;  
mpc.sc.nmin = 0;
mpc.sc.cand   = 2:nb;
mpc.sc.Qstep  = 0.15;                         % MVAr per step
mpc.sc.nsteps = 8;                            % steps per bank
mpc.sc.c_fix  = 3e3;                          % $ per bank
mpc.sc.c_var  = 25e3;                         % $ per MVAr

%% ---- OLTC of the substation transformer ------------------------------------------------
mpc.oltc.tap_min = -8;  mpc.oltc.tap_max = 8;
mpc.oltc.step    = 0.0125;                    % p.u. per tap step
mpc.oltc.c_op    = 5;                         % $ per tap step
mpc.oltc.max_ops = Inf;                       % tap steps per day
mpc.oltc.ramp    = 8;                         % max. tap change per hour (steps)

%% ---- TX: in-line tap-changing transformers (step voltage regulators) --------------------
mpc.tx.nmax    = 1;  mpc.tx.nmin = 0;
mpc.tx.cand    = setdiff(find(mpc.branch(:, 11) > 0)', mpc.pcc.branch);   % candidate branches
mpc.tx.tap_min = -16;  mpc.tx.tap_max = 16;   % 32-step regulator, +/-10 %
mpc.tx.step    = 0.00625;
mpc.tx.c_fix   = 40e3;                        % $ per regulator
mpc.tx.c_op    = 5;                           % $ per tap step
mpc.tx.max_ops = Inf;
mpc.tx.ramp    = 16;                          % max. tap change per hour (steps)

%% ---- D-STATCOM -----------------------------------------------------------------------
mpc.dstatcom.nmax = 1;  mpc.dstatcom.nmin = 0;
mpc.dstatcom.cand = 2:nb;
mpc.dstatcom.Qmax = 1.0;                      % MVAr per unit
mpc.dstatcom.mode = 'pq';                     % 'pq': Q(t) set-point | 'pv': V control
mpc.dstatcom.vset = 1.0;                      % voltage set-point in 'pv' mode (p.u.)
mpc.dstatcom.c_var = 50e3;                    % $ per MVAr

%% ---- UPFC ------------------------------------------------------------------------------
mpc.upfc.nmax    = 1;  mpc.upfc.nmin = 0;
mpc.upfc.cand    = setdiff(find(mpc.branch(:, 11) > 0)', mpc.pcc.branch);
mpc.upfc.r_max   = 0.08;                      % series voltage / sending voltage (p.u.)
mpc.upfc.Qsh_max = 0.5;                       % shunt converter (MVAr)
mpc.upfc.c_curve = [0.0003 -0.2691 188.22];   % $/kVAr = a*S^2 + b*S + c, S in MVA

%% ---- EV charging stations ---------------------------------------------------------------
mpc.evcs.nmax    = 5;  
mpc.evcs.nmin = 2;
mpc.evcs.cand    = 2:nb;
mpc.evcs.nch_max = 20;                        % chargers per station
mpc.evcs.rates   = [7.4 22 50 150];           % kW per charger (types)
mpc.evcs.c_ch    = [1.5e3 5e3 30e3 90e3];     % $ per charger
mpc.evcs.c_site  = 25e3;                      % $ per station
mpc.evcs.N       = 200;                       % EVs per day
mpc.evcs.E       = 18;                        % kWh per charging session
mpc.evcs.eta     = 0.92;  mpc.evcs.pf = 1.0;
mpc.evcs.rho_max = 0.95;                      % max. utilisation of the chargers
mpc.evcs.w_pen   = 60;                        % waiting penalty of a lost EV (min)
mpc.evcs.tariff  = 0.25;                      % $/kWh
mpc.evcs.serve_min = 0;                       % min. served share of the EV energy
mpc.evcs.cap_min = 0.3;                       % min. hourly charging-power cap (p.u.)

%% ---- battery swapping stations -----------------------------------------------------------
mpc.bss.nmax    = 5;  
mpc.bss.nmin = 2;
mpc.bss.cand    = 2:nb;
mpc.bss.nch_max = 20;
mpc.bss.rates   = [11 22 44 60];              % kW per charger
mpc.bss.c_ch    = [2e3 4e3 9e3 12e3];
mpc.bss.c_site  = 150e3;
mpc.bss.bpc     = 2;                          % batteries per charger
mpc.bss.Ebat    = 40;                         % kWh per battery
mpc.bss.c_bat   = 130;                        % $ per kWh
mpc.bss.N       = 120;                        % swaps per day
mpc.bss.soc_arr = 0.2;  mpc.bss.soc_full = 0.9;  mpc.bss.soc0 = 0.7;
mpc.bss.eta     = 0.95;  mpc.bss.crate = 1.0;  mpc.bss.b2g = 1.0;  mpc.bss.pf = 1.0;
mpc.bss.hres    = 2;                          % swap reserve look-ahead (h)
mpc.bss.control = 'schedule';                 % 'schedule': u(t) | 'rule': thresholds(t)
mpc.bss.signal  = 'load';                     % 'rule' signal: 'load' or 'price'
mpc.bss.ach     = [0.40 0.95];                % bounds of alpha_ch(t)
mpc.bss.adis    = [0.60 1.00];                % bounds of alpha_dis(t)
mpc.bss.gap     = 0.05;                       % alpha_dis - alpha_ch >= gap
mpc.bss.tariff  = 0.25;                       % $/kWh swapped
mpc.bss.gf      = 1;                          % grid-forming in islands
mpc.bss.deg = struct('model', 'rainflow', 'N100', 3000, 'kp', 1.5, ...
                     'dod_ref', 0.8, 'cal_life', 15, 'swap', false);

%% ---- demand response management ------------------------------------------------------------
mpc.drm.part        = zeros(nb, 1);
mpc.drm.part(loadBus) = 0.3;                  % enrolled share of each bus load
mpc.drm.cand        = loadBus(:)';            % buses that may enrol
mpc.drm.site_select = false;                  % true: enrolment z_i is a decision
mpc.drm.rho_max     = 0.3;                    % max. shifted share per hour
mpc.drm.chi_max     = 0.2;                    % max. curtailed share per hour
mpc.drm.il_frac     = 0.1;                    % daily curtailed energy / enrolled energy
mpc.drm.c_shift     = 15;                     % $ per MWh shifted (reduction incentive)
mpc.drm.c_il        = 120;                    % $ per MWh curtailed
mpc.drm.c_enrol     = 2;                      % $ per enrolled bus and day

%% ---- reliability, economy, limits --------------------------------------------------------
mpc.rel.lambda = 0.25 * ones(nl, 1);          % failures per year
mpc.rel.r_rep  = 5 * ones(nl, 1);             % repair time (h)
mpc.rel.Ncust  = round(Pd * 1000 / 10);       % customers per bus
mpc.econ.ir = 0.08;  mpc.econ.life = 20;      % interest rate, project life (yr)
mpc.econ.sell_frac = 0.5;                     % export price / import price
mpc.econ.ef_grid   = 500;                     % kg CO2 per MWh imported
mpc.limits.Vmin = 0.90;  mpc.limits.Vmax = 1.05;
mpc.limits.zip_np = 0;   mpc.limits.zip_nq = 0;   % load voltage exponents
mpc.limits.fbs_tol = 1e-8;  mpc.limits.fbs_maxit = 100;
end

% ======================================================================================
function mpc = case33_embedded()
%CASE33_EMBEDDED  Baran & Wu (1989) 33-bus feeder, identical to MATPOWER case33bw.
mpc.version = '2';
mpc.baseMVA = 10;
%   bus  type  Pd(kW)  Qd(kVAr)  Gs  Bs  area  Vm  Va  baseKV  zone  Vmax  Vmin
bd = [ 1 0 0;  2 100 60;  3 90 40;  4 120 80;  5 60 30;  6 60 20;  7 200 100;
       8 200 100;  9 60 20; 10 60 20; 11 45 30; 12 60 35; 13 60 35; 14 120 80;
      15 60 10; 16 60 20; 17 60 20; 18 90 40; 19 90 40; 20 90 40; 21 90 40;
      22 90 40; 23 90 50; 24 420 200; 25 420 200; 26 60 25; 27 60 25; 28 60 20;
      29 120 70; 30 200 600; 31 150 70; 32 210 100; 33 60 40];
nb  = size(bd, 1);
mpc.bus = [bd(:, 1), ones(nb, 1), bd(:, 2:3) / 1000, zeros(nb, 2), ones(nb, 1), ...
           ones(nb, 1), zeros(nb, 1), 12.66 * ones(nb, 1), ones(nb, 1), ...
           1.1 * ones(nb, 1), 0.9 * ones(nb, 1)];
mpc.bus(1, [2 12 13]) = [3 1 1];
%   gen: bus Pg Qg Qmax Qmin Vg mBase status Pmax Pmin (+ 11 zero columns)
mpc.gen = [1 0 0 10 -10 1 100 1 10 0 zeros(1, 11)];
%   branch: from to R(ohm) X(ohm); branches 33-37 = tie lines (open)
br = [ 1  2 0.0922 0.0470;  2  3 0.4930 0.2511;  3  4 0.3660 0.1864;
       4  5 0.3811 0.1941;  5  6 0.8190 0.7070;  6  7 0.1872 0.6188;
       7  8 0.7114 0.2351;  8  9 1.0300 0.7400;  9 10 1.0440 0.7400;
      10 11 0.1966 0.0650; 11 12 0.3744 0.1238; 12 13 1.4680 1.1550;
      13 14 0.5416 0.7129; 14 15 0.5910 0.5260; 15 16 0.7463 0.5450;
      16 17 1.2890 1.7210; 17 18 0.7320 0.5740;  2 19 0.1640 0.1565;
      19 20 1.5042 1.3554; 20 21 0.4095 0.4784; 21 22 0.7089 0.9373;
       3 23 0.4512 0.3083; 23 24 0.8980 0.7091; 24 25 0.8960 0.7011;
       6 26 0.2030 0.1034; 26 27 0.2842 0.1447; 27 28 1.0590 0.9337;
      28 29 0.8042 0.7006; 29 30 0.5075 0.2585; 30 31 0.9744 0.9630;
      31 32 0.3105 0.3619; 32 33 0.3410 0.5302;
      21  8 2.0    2.0   ;  9 15 2.0    2.0   ; 12 22 2.0    2.0   ;
      18 33 0.5    0.5   ; 25 29 0.5    0.5   ];
nl  = size(br, 1);
Zb  = 12.66^2 / mpc.baseMVA;                 % ohm
st  = [ones(32, 1); zeros(5, 1)];
mpc.branch = [br(:, 1:2), br(:, 3:4) / Zb, zeros(nl, 5), zeros(nl, 1), st, ...
              -360 * ones(nl, 1), 360 * ones(nl, 1)];
mpc.gencost = [2 0 0 3 0 20 0];
end
