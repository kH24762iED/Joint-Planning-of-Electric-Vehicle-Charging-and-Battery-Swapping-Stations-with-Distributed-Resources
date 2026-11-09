%% =====================================================================
%  Multi-period planning and operation of the IEEE 33-bus distribution
%  network (grid-connected and islanded) with any combination of
%     'ONR' 'DG' 'SC' 'TAP' 'TX' 'DSTATCOM' 'UPFC' 'EVCS' 'BSS' 'DRM' 'ISL'
%
%  Basic objectives (P.allObj, always computed):
%     'SWCOST' 'INVCOST' 'COST' 'DEGCOST' 'OPCOST' 'TAC' 'ENS' 'SAIDI' 'SAIFI'
%     'LOSS' 'VD' 'VSI' 'PEAK' 'LVAR' 'PAR' 'EMISSION' 'EVUNMET' 'EVWAIT'
%     'DRCOST' 'DISCOMF' 'LSHED' 'WLSHED' 'FDEV' 'SAI' 'SAD' 'EVCSUN' 'BSSUN'
%  User objective functions (opts.objectives), see RC3_objlib / RC3_objdef:
%     f1 = TAC
%     f2 = 0.3*LOSS + 0.3*VD + 0.4*VSI
%     f3 = 0.4*EMISSION + 0.2*EVCSUN + 0.15*EVWAIT + 0.15*BSSUN + 0.1*DEGCOST
%     f4 = 0.3*ENS + 0.35*SAIDI + 0.35*SAIFI
%     f5 = 0.3*WLSHED + 0.15*FDEV + 0.15*SAD + 0.15*DRCOST + 0.1*DISCOMF + 0.15*PAR
%     (each basic objective in a weighted sum is divided by its reference
%     value P.objRef: base-case value, or default scale if that is zero)
%     or any expression, e.g. 'myF = 0.5*LOSS + 0.5*ENS', or a function
%     handle @(V, Vn) ... of the raw (V) and normalised (Vn) values.
%  m = 1: single objective  sum_k w_k f_k / f_k(base)  (opts.weights)
%  m = numel(objectives) = 2, 3, 4, 5: multi-objective (Pareto set)
%
%  Data     : MATPOWER case + device fields (RC3_case33)
%  Model    : MP-Opt-Model (RC3_om_build): variables, linear and nonlinear
%             constraints, one nonlinear cost per user objective
%  Solver   : PlatEMO 4.x (GA / NSGA-II / NSGA-III) or built-in RC3_nsga2;
%             both use the same initial population (RC3_init_pop), repair
%             (RC3_repair) and evaluation (RC3_evaluate) - see RC3_check_solvers
%  Check    : MATPOWER runpf of the best solution (RC3_validate_solution)
%% =====================================================================
clear; clc; close all;
addpath(genpath(cd));

%% ---- 1) data: MATPOWER case with all device parameters -------------------------
mpc = RC3_case33();                    % edit mpc.<device> fields here, e.g.
mpc.dg.nmax  = 3;                     % up to 3 DGs (0..3 are installed)
mpc.dg.allowed = [1 3];               % PV and microturbine
mpc.island.outage = [18 21];          % main-grid outage window (island_mode 'outage')

%% ---- 2) study settings -------------------------------------------------------------
T                = 24;                % 24 = hourly day, 1 = peak-hour snapshot
opts.mpc         = mpc;
opts.devices     = {'ONR','DG','SC','TAP','TX','DSTATCOM','UPFC','EVCS','BSS','DRM','ISL'};
opts.objectives  = {'f1', 'f2', 'f3', 'f4', 'f5'};   % library (RC3_objlib)
opts.m           = 5;                 % 5 objectives -> Pareto set (NSGA-III)
% opts.m       = 1;  opts.weights = [0.2 0.2 0.2 0.2 0.2];   % weighted single objective
% opts.objectives = {'f1', 'tech = 0.5*LOSS + 0.5*VD', ...     % own definitions
%                    struct('name', 'rel', 'expr', '0.3*ENS + 0.7*SAIDI', 'label', 'Reliability')};
% opts.objref  = struct('LOSS', 202.68);                     % own reference values
opts.island_mode = 'outage';          % 'off' | 'outage' | 'optimal' | 'islanded'
opts.encoding    = 'loop';            % ONR encoding: 'loop' or 'binary'

%% ---- 3) solver ---------------------------------------------------------------------
opts.solver   = 'auto';               % 'platemo' | 'builtin' | 'auto'
opts.N        = 100;                  % population (NSGA-III adjusts it for m = 5)
opts.maxFE    = 30000;                % evaluations (raise for T = 24 with many devices)
opts.seed     = 1;                    % same initial population for both solvers
opts.validate = true;                 % MATPOWER runpf check of the best solution
opts.refine   = false;                % m = 1 only: MP-Opt-Model NLP polish

%% ---- 4) solve, report, plot, save ----------------------------------------------------
R    = RC3_solve(T, opts);
RC3_om_eval(R.P, 'display');           % MP-Opt-Model formulation summary
Tb   = RC3_report(R);                  % Tb.objectives, Tb.library, Tb.objterms, ...
figs = RC3_plot_results(R);
RC3_save_results(R, [], Tb, figs);
%% ----- 5) RC experiments
clear; clc; close all;
addpath(genpath(cd));
rng(1);
study.T=24;
study.N=100;
study.maxFE=600000;
study.runs=1;
study.solver='auto';
study.mpc=RC3_case33();
% study.folder=()
E = RC3_experiments('L_M', study);
test=1;