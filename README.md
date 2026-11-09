# IEEE 33-bus planning and operation model (grid-connected and islanded) — version 3.0 (RC3_ files)

All MATLAB files and functions carry the prefix `RC3_` (e.g. `RC3_main.m`, `RC3_solve.m`, class
`RC3_PlatEMO`), so this version can be on the MATLAB path together with the earlier `RC_` and `RC2_` files.

Devices: ONR, DG (PV/WT/MT), SC, OLTC (`TAP`), in-line transformers (`TX`), D-STATCOM, UPFC,
EVCS, BSS, demand response (`DRM`), islanding decisions (`ISL`).
Data in a MATPOWER case (`RC3_case33.m`), model in MP-Opt-Model (`RC3_om_build.m`), solved with
PlatEMO (GA / NSGA-II / NSGA-III) or the built-in `RC3_nsga2.m`. The formulation, notation and
experimental setup are in `RC3_formulation.docx` (version 3.0).

## What changed in 3.0
- All files, functions and the PlatEMO problem class carry the prefix `RC3_` (`RC3_main.m`,
  `RC3_solve.m`, class `RC3_PlatEMO`, ...).
- Group L of `RC3_experiments.m`: comparison of recent large-scale single- and multi-objective EAs
  of PlatEMO on the all-device problem (hypervolume, best objective, Wilcoxon rank-sum test).
- `opts.alg` accepts a function handle, a name or a cell with algorithm parameters.


## Requirements
- MATLAB R2019b or later
- MATPOWER 8.1 (includes MP-Opt-Model 5.0 and MIPS); MATPOWER 7.1+ works with the legacy `opt_model`
- PlatEMO 4.x (optional; without it the built-in solver is used)
- Copy the folder "1.RC3_v1.m" into the Problem folder of Platemo


## Quick start
```matlab
addpath(genpath('<matpower folder>'));  addpath(genpath('<PlatEMO folder>'));  addpath('<RC3 folder>');
RC3_validate          % load-flow checks (published 33-bus results, MATPOWER runpf), f1..f5 of a test solution
C = RC3_check_solvers % PlatEMO and built-in evaluate identically (differences should be 0)
RC3_main              % one case: edit devices, objectives, m, T, island mode at the top
E = RC3_experiments('all', struct('T', 24, 'N', 100, 'maxFE', 30000, 'runs', 10));
L = RC3_experiments({'L'}, struct('T', 24, 'N', 100, 'maxFE', 60000, 'runs', 10));  % algorithms
```

## Objective functions (`RC3_objlib.m`)
| name | definition (each basic objective divided by its reference value `P.objRef`) |
|---|---|
| f1 | `TAC` (USD/yr, raw) |
| f2 | `0.3*LOSS + 0.3*VD + 0.4*VSI` |
| f3 | `0.4*EMISSION + 0.2*EVCSUN + 0.15*EVWAIT + 0.15*BSSUN + 0.1*DEGCOST` |
| f4 | `0.3*ENS + 0.35*SAIDI + 0.35*SAIFI` |
| f5 | `0.3*WLSHED + 0.15*FDEV + 0.15*SAD + 0.15*DRCOST + 0.1*DISCOMF + 0.15*PAR` |

Basic objectives: `SWCOST INVCOST COST DEGCOST OPCOST TAC ENS SAIDI SAIFI LOSS VD VSI PEAK LVAR PAR
EMISSION EVUNMET EVWAIT DRCOST DISCOMF LSHED WLSHED FDEV SAI SAD EVCSUN BSSUN`.

## Main options (`opts`)
| field | values |
|---|---|
| `devices` | any of `'ONR' 'DG' 'SC' 'TAP' 'TX' 'DSTATCOM' 'UPFC' 'EVCS' 'BSS' 'DRM' 'ISL'` |
| `objectives` | `{'f1','f2',...}`, basic names, expressions `'name = 0.5*LOSS + 0.5*VD'`, structs, handles |
| `objref` | struct of reference values overriding `P.objRef`, e.g. `struct('LOSS', 202.68)` |
| `m`, `weights` | 1 (weighted single objective) or `numel(objectives)` |
| `island_mode` | `'off'` `'outage'` `'optimal'` `'islanded'` |
| `T` (argument) | 24 (hourly) or 1 (peak snapshot) |
| `mpc`, `set` | case struct / `{'dg.nmax', 4; ...}` overrides |
| `solver`, `alg`, `N`, `maxFE`, `seed` | `'auto'`/`'platemo'`/`'builtin'`, PlatEMO algorithm (`@LMOCSO`, `'LMOCSO'`, `{@FDV, 0.8, 0.4, 5}`), population, budget, seed |
| `validate`, `refine` | MATPOWER check of the best solution; MIPS polish (m = 1) |

Results are written to `results/<case>_<date>/` (MAT, Excel/CSV tables incl. `objectives`,
`library`, `objterms`, `check`; PNG/FIG figures incl. `objectives`; summary).
