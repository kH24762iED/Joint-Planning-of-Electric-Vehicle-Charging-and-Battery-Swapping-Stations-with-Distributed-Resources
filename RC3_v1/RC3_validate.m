%% RC3_validate: checks of the load-flow models against published results and
%  against MATPOWER (runpf, Newton-Raphson) for hand-made solutions with every
%  device, in grid-connected and islanded operation.
clear; clc;
hasMP = exist('runpf', 'file') == 2;

%% 1) IEEE 33-bus: base case and minimum-loss configuration (Baran & Wu 1989)
P = RC3_setup(1, struct('devices', {{'ONR'}}, 'encoding', 'binary', 'objectives', {{'LOSS', 'VD'}}));
ref = {'Base case, open 33-37', [33 34 35 36 37], 202.68; ...
       'Minimum loss, open 7 9 14 32 37', [7 9 14 32 37], 139.55};
for k = 1:size(ref, 1)
    x = RC3_set(P, P.xbase, 'ONR_sw', 1);  sw = ones(P.nl, 1);  sw(ref{k, 2}) = 0;
    x = RC3_set(P, x, 'ONR_sw', sw);
    tic;  [~, ~, o] = RC3_objectives(x, P);  ms = 1e3 * toc;
    [vmin, ib] = min(o.Vm);
    fprintf('%-32s loss %8.3f kW (published %.2f)  Vmin %.4f p.u. at bus %2d  %d it  %.2f ms\n', ...
            ref{k, 1}, o.LOSS, ref{k, 3}, vmin, ib, o.fbs_it, ms);
    if hasMP, disp(RC3_validate_solution(P, o)); end
end

%% 2) every device in one grid-connected snapshot (T = 1)
dev = {'ONR','DG','SC','TAP','TX','DSTATCOM','UPFC','EVCS','BSS'};
P = RC3_setup(1, struct('devices', {dev}, 'objectives', {{'LOSS', 'VD'}}));
x = P.xbase;
x = RC3_set(P, x, 'DG_y', 1, 1);   x = RC3_set(P, x, 'DG_bus', 14, 1);  x = RC3_set(P, x, 'DG_type', 'MT', 1);
x = RC3_set(P, x, 'DG_P', 0.6, 1, 1);
x = RC3_set(P, x, 'SC_y', 1, 1);   x = RC3_set(P, x, 'SC_bus', 30, 1);  x = RC3_set(P, x, 'SC_step', 4, 1, 1);
x = RC3_set(P, x, 'TAP_oltc', 2);
x = RC3_set(P, x, 'TX_y', 1, 1);   x = RC3_set(P, x, 'TX_br', 6, 1);    x = RC3_set(P, x, 'TX_tap', 4, 1, 1);
x = RC3_set(P, x, 'ST_y', 1, 1);   x = RC3_set(P, x, 'ST_bus', 25, 1);  x = RC3_set(P, x, 'ST_Q', 0.3, 1, 1);
x = RC3_set(P, x, 'UP_y', 1, 1);   x = RC3_set(P, x, 'UP_br', 26, 1);
x = RC3_set(P, x, 'UP_r', 0.03, 1, 1);  x = RC3_set(P, x, 'UP_th', 1.0, 1, 1);  x = RC3_set(P, x, 'UP_Qsh', 0.1, 1, 1);
x = RC3_set(P, x, 'EV_y', 1, 1);   x = RC3_set(P, x, 'EV_bus', 20, 1);
x = RC3_set(P, x, 'BSS_y', 1, 1);  x = RC3_set(P, x, 'BSS_bus', 31, 1);
x = RC3_repair(x, P);
[F, G, o] = RC3_objectives(x, P);
fprintf('\nAll devices, grid-connected: loss %.3f kW, Vmin %.4f, feasible %d, G = %s\n', ...
        o.LOSS, min(o.Vmin), o.feasible, mat2str(G, 3));
if hasMP, disp(RC3_validate_solution(P, o)); end

%% 3) islanded network: two microgrids (MT at bus 18, BSS at bus 30 / MT at bus 24)
P = RC3_setup(24, struct('devices', {{'DG','BSS','ISL'}}, 'objectives', {{'LSHED', 'FDEV'}}, ...
                        'island_mode', 'outage'));
x = P.xbase;
for u = 1:2
    x = RC3_set(P, x, 'DG_y', 1, u);  x = RC3_set(P, x, 'DG_type', 'MT', u);
    x = RC3_set(P, x, 'DG_P', 0.9, u, 1:P.T);
end
x = RC3_set(P, x, 'DG_bus', 18, 1);  x = RC3_set(P, x, 'DG_bus', 24, 2);
x = RC3_set(P, x, 'BSS_y', 1, 1);    x = RC3_set(P, x, 'BSS_bus', 30, 1);  x = RC3_set(P, x, 'BSS_nch', 10, 1);
x = RC3_set(P, x, 'ISL_cut', 3, 1, 1:P.T);          % open branch 3-4 while islanded
x = RC3_repair(x, P);
[F, G, o] = RC3_objectives(x, P);
hI = find(o.ug < 0.5);
fprintf('\nIslanded slots %s: islands %s, df %s Hz, shed %.1f kWh, feasible %d\n', mat2str(hI), ...
        mat2str(o.nIsl(hI)), mat2str(o.df(hI), 3), o.LSHED, o.feasible);
if hasMP, disp(RC3_validate_solution(P, o, hI)); end

%% 4) MP-Opt-Model formulation of the last model
RC3_om_eval(P, 'display');
[x0, xl, xu, vt] = RC3_om_eval(P, 'var');
fprintf('Variables: %d (C %d, I %d, B %d)\n', numel(x0), nnz(vt == 'C'), nnz(vt == 'I'), nnz(vt == 'B'));

%% 5) library objectives f1..f5 of the islanded test solution and their terms
fprintf('\nLibrary objectives (base -> test solution):\n');
[~, ~, ob] = RC3_objectives(P.xbase, P);
for k = 1:numel(P.objLib)
    fprintf('  %-3s = %-72s %12.5g -> %12.5g\n', P.objLib(k).name, P.objLib(k).expr, ob.Flib(k), o.Flib(k));
end
fprintf('Reference values used for normalisation (P.objRef):\n');
disp(array2table(P.objRef, 'VariableNames', P.allObj));
