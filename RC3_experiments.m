function E = RC3_experiments(groups, study)
%RC3_EXPERIMENTS  Experimental study with the weighted-sum objective functions
%of RC3_objlib:
%   f1 = TAC
%   f2 = 0.3*LOSS + 0.3*VD + 0.4*VSI                                (technical)
%   f3 = 0.4*EMISSION + 0.2*EVCSUN + 0.15*EVWAIT + 0.15*BSSUN + 0.1*DEGCOST
%                                                     (environment + e-mobility)
%   f4 = 0.3*ENS + 0.35*SAIDI + 0.35*SAIFI                          (reliability)
%   f5 = 0.3*WLSHED + 0.15*FDEV + 0.15*SAD + 0.15*DRCOST + 0.1*DISCOMF + 0.15*PAR
%                                                 (islanding + demand response)
%Each basic objective in f2..f5 is divided by its reference value (base
%case, or default scale), so f2..f5 are dimensionless and equal to the sum
%of the weights of their active terms in the base case.
%
%   E = RC3_experiments()                    all groups, default budget
%   E = RC3_experiments({'A','W'})           selected groups
%   E = RC3_experiments('all', study)        study options:
%       .T (24) .N (60) .maxFE (12000) .runs (1) .seed0 (1) .solver ('auto')
%       .mpc (RC3_case33()) .folder (results/experiments_<date>)
%       .objsets  objective sets of group D for m = 1..5
%                 (default {f1..f5}, {f1,f2}, {f1,f2,f4}, {f1,f2,f3,f4}, {f1..f5})
%       .weights  weight scenarios of group W (rows, 5 columns)
%
%   Group A  each device alone, grid-connected, objectives relevant to it
%   Group B  combinations of devices, grid-connected
%   Group C  islanded operation: outage, microgrid formation, intentional
%            (optimal) islanding, standalone network
%   Group D  all devices, outage 18-21 h: m = 1 (equal-weight sum of f1..f5),
%            m = 2, 3, 4, 5 (Pareto sets of f1..f5 subsets)
%   Group W  all devices, m = 1: weight scenarios (economy, technical,
%            environment/e-mobility, reliability, resilience/DR)
%   Group E  peak-hour snapshot (T = 1) versus 24-h profile (T = 24)
%   Group L  all devices, outage: comparison of recent large-scale and
%            constrained EAs of PlatEMO on the same problems (needs PlatEMO)
%            L-S  m = 1 (equal-weight sum of f1..f5): study.algs_single
%                 (default GA, CSO, IMODE, SHADE, AESSPSO)
%            L-M3 m = 3 (f1, f2, f4) and L-M5 m = 5 (f1..f5): study.algs_multi
%                 (default NSGAII, NSGAIII, LMOCSO, CMOCSO, IMTCMO_BS, DSSEA,
%                 LCMEA, ILCMO, PRCEA); the first algorithm of each list is
%                 the reference of the Wilcoxon rank-sum test
%            Indicators: normalised hypervolume (m > 1) or best objective
%            (m = 1) over study.runs runs, feasible runs, time, convergence
%
%   For every case the five library objectives f1..f5 and all basic
%   objectives of the base case and of the best solution are tabulated, so
%   cases optimised for different objectives can be compared on all of them.
%   E.cases  case definitions     E.table  best solutions (f1..f5 + indicators)
%   E.base   base cases           E.rel    f_k(best) / f_k(base)
%   E.stats  mean / std / best of f1..f5 over the runs of each case
%   E.L      group L: algorithm comparison table (sheet 'algorithms')
%   'all' runs groups A-E and W; group L is run only when requested,
%   e.g. RC3_experiments({'L'}, struct('runs', 10, 'maxFE', 60000)).
if nargin < 1 || isempty(groups) || (ischar(groups) && strcmpi(groups, 'all'))
    groups = {'A', 'B', 'C', 'D', 'W', 'E'};
end
groups = upper(cellstr(groups));
if nargin < 2, study = struct(); end
T      = gf(study, 'T', 24);       
N = gf(study, 'N', 60);   
maxFE = gf(study, 'maxFE', 12000);
runs   = gf(study, 'runs', 1);     
seed0 = gf(study, 'seed0', 1);
solver = gf(study, 'solver', 'auto');
mpc0   = gf(study, 'mpc', RC3_case33());
folder = gf(study, 'folder', fullfile('results', ['experiments_' datestr(now, 'yyyymmdd_HHMMSS')])); %#ok<TNOW1,DATST>
if ~exist(folder, 'dir'), mkdir(folder); end
F5 = {'f1', 'f2', 'f3', 'f4', 'f5'};
objsets = gf(study, 'objsets', {F5, {'f1','f2'}, {'f1','f2','f4'}, {'f1','f2','f3','f4'}, F5});
Wsc = gf(study, 'weights', [0.6 0.1 0.1 0.1 0.1;  0.1 0.6 0.1 0.1 0.1;  0.1 0.1 0.6 0.1 0.1; ...
                            0.1 0.1 0.1 0.6 0.1;  0.1 0.1 0.1 0.1 0.6]);
Wname = {'economy', 'technical', 'environment and e-mobility', 'reliability', 'resilience and DR'};

%% ---- case matrix ---------------------------------------------------------------------
C = struct('id', {}, 'group', {}, 'devices', {}, 'objectives', {}, 'm', {}, 'weights', {}, ...
           'island', {}, 'T', {}, 'note', {}, 'alg', {}, 'prob', {});
net = {'f1', 'f2', 'f4'};
if ismember('A', groups)
    A = {'ONR', net;  'DG', net;  'SC', {'f1','f2'};  'TAP', {'f1','f2'};  'TX', {'f1','f2'}; ...
         'DSTATCOM', {'f1','f2'};  'UPFC', {'f1','f2'};  'EVCS', {'f1','f2','f3'}; ...
         'BSS', {'f1','f3','f5'};  'DRM', {'f1','f2','f5'}};
    for k = 1:size(A, 1)
        C(end+1) = mk(sprintf('A%02d', k), 'A', A(k, 1), A{k, 2}, [], 'off', T, ['only ' A{k, 1}]); %#ok<AGROW>
    end
end
if ismember('B', groups)
    Bc = {{'ONR','DG'}, net, 'ONR + DG'; ...
          {'ONR','DG','SC'}, net, 'ONR + DG + SC'; ...
          {'DG','SC','TAP','TX'}, {'f1','f2'}, 'voltage control: DG + SC + OLTC + TX'; ...
          {'DG','DSTATCOM','UPFC'}, {'f1','f2'}, 'FACTS: DG + D-STATCOM + UPFC'; ...
          {'EVCS','BSS'}, {'f1','f2','f3'}, 'e-mobility: EVCS + BSS'; ...
          {'EVCS','BSS','DRM'}, {'f1','f3','f5'}, 'e-mobility + DR'; ...
          {'ONR','DG','SC','TAP','TX','DSTATCOM','UPFC'}, net, 'all network devices'; ...
          {'ONR','DG','SC','TAP','TX','DSTATCOM','UPFC','EVCS','BSS','DRM'}, {'f1','f2','f3','f4'}, 'all devices'};
    for k = 1:size(Bc, 1)
        C(end+1) = mk(sprintf('B%02d', k), 'B', Bc{k, 1}, Bc{k, 2}, [], 'off', T, Bc{k, 3}); %#ok<AGROW>
    end
end
if ismember('C', groups)
    Ic = {{'DG','BSS'}, {'f1','f5'}, 'outage', 'outage 18-21 h, one island'; ...
          {'ONR','DG','BSS','DRM','ISL'}, {'f1','f2','f5'}, 'outage', 'outage, optimal microgrid formation'; ...
          {'DG','BSS','ISL'}, {'f1','f4','f5'}, 'optimal', 'intentional (optimal) islanding'; ...
          {'DG','SC','BSS','DRM','ISL'}, {'f1','f2','f5'}, 'islanded', 'standalone islanded network'};
    for k = 1:size(Ic, 1)
        C(end+1) = mk(sprintf('C%02d', k), 'C', Ic{k, 1}, Ic{k, 2}, [], Ic{k, 3}, T, Ic{k, 4}); %#ok<AGROW>
    end
end
allD = {'ONR','DG','SC','TAP','TX','DSTATCOM','UPFC','EVCS','BSS','DRM','ISL'};
if ismember('D', groups)
    for m = 1:numel(objsets)
        os = objsets{m};
        if m == 1
            C(end+1) = mk('D01', 'D', allD, os, ones(1, numel(os)) / numel(os), 'outage', T, ...
                          sprintf('m = 1: equal-weight sum of %s', strjoin(os, ', '))); %#ok<AGROW>
        else
            C(end+1) = mk(sprintf('D%02d', m), 'D', allD, os, [], 'outage', T, ...
                          sprintf('m = %d: %s', numel(os), strjoin(os, ', '))); %#ok<AGROW>
        end
    end
end
if ismember('W', groups)
    for k = 1:size(Wsc, 1)
        nm = sprintf('scenario %d', k);  
        if k <= numel(Wname)
           nm = Wname{k}; 
        end
        C(end+1) = mk(sprintf('W%02d', k), 'W', allD, F5, Wsc(k, :), 'outage', T, ...
                      sprintf('m = 1, %s weights %s', nm, mat2str(Wsc(k, :)))); %#ok<AGROW>
    end
end
if ismember('E', groups)
    dE = {'ONR','DG','SC','TAP','DSTATCOM'};
    C(end+1) = mk('E01', 'E', dE, {'f1','f2'}, [], 'off', 1,  'peak-hour snapshot');
    C(end+1) = mk('E02', 'E', dE, {'f1','f2'}, [], 'off', 24, '24-h profile');
end

if ismember('L', groups)
    if exist('platemo', 'file') ~= 2
        error('RC3_experiments: group L needs PlatEMO on the MATLAB path.');
    end
    algS = cellstr(gf(study, 'algs_single', {'GA', 'CSO', 'IMODE', 'SHADE', 'AESSPSO'}));
    algM = cellstr(gf(study, 'algs_multi', {'NSGAII', 'NSGAIII', 'LMOCSO', 'CMOCSO', ...
                   'IMTCMO_BS', 'DSSEA', 'LCMEA', 'ILCMO', 'PRCEA'}));
    Lp = {'L-S', F5, ones(1, 5) / 5, algS;  'L-M3', {'f1','f2','f4'}, [], algM;  'L-M5', F5, [], algM};
    k = 0;
    for p = 1:size(Lp, 1)
        for a = 1:numel(Lp{p, 4})
            k = k + 1;
            c = mk(sprintf('L%02d', k), 'L', allD, Lp{p, 2}, Lp{p, 3}, 'outage', T, ...
                   sprintf('%s %s', Lp{p, 1}, Lp{p, 4}{a}));
            c.alg = Lp{p, 4}{a};  c.prob = Lp{p, 1};
            C(end+1) = c; %#ok<AGROW>
        end
    end
end
%% Selection of best Algorithm and Create our method for the proposed problem
if ismember('L_M', groups)
   allD = {'ONR','DG','SC','TAP','TX','DSTATCOM','EVCS','BSS','ISL'};
   outage= {'off'; 'outage' ; 'optimal' ; 'islanded'};
   if exist('platemo', 'file') ~= 2
      error('RC3_experiments: group L needs PlatEMO on the MATLAB path.');
   end
   algS = cellstr(gf(study, 'algs_single', {'GA', 'CSO', 'IMODE', 'SHADE', 'AESSPSO', 'MASSRMADE'}));
   algM = cellstr(gf(study, 'algs_multi', {'NSGAII', 'LMOCSO', 'CMOCSO', ...
      'IMTCMO_BS', 'DSSEA', 'LCMEA', 'ILCMO', 'PRCEA','MASSRCMODRL', 'MASSRMAMOEAD'}));
   Lp = {'M_L-S', {'f1','f3', 'f5'}, ones(1, 3) / 3, algS;
              'M_L-M3', {'f1','f3', 'f5'}, [], algM;};
   k = 0;
   for p = 1:size(Lp, 1)
      for a = 1:numel(Lp{p, 4})
         for i=1:numel(outage)
            k = k + 1;
            c = mk(sprintf('M_L%02d', k), 'M_L', allD, Lp{p, 2}, Lp{p, 3}, outage{i}, T, ...
               sprintf('%s %s', Lp{p, 1}, Lp{p, 4}{a}));
            c.alg = Lp{p, 4}{a};  c.prob = Lp{p, 1};
            C(end+1) = c; %#ok<AGROW>
         end
      end
   end
end



%% ---- solve every case ---------------------------------------------------------------------
ind = {'TAC','INVCOST','OPCOST','LOSS','VD','VSI','ENS','SAIDI','SAIFI','PEAK','PAR','EMISSION', ...
       'EVCSUN','BSSUN','EVWAIT','DEGCOST','DRCOST','DISCOMF','LSHED','WLSHED','FDEV','SAI','SAD'};
nC = numel(C);  
nL = 5;
E.cases = C;  
E.R = cell(nC, 1);  
E.folder = folder;
fb = nan(nC, nL);  
fbase = fb;  
vals = nan(nC, numel(ind));  
base = vals;
extra = nan(nC, 7); 
rows = cell(nC, 1);  
devs = cell(nC, 1);  
stats = cell(nC, 1);
fronts = cell(nC, 1); 
sbest = cell(nC, 1);  
tRun = cell(nC, 1);  
okRun = cell(nC, 1);
conv = cell(nC, 1);
if exist('RC3_exp_Workspace.mat')==2
   load('RC3_exp_Workspace.mat')
end
for c = 1:nC
    opts = struct('mpc', mpc0, 'devices', {C(c).devices}, 'objectives', {C(c).objectives}, ...
                  'island_mode', C(c).island, 'N', N, 'maxFE', maxFE, 'solver', solver, ...
                  'validate', true, 'verbose', true);
    if isempty(C(c).weights)
        opts.m = numel(C(c).objectives);
    else
        opts.m = 1;
        opts.weights = C(c).weights;
    end
    if ~isempty(C(c).alg)
       opts.alg = C(c).alg;  
       opts.solver = 'platemo';
    end
    Fr = nan(runs, nL);  nf = zeros(runs, 1);
    fronts{c} = cell(runs, 1);  sbest{c} = nan(runs, 1);  tRun{c} = nan(runs, 1);
    okRun{c} = false(runs, 1);
    for r = 1:runs
        opts.seed = seed0 + r - 1;
        R = RC3_solve(C(c).T, opts);
        Fr(r, :) = R.Flib_best;  nf(r) = R.nfeas;
        fronts{c}{r} = R.Obj;  sbest{c}(r) = R.Obj(R.ibest, 1);
        tRun{c}(r) = R.time;  okRun{c}(r) = R.nfeas > 0;
        if r == 1
            R1 = R;
            H = R.hist;  fz = H.F(:, 1);  fz(H.CV > 0) = Inf;   % best feasible value so far
            conv{c} = cummin(fz);
        end
    end
    R = R1;
    sub  = fullfile(folder, C(c).id);
    Tb   = RC3_report(R, false);
    figs = RC3_plot_results(R, 'off');
    RC3_save_results(R, sub, Tb, figs);
    close(figs);
    E.R{c} = nohandle(R);
    [~, ii] = ismember(ind, R.P.allObj);
    fb(c, :) = R.Flib_best;  fbase(c, :) = R.Flib_base;
    vals(c, :) = R.out_best.Fall(ii);  base(c, :) = R.out_base.Fall(ii);
    extra(c, :) = [min(R.out_best.Vmin), R.out_best.VSImin, R.P.D, R.time, ...
                   double(R.out_best.feasible), size(R.Obj, 1), mean(nf)];
    rows{c} = sprintf('%s %s', C(c).id, C(c).note);
    dl = cellfun(@(a, b2, s) sprintf('%s@%s %g', a, s, b2), Tb.devices.Device, ...
                 num2cell(Tb.devices.Site), Tb.devices.SiteType, 'UniformOutput', false);
    devs{c} = strjoin(dl(:)', ', ');
    stats{c} = array2table([mean(Fr, 1); std(Fr, 0, 1); min(Fr, [], 1)], ...
                           'VariableNames', F5, 'RowNames', {'mean', 'std', 'best'});
    save('RC3_exp_Workspace')
end

%% ---- comparison tables -------------------------------------------------------------------------
rn = matlab.lang.makeUniqueStrings(rows);
vn = [F5, ind, {'Vmin_pu', 'VSImin', 'D', 'Time_s', 'Feasible', 'nPareto', 'FeasibleFinalPop'}];
E.table = array2table([fb, vals, extra], 'VariableNames', vn, 'RowNames', rn);
E.table.Optimised = arrayfun(@(q) objtext(q), C, 'UniformOutput', false)';
E.table.Devices = devs;
E.table.Island = {C.island}';
E.base = array2table([fbase, base], 'VariableNames', [F5, ind], 'RowNames', rn);
rel = fb ./ fbase;  rel(abs(fbase) < 1e-12) = NaN;
E.rel = array2table(rel, 'VariableNames', F5, 'RowNames', rn);
chg = 100 * ([fb, vals] - [fbase, base]) ./ max(abs([fbase, base]), 1e-12);
chg(abs([fbase, base]) < 1e-12) = NaN;
E.change = array2table(chg, 'VariableNames', [F5, ind], 'RowNames', rn);
E.stats = stats;
xl = fullfile(folder, 'comparison.xlsx');
if exist(xl, 'file'), delete(xl); end
writetable(E.table,  xl, 'Sheet', 'best',       'WriteRowNames', true);
writetable(E.base,   xl, 'Sheet', 'base',       'WriteRowNames', true);
writetable(E.rel,    xl, 'Sheet', 'f_relative', 'WriteRowNames', true);
writetable(E.change, xl, 'Sheet', 'change_pct', 'WriteRowNames', true);
lib = RC3_objlib();
writetable(struct2table(lib), xl, 'Sheet', 'objective_definitions');
writetable(E.table, fullfile(folder, 'comparison.csv'), 'WriteRowNames', true);
for c = 1:nC
    writetable(stats{c}, xl, 'Sheet', ['stats_' C(c).id], 'WriteRowNames', true);
end

%% ---- group L: algorithm comparison -------------------------------------------------------------
kL = find(strcmp({C.group}, 'L'));
if ~isempty(kL)
    [E.L, figL] = algcompare(C(kL), fronts(kL), sbest(kL), tRun(kL), okRun(kL), conv(kL));
    writetable(E.L, xl, 'Sheet', 'algorithms', 'WriteRowNames', false);
    writetable(E.L, fullfile(folder, 'algorithms.csv'));
    for k = 1:numel(figL)
        f = fullfile(folder, get(figL(k), 'Name'));
        print(figL(k), f, '-dpng', '-r300');  savefig(figL(k), [f '.fig']);
    end
    close(figL);
end

%% ---- comparison figures -------------------------------------------------------------------------
figs = gobjects(0);
for g = unique({C.group}, 'stable')
    k = find(strcmp({C.group}, g{1}));
    figs(end+1) = figure('Name', ['compare_' g{1}], 'Color', 'w', 'Visible', 'off'); %#ok<AGROW>
    Y = rel(k, :);  if numel(k) == 1, Y = [Y; nan(1, nL)]; end   %#ok<AGROW>
    bar(Y);  grid on;  xlim([0.5, numel(k) + 0.5]);
    set(gca, 'XTick', 1:numel(k), 'XTickLabel', {C(k).id});
    yline(1, 'k--');  ylabel('f_k(best) / f_k(base)');
    legend([F5, {'base'}], 'Location', 'eastoutside');
    title(sprintf('Group %s: library objectives relative to the base case (lower is better)', g{1}));
end
kw = find(strcmp({C.group}, 'W'));
if numel(kw) > 1
    figs(end+1) = figure('Name', 'compare_W_tradeoff', 'Color', 'w', 'Visible', 'off');
    plot(rel(kw, :)', '-o', 'LineWidth', 1.5);  grid on;
    set(gca, 'XTick', 1:nL, 'XTickLabel', F5);  ylabel('f_k(best) / f_k(base)');
    legend({C(kw).id}, 'Location', 'best');  title('Weight scenarios (m = 1): trade-off between f1..f5');
end
figs(end+1) = figure('Name', 'compare_heatmap', 'Color', 'w', 'Visible', 'off');
Z = chg;  Z(~isfinite(Z)) = 0;
imagesc(max(min(Z, 100), -100));  colormap(gca, redblue(64));  colorbar;
set(gca, 'XTick', 1:numel(vn) - 7, 'XTickLabel', [F5, ind], 'XTickLabelRotation', 45, ...
         'YTick', 1:nC, 'YTickLabel', {C.id}, 'FontSize', 7);
title('Change against the base case (%, blue = improvement)');
for k = 1:numel(figs)
    f = fullfile(folder, get(figs(k), 'Name'));
    print(figs(k), f, '-dpng', '-r300');  savefig(figs(k), [f '.fig']);
end
close(figs);
save(fullfile(folder, 'experiments.mat'), 'E', '-v7.3');
fprintf('\nExperimental study saved in %s\n', folder);
disp(E.table(:, [F5, {'Vmin_pu', 'Feasible', 'Time_s'}]));
end

% =========================================================================================
function c = mk(id, group, devices, objectives, weights, island, T, note)
c = struct('id', id, 'group', group, 'devices', {devices}, 'objectives', {objectives}, ...
           'm', [], 'weights', weights, 'island', island, 'T', T, 'note', note, ...
           'alg', '', 'prob', '');
end

function [Tb, figs] = algcompare(C, fronts, sbest, tRun, okRun, conv)
% Group L: normalised hypervolume (m > 1) or best objective value (m = 1) of
% every algorithm over the runs, Wilcoxon rank-sum test against the first
% algorithm of each problem (+ better, - worse, = no significant difference).
probs = unique({C.prob}, 'stable');
nC = numel(C);  score = cell(nC, 1);  kind = cell(nC, 1);
for p = 1:numel(probs)
    k = find(strcmp({C.prob}, probs{p}));
    if isempty(C(k(1)).weights)                     % m > 1: hypervolume
        Y = [];
        for c = k
            for r = 1:numel(fronts{c})
                if okRun{c}(r), Y = [Y; fronts{c}{r}]; end %#ok<AGROW>
            end
        end
        if isempty(Y), lo = 0; hi = 1; else, Y = ndfilter(Y); lo = min(Y, [], 1); hi = max(Y, [], 1); end
        sc = hi - lo;  sc(sc < 1e-12) = 1;
        for c = k
            score{c} = zeros(numel(fronts{c}), 1);
            for r = 1:numel(fronts{c})
                if okRun{c}(r)
                    score{c}(r) = hvmc((fronts{c}{r} - lo) ./ sc, 1.1 * ones(1, numel(lo)), 1e5);
                end
            end
            kind{c} = 'HV (max)';
        end
    else                                            % m = 1: best objective value
        for c = k
            score{c} = sbest{c};  score{c}(~okRun{c}) = Inf;  kind{c} = 'F (min)';
        end
    end
end
mu = cellfun(@(v) mean(v(isfinite(v))), score);  sd = cellfun(@(v) std(v(isfinite(v))), score);
bs = zeros(nC, 1);  test = repmat({''}, nC, 1);  pval = nan(nC, 1);
for c = 1:nC
    k = find(strcmp({C.prob}, C(c).prob));  ref = k(1);
    hv = strcmp(kind{c}, 'HV (max)');
    if hv, bs(c) = max(score{c}); else, bs(c) = min(score{c}); end
    if c ~= ref && exist('ranksum', 'file') == 2
        a = score{c};  b = score{ref};
        if ~hv, a(~isfinite(a)) = realmax; b(~isfinite(b)) = realmax; end
        pval(c) = ranksum(a, b);
        better = (hv && median(a) > median(b)) || (~hv && median(a) < median(b));
        if pval(c) >= 0.05, test{c} = '='; elseif better, test{c} = '+'; else, test{c} = '-'; end
    end
end
Tb = table({C.id}', {C.prob}', {C.alg}', kind, mu, sd, bs, pval, test, ...
           cellfun(@(v) sum(v), okRun), cellfun(@(v) numel(v), okRun), cellfun(@mean, tRun), ...
           'VariableNames', {'Case', 'Problem', 'Algorithm', 'Indicator', 'Mean', 'Std', 'Best', ...
           'p_ranksum', 'Test', 'FeasibleRuns', 'Runs', 'MeanTime_s'});
figs = gobjects(0);
for p = 1:numel(probs)
    k = find(strcmp({C.prob}, probs{p}));
    figs(end+1) = figure('Name', ['algorithms_' strrep(probs{p}, '-', '')], 'Color', 'w', ...
                         'Visible', 'off'); %#ok<AGROW>
    bar(mu(k));  hold on;  errorbar(1:numel(k), mu(k), sd(k), 'k.', 'LineWidth', 1);  grid on;
    set(gca, 'XTick', 1:numel(k), 'XTickLabel', {C(k).alg}, 'XTickLabelRotation', 30);
    ylabel(kind{k(1)});  title(sprintf('%s: mean and standard deviation over the runs', probs{p}));
    if ~isempty(C(k(1)).weights)                    % convergence of the single-objective runs
        figs(end+1) = figure('Name', ['convergence_' strrep(probs{p}, '-', '')], 'Color', 'w', ...
                             'Visible', 'off'); %#ok<AGROW>
        hold on;
        for c = k, plot(conv{c}, 'LineWidth', 1.3); end
        grid on;  xlabel('function evaluations');  ylabel('best feasible objective');
        legend({C(k).alg}, 'Location', 'northeast');  title([probs{p} ': convergence (run 1)']);
    end
end
end

function Y = ndfilter(Y)
Y = unique(Y, 'rows');  nd = true(size(Y, 1), 1);
for i = 1:size(Y, 1)
    nd(i) = ~any(all(Y <= Y(i, :), 2) & any(Y < Y(i, :), 2));
end
Y = Y(nd, :);
end

function h = hvmc(Y, ref, ns)
% Monte Carlo hypervolume of the normalised points Y (minimisation) w.r.t. ref
Y = Y(all(Y < ref, 2), :);
if isempty(Y), h = 0; return; end
st = rng;  rng(0, 'twister');                 % same samples for every algorithm of a problem
lo = zeros(1, numel(ref));  S = rand(ns, numel(ref)) .* ref;   % normalised values are >= 0
dom = false(ns, 1);
for i = 1:size(Y, 1), dom = dom | all(S >= Y(i, :), 2); end
rng(st);
h = prod(ref - lo) * mean(dom);
end

function s = objtext(q)
if isempty(q.weights)
    s = sprintf('min [%s]', strjoin(q.objectives, ', '));
else
    t = arrayfun(@(w, o) sprintf('%.2g*%s', w, o{1}), q.weights, q.objectives, 'UniformOutput', false);
    s = ['min ' strjoin(t, ' + ') ' (normalised)'];
end
end

function v = gf(s, name, default)
if isfield(s, name) && ~isempty(s.(name))
   v = s.(name);
else
   v = default;
end
end

function R = nohandle(R)
if isfield(R, 'P') && isfield(R.P, 'om'), R.P = rmfield(R.P, 'om'); end
end

function c = redblue(n)
% blue (negative = improvement) - white - red colormap
h = floor(n / 2);  u = linspace(0, 1, h)';
c = [[u u ones(h, 1)]; [ones(n - h, 1) flipud(linspace(0, 1, n - h)') flipud(linspace(0, 1, n - h)')]];
end
