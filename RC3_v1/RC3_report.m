function Tb = RC3_report(R, verbose)
%RC3_REPORT  Result tables of a solved case (all also written to Excel by
%RC3_save_results) and a printed summary.
%
%   Tb = RC3_report(R)            R from RC3_solve
%   Tb = RC3_report(R, false)     tables only, nothing printed
%
%   Tb.summary      all indicators: base case, best solution, change (%)
%   Tb.objectives   user objectives (opts.objectives): definition, base, best,
%                   change, weight and normalisation for m = 1
%   Tb.library      library objectives f1..f5 (RC3_objlib) of base and best
%   Tb.objterms     every term of every weighted sum: weight, reference value,
%                   raw, normalised and weighted values of base and best
%   Tb.devices      installed units: site (bus / branch), technology, rating, cost
%   Tb.hourly       slot results: grid import, losses, Vmin, VSI, PCC, islands,
%                   frequency, load shed, DR, EV, BSS, OLTC
%   Tb.dg Tb.sc Tb.tx Tb.dstatcom Tb.upfc Tb.ev Tb.bss   set-points per slot
%   Tb.switches     open branches per slot (ONR, PCC, microgrid cuts)
%   Tb.voltage      bus voltages per slot (best),  Tb.voltage_base
%   Tb.pareto       non-dominated set (selected objectives, fuzzy score, all objectives)
%   Tb.constraints  constraint values of base and best
%   Tb.variables    MP-Opt-Model variable sets (type, size, bounds, time dependence)
%   Tb.lincons      MP-Opt-Model linear constraint sets
%   Tb.validation   MATPOWER (runpf) validation of the best solution
if nargin < 2, verbose = true; end
P = R.P;  b = R.out_base;  p = R.out_best;  S = p.S;  T = P.T;
slot = arrayfun(@(t) sprintf('t%02d', t), 1:T, 'UniformOutput', false);

%% ---- summary of all indicators -----------------------------------------------------
names = [P.allObj, {'EGCOST','FUEL','SHEDCOST','REVENUE','VMIN','VSIMIN','NTRANS'}];
desc  = [cellfun(@(d, u) sprintf('%s (%s)', d, u), P.objDesc, P.objUnit, 'UniformOutput', false), ...
         {'grid energy cost (USD/day)', 'DG fuel cost (USD/day)', 'lost-load cost (USD/day)', ...
          'EV/BSS revenue (USD/day)', 'minimum voltage (p.u.)', 'minimum VSI (-)', 'PCC transitions (-)'}];
xv = @(o) [o.Fall, o.EGCOST, o.FUEL, o.SHEDCOST, o.REVENUE, min(o.Vmin), o.VSImin, o.nTrans];
Base = xv(b)';  Best = xv(p)';
Change_pct = 100 * (Best - Base) ./ abs(Base);
Change_pct(abs(Base) < 1e-12) = NaN;
Description = desc(:);
Ref = [P.objRef, nan(1, 7)]';
Tb.summary = table(Description, Base, Best, Change_pct, Ref, 'RowNames', names(:));

%% ---- user objective functions and library f1..f5 ---------------------------------------
Base = R.F_base(:);  Best = R.F_best(:);
Change_pct = 100 * (Best - Base) ./ abs(Base);  Change_pct(abs(Base) < 1e-12) = NaN;
Weight = P.weights(:);  Fref = P.Fref(:);
Definition = {P.objdef.expr}';  Label = P.objLabels(:);
Tb.objectives = table(Label, Definition, Base, Best, Change_pct, Weight, Fref, 'RowNames', P.objNames(:));
Base = R.Flib_base(:);  Best = R.Flib_best(:);
Change_pct = 100 * (Best - Base) ./ abs(Base);  Change_pct(abs(Base) < 1e-12) = NaN;
Definition = {P.objLib.expr}';  Label = {P.objLib.label}';
Tb.library = table(Label, Definition, Base, Best, Change_pct, 'RowNames', {P.objLib.name}');
Tb.objterms = term_table([b.terms(:); b.libTerms(:)], [p.terms(:); p.libTerms(:)]);

%% ---- installed devices -----------------------------------------------------------------
dv = cell(0, 7);
for k = 1:numel(S.dg.bus)
    dv(end+1, :) = {'DG', S.dg.unit(k), S.dg.bus(k), 'bus', P.dg.types{S.dg.type(k)}, S.dg.cap(k), 'MW'}; %#ok<AGROW>
end
for k = 1:numel(S.sc.bus)
    dv(end+1, :) = {'SC', S.sc.unit(k), S.sc.bus(k), 'bus', sprintf('%d x %.3g MVAr', P.sc.nsteps, P.sc.Qstep), S.sc.cap(k), 'MVAr'}; %#ok<AGROW>
end
if P.use.TAP
    dv(end+1, :) = {'OLTC', 1, P.pcc, 'branch', sprintf('+/-%d x %.4g', P.oltc.tap_max, P.oltc.step), NaN, '-'};
end
for k = 1:numel(S.tx.br)
    dv(end+1, :) = {'TX', S.tx.unit(k), S.tx.br(k), 'branch', sprintf('%d-%d', P.fb(S.tx.br(k)), P.tb(S.tx.br(k))), ...
                    max(abs(S.tx.a(k, :) - 1)) * 100, '% reg.'}; %#ok<AGROW>
end
for k = 1:numel(S.st.bus)
    dv(end+1, :) = {'DSTATCOM', S.st.unit(k), S.st.bus(k), 'bus', P.st.mode, max(abs(p.Qst(k, :))), 'MVAr'}; %#ok<AGROW>
end
for k = 1:numel(S.up.br)
    dv(end+1, :) = {'UPFC', S.up.unit(k), S.up.br(k), 'branch', sprintf('%d-%d', P.fb(S.up.br(k)), P.tb(S.up.br(k))), ...
                    p.Sup(k), 'MVA'}; %#ok<AGROW>
end
for k = 1:numel(S.ev.bus)
    dv(end+1, :) = {'EVCS', S.ev.unit(k), S.ev.bus(k), 'bus', sprintf('%d x %g kW', S.ev.n(k), S.ev.Pch(k)), ...
                    S.ev.n(k) * S.ev.Pch(k), 'kW'}; %#ok<AGROW>
end
for k = 1:numel(S.bss.bus)
    dv(end+1, :) = {'BSS', S.bss.unit(k), S.bss.bus(k), 'bus', sprintf('%d x %g kW, %d batt.', S.bss.n(k), ...
                    S.bss.Pch(k), S.bss.Nbat(k)), S.bss.Ecap(k), 'kWh'}; %#ok<AGROW>
end
if P.use.DRM
    dv(end+1, :) = {'DRM', 1, sum(S.dr.z), 'buses', 'enrolled buses', sum(S.dr.kappa .* P.Pd), 'MW enrolled'};
end
if isempty(dv), dv = {'none', 0, 0, '-', '-', 0, '-'}; end
Tb.devices = cell2table(dv, 'VariableNames', {'Device','Unit','Site','SiteType','Technology','Rating','RatingUnit'});
Tb.investment = struct2table(p.INV);

%% ---- slot results -----------------------------------------------------------------------
Hour = P.hours(:);  LoadMult = P.mult(:) .* ones(T, 1);  Price = P.price(:);
PCC = p.ug(:);  Psub_base = b.Psub(:);  Psub = p.Psub(:);  Loss_base = b.loss_kW(:);  Loss = p.loss_kW(:);
Vmin_base = b.Vmin(:);  Vmin = p.Vmin(:);  VSImin = p.VSI_t(:);  Islands = p.nIsl(:);
Freq_dev_Hz = p.df(:);  Shed_kW = 1000 * sum(p.shed, 1)';  SA = p.SA(:);
OLTC = S.oltc.tap(:);  DR_shift = S.dr.shift(:);  DR_curt = S.dr.il(:);
EV_kW = sum(p.ev.P, 1)';  BSS_kW = sum(p.bss.P, 1)';  DG_MW = sum(p.Pdg, 1)';
if isempty(EV_kW), EV_kW = zeros(T, 1); end
if isempty(BSS_kW), BSS_kW = zeros(T, 1); end
if isempty(DG_MW), DG_MW = zeros(T, 1); end
Tb.hourly = table(Hour, LoadMult, Price, PCC, Psub_base, Psub, Loss_base, Loss, Vmin_base, Vmin, ...
                  VSImin, Islands, Freq_dev_Hz, Shed_kW, SA, OLTC, DR_shift, DR_curt, EV_kW, BSS_kW, DG_MW, ...
                  'RowNames', slot(:));

%% ---- set-points per slot --------------------------------------------------------------------
Tb.dg = sched(p.Pdg, arrayfun(@(k) sprintf('DG%d bus %d %s (MW)', S.dg.unit(k), S.dg.bus(k), ...
              P.dg.types{S.dg.type(k)}), 1:numel(S.dg.bus), 'UniformOutput', false), slot);
Tb.sc = sched(S.sc.Q, arrayfun(@(k) sprintf('SC%d bus %d (MVAr)', S.sc.unit(k), S.sc.bus(k)), ...
              1:numel(S.sc.bus), 'UniformOutput', false), slot);
Tb.tx = sched([S.oltc.tap; S.tx.tap], [{'OLTC tap'}, arrayfun(@(k) sprintf('TX%d branch %d tap', ...
              S.tx.unit(k), S.tx.br(k)), 1:numel(S.tx.br), 'UniformOutput', false)], slot);
Tb.dstatcom = sched(p.Qst, arrayfun(@(k) sprintf('ST%d bus %d Q (MVAr)', S.st.unit(k), S.st.bus(k)), ...
              1:numel(S.st.bus), 'UniformOutput', false), slot);
nu = numel(S.up.br);  rows = cell(1, 4 * nu);  M = zeros(4 * nu, T);
for k = 1:nu
    M(4*k-3:4*k, :) = [S.up.r(k, :); S.up.th(k, :); S.up.Qsh(k, :); p.Sse(k, :)];
    rows(4*k-3:4*k) = {sprintf('UPFC%d br %d r (p.u.)', S.up.unit(k), S.up.br(k)), ...
                       sprintf('UPFC%d br %d theta (rad)', S.up.unit(k), S.up.br(k)), ...
                       sprintf('UPFC%d br %d Qsh (MVAr)', S.up.unit(k), S.up.br(k)), ...
                       sprintf('UPFC%d br %d Sse (MVA)', S.up.unit(k), S.up.br(k))};
end
Tb.upfc = sched(M, rows, slot);
ne = numel(S.ev.bus);  rows = cell(1, 3 * ne);  M = zeros(3 * ne, T);
for k = 1:ne
    M(3*k-2:3*k, :) = [p.ev.P(k, :); S.ev.capf(k, :); p.ev.Wq(k, :)];
    rows(3*k-2:3*k) = {sprintf('EVCS%d bus %d P (kW)', S.ev.unit(k), S.ev.bus(k)), ...
                       sprintf('EVCS%d bus %d cap (p.u.)', S.ev.unit(k), S.ev.bus(k)), ...
                       sprintf('EVCS%d bus %d wait (min)', S.ev.unit(k), S.ev.bus(k))};
end
Tb.ev = sched(M, rows, slot);
nbs = numel(S.bss.bus);  rows = cell(1, 3 * nbs);  M = zeros(3 * nbs, T);
for k = 1:nbs
    M(3*k-2:3*k, :) = [p.bss.P(k, :); 100 * p.bss.soc(k, 2:end); p.bss.swap(k, :)];
    rows(3*k-2:3*k) = {sprintf('BSS%d bus %d P (kW)', S.bss.unit(k), S.bss.bus(k)), ...
                       sprintf('BSS%d bus %d SOC end (%%)', S.bss.unit(k), S.bss.bus(k)), ...
                       sprintf('BSS%d bus %d swaps', S.bss.unit(k), S.bss.bus(k))};
end
Tb.bss = sched(M, rows, slot);
Open = cellfun(@mat2str, p.open(:), 'UniformOutput', false);
Cuts = arrayfun(@(t) mat2str(nonzeros(S.cut(:, t))'), (1:T)', 'UniformOutput', false);
Roots = arrayfun(@(t) mat2str([p.roots{t}.bus]), (1:T)', 'UniformOutput', false);
Tb.switches = table(Open, Cuts, Roots, 'RowNames', slot(:));
busN = arrayfun(@(i) sprintf('bus%02d', i), 1:P.nb, 'UniformOutput', false);
Tb.voltage      = sched(p.Vm, busN, slot);
Tb.voltage_base = sched(b.Vm, busN, slot);

%% ---- Pareto set, constraints, model -------------------------------------------------------
if P.m > 1
    Fz = R.score;  if isempty(Fz), Fz = nan(size(R.Obj, 1), 1); end
    Tp = array2table(R.Obj, 'VariableNames', P.objNames);
    Tp.FuzzyScore = Fz;  Tp.Best = (1:size(R.Obj, 1))' == R.ibest;
else
    Tp = table(R.Obj, 'VariableNames', {'WeightedObjective'});
    Tp = [Tp, array2table(R.Fuser_pareto, 'VariableNames', P.objNames)];
    Tp.Best = (1:size(R.Obj, 1))' == R.ibest;
end
Tlib = array2table(R.Flib_pareto, 'VariableNames', strcat('lib_', {P.objLib.name}));
Tall = array2table(R.Fall_pareto, 'VariableNames', strcat('all_', P.allObj));
Tb.pareto = [Tp, Tlib, Tall];
Tb.constraints = table(b.G(:), p.G(:), 'VariableNames', {'Base', 'Best'}, 'RowNames', P.conNl(:));
blk = P.blocks;
tv  = arrayfun(@(q) q.dims(2) == T && T > 1, blk)';
vtc = 'CI-B';
vt  = arrayfun(@(q) vtc(q.enc(1)), blk)';
Tb.variables = table({blk.name}', cellstr(vt), arrayfun(@(q) numel(q.idx), blk)', ...
                     arrayfun(@(q) min(q.lo), blk)', arrayfun(@(q) max(q.up), blk)', tv, {blk.desc}', ...
                     'VariableNames', {'Name','Type','Size','Lower','Upper','TimeVarying','Description'});
if isempty(P.linsets)
    Tb.lincons = table({'none'}, 0, {''}, 'VariableNames', {'Name','Rows','Description'});
else
    Tb.lincons = table({P.linsets.name}', arrayfun(@(c) size(c.A, 1), P.linsets)', {P.linsets.desc}', ...
                       arrayfun(@(c) c.repaired, P.linsets)', 'VariableNames', {'Name','Rows','Description','Repaired'});
end
if ~isempty(R.val), Tb.validation = R.val; end
Tb.check = table({R.solver}, R.check.maxdObj, R.check.maxdCon, R.nfeas, ...
                 'VariableNames', {'Solver', 'MaxRelDiffObj', 'MaxRelDiffCon', 'FeasibleInFinalPop'});

%% ---- print -------------------------------------------------------------------------------
if verbose
    fprintf('\n=== %s | objectives: %s (m = %d) | T = %d | island mode: %s | MP-Opt-Model: %s ===\n', ...
            strjoin(P.devices, ' + '), strjoin(P.objNames, ', '), P.m, T, P.isl.mode, P.om_api);
    disp(Tb.objectives(:, 2:end));
    disp(Tb.library(:, 2:end));
    disp(Tb.summary(:, 2:4));
    disp(Tb.devices);
    if ~p.feasible, warning('Best solution violates constraints: G = %s', mat2str(p.G, 3)); end
    fprintf('Open switches (slot 1): %s\n', mat2str(p.open{1}));
    if ~isempty(R.val)
        fprintf('MATPOWER check: max |dV| = %.2e p.u., max |dLoss| = %.3g kW\n', ...
                max(R.val.MaxDeltaV_pu), max(abs(R.val.Loss_FBS_kW - R.val.Loss_MP_kW)));
    end
end
end

function Tt = term_table(Tbase, Tbest)
% one row per term of every weighted-sum objective (base and best solution)
c = cell(0, 10);
for k = 1:numel(Tbase)
    a = Tbase(k);  z = Tbest(k);
    for j = 1:numel(a.terms)
        c(end+1, :) = {a.name, a.terms{j}, a.w(j), a.ref(j), a.raw(j), z.raw(j), ...
                       a.z(j), z.z(j), a.contrib(j), z.contrib(j)}; %#ok<AGROW>
    end
end
vn = {'Objective','Term','Weight','Reference','RawBase','RawBest','NormBase','NormBest', ...
      'ContribBase','ContribBest'};
if isempty(c), c = {'none', '-', 0, 1, 0, 0, 0, 0, 0, 0}; end
Tt = cell2table(c, 'VariableNames', vn);
end

function Tt = sched(M, rows, slot)
% schedule matrix -> table with one column per slot
if isempty(M) || isempty(rows)
    Tt = array2table(nan(1, numel(slot)), 'VariableNames', slot, 'RowNames', {'none'});
else
    Tt = array2table(M, 'VariableNames', slot, 'RowNames', rows(:));
end
end
