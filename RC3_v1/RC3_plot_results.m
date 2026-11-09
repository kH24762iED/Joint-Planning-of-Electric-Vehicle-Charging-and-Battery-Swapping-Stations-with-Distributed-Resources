function figs = RC3_plot_results(R, show)
%RC3_PLOT_RESULTS  Figures of a solved case for comparison and reporting.
%
%   figs = RC3_plot_results(R)          R from RC3_solve
%   figs = RC3_plot_results(R, 'off')   figures created invisible (batch runs)
%
%   pareto       Pareto front (2-D, 3-D) or parallel coordinates (m >= 4)
%   convergence  best feasible objective(s) vs function evaluations
%   base_vs_best key indicators of the base case and the best solution
%   objectives   library objectives f1..f5 (best / base) and the weighted
%                contributions of every term of f2..f5, base vs best
%   voltage      voltage profile in the worst slot (+ bus-slot heat maps)
%   network      single-line graph with sites of all installed devices
%   hourly       grid import, losses, minimum voltage per slot        (T > 1)
%   devices      DG, SC, D-STATCOM, UPFC set-points and tap schedules (T > 1)
%   switching    open branches per slot                               (T > 1)
%   ev_bss       EV charging and BSS operation                        (T > 1)
%   drm          load curve before/after demand response             (T > 1)
%   islanding    PCC status, islands, frequency deviation, load shed
if nargin < 2, show = 'on'; end
P = R.P;  b = R.out_base;  p = R.out_best;  S = p.S;  T = P.T;
cB = [0.55 0.55 0.55];  cP = [0.15 0.45 0.80];  cS = [0.85 0.33 0.10];  cG = [0.20 0.60 0.30];
figs = gobjects(0);
nf = @(name) figure('Name', name, 'Color', 'w', 'Visible', show);
ttl = sprintf('%s', strjoin(P.devices, ' + '));

%% (1) Pareto front
if P.m > 1 && size(R.Obj, 1) >= 1
    figs(end+1) = nf('pareto');
    if P.M == 2
        plot(R.Obj(:, 1), R.Obj(:, 2), '-o', 'Color', cP, 'MarkerFaceColor', cP, 'LineWidth', 1.2); hold on;
        plot(R.F_base(1), R.F_base(2), 's', 'MarkerSize', 10, 'MarkerFaceColor', cB, 'MarkerEdgeColor', 'k');
        plot(R.Obj(R.ibest, 1), R.Obj(R.ibest, 2), 'p', 'MarkerSize', 16, 'MarkerFaceColor', cS, 'MarkerEdgeColor', 'k');
        xlabel(P.objNames{1}); ylabel(P.objNames{2});
        legend('Pareto front', 'Base case', 'Best compromise', 'Location', 'best');
    elseif P.M == 3
        scatter3(R.Obj(:, 1), R.Obj(:, 2), R.Obj(:, 3), 30, cP, 'filled'); hold on;
        scatter3(R.F_base(1), R.F_base(2), R.F_base(3), 90, cB, 's', 'filled');
        scatter3(R.Obj(R.ibest, 1), R.Obj(R.ibest, 2), R.Obj(R.ibest, 3), 160, cS, 'p', 'filled');
        xlabel(P.objNames{1}); ylabel(P.objNames{2}); zlabel(P.objNames{3}); view(45, 25);
        legend('Pareto front', 'Base case', 'Best compromise', 'Location', 'best');
    else
        F  = [R.Obj; R.F_base];  lo = min(F, [], 1);  rg = max(F, [], 1) - lo;  rg(rg == 0) = 1;
        Fn = (F - lo) ./ rg;  np = size(R.Obj, 1);
        plot(1:P.M, Fn(1:np, :)', '-', 'Color', [cP 0.35]); hold on;
        plot(1:P.M, Fn(R.ibest, :), '-p', 'Color', cS, 'LineWidth', 2.5, 'MarkerFaceColor', cS);
        plot(1:P.M, Fn(end, :), '--s', 'Color', 'k', 'LineWidth', 1.5);
        set(gca, 'XTick', 1:P.M, 'XTickLabel', P.objNames);  ylabel('normalised value');
    end
    grid on;  title(sprintf('Pareto front (%s)', ttl));
end

%% (2) convergence
H = R.hist;
if ~isempty(H) && isfield(H, 'F') && ~isempty(H.F)
    figs(end+1) = nf('convergence');
    nO = size(H.F, 2);  feas = H.CV <= 1e-9;
    for k = 1:nO
        subplot(nO, 1, k);
        f = H.F(:, k);  f(~feas) = Inf;  bsf = cummin(f);  bsf(~isfinite(bsf)) = NaN;
        plot(1:numel(bsf), bsf, '-', 'Color', cP, 'LineWidth', 1.5); grid on;
        if P.m == 1, ylabel('weighted objective'); else, ylabel(P.objNames{k}); end
        if k == 1, title(sprintf('Best feasible value vs evaluations (%s)', R.alg)); end
    end
    xlabel('function evaluations');
end

%% (3) base vs best
figs(end+1) = nf('base_vs_best');
key = {'TAC','INVCOST','LOSS','VD','VSI','ENS','SAIDI','PEAK','EMISSION','LSHED','EVUNMET','OPCOST'};
lab = {'TAC ($/yr)','Investment ($)','Loss (kW)','Voltage deviation','1/min VSI','ENS (kWh/yr)', ...
       'SAIDI (h/cust.yr)','Peak import (MW)','CO_2 (kg/day)','Load shed (kWh/day)','EV unserved (kWh/day)', ...
       'Operating cost ($/day)'};
for k = 1:numel(key)
    i = find(strcmp(P.allObj, key{k}));
    v = [b.Fall(i), p.Fall(i)];
    subplot(3, 4, k);
    h = bar(v, 0.6, 'FaceColor', 'flat');  h.CData = [cB; cP];
    set(gca, 'XTickLabel', {'Base', 'Best'});  title(lab{k});  grid on;
    text(1:2, v, compose('%.4g', v), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 7);
end
sgt(sprintf('Base case vs best solution (%s)', ttl));

%% (3b) library objectives f1..f5 and their weighted terms
if ~isempty(p.Flib)
    figs(end+1) = nf('objectives');
    nL = numel(P.objLib);
    subplot(1, 2, 1);
    r = p.Flib ./ max(abs(b.Flib), 1e-12);  r(abs(b.Flib) < 1e-12) = NaN;
    h = bar(1:nL, [ones(1, nL); r]', 0.8);  h(1).FaceColor = cB;  h(2).FaceColor = cP;
    set(gca, 'XTick', 1:nL, 'XTickLabel', {P.objLib.name});  grid on;
    ylabel('value / base-case value');  legend('Base', 'Best', 'Location', 'best');
    title('Library objectives (relative to the base case)');
    subplot(1, 2, 2);
    lin = find(arrayfun(@(d) strcmp(d.type, 'linear'), P.objLib));
    C = [];  lab = {};  names = {};
    for k = lin(:)'
        tt = b.libTerms(k).terms;
        names = [names, tt(~ismember(tt, names))]; %#ok<AGROW>
    end
    for k = lin(:)'
        tb = b.libTerms(k);  tp = p.libTerms(k);
        cb = zeros(1, numel(names));  cp = cb;
        [~, j] = ismember(tb.terms, names);  cb(j) = tb.contrib;  cp(j) = tp.contrib;
        C = [C; cb; cp]; %#ok<AGROW>
        lab = [lab, {[P.objLib(k).name ' base'], [P.objLib(k).name ' best']}]; %#ok<AGROW>
    end
    if ~isempty(C)
        bar(C, 'stacked');  set(gca, 'XTick', 1:numel(lab), 'XTickLabel', lab, 'XTickLabelRotation', 45);
        legend(names, 'Location', 'eastoutside', 'FontSize', 7);  grid on;
        ylabel('w_j F_j / F_j^{ref}');  title('Weighted terms of the composite objectives');
    end
end

%% (4) voltage profile and heat maps
[~, tp] = min(min([b.Vmin; p.Vmin], [], 1));
figs(end+1) = nf('voltage');
if T > 1, subplot(1, 3, 1); end
plot(1:P.nb, b.Vm(:, tp), '-s', 'Color', cB, 'LineWidth', 1.3); hold on;
plot(1:P.nb, p.Vm(:, tp), '-o', 'Color', cP, 'LineWidth', 1.3);
yline(P.Vmin, 'r--');  yline(P.Vmax, 'r--');
xlabel('Bus'); ylabel('Voltage (p.u.)'); grid on; xlim([1 P.nb]);
legend('Base', 'Best', 'Location', 'southwest');  title(sprintf('Slot %d (lowest voltage)', tp));
if T > 1
    subplot(1, 3, 2); imagesc(1:T, 1:P.nb, b.Vm); axis xy; colorbar; caxis([0.88 1.06]);
    xlabel('Slot'); ylabel('Bus'); title('Base case |V|');
    subplot(1, 3, 3); imagesc(1:T, 1:P.nb, p.Vm); axis xy; colorbar; caxis([0.88 1.06]);
    xlabel('Slot'); ylabel('Bus'); title('Best solution |V|');
end

%% (5) network with device sites (slot of peak load)
try
    [~, tk] = max(P.mult);
    figs(end+1) = nf('network');
    G  = graph(P.fb, P.tb);
    hG = plot(G, 'Layout', 'layered', 'Sources', 1, 'NodeColor', [0.3 0.3 0.3], 'EdgeColor', [0.5 0.5 0.5], ...
              'LineWidth', 1.2, 'MarkerSize', 4);
    eo = find(~p.X(:, tk));  ei = findedge(G, P.fb(eo), P.tb(eo));
    highlight(hG, 'Edges', ei, 'LineStyle', '--', 'EdgeColor', [0.8 0.2 0.2]);
    hold on;  xd = hG.XData;  yd = hG.YData;
    mark = {S.dg.bus, 'DG', [0.95 0.65 0.1];  S.sc.bus, 'SC', [0.2 0.6 0.9];  S.st.bus, 'ST', [0.6 0.3 0.8]; ...
            S.ev.bus, 'EV', [0.1 0.7 0.4];  S.bss.bus, 'BSS', [0.85 0.2 0.5]};
    for k = 1:size(mark, 1)
        bb = mark{k, 1};
        if isempty(bb), continue; end
        scatter(xd(bb), yd(bb), 90, mark{k, 3}, 'filled', 'MarkerEdgeColor', 'k');
        text(xd(bb) + 0.08, yd(bb), mark{k, 2}, 'FontSize', 8, 'FontWeight', 'bold');
    end
    brs = [S.tx.br; S.up.br];  bl = [repmat({'TX'}, numel(S.tx.br), 1); repmat({'UPFC'}, numel(S.up.br), 1)];
    for k = 1:numel(brs)
        xm = mean(xd([P.fb(brs(k)) P.tb(brs(k))]));  ym = mean(yd([P.fb(brs(k)) P.tb(brs(k))]));
        plot(xm, ym, 'd', 'MarkerSize', 10, 'MarkerFaceColor', [0.95 0.9 0.2], 'MarkerEdgeColor', 'k');
        text(xm + 0.08, ym, bl{k}, 'FontSize', 8, 'FontWeight', 'bold');
    end
    axis off;  title(sprintf('Network and device sites (open branches dashed, slot %d)', tk));
catch err
    warning('RC3_plot_results:network', 'network plot skipped: %s', err.message);
end

if T > 1
    hr = 1:T;
    %% (6) slot results
    figs(end+1) = nf('hourly');
    subplot(3, 1, 1);
    plot(hr, b.Psub, '-s', 'Color', cB, 'LineWidth', 1.3); hold on;
    plot(hr, p.Psub, '-o', 'Color', cP, 'LineWidth', 1.3);
    plot(hr, sum(p.Pl, 1), ':', 'Color', cS, 'LineWidth', 1.3);
    ylabel('MW'); grid on; legend('Import base', 'Import best', 'Load (after DR)', 'Location', 'northwest');
    title('Grid import');
    subplot(3, 1, 2); bar(hr, [b.loss_kW(:) p.loss_kW(:)]); grid on; ylabel('kW');
    legend('Base', 'Best'); title('Power loss');
    subplot(3, 1, 3);
    plot(hr, b.Vmin, '-s', 'Color', cB, 'LineWidth', 1.3); hold on;
    plot(hr, p.Vmin, '-o', 'Color', cP, 'LineWidth', 1.3); yline(P.Vmin, 'r--');
    grid on; xlabel('Slot'); ylabel('p.u.'); legend('Base', 'Best'); title('Minimum voltage');

    %% (7) device set-points
    rows = {};
    if ~isempty(S.dg.bus), rows{end+1} = 'dg'; end
    if ~isempty(S.sc.bus), rows{end+1} = 'sc'; end
    if ~isempty(S.st.bus), rows{end+1} = 'st'; end
    if P.use.TAP || ~isempty(S.tx.br), rows{end+1} = 'tap'; end
    if ~isempty(S.up.br), rows{end+1} = 'up'; end
    if ~isempty(rows)
        figs(end+1) = nf('devices');
        for k = 1:numel(rows)
            subplot(numel(rows), 1, k);
            switch rows{k}
                case 'dg'
                    area(hr, p.Pdg'); grid on; ylabel('MW'); title('DG output');
                    legend(arrayfun(@(i) sprintf('DG%d %s bus %d', S.dg.unit(i), P.dg.types{S.dg.type(i)}, ...
                           S.dg.bus(i)), 1:numel(S.dg.bus), 'UniformOutput', false), 'Location', 'northwest');
                case 'sc'
                    stairs(hr, S.sc.Q', 'LineWidth', 1.3); grid on; ylabel('MVAr'); title('Capacitor banks');
                case 'st'
                    stairs(hr, p.Qst', 'LineWidth', 1.3); grid on; ylabel('MVAr'); title('D-STATCOM reactive power');
                case 'tap'
                    stairs(hr, [S.oltc.tap; S.tx.tap]', 'LineWidth', 1.3); grid on; ylabel('tap');
                    legend([{'OLTC'}, arrayfun(@(i) sprintf('TX br %d', S.tx.br(i)), 1:numel(S.tx.br), ...
                           'UniformOutput', false)], 'Location', 'best'); title('Tap positions');
                case 'up'
                    yyaxis left;  plot(hr, S.up.r', '-o', 'LineWidth', 1.2); ylabel('r (p.u.)');
                    yyaxis right; plot(hr, p.Sse', '-s', 'LineWidth', 1.2); ylabel('S_{se} (MVA)');
                    grid on; title('UPFC series voltage and power');
            end
        end
        xlabel('Slot');
    end

    %% (8) switching
    if P.use.ONR || P.use.ISL
        figs(end+1) = nf('switching');
        imagesc(hr, 1:P.nl, ~p.X); colormap(gca, [1 1 1; cS]); axis xy;
        xlabel('Slot'); ylabel('Branch'); title('Open branches (coloured): ONR, PCC and microgrid cuts');
    end

    %% (9) EV and BSS
    if P.use.EVCS || P.use.BSS
        figs(end+1) = nf('ev_bss');
        subplot(3, 1, 1);
        plot(hr, b.Psub, '-s', 'Color', cB, 'LineWidth', 1.3); hold on;
        plot(hr, p.Psub, '-o', 'Color', cP, 'LineWidth', 1.3);
        grid on; ylabel('MW'); legend('Base', 'Best'); title('Grid import');
        subplot(3, 1, 2);
        evp = sum(p.ev.P, 1);  bch = sum(p.bss.Pch, 1);  bdi = sum(p.bss.Pdis, 1);
        if isempty(evp), evp = zeros(1, T); end
        if isempty(bch), bch = zeros(1, T); bdi = zeros(1, T); end
        bar(hr, [evp; bch]', 'stacked'); hold on; bar(hr, -bdi, 'FaceColor', cS);
        grid on; ylabel('kW'); legend('EV charging', 'BSS charging', 'BSS discharging', 'Location', 'northwest');
        title('EV and BSS power (battery-to-grid negative)');
        subplot(3, 1, 3);
        if ~isempty(S.bss.n), plot(0:T, 100 * p.bss.soc', 'LineWidth', 1.3); end
        grid on; xlabel('Slot'); ylabel('SOC (%)'); xlim([0 T]); title('BSS battery pool state of charge');
    end

    %% (10) demand response
    if P.use.DRM
        figs(end+1) = nf('drm');
        subplot(2, 1, 1);
        plot(hr, P.Ptot * P.mult, '-s', 'Color', cB, 'LineWidth', 1.3); hold on;
        plot(hr, sum(p.Pl, 1), '-o', 'Color', cP, 'LineWidth', 1.3);
        grid on; ylabel('MW'); legend('Load without DR', 'Load with DR'); title('Demand response');
        subplot(2, 1, 2);
        bar(hr, [S.dr.shift; -S.dr.il]', 'grouped'); grid on; xlabel('Slot');
        ylabel('share of enrolled load'); legend('shift \rho(t)', 'curtailment -\chi(t)');
    end
end

%% (11) islanding
if P.isl.enabled
    figs(end+1) = nf('islanding');
    hr = 1:T;
    subplot(4, 1, 1); stairs(hr, p.ug, 'LineWidth', 1.5, 'Color', 'k'); hold on;
    stairs(hr, p.nIsl, 'LineWidth', 1.5, 'Color', cP); ylim([-0.2 max(2, max(p.nIsl) + 0.5)]);
    grid on; legend('PCC closed', 'energized islands'); title('Operating mode');
    subplot(4, 1, 2); bar(hr, p.df, 'FaceColor', cP); yline(P.isl.f_lim, 'r--');
    grid on; ylabel('Hz'); title('Frequency deviation of the islands');
    subplot(4, 1, 3); bar(hr, 1000 * [sum(b.shed, 1); sum(p.shed, 1)]'); grid on; ylabel('kW');
    legend('Base', 'Best'); title('Load shed');
    subplot(4, 1, 4); bar(hr, p.SA, 'FaceColor', cG); yline(1, 'r--'); grid on;
    xlabel('Slot'); ylabel('-'); title('Supply adequacy (GF + renewable / demand)');
end
end

function sgt(s)
try, sgtitle(s); catch, end %#ok<NOCOMMA>
end
