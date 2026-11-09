function folder = RC3_save_results(R, folder, Tb, figs)
%RC3_SAVE_RESULTS  Save a solved case for study and reporting:
%   <folder>/results.mat      R (without the MP-Opt-Model object and function
%                             handles) and the result tables Tb
%   <folder>/tables.xlsx      one sheet per table of RC3_report
%   <folder>/csv/*.csv        the same tables as CSV files
%   <folder>/fig/*.png|.fig   all figures of RC3_plot_results (300 dpi)
%   <folder>/summary.txt      settings, objectives, devices, run time
%
%   folder = RC3_save_results(R)                     results/<case>_<time>
%   folder = RC3_save_results(R, folder, Tb, figs)
if nargin < 2 || isempty(folder)
    tag = regexprep(strjoin(R.P.devices, '_'), '[^A-Za-z0-9_]', '');
    folder = fullfile('results', sprintf('%s_T%d_m%d_%s', tag, R.P.T, R.P.m, datestr(now, 'yyyymmdd_HHMMSS'))); %#ok<TNOW1,DATST>
end
if nargin < 3 || isempty(Tb), Tb = RC3_report(R, false); end
if nargin < 4, figs = []; end
if ~exist(folder, 'dir'), mkdir(folder); end
if ~exist(fullfile(folder, 'csv'), 'dir'), mkdir(fullfile(folder, 'csv')); end
if ~exist(fullfile(folder, 'fig'), 'dir'), mkdir(fullfile(folder, 'fig')); end

%% MAT file (no handles)
Rs = R;  Rs.P = rmfield(Rs.P, intersect(fieldnames(Rs.P), {'om'}));
save(fullfile(folder, 'results.mat'), 'Rs', 'Tb', '-v7.3');

%% tables
xl = fullfile(folder, 'tables.xlsx');
if exist(xl, 'file'), delete(xl); end
fn = fieldnames(Tb);
for k = 1:numel(fn)
    Tk = Tb.(fn{k});
    if ~istable(Tk), continue; end
    sh = fn{k};  if numel(sh) > 31, sh = sh(1:31); end
    hasRows = ~isempty(Tk.Properties.RowNames);
    try
        writetable(Tk, xl, 'Sheet', sh, 'WriteRowNames', hasRows);
    catch err
        warning('RC3_save_results:xlsx', 'sheet %s not written: %s', sh, err.message);
    end
    writetable(Tk, fullfile(folder, 'csv', [fn{k} '.csv']), 'WriteRowNames', hasRows);
end

%% figures
for k = 1:numel(figs)
    if ~isgraphics(figs(k)), continue; end
    nm = get(figs(k), 'Name');  if isempty(nm), nm = sprintf('figure%d', k); end
    f  = fullfile(folder, 'fig', nm);
    try
        set(figs(k), 'PaperPositionMode', 'auto');
        print(figs(k), f, '-dpng', '-r300');
        savefig(figs(k), [f '.fig']);
    catch err
        warning('RC3_save_results:fig', 'figure %s not saved: %s', nm, err.message);
    end
end

%% text summary
fid = fopen(fullfile(folder, 'summary.txt'), 'w');
P = R.P;
fprintf(fid, 'Devices      : %s\n', strjoin(P.devices, ' + '));
fprintf(fid, 'Objectives   : %s (m = %d)\n', strjoin(P.objNames, ', '), P.m);
for k = 1:P.M
    fprintf(fid, '   %-8s = %s   [%s]\n', P.objNames{k}, P.objdef(k).expr, P.objLabels{k});
end
fprintf(fid, 'Time slots   : T = %d (dt = %.3g h)\n', P.T, P.dt);
fprintf(fid, 'Island mode  : %s\n', P.isl.mode);
fprintf(fid, 'Variables    : D = %d   constraints: %d (nonlinear %d)\n', P.D, P.nCon, P.nConNl);
fprintf(fid, 'MP-Opt-Model : %s\n', P.om_api);
fprintf(fid, 'Solver       : %s / %s, N = %d, maxFE = %d, time %.1f s\n', R.solver, R.alg, R.N, R.maxFE, R.time);
fprintf(fid, '\n%-10s %16s %16s %10s\n', 'Objective', 'Base', 'Best', 'Change %');
for k = 1:P.M
    fprintf(fid, '%-10s %16.6g %16.6g %10.2f\n', P.objNames{k}, R.F_base(k), R.F_best(k), ...
            100 * (R.F_best(k) - R.F_base(k)) / max(abs(R.F_base(k)), 1e-12));
end
fprintf(fid, '\nLibrary objectives (RC3_objlib):\n');
for k = 1:numel(P.objLib)
    fprintf(fid, '%-6s %-75s base %12.6g  best %12.6g\n', P.objLib(k).name, P.objLib(k).expr, ...
            R.Flib_base(k), R.Flib_best(k));
end
fprintf(fid, '\nSolver check: max rel. difference of re-evaluated objectives %.3g, constraints %.3g; %d feasible in the final population\n', ...
        R.check.maxdObj, R.check.maxdCon, R.nfeas);
fprintf(fid, '\nFeasible best solution: %d   constraint values: %s\n', R.out_best.feasible, mat2str(R.out_best.G, 4));
fprintf(fid, '\nInstalled devices:\n');
D = Tb.devices;
for k = 1:height(D)
    fprintf(fid, '  %-9s unit %d at %s %d  %-28s %10.4g %s\n', D.Device{k}, D.Unit(k), D.SiteType{k}, ...
            D.Site(k), D.Technology{k}, D.Rating(k), D.RatingUnit{k});
end
fclose(fid);
fprintf('Results saved in %s\n', folder);
end
