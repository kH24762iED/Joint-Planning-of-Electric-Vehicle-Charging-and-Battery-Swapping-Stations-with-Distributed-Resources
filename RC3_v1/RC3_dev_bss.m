function [b, out2] = RC3_dev_bss(action, b, varargin)
%RC3_DEV_BSS  Battery swapping stations: the battery pool of each station is
%one equivalent BESS whose energy is updated slot by slot.
%
%   b        = RC3_dev_bss('init', [], S, P, Lnet)    state at the start of the day
%   b        = RC3_dev_bss('swap', b, t)              1) serve the swap requests
%   [b, Pkw] = RC3_dev_bss('grid', b, t)              grid-connected request (kW,
%                                                    + charge / - battery-to-grid)
%   [b, cap] = RC3_dev_bss('island', b, t)            GF capability: cap.Pdis, cap.Pch (kW)
%   b        = RC3_dev_bss('commit', b, t, Pkw)       apply the slot power (kW)
%   b        = RC3_dev_bss('final', b)                SOC, energies, degradation
%
%Grid-connected request (Sarker et al. 2015; Zheng et al. 2014):
%   2  must-charge to keep a reserve for the swaps of the next hres hours
%      and to be able to return to the initial energy at the end of the day
%   3  'schedule': charge u(t)*Pch (u > 0) or discharge -u(t)*Pdis (u < 0);
%      'rule'    : valley filling while signal(t) <= alpha_ch(t), peak
%                  shaving while signal(t) >= alpha_dis(t)
%   4  limits: charger power, pool energy between soc_arr and soc_full
%Degradation (Xu et al. 2018): rainflow cycle counting of the pool SOC with
%a Woehler curve N(DoD) = N100 * DoD^-kp, plus calendar ageing.
out2 = [];
switch action
    case 'init'
        [S, P, Lnet] = varargin{:};
        n = numel(S.bss.n);  T = P.T;  dt = P.dt;  z = zeros(n, T);
        b = struct('n', n, 'T', T, 'dt', dt, 'P', z, 'Pch', z, 'Pdis', z, 'swap', z, ...
                   'E', zeros(n, T + 1), 'soc', zeros(n, T + 1), 'deg_k', zeros(n, 1), ...
                   'deg', 0, 'unmet', 0, 'Eswap', 0, 'Edem', 0, 'g_soc', 0, 'use', P.use.BSS);
        if ~P.use.BSS, return; end
        q = P.bss;
        b.eta  = q.eta;  b.esw = q.Ebat * (q.soc_full - q.soc_arr);      % kWh per swap
        b.Ecap = S.bss.Ecap(:);
        b.Elow = q.soc_arr * b.Ecap;  b.Ehigh = q.soc_full * b.Ecap;  b.E0 = q.soc0 * b.Ecap;
        b.Pchm = S.bss.n(:) .* min(S.bss.Pch(:), q.crate * q.Ebat);     % kW
        b.Pdism = q.b2g * b.Pchm;
        b.dem  = q.N * P.bss_arr;                                       % swaps per hour
        b.shD  = S.bss.Nbat(:) / max(sum(S.bss.Nbat), eps);
        b.shP  = b.Pchm / max(sum(b.Pchm), eps);
        b.rule = strcmp(q.control, 'rule');  b.price = strcmp(q.signal, 'price');
        Lpk    = max(max(abs(Lnet)), eps);
        if b.price, b.sig = P.price / max(P.price); else, b.sig = Lnet / Lpk; end
        b.Lpk  = 1000 * Lpk;                                            % kW
        b.ach  = S.bss.ach;  b.adis = S.bss.adis;  b.u = S.bss.u;
        b.m    = max(1, round(q.hres / dt));
        b.Ecur = b.E0;  b.E(:, 1) = b.E0;
        b.q    = q;
    case 'swap'
        t = varargin{1};
        if ~b.use, return; end
        req = b.dem(t) * b.dt * b.shD;
        ns  = min(req, max(b.Ecur - b.Elow, 0) / b.esw);
        b.Ecur = b.Ecur - ns * b.esw;  b.swap(:, t) = ns;
    case 'grid'
        t = varargin{1};  out2 = zeros(b.n, 1);
        if ~b.use || b.n == 0, return; end
        T = b.T;  dt = b.dt;  E = b.Ecur;  eta = b.eta;
        nxt  = mod(t + (0:b.m-1), T) + 1;
        Eres = min(b.Elow + b.esw * dt * sum(b.dem(nxt)) * b.shD, b.Ehigh);
        Eres = max(Eres, b.E0 - eta * b.Pchm * dt * (T - t));
        pmust = max(Eres - E, 0) / (eta * dt);
        if b.rule
            if b.price
                pval = b.Pchm  .* (b.sig(t) <= b.ach(:, t));
                pdr  = b.Pdism .* (b.sig(t) >= b.adis(:, t));
            else
                pval = max(b.ach(:, t) - b.sig(t), 0) * b.Lpk .* b.shP;
                pdr  = max(b.sig(t) - b.adis(:, t), 0) * b.Lpk .* b.shP;
            end
        else
            pval = max(b.u(:, t), 0) .* b.Pchm;
            pdr  = max(-b.u(:, t), 0) .* b.Pdism;
        end
        pd = min([pdr, b.Pdism, max(E - Eres, 0) * eta / dt], [], 2);
        pc = max(pmust, pval);  pc(pd > 0) = 0;
        pc = min([pc, b.Pchm, max(b.Ehigh - E, 0) / (eta * dt)], [], 2);
        out2 = pc - pd;
    case 'island'
        t = varargin{1};  %#ok<NASGU>
        out2 = struct('Pdis', zeros(b.n, 1), 'Pch', zeros(b.n, 1));
        if ~b.use || b.n == 0, return; end
        out2.Pdis = min(b.Pdism, max(b.Ecur - b.Elow, 0) * b.eta / b.dt);
        out2.Pch  = min(b.Pchm,  max(b.Ehigh - b.Ecur, 0) / (b.eta * b.dt));
    case 'commit'
        [t, Pkw] = varargin{:};
        if ~b.use || b.n == 0, return; end
        pc = max(Pkw(:), 0);  pd = max(-Pkw(:), 0);
        b.Ecur = b.Ecur + b.eta * pc * b.dt - pd * b.dt / b.eta;
        b.Pch(:, t) = pc;  b.Pdis(:, t) = pd;  b.P(:, t) = pc - pd;  b.E(:, t + 1) = b.Ecur;
    case 'final'
        if ~b.use, return; end
        q = b.q;
        b.soc   = b.E ./ max(b.Ecap, eps);
        b.Eswap = b.esw * sum(b.swap(:));                               % kWh/day swapped
        b.Edem  = b.esw * b.dt * sum(b.dem);                            % kWh/day requested
        b.unmet = b.Edem - b.Eswap;
        b.g_soc = sum(max(b.E0 - b.E(:, end) - 1, 0)) / 1000;           % MWh (1 kWh tolerance)
        for k = find(b.Ecap' > 0)
            if strcmp(q.deg.model, 'linear')                            % throughput model
                thr = b.dt * (b.eta * sum(b.Pch(k, :)) + sum(b.Pdis(k, :)) / b.eta);
                dmg = thr / (2 * b.Ecap(k) * q.deg.dod_ref * q.deg.N100 * q.deg.dod_ref^(-q.deg.kp));
            else                                                        % rainflow + Woehler
                [dod, cnt] = RC3_rainflow(b.soc(k, :));
                dmg = sum(cnt .* dod.^q.deg.kp) / q.deg.N100;
            end
            if q.deg.swap                                               % mobility cycles
                dmg = dmg + sum(b.swap(k, :)) * (q.soc_full - q.soc_arr)^q.deg.kp ...
                            / q.deg.N100 * q.Ebat / b.Ecap(k);
            end
            dmg = dmg + 1 / (365 * q.deg.cal_life);                     % calendar ageing
            b.deg_k(k) = q.c_bat * b.Ecap(k) * dmg;
        end
        b.deg = sum(b.deg_k);
    otherwise
        error('RC3_dev_bss: unknown action %s', action);
end
end
