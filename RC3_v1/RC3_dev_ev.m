function ev = RC3_dev_ev(S, P)
%RC3_DEV_EV  EV charging stations as M/M/c queues (Bayram et al. 2013).
%The daily EV arrivals follow P.ev_arr and are shared among the stations
%in proportion to their charging capacity. In slot t the charger power is
%limited to capf(t) * rated power (smart-charging cap, a decision
%variable), so the service rate is mu = eta * Pch * capf / E. A station
%serves at most rho_max of its service capacity; the rest of the demand is
%lost (unserved). EVs that are shed in islanded slots are added by
%RC3_objectives.
%   ev.P (n x T, kW)  ev.lam  ev.Wq (min)  ev.unmet / Edel / Edem (kWh/day)
%   ev.wait mean waiting time over all EVs, lost EVs included (min)
T = P.T;  n = numel(S.ev.n);
ev = struct('P', zeros(n, T), 'lam', zeros(n, T), 'Wq', zeros(n, T), ...
            'unmet', 0, 'Edel', 0, 'Edem', 0, 'wait', 0, 'lost', zeros(1, T));
if ~P.use.EVCS, return; end
lamTot = P.ev.N * P.ev_arr;                       % EV arrivals (1/h), 1 x T
num    = 0;
for t = 1:T
    if n > 0
        Pc   = S.ev.Pch .* S.ev.capf(:, t);       % kW per charger in slot t
        cap  = S.ev.n .* Pc;                      % station capacity (kW)
        mu   = P.ev.eta * Pc / P.ev.E;            % service rate of one charger (1/h)
        lam  = lamTot(t) * cap / max(sum(cap), eps);
        lamS = min(lam, P.ev.rho_max * S.ev.n .* mu);
        for k = find(lamS > 0)'
            ev.Wq(k, t) = 60 * erlangC(S.ev.n(k), lamS(k) / mu(k)) / (S.ev.n(k) * mu(k) - lamS(k));
        end
    else
        lamS = zeros(0, 1);
    end
    lost = lamTot(t) - sum(lamS);
    num  = num + sum(lamS .* ev.Wq(:, t)) + lost * P.ev.w_pen;
    ev.lam(:, t) = lamS;
    ev.P(:, t)   = lamS * P.ev.E / P.ev.eta;      % average grid power (kW)
    ev.lost(t)   = lost;
    ev.unmet = ev.unmet + lost * P.ev.E * P.dt;   % kWh/day not delivered
    ev.Edel  = ev.Edel + sum(lamS) * P.ev.E * P.dt;
    ev.Edem  = ev.Edem + lamTot(t) * P.ev.E * P.dt;
end
ev.wait = num / max(sum(lamTot), eps);
end

function C = erlangC(c, a)
% Erlang-C probability of waiting for an M/M/c queue (a = lambda/mu < c)
B = 1;
for k = 1:c, B = a * B / (k + a * B); end         % Erlang-B recursion
C = B / (1 - (a / c) * (1 - B));
end
