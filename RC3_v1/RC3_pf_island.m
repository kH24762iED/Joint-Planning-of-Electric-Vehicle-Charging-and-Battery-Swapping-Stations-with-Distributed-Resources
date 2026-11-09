function [r, il] = RC3_pf_island(net, gf, L, Ysh, pv, up, P, Vinit)
%RC3_PF_ISLAND  Load flow of one islanded radial component (microgrid) with
%droop-controlled grid-forming (GF) units, after Hameed et al. (2019) and
%Mumtaz et al. (2016):
%   P-f droop  P_g = (f0 - f)/m_g,  m_g = df_max / Pmax_g   (sharing by rating)
%   Q-V droop  Q_g = (V_nom - |V_g|)/n_g, n_g = dV_max / Qmax_g
%The GF unit with the largest rating is the reference (root, node 1): its
%voltage follows its Q-V droop and it balances the island; the other GF
%units are P/Q injections updated in an outer loop until the losses, the
%root voltage and the droop reactive powers converge.
%If the GF capability cannot cover the island demand, load is shed in the
%order EV charging -> loads of lowest priority weight (Chen et al. 2016);
%excess renewable output is curtailed when the GF units cannot absorb it.
%
%   gf : .node (local, gf.node(1) = 1) .Pmax .Pabs (absorbing) .Qmax  (p.u.)
%   L  : local vectors (p.u., load convention at 1 p.u.):
%        .Pl .Ql loads   .Pev .Qev EV charging   .Pres .Qres renewable
%        injections (positive = generation)   .Pfix .Qfix other net load
%        .w priority weights
%   il : .phi loading of the GF units (-1..1) .df frequency deviation (Hz)
%        .shed .shedEV (p.u., per node) .fracEV .fracCurt .Pgf .Qgf (p.u.)
%        .SA supply adequacy (GF + renewable capability / demand) .Vroot .ok
n   = net.n;  ng = numel(gf.node);
isl = P.isl;  lim = P.lim;
PmaxI = sum(gf.Pmax);  PabsI = sum(gf.Pabs);
nq  = isl.dV_max ./ max(gf.Qmax, 1e-9);
Qg  = zeros(ng, 1);  Vroot = isl.V_nom;  loss = 0;  V = [];
if nargin > 7, V = Vinit; end
dem0  = sum(L.Pl) + sum(L.Pev) + sum(L.Pfix);         % demand before shedding
resav = sum(L.Pres);
il.SA = (PmaxI + resav) / max(dem0, 1e-9);
[~, ord] = sort(L.w, 'ascend');
il.ok = false;
for outer = 1:30
    D = dem0 - resav + loss;                           % demand on the GF units
    shed = zeros(n, 1);  fracEV = 0;  fracC = 0;
    if D > PmaxI                                       % shed: EV first, then priority
        need = D - PmaxI;  sev = sum(L.Pev);
        if sev > 0, fracEV = min(1, need / sev);  need = need - fracEV * sev; end
        for i = ord(:)'
            if need <= 0, break; end
            s = min(L.Pl(i), need);  shed(i) = s;  need = need - s;
        end
    elseif D < -PabsI && resav > 0                     % curtail renewables
        fracC = min(1, (-PabsI - D) / resav);
    end
    Dn  = D - sum(shed) - fracEV * sum(L.Pev) + fracC * resav;
    if Dn >= 0, phi = Dn / max(PmaxI, 1e-12); else, phi = Dn / max(PabsI, 1e-12); end
    phi = min(max(phi, -1), 1);
    if phi >= 0, Pg = phi * gf.Pmax; else, Pg = phi * gf.Pabs; end
    ks  = ones(n, 1);  ks(L.Pl > 0) = 1 - shed(L.Pl > 0) ./ L.Pl(L.Pl > 0);
    S0  = complex(L.Pl .* ks + (1 - fracEV) * L.Pev - (1 - fracC) * L.Pres + L.Pfix, ...
                  L.Ql .* ks + (1 - fracEV) * L.Qev - (1 - fracC) * L.Qres + L.Qfix);
    if ng > 1                                          % non-reference GF units
        S0 = S0 - complex(accumarray(gf.node(2:end), Pg(2:end), [n 1]), ...
                          accumarray(gf.node(2:end), Qg(2:end), [n 1]));
    end
    Vr = min(max(isl.V_nom - nq(1) * Qg(1), isl.V_nom - isl.dV_max), isl.V_nom + isl.dV_max);
    r  = RC3_pf_fbs(net, Vr, S0, Ysh, lim.zip_np, lim.zip_nq, pv, up, lim.fbs_tol, lim.fbs_maxit, V);
    if ~r.success, break; end
    V = r.Vt;
    Qn = Qg;  Qn(1) = imag(r.Sroot);                   % reference unit balances Q
    if ng > 1
        Vg = abs(r.V(gf.node(2:end)));
        Qn(2:end) = min(max((isl.V_nom - Vg) ./ nq(2:end), -gf.Qmax(2:end)), gf.Qmax(2:end));
    end
    dl = abs(real(r.loss) - loss);  dq = max(abs(Qn - Qg));  dv = abs(Vr - Vroot);
    loss = real(r.loss);  Qg = 0.5 * Qg + 0.5 * Qn;  Qg(1) = Qn(1);  Vroot = Vr;
    if outer > 1 && dl < 1e-7 && dq < 1e-6 && dv < 1e-7, il.ok = true; break; end
end
if ~il.ok && r.success, il.ok = true; end              % accept the last converged sweep
Pg(1) = real(r.Sroot);
il.phi = phi;  il.df = phi * isl.df_max;  il.shed = shed;  il.fracEV = fracEV;
il.fracCurt = fracC;  il.Pgf = Pg;  il.Qgf = Qg;  il.Vroot = Vroot;  il.it = outer;
end
