function varargout = RC3_om_eval(P, what, x, name)
%RC3_OM_EVAL  Access to the MP-Opt-Model object of RC3_setup / RC3_om_build.
%
%   [A, l, u]       = RC3_om_eval(P, 'lin')          linear constraints
%   [x0, xl, xu, vt] = RC3_om_eval(P, 'var')         variables
%   f = RC3_om_eval(P, 'cost', x, 'f_LOSS')          one named cost
%   f = RC3_om_eval(P, 'costall', x)                 sum of all costs (m = 1)
%   g = RC3_om_eval(P, 'nln', x)                     nonlinear constraints 'net'
%   RC3_om_eval(P, 'display')                        print the model
om = P.om;
if nargin > 2, x = reshape(x, [], 1); end
switch what
    case 'lin'
        if isempty(P.linsets)
            varargout = {sparse(0, P.D), zeros(0, 1), zeros(0, 1)};  return
        end
        switch P.om_api
            case 'mp',     [A, l, u] = om.lin.params(om.var);
            case 'legacy', [A, l, u] = om.params_lin_constraint();
            otherwise
                A = vertcat(om.lin.A);  l = vertcat(om.lin.l);  u = vertcat(om.lin.u);
        end
        varargout = {A, full(l(:)), full(u(:))};
    case 'var'
        switch P.om_api
            case 'mp',     [x0, xl, xu, vt] = om.var.params();
            case 'legacy', [x0, xl, xu, vt] = om.params_var();
            otherwise
                x0 = vertcat(om.var.v0);  xl = vertcat(om.var.vl);
                xu = vertcat(om.var.vu);  vt = [om.var.vt];
        end
        varargout = {x0, xl, xu, vt};
    case 'cost'
        switch P.om_api
            case 'mp',     f = om.nlc.eval(om.var, x, name);
            case 'legacy', f = om.eval_nln_cost(x, name);
            otherwise,     h = om.nlc(strcmp({om.nlc.name}, name)).fcn;  f = h(x);
        end
        varargout = {f};
    case 'costall'
        switch P.om_api
            case 'mp',     f = om.nlc.eval(om.var, x);
            case 'legacy', f = om.eval_nln_cost(x);
            otherwise
                f = 0;
                for k = 1:numel(om.nlc), h = om.nlc(k).fcn;  f = f + h(x); end
        end
        varargout = {f};
    case 'nln'
        switch P.om_api
            case 'mp',     g = om.nli.eval(om.var, x, 'net');
            case 'legacy', g = om.eval_nln_constraint(x, 0, 'net');
            otherwise,     h = om.nli(1).fcn;  g = h(x);
        end
        varargout = {g};
    case 'display'
        if any(strcmp(P.om_api, {'mp', 'legacy'}))
            om.display();
        else
            fprintf('Internal model: %d variable sets, %d linear sets, %d costs\n', ...
                    numel(om.var), numel(om.lin), numel(om.nlc));
        end
    otherwise
        error('RC3_om_eval: unknown request ''%s''.', what);
end
end
