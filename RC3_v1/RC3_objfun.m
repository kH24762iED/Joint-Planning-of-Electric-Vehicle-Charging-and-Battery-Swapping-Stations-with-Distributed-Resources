function [F, T] = RC3_objfun(Fall, P, D)
%RC3_OBJFUN  Values of the user objective functions from the basic objectives.
%
%   F      = RC3_objfun(Fall, P)       objectives P.objdef (opts.objectives)
%   [F, T] = RC3_objfun(Fall, P, D)    objectives D (e.g. P.objLib) and the
%                                     breakdown T of every weighted sum
%
%   Fall : 1 x numel(P.allObj), basic objectives from RC3_objectives
%   basic  : F_k = F_j                               (raw value)
%   linear : F_k = sum_j w_j * F_j / ref_j           (norm = true)
%            F_k = sum_j w_j * F_j                   (norm = false)
%            ref_j = P.objRef(j); F_k = P.bigM if any term is P.bigM
%   handle : F_k = fn(V, Vn), V.(name) = F_j, Vn.(name) = F_j / ref_j
%
%   T(k): .name .terms .w .raw .ref .z (normalised) .contrib (w .* z)
if nargin < 3, D = P.objdef; end
n = numel(D);  F = zeros(1, n);
if nargout > 1
    T = repmat(struct('name', '', 'terms', {{}}, 'w', [], 'raw', [], 'ref', [], 'z', [], 'contrib', []), n, 1);
end
for k = 1:n
    d = D(k);
    raw = [];  ref = [];  z = [];  c = [];
    switch d.type
        case 'basic'
            raw = Fall(d.idx);  ref = 1;  z = raw;  c = raw;
            F(k) = raw;
        case 'linear'
            raw = Fall(d.idx);
            if d.norm, ref = P.objRef(d.idx); else, ref = ones(size(raw)); end
            z = raw ./ ref;  c = d.w .* z;
            if any(raw >= P.bigM), F(k) = P.bigM; else, F(k) = sum(c); end
        case 'handle'
            V  = cell2struct(num2cell(Fall(:)), P.allObj(:), 1);
            Vn = cell2struct(num2cell(Fall(:) ./ P.objRef(:)), P.allObj(:), 1);
            F(k) = d.fn(V, Vn);
            if ~isfinite(F(k)), F(k) = P.bigM; end
    end
    if nargout > 1
        T(k) = struct('name', d.name, 'terms', {P.allObj(d.idx)}, 'w', d.w, 'raw', raw, ...
                      'ref', ref, 'z', z, 'contrib', c);
    end
end
end
