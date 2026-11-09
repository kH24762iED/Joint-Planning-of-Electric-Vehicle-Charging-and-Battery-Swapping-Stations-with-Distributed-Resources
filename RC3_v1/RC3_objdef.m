function D = RC3_objdef(spec, P)
%RC3_OBJDEF  Parse the objective specification (opts.objectives) into the
%user objective functions evaluated by RC3_objfun.
%
%   D = RC3_objdef(spec, P)      spec: cell array, one entry per objective
%
%   Entry                                   Meaning
%   'TAC'                                   basic objective, raw value (its unit)
%   'f2'                                    function of the library RC3_objlib
%   '0.3*LOSS + 0.3*VD + 0.4*VSI'           sum_j w_j F_j / ref_j (normalised)
%   'tech = 0.5*LOSS + 0.5*VD'              the same, named 'tech'
%   struct('name','f6', 'expr','0.5*LOSS + 0.5*ENS', 'norm',true, 'label','...')
%                                           norm = false: raw weighted sum
%   @(V, Vn) Vn.LOSS^2 + Vn.VD              any function of the raw values V and
%   struct('name','f7', 'fn', @(V,Vn) ..)   normalised values Vn (structs with one
%                                           field per basic objective)
%
%   D(k): .name .label .expr .type ('basic' | 'linear' | 'handle')
%         .idx (indices in P.allObj) .w (weights) .norm (logical) .fn
if ~iscell(spec), spec = {spec}; end
lib = RC3_objlib();
D = repmat(struct('name', '', 'label', '', 'expr', '', 'type', '', 'idx', [], 'w', [], ...
                  'norm', false, 'fn', []), 0, 1);
for k = 1:numel(spec)
    s = spec{k};
    d = struct('name', sprintf('F%d', k), 'label', '', 'expr', '', 'type', 'linear', ...
               'idx', [], 'w', [], 'norm', true, 'fn', []);
    normSet = false;
    if isa(s, 'function_handle')
        d.type = 'handle';  d.fn = s;  d.expr = func2str(s);  d.norm = false;
    elseif isstruct(s)
        if isfield(s, 'name') && ~isempty(s.name), d.name = s.name; end
        if isfield(s, 'label') && ~isempty(s.label), d.label = s.label; end
        if isfield(s, 'fn') && ~isempty(s.fn)
            d.type = 'handle';  d.fn = s.fn;  d.expr = func2str(s.fn);  d.norm = false;
        else
            assert(isfield(s, 'expr') && ~isempty(s.expr), 'RC3_objdef: struct entry %d needs expr or fn.', k);
            d.expr = s.expr;
            if isfield(s, 'norm') && ~isempty(s.norm), d.norm = logical(s.norm); normSet = true; end
        end
    elseif ischar(s) || isstring(s)
        s  = strtrim(char(s));
        il = find(strcmpi({lib.name}, s), 1);
        eq = regexp(s, '^([A-Za-z]\w*)\s*=\s*(.+)$', 'tokens', 'once');
        if ~isempty(il)
            d.name = lib(il).name;  d.label = lib(il).label;  d.expr = lib(il).expr;
        elseif ~isempty(eq)
            d.name = eq{1};  d.expr = eq{2};
        elseif any(strcmpi(P.allObj, s))
            d.name = upper(s);  d.expr = upper(s);
        else
            d.expr = s;
        end
    else
        error('RC3_objdef: entry %d must be a name, an expression, a struct or a function handle.', k);
    end
    if ~strcmp(d.type, 'handle')
        [d.idx, d.w] = parse_expr(d.expr, P.allObj);
        if isscalar(d.idx) && d.w == 1 && ~normSet
            d.type = 'basic';  d.norm = false;               % one basic objective: raw value
        end
    end
    if isempty(d.label)
        if strcmp(d.type, 'basic')
            d.label = sprintf('%s (%s)', P.objDesc{d.idx}, P.objUnit{d.idx});
        else
            d.label = d.expr;
        end
    end
    D(end+1, 1) = d; %#ok<AGROW>
end
nm = matlab.lang.makeValidName({D.name});
nm = matlab.lang.makeUniqueStrings(nm);
for k = 1:numel(D), D(k).name = nm{k}; end
end

% =========================================================================================
function [idx, w] = parse_expr(expr, names)
% linear combination  [+-][coef*]NAME [+-][coef*]NAME ...
s = regexprep(char(expr), '\s+', '');
[tok, mt] = regexp(s, '([+-]?)((?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?\*|)([A-Za-z]\w*)', 'tokens', 'match');
if isempty(tok) || ~strcmp([mt{:}], s)
    error('RC3_objdef:parse', ['Cannot parse the objective "%s". Write a weighted sum such as ' ...
          '''0.3*LOSS + 0.3*VD + 0.4*VSI''.'], expr);
end
idx = zeros(1, numel(tok));  w = ones(1, numel(tok));
for j = 1:numel(tok)
    c = tok{j}{2};
    if ~isempty(c), w(j) = str2double(c(1:end-1)); end
    if strcmp(tok{j}{1}, '-'), w(j) = -w(j); end
    i = find(strcmpi(names, tok{j}{3}), 1);
    if isempty(i)
        error('RC3_objdef:name', 'Unknown objective "%s" in "%s". Basic objectives: %s', ...
              tok{j}{3}, expr, strjoin(names, ', '));
    end
    idx(j) = i;
end
[u, ~, ic] = unique(idx, 'stable');                         % merge repeated names
w   = accumarray(ic(:), w(:))';
idx = u;
end
