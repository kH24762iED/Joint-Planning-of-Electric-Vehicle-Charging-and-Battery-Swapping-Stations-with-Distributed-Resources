function x = RC3_set(P, x, name, value, k, t)
%RC3_SET  Set a decision-variable block of x by its name (see P.blocks or
%Tb.variables). Useful to build hand-made solutions for tests or to seed
%the optimizer.
%
%   x = RC3_set(P, x, 'DG_y', [1 1 0])         whole block (numel = block size)
%   x = RC3_set(P, x, 'DG_P', 0.6, 1, 1:P.T)   unit (row) k, slots (columns) t
%   x = RC3_set(P, x, 'DG_bus', 14, 1)         sites are given as bus / branch
%   x = RC3_set(P, x, 'TX_br', 6, 1)           numbers and converted to the
%   x = RC3_set(P, x, 'DG_type', 'MT', 1)      candidate index; DG types by name
%   x = RC3_base_vector(P) gives the base case to start from.
if isempty(x), x = P.xbase; end
assert(isfield(P.bidx, name), 'RC3_set: block %s is not part of this model.', name);
b = P.blocks(strcmp({P.blocks.name}, name));
cand = struct('DG_bus', 'dg', 'SC_bus', 'sc', 'ST_bus', 'st', 'EV_bus', 'ev', ...
              'BSS_bus', 'bss', 'TX_br', 'tx', 'UP_br', 'up');
if isfield(cand, name)
    lst = P.(cand.(name)).cand;
    [ok, value] = ismember(value, lst);
    assert(all(ok), 'RC3_set: %s is not a candidate site of %s.', mat2str(value(~ok)), name);
end
if strcmp(name, 'DG_type') && (ischar(value) || iscell(value))
    [~, ty] = ismember(upper(cellstr(value)), upper(P.dg.types));
    [ok, value] = ismember(ty, P.dg.allowed);
    assert(all(ok), 'RC3_set: DG type not allowed (P.dg.allowed).');
end
M = reshape(x(b.idx), b.dims);
if nargin < 5
    M(:) = value;
else
    if nargin < 6, t = 1:b.dims(2); end
    M(k, t) = value;
end
x(b.idx) = M(:)';
end
