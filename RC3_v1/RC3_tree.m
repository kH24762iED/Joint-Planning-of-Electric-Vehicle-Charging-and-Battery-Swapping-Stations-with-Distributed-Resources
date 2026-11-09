function [par, pbr, order, isRadial, nReach] = RC3_tree(fb, tb, closed, nb, root)
%RC3_TREE  Breadth-first search of the network formed by the closed branches,
%starting from the root bus (default: the substation, bus 1).
%   par(i)  : parent bus of bus i          pbr(i) : branch feeding bus i
%   order   : buses in search order        nReach : buses reached from root
%   isRadial: true if the closed branches form a spanning tree, i.e.
%             exactly nb-1 closed branches AND all buses connected.
if nargin < 5 || isempty(root), root = 1; end
closed = logical(closed(:));
adj = cell(nb, 1);
for e = find(closed)'
    adj{fb(e)}(end+1) = e;
    adj{tb(e)}(end+1) = e;
end
par = zeros(nb, 1);  pbr = zeros(nb, 1);  order = zeros(nb, 1);
seen = false(nb, 1);  seen(root) = true;  order(1) = root;
head = 1;  tail = 1;
while head <= tail
    u = order(head);  head = head + 1;
    for e = adj{u}
        v = fb(e) + tb(e) - u;                       % other end of branch e
        if ~seen(v)
            seen(v) = true;  par(v) = u;  pbr(v) = e;
            tail = tail + 1;  order(tail) = v;
        end
    end
end
nReach   = tail;
order    = order(1:tail);
isRadial = (nReach == nb) && (nnz(closed) == nb - 1);
end
