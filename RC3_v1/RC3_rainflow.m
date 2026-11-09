function [dod, cnt] = RC3_rainflow(s)
%RC3_RAINFLOW  Cycle depths and counts (1 = full, 0.5 = half cycle) of the
%signal s by the ASTM E1049 three-point rainflow algorithm (no toolbox needed).
s   = s(:)';  s = s([true, abs(diff(s)) > 1e-12]);        % drop flat steps
dod = zeros(1, 0);  cnt = zeros(1, 0);
if numel(s) < 2, return; end
d   = diff(s);
tp  = s([1, find(d(1:end-1) .* d(2:end) < 0) + 1, numel(s)]);   % turning points
st  = zeros(1, 0);
for p = tp
    st(end+1) = p;                                              %#ok<AGROW>
    while numel(st) >= 3
        X = abs(st(end) - st(end-1));  Y = abs(st(end-1) - st(end-2));
        if X < Y, break; end
        dod(end+1) = Y;                                         %#ok<AGROW>
        if numel(st) == 3                                       % range contains the start
            cnt(end+1) = 0.5;  st(1) = [];                      %#ok<AGROW>
        else
            cnt(end+1) = 1;    st(end-2:end-1) = [];            %#ok<AGROW>
        end
    end
end
for i = 1:numel(st) - 1                                         % residue: half cycles
    dod(end+1) = abs(st(i+1) - st(i));  cnt(end+1) = 0.5;       %#ok<AGROW>
end
end
