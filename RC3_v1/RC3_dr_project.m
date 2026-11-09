function s = RC3_dr_project(s, w, lo, hi)
%RC3_DR_PROJECT  Load-shift shares s (1 x T) made energy-neutral, w*s' = 0,
%inside the box lo <= s <= hi, by alternating projections onto the
%hyperplane and the box. A vector that already satisfies the balance is
%returned unchanged (the repair is idempotent); if no balanced point is
%found the shift is set to zero.
w   = reshape(w, 1, []);  s = reshape(s, 1, []);
tol = 1e-10 * sum(abs(w));
if abs(s * w') <= tol, return; end
ww = w * w';
for it = 1:100
    s = s - w * ((s * w') / ww);
    s = min(max(s, lo), hi);
    if abs(s * w') <= tol, return; end
end
s(:) = 0;
end
