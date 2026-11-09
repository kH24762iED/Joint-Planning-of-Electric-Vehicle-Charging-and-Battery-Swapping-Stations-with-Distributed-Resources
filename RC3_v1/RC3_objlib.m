function L = RC3_objlib()
%RC3_OBJLIB  Library of user objective functions built as weighted sums of
%the basic objectives (P.allObj). Every basic objective in an expression
%is normalised by its reference value P.objRef (base-case value, or the
%default scale when the base value is zero), so the terms are
%dimensionless and the base case of each index f2..f5 equals the sum of
%the weights of its active terms (1 when all terms are active).
%
%   Use the names in opts.objectives, e.g. {'f1','f2','f4'}, or write any
%   other expression directly, e.g. '0.5*LOSS + 0.5*ENS'. Edit the weights
%   below to change the library; a single basic name ('TAC') keeps its raw
%   value and unit.
%
%   f1  economy        : TAC
%   f2  technical      : losses, voltage deviation, voltage stability
%   f3  environment and e-mobility : CO2, EV charging and swap service,
%                        EV waiting time, battery degradation
%   f4  reliability    : ENS, SAIDI, SAIFI
%   f5  islanding and demand response : priority-weighted load shed,
%                        frequency deviation, supply-adequacy deficit, DR cost,
%                        DR discomfort, peak-to-average ratio
L = struct('name', {}, 'expr', {}, 'label', {});
L(end+1) = struct('name', 'f1', 'expr', 'TAC', ...
                  'label', 'Total annual cost (USD/yr)');
L(end+1) = struct('name', 'f2', 'expr', '0.3*LOSS + 0.3*VD + 0.4*VSI', ...
                  'label', 'Technical index (p.u.)');
L(end+1) = struct('name', 'f3', 'expr', '0.4*EMISSION + 0.2*EVCSUN + 0.15*EVWAIT + 0.15*BSSUN + 0.1*DEGCOST', ...
                  'label', 'Environmental and e-mobility index (p.u.)');
L(end+1) = struct('name', 'f4', 'expr', '0.3*ENS + 0.35*SAIDI + 0.35*SAIFI', ...
                  'label', 'Reliability index (p.u.)');
L(end+1) = struct('name', 'f5', 'expr', '0.3*WLSHED + 0.15*FDEV + 0.15*SAD + 0.15*DRCOST + 0.1*DISCOMF + 0.15*PAR', ...
                  'label', 'Islanding and demand-response index (p.u.)');
end
