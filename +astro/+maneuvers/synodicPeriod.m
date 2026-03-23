function Tsyn = synodicPeriod(T1, T2)
%SYNODICPERIOD Computes synodic period between two orbital periods.
%
% INPUTS
%   T1, T2 : orbital periods [s]
%
% OUTPUT
%   Tsyn   : synodic period [s]
%
% RELATION
%   1/Tsyn = |1/T1 - 1/T2|

den = abs(1./T1 - 1./T2);

if any(den == 0)
    error('Synodic period is undefined when the two periods are equal.');
end

Tsyn = 1 ./ den;
end