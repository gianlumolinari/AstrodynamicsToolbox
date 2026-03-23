function dates = unpackDatesFromDecisionVector(t0UTC, x)
%UNPACKDATESFROMDECISIONVECTOR Build encounter dates from departure offset + TOFs.
%
% INPUTS
%   t0UTC : reference UTC string
%   x     : [departureOffsetDays; tof1; tof2; ...]
%
% OUTPUT
%   dates : cell array of UTC strings

t0 = datetime(t0UTC, 'InputFormat', 'yyyy-MM-dd HH:mm:ss');

depOffsetDays = x(1);
tofDays = x(2:end);

tk = t0 + days(depOffsetDays);

dates = cell(numel(tofDays) + 1, 1);
dates{1} = datestr(tk, 'yyyy-mm-dd HH:MM:SS');

for k = 1:numel(tofDays)
    tk = tk + days(tofDays(k));
    dates{k+1} = datestr(tk, 'yyyy-mm-dd HH:MM:SS');
end
end