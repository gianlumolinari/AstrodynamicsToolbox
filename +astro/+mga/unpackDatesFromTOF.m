function dates = unpackDatesFromTOF(t0UTC, tofDays)
%UNPACKDATESFROMTOF Build encounter date strings from departure epoch and leg TOFs.
%
% INPUTS
%   t0UTC   : departure UTC string, e.g. '2026-01-01 00:00:00'
%   tofDays : vector of leg times of flight [days]
%
% OUTPUT
%   dates   : cell array of UTC strings, one per encounter

t0 = datetime(t0UTC, 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
nLegs = numel(tofDays);

dates = cell(nLegs + 1, 1);
dates{1} = datestr(t0, 'yyyy-mm-dd HH:MM:SS');

tk = t0;
for k = 1:nLegs
    tk = tk + days(tofDays(k));
    dates{k+1} = datestr(tk, 'yyyy-mm-dd HH:MM:SS');
end
end