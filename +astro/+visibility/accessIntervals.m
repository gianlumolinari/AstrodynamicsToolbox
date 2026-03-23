function intervals = accessIntervals(t, visibleMask)
%ACCESSINTERVALS Extracts contiguous visibility intervals from a boolean mask.
%
% INPUTS
%   t           : time vector [Nx1] numeric or datetime
%   visibleMask : logical visibility vector [Nx1]
%
% OUTPUT
%   intervals : struct array with fields
%       .start
%       .stop
%       .duration
%       .iStart
%       .iStop

t = t(:);
visibleMask = logical(visibleMask(:));

if numel(t) ~= numel(visibleMask)
    error('t and visibleMask must have the same length.');
end

d = diff([false; visibleMask; false]);
iStart = find(d == 1);
iStop = find(d == -1) - 1;

nInt = numel(iStart);
intervals = struct('start', {}, 'stop', {}, 'duration', {}, 'iStart', {}, 'iStop', {});

for k = 1:nInt
    intervals(k).start = t(iStart(k));
    intervals(k).stop = t(iStop(k));
    intervals(k).duration = t(iStop(k)) - t(iStart(k));
    intervals(k).iStart = iStart(k);
    intervals(k).iStop = iStop(k);
end
end