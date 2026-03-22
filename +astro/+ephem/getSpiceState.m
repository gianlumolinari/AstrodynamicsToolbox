function state = getSpiceState(target, epochUTC, observer, frame, abcorr)
%GETSPICESTATE Returns a Cartesian state from SPICE.
%
% INPUTS
%   target   : SPICE target, e.g. 'EARTH BARYCENTER', 'MARS BARYCENTER', 'SUN'
%   epochUTC : UTC epoch string, e.g. '2026-11-12 00:00:00'
%   observer : observer, default 'SUN'
%   frame    : reference frame, default 'J2000'
%   abcorr   : aberration correction, default 'NONE'
%
% OUTPUT
%   state : struct with fields
%       .target
%       .observer
%       .epochUTC
%       .et
%       .r
%       .v
%
% Notes:
%   For interplanetary planetary states, barycenters are often the safest
%   choice with generic DE ephemeris kernels.

if nargin < 3 || isempty(observer)
    observer = 'SUN';
end
if nargin < 4 || isempty(frame)
    frame = 'J2000';
end
if nargin < 5 || isempty(abcorr)
    abcorr = 'NONE';
end

target = char(string(target));
observer = char(string(observer));
frame = char(string(frame));
abcorr = char(string(abcorr));

% ----------------------------------------------------------
% Friendly mapping to robust SPICE targets
% ----------------------------------------------------------
switch upper(strtrim(target))
    case 'MERCURY'
        target = 'MERCURY BARYCENTER';
    case 'VENUS'
        target = 'VENUS BARYCENTER';
    case 'EARTH'
        target = 'EARTH BARYCENTER';
    case 'MARS'
        target = 'MARS BARYCENTER';
    case 'JUPITER'
        target = 'JUPITER BARYCENTER';
    case 'SATURN'
        target = 'SATURN BARYCENTER';
    case 'URANUS'
        target = 'URANUS BARYCENTER';
    case 'NEPTUNE'
        target = 'NEPTUNE BARYCENTER';
end

switch upper(strtrim(observer))
    case 'MERCURY'
        observer = 'MERCURY BARYCENTER';
    case 'VENUS'
        observer = 'VENUS BARYCENTER';
    case 'EARTH'
        observer = 'EARTH BARYCENTER';
    case 'MARS'
        observer = 'MARS BARYCENTER';
    case 'JUPITER'
        observer = 'JUPITER BARYCENTER';
    case 'SATURN'
        observer = 'SATURN BARYCENTER';
    case 'URANUS'
        observer = 'URANUS BARYCENTER';
    case 'NEPTUNE'
        observer = 'NEPTUNE BARYCENTER';
end

et = cspice_str2et(epochUTC);
sv = cspice_spkezr(target, et, frame, abcorr, observer);

state.target = target;
state.observer = observer;
state.epochUTC = epochUTC;
state.et = et;
state.r = sv(1:3);
state.v = sv(4:6);
end