function data = queryJPLPeriodicOrbits(varargin)
%QUERYJPLPERIODICORBITS Query JPL Three-Body Periodic Orbits API.
%
% USAGE EXAMPLES
%   data = astro.cr3bp.queryJPLPeriodicOrbits( ...
%       'sys','earth-moon','family','halo','libr',1,'branch','N');
%
%   data = astro.cr3bp.queryJPLPeriodicOrbits( ...
%       'sys','earth-moon','family','lyapunov','libr',1, ...
%       'periodmax', 3.0, 'periodunits','TU');
%
% NOTES
%   JPL API docs:
%   https://ssd-api.jpl.nasa.gov/doc/periodic_orbits.html
%
%   Required parameters always include:
%       sys, family
%   Depending on family, libr and/or branch may also be required.

p = inputParser;
addParameter(p, 'sys', '', @(x)ischar(x)||isstring(x));
addParameter(p, 'family', '', @(x)ischar(x)||isstring(x));
addParameter(p, 'libr', [], @(x)isempty(x)||isnumeric(x));
addParameter(p, 'branch', '', @(x)ischar(x)||isstring(x));
addParameter(p, 'periodmin', [], @(x)isempty(x)||isnumeric(x));
addParameter(p, 'periodmax', [], @(x)isempty(x)||isnumeric(x));
addParameter(p, 'periodunits', '', @(x)ischar(x)||isstring(x));
addParameter(p, 'jacobimin', [], @(x)isempty(x)||isnumeric(x));
addParameter(p, 'jacobimax', [], @(x)isempty(x)||isnumeric(x));
addParameter(p, 'stabmin', [], @(x)isempty(x)||isnumeric(x));
addParameter(p, 'stabmax', [], @(x)isempty(x)||isnumeric(x));
parse(p, varargin{:});

R = p.Results;

if strlength(string(R.sys)) == 0 || strlength(string(R.family)) == 0
    error('JPL periodic orbit query requires at least sys and family.');
end

baseURL = 'https://ssd-api.jpl.nasa.gov/periodic_orbits.api';

q = {};
q(end+1:end+2) = {'sys', char(string(R.sys))};
q(end+1:end+2) = {'family', char(string(R.family))};

if ~isempty(R.libr)
    q(end+1:end+2) = {'libr', num2str(R.libr)};
end
if strlength(string(R.branch)) > 0
    q(end+1:end+2) = {'branch', char(string(R.branch))};
end
if ~isempty(R.periodmin)
    q(end+1:end+2) = {'periodmin', num2str(R.periodmin, '%.16g')};
end
if ~isempty(R.periodmax)
    q(end+1:end+2) = {'periodmax', num2str(R.periodmax, '%.16g')};
end
if strlength(string(R.periodunits)) > 0
    q(end+1:end+2) = {'periodunits', char(string(R.periodunits))};
end
if ~isempty(R.jacobimin)
    q(end+1:end+2) = {'jacobimin', num2str(R.jacobimin, '%.16g')};
end
if ~isempty(R.jacobimax)
    q(end+1:end+2) = {'jacobimax', num2str(R.jacobimax, '%.16g')};
end
if ~isempty(R.stabmin)
    q(end+1:end+2) = {'stabmin', num2str(R.stabmin, '%.16g')};
end
if ~isempty(R.stabmax)
    q(end+1:end+2) = {'stabmax', num2str(R.stabmax, '%.16g')};
end

opts = weboptions('Timeout', 30);
data = webread(baseURL, q{:}, opts);
end