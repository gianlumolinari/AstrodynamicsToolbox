function orbit = parseJPLPeriodicOrbit(data, idx)
%PARSEJPLPERIODICORBIT Parse one orbit row from JPL periodic-orbit API output.
%
% INPUTS
%   data : struct returned by astro.cr3bp.queryJPLPeriodicOrbits
%   idx  : row index into data.data
%
% OUTPUT
%   orbit : struct with fields
%       .state      [6x1]
%       .jacobi
%       .period
%       .stability
%       .mu
%       .lunit_km
%       .tunit_s
%       .systemName
%       .family
%       .branch
%       .librationPoint
%
% NOTES
%   The API provides the field names in data.fields and the rows in data.data.

if nargin < 2 || isempty(idx)
    idx = 1;
end

if ~isfield(data, 'data') || isempty(data.data)
    error('No orbit rows available in JPL periodic-orbit response.');
end

if idx < 1 || idx > numel(data.data)
    error('Requested index exceeds number of returned orbit rows.');
end

row = data.data{idx};
fields = data.fields;

orbit = struct();
orbit.state = zeros(6,1);

for k = 1:numel(fields)
    fname = lower(string(fields{k}));
    val = str2double(row{k});

    switch fname
        case "x"
            orbit.state(1) = val;
        case "y"
            orbit.state(2) = val;
        case "z"
            orbit.state(3) = val;
        case "vx"
            orbit.state(4) = val;
        case "vy"
            orbit.state(5) = val;
        case "vz"
            orbit.state(6) = val;
        case "jacobi"
            orbit.jacobi = val;
        case "period"
            orbit.period = val;
        case "stability"
            orbit.stability = val;
    end
end

orbit.mu = str2double(data.system.mass_ratio);
orbit.lunit_km = str2double(data.system.lunit);
orbit.tunit_s = str2double(data.system.tunit);
orbit.systemName = data.system.name;
orbit.family = data.family;

if isfield(data, 'branch')
    orbit.branch = data.branch;
else
    orbit.branch = '';
end

if isfield(data, 'libration_point')
    orbit.librationPoint = data.libration_point;
else
    orbit.librationPoint = '';
end
end