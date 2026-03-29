function orbitData = loadValidatedOrbit(name)
%LOADVALIDATEDORBIT Load a validated periodic orbit from disk.
%
% INPUT
%   name : string
%
% Supported names:
%   'planar_lyapunov_L1'
%   'halo_L1_north'
%   'halo_L2_south'
%
% OUTPUT
%   orbitData : struct loaded from .mat file

switch lower(string(name))
    case "planar_lyapunov_l1"
        file = 'data/validated/cr3bp/earth_moon/planar_lyapunov_L1_validated.mat';

    case "halo_l1_north"
        file = 'data/validated/cr3bp/earth_moon/halo_L1_north_validated.mat';

    case "halo_l2_south"
        file = 'data/validated/cr3bp/earth_moon/halo_L2_south_validated.mat';

    otherwise
        error('Unknown validated orbit name: %s', name);
end

if ~exist(file, 'file')
    error('Validated orbit file not found: %s', file);
end

orbitData = load(file);
end