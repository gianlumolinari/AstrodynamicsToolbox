function file = getValidatedOrbitPath(name)
%GETVALIDATEDORBITPATH Return the file path for a validated orbit key.
%
% INPUT
%   name : string
%
% Examples:
%   'planar_lyapunov_L1'
%   'planar_lyapunov_L1_small'
%   'halo_L1_north'
%   'halo_L1_north_medium'

switch lower(string(name))
    case "planar_lyapunov_l1"
        file = 'data/validated/cr3bp/earth_moon/planar_lyapunov_L1_validated.mat';

    case "planar_lyapunov_l1_small"
        file = 'data/validated/cr3bp/earth_moon/planar_lyapunov_L1_small.mat';

    case "planar_lyapunov_l1_medium"
        file = 'data/validated/cr3bp/earth_moon/planar_lyapunov_L1_medium.mat';

    case "planar_lyapunov_l1_large"
        file = 'data/validated/cr3bp/earth_moon/planar_lyapunov_L1_large.mat';

    case "halo_l1_north"
        file = 'data/validated/cr3bp/earth_moon/halo_L1_north_validated.mat';

    case "halo_l1_north_small"
        file = 'data/validated/cr3bp/earth_moon/halo_L1_north_small.mat';

    case "halo_l1_north_medium"
        file = 'data/validated/cr3bp/earth_moon/halo_L1_north_medium.mat';

    case "halo_l1_north_large"
        file = 'data/validated/cr3bp/earth_moon/halo_L1_north_large.mat';

    case "halo_l2_south"
        file = 'data/validated/cr3bp/earth_moon/halo_L2_south_validated.mat';

    otherwise
        error('Unknown validated orbit name: %s', name);
end
end