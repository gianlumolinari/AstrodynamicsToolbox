function key = makeOrbitKey(familyType, libr, branch, qualifier)
%MAKEORBITKEY Build a canonical orbit key string.
%
% Examples:
%   makeOrbitKey('planar_lyapunov',1,'','validated')
%     -> 'planar_lyapunov_L1_validated'
%
%   makeOrbitKey('halo',1,'north','medium')
%     -> 'halo_L1_north_medium'

if nargin < 3 || isempty(branch)
    branch = '';
end
if nargin < 4 || isempty(qualifier)
    qualifier = '';
end

base = sprintf('%s_L%d', lower(string(familyType)), libr);

if strlength(string(branch)) > 0
    base = sprintf('%s_%s', base, lower(string(branch)));
end

if strlength(string(qualifier)) > 0
    key = sprintf('%s_%s', base, lower(string(qualifier)));
else
    key = base;
end
end