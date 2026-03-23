function M = monodromyMatrix(x0, T, mu)
%MONODROMYMATRIX Compute CR3BP monodromy matrix over one period.
%
% INPUTS
%   x0 : periodic-orbit initial state [6x1]
%   T  : orbit period
%   mu : CR3BP mass parameter
%
% OUTPUT
%   M  : monodromy matrix [6x6]

x0 = x0(:);

prop = astro.cr3bp.propagateWithSTM( ...
    x0, [0 T], mu, odeset('RelTol',1e-12,'AbsTol',1e-12));

M = prop.PhiFinal;
end