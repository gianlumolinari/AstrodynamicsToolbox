function L = lagrangePoints(mu)
%LAGRANGEPOINTS Computes the five Lagrange points in the nondimensional CR3BP.
%
% INPUT
%   mu : mass parameter
%
% OUTPUT
%   L : struct with fields L1, L2, L3, L4, L5 (each 3x1)

opts = optimset('TolX',1e-14,'Display','off');

f = @(x) x ...
    - (1-mu)*(x + mu)/abs(x + mu)^3 ...
    - mu*(x - 1 + mu)/abs(x - 1 + mu)^3;

xL1 = fzero(f, [0.5, 1-mu-1e-6], opts);
xL2 = fzero(f, [1-mu+1e-6, 1.5], opts);
xL3 = fzero(f, [-1.5, -mu-1e-6], opts);

L.L1 = [xL1; 0; 0];
L.L2 = [xL2; 0; 0];
L.L3 = [xL3; 0; 0];
L.L4 = [0.5 - mu;  sqrt(3)/2; 0];
L.L5 = [0.5 - mu; -sqrt(3)/2; 0];
end