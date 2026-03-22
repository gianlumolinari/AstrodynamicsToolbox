function c3 = C3(vInf)
%C3 Computes launch characteristic energy C3.
%
% INPUT
%   vInf : hyperbolic excess speed or vector [km/s]
%
% OUTPUT
%   c3   : characteristic energy [km^2/s^2]

if isvector(vInf) && numel(vInf) == 3
    c3 = dot(vInf, vInf);
else
    c3 = vInf.^2;
end
end