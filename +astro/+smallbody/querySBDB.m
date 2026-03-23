function data = querySBDB(target)
%QUERYSBDB Query JPL SBDB API for a single small body.
%
% INPUT
%   target : designation, number, or name string
%
% OUTPUT
%   data : decoded JSON struct returned by SBDB API
%
% NOTES
%   Requests physical parameters too (phys-par=1), so H can be returned.

if nargin < 1 || isempty(target)
    error('A target designation/name must be provided.');
end

baseURL = 'https://ssd-api.jpl.nasa.gov/sbdb.api';
query = ['?sstr=' urlencode(char(string(target))) ...
         '&phys-par=1' ...
         '&full-prec=1'];

opts = weboptions('Timeout', 30);
data = webread([baseURL query], opts);
end