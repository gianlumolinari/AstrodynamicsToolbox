function state = getHorizonsState(command, epochUTC, center)
%GETHORIZONSSTATE Query JPL Horizons API for one Cartesian state.
%
% INPUTS
%   command  : Horizons target identifier as char/string
%   epochUTC : UTC epoch string, e.g. '2026-11-12 00:00:00'
%   center   : Horizons center code, e.g. '500@10'
%
% OUTPUT
%   state : struct with fields
%       .target
%       .center
%       .epochUTC
%       .r
%       .v
%       .rawText

if nargin < 3
    center = '500@10';
end

command  = char(string(command));
epochUTC = char(string(epochUTC));
center   = char(string(center));

baseURL = 'https://ssd.jpl.nasa.gov/api/horizons.api';

% Use TLIST for a single discrete epoch
query = [ ...
    '?format=text' ...
    '&COMMAND='     urlencode(['''' command '''']) ...
    '&OBJ_DATA='    urlencode('''YES''') ...
    '&MAKE_EPHEM='  urlencode('''YES''') ...
    '&EPHEM_TYPE='  urlencode('''VECTORS''') ...
    '&CENTER='      urlencode(['''' center '''']) ...
    '&TLIST='       urlencode(['''' epochUTC '''']) ...
    '&TLIST_TYPE='  urlencode('''CAL''') ...
    '&REF_PLANE='   urlencode('''ECLIPTIC''') ...
    '&REF_SYSTEM='  urlencode('''J2000''') ...
    '&VEC_TABLE='   urlencode('''2''') ...
    ];

fullURL = [baseURL query];
txt = webread(fullURL);

soeIdx = strfind(txt, '$$SOE');
eoeIdx = strfind(txt, '$$EOE');

if isempty(soeIdx) || isempty(eoeIdx)
    error('Could not find $$SOE/$$EOE in Horizons response.\nResponse was:\n%s', txt(1:min(end,3000)));
end

block = txt(soeIdx(1):eoeIdx(1));
lines = regexp(block, '\r\n|\n|\r', 'split');

xyzLine = '';
velLine = '';

for k = 1:numel(lines)
    L = strtrim(lines{k});
    if contains(L, 'X =') && contains(L, 'Y =') && contains(L, 'Z =')
        xyzLine = L;
    end
    if contains(L, 'VX=') && contains(L, 'VY=') && contains(L, 'VZ=')
        velLine = L;
    end
end

if isempty(xyzLine) || isempty(velLine)
    error('Could not parse state lines from Horizons response.\nBlock was:\n%s', block);
end

r = parseXYZLine(xyzLine);
v = parseVLine(velLine);

state.target   = command;
state.center   = center;
state.epochUTC = epochUTC;
state.r        = r;
state.v        = v;
state.rawText  = txt;
end

function r = parseXYZLine(line)
tok = regexp(line, ...
    'X\s*=\s*([-\+\d\.Ee]+)\s*Y\s*=\s*([-\+\d\.Ee]+)\s*Z\s*=\s*([-\+\d\.Ee]+)', ...
    'tokens', 'once');

if isempty(tok)
    error('Failed to parse XYZ line.');
end

r = [str2double(tok{1}); str2double(tok{2}); str2double(tok{3})];
end

function v = parseVLine(line)
tok = regexp(line, ...
    'VX\s*=\s*([-\+\d\.Ee]+)\s*VY\s*=\s*([-\+\d\.Ee]+)\s*VZ\s*=\s*([-\+\d\.Ee]+)', ...
    'tokens', 'once');

if isempty(tok)
    error('Failed to parse velocity line.');
end

v = [str2double(tok{1}); str2double(tok{2}); str2double(tok{3})];
end