function state = getState(backend, target, epochUTC, varargin)
%GETSTATE Unified ephemeris interface.
%
% USAGE
%   state = astro.ephem.getState('spice', target, epochUTC)
%   state = astro.ephem.getState('horizons', target, epochUTC)
%
% OPTIONAL INPUTS
%   For SPICE:
%       observer, frame, abcorr
%   For Horizons:
%       center
%
% NOTES
%   - target may be simple names like 'earth', 'mars', etc.
%   - for SPICE, common planet names are mapped inside getSpiceState
%   - for Horizons, numeric IDs are used internally

backend = lower(string(backend));
target = lower(string(target));
epochUTC = char(string(epochUTC));

switch backend
    case "spice"
        observer = 'SUN';
        frame = 'J2000';
        abcorr = 'NONE';

        if numel(varargin) >= 1 && ~isempty(varargin{1}), observer = varargin{1}; end
        if numel(varargin) >= 2 && ~isempty(varargin{2}), frame    = varargin{2}; end
        if numel(varargin) >= 3 && ~isempty(varargin{3}), abcorr   = varargin{3}; end

        spiceTarget = localTargetToSpice(target);
        state = astro.ephem.getSpiceState(spiceTarget, epochUTC, observer, frame, abcorr);

    case "horizons"
        center = '500@10';
        if numel(varargin) >= 1 && ~isempty(varargin{1}), center = varargin{1}; end

        horizonsID = localTargetToHorizons(target);
        state = astro.ephem.getHorizonsState(horizonsID, epochUTC, center);

    otherwise
        error('Unsupported ephemeris backend: %s', backend);
end

end

function spiceTarget = localTargetToSpice(target)
switch target
    case "sun"
        spiceTarget = 'SUN';
    case "mercury"
        spiceTarget = 'MERCURY';
    case "venus"
        spiceTarget = 'VENUS';
    case "earth"
        spiceTarget = 'EARTH';
    case "moon"
        spiceTarget = 'MOON';
    case "mars"
        spiceTarget = 'MARS';
    case "jupiter"
        spiceTarget = 'JUPITER';
    case "saturn"
        spiceTarget = 'SATURN';
    case "uranus"
        spiceTarget = 'URANUS';
    case "neptune"
        spiceTarget = 'NEPTUNE';
    otherwise
        error('Unsupported target for SPICE backend: %s', target);
end
end

function horizonsID = localTargetToHorizons(target)
switch target
    case "sun"
        horizonsID = '10';
    case "mercury"
        horizonsID = '199';
    case "venus"
        horizonsID = '299';
    case "earth"
        horizonsID = '399';
    case "moon"
        horizonsID = '301';
    case "mars"
        horizonsID = '499';
    case "jupiter"
        horizonsID = '599';
    case "saturn"
        horizonsID = '699';
    case "uranus"
        horizonsID = '799';
    case "neptune"
        horizonsID = '899';
    otherwise
        error('Unsupported target for Horizons backend: %s', target);
end
end