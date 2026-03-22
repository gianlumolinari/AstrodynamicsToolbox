function body = getBody(name)
%GETBODY Returns a body struct by name.
%
% INPUT
%   name : string or char
%
% OUTPUT
%   body : struct returned by the corresponding body function
%
% Supported names:
%   'sun', 'mercury', 'venus', 'earth', 'moon',
%   'mars', 'jupiter', 'saturn', 'uranus', 'neptune'

name = lower(string(name));

switch name
    case "sun"
        body = astro.bodies.sun();
    case "mercury"
        body = astro.bodies.mercury();
    case "venus"
        body = astro.bodies.venus();
    case "earth"
        body = astro.bodies.earth();
    case "moon"
        body = astro.bodies.moon();
    case "mars"
        body = astro.bodies.mars();
    case "jupiter"
        body = astro.bodies.jupiter();
    case "saturn"
        body = astro.bodies.saturn();
    case "uranus"
        body = astro.bodies.uranus();
    case "neptune"
        body = astro.bodies.neptune();
    otherwise
        error('Unsupported body name: %s', name);
end
end