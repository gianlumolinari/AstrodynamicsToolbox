function results = screenNEOTargets(targetList)
%SCREENNEOTARGETS Query and summarize a list of candidate NEAs.

if nargin < 1 || isempty(targetList)
    error('targetList must be a non-empty cell array of names/designations.');
end

template = struct( ...
    'query',    '', ...
    'fullname', '', ...
    'spkid',    '', ...
    'H',        NaN, ...
    'pha',      false, ...
    'neo',      false, ...
    'a_au',     NaN, ...
    'e',        NaN, ...
    'i_deg',    NaN, ...
    'q_au',     NaN, ...
    'Q_au',     NaN, ...
    'moid_au',  NaN, ...
    'period_d', NaN);

results = repmat(template, numel(targetList), 1);

for k = 1:numel(targetList)
    target = targetList{k};
    entry = template;
    entry.query = char(string(target));

    try
        sb = astro.smallbody.querySBDB(target);

        if isfield(sb, 'object') && isfield(sb.object, 'fullname')
            entry.fullname = char(string(sb.object.fullname));
        end

        if isfield(sb, 'object') && isfield(sb.object, 'spkid')
            entry.spkid = char(string(sb.object.spkid));
        end

        if isfield(sb, 'phys_par') && isfield(sb.phys_par, 'H')
            entry.H = localToDouble(sb.phys_par.H);
        end

        if isfield(sb, 'object') && isfield(sb.object, 'pha')
            entry.pha = localToLogical(sb.object.pha);
        end

        if isfield(sb, 'object') && isfield(sb.object, 'neo')
            entry.neo = localToLogical(sb.object.neo);
        end

        elems = [];
        if isfield(sb, 'orbit') && isfield(sb.orbit, 'elements')
            elems = sb.orbit.elements;
        end

        entry.a_au     = localReadElement(elems, 'a');
        entry.e        = localReadElement(elems, 'e');
        entry.i_deg    = localReadElement(elems, 'i');
        entry.q_au     = localReadElement(elems, 'q');
        entry.Q_au     = localReadElement(elems, 'ad');
        entry.period_d = localReadElement(elems, 'per');

        if isfield(sb, 'orbit') && isfield(sb.orbit, 'moid')
            entry.moid_au = localToDouble(sb.orbit.moid);
        end

    catch
        % Keep defaults if query/parsing fails
    end

    results(k) = entry;
end

end

function val = localReadElement(elems, name)
val = NaN;

if isempty(elems)
    return
end

if iscell(elems)
    for ii = 1:numel(elems)
        el = elems{ii};
        if isstruct(el) && isfield(el, 'name') && strcmpi(el.name, name)
            if isfield(el, 'value')
                val = localToDouble(el.value);
                return
            end
        end
    end
elseif isstruct(elems)
    for ii = 1:numel(elems)
        el = elems(ii);
        if isfield(el, 'name') && strcmpi(el.name, name)
            if isfield(el, 'value')
                val = localToDouble(el.value);
                return
            end
        end
    end
end
end

function x = localToDouble(v)
if isnumeric(v)
    x = double(v);
elseif isstring(v) || ischar(v)
    x = str2double(v);
else
    x = NaN;
end
end

function tf = localToLogical(v)
if islogical(v)
    tf = v;
elseif isnumeric(v)
    tf = logical(v);
elseif ischar(v) || isstring(v)
    s = lower(char(string(v)));
    tf = any(strcmp(s, {'y','yes','true','1'}));
else
    tf = false;
end
end