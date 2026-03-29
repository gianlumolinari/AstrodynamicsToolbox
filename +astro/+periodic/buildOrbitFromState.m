function orbit = buildOrbitFromState(state0, T, mu, meta, opts)
%BUILDORBITFROMSTATE Backward-compatible wrapper for buildOrbitStruct.

if nargin < 4 || isempty(meta)
    meta = struct();
end
if nargin < 5 || isempty(opts)
    opts = struct();
end

orbit = astro.periodic.buildOrbitStruct(state0, T, mu, meta, opts);
end