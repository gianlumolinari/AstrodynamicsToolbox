function atlasOut = sortTubeStrands(atlas, mode)
%SORTTUBESTRANDS Sort manifold strands for cleaner plotting.

    if nargin < 2 || isempty(mode)
        mode = 'signphase';
    end

    branches = atlas.branches;

    switch lower(mode)
        case 'phase'
            key = [[branches.phaseIndex].'];
            [~,I] = sort(key);

        case 'signphase'
            key = [[branches.sign].', [branches.phaseIndex].'];
            [~,I] = sortrows(key,[1 2]);

        otherwise
            error('mode must be ''phase'' or ''signphase''.');
    end

    atlasOut = atlas;
    atlasOut.branches = branches(I);

    if isfield(atlas,'seeds') && numel(atlas.seeds) == numel(branches)
        atlasOut.seeds = atlas.seeds(I);
    end
end