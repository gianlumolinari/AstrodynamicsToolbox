function startup
    root = fileparts(mfilename('fullpath'));
    addpath(root);

    % Add main toolbox packages / utilities if needed
    addpath(genpath(fullfile(root, 'dev')));

    % Add MICE wrapper functions
    miceSrcPath = fullfile(root, 'external', 'mice', 'src', 'mice');
    if exist(miceSrcPath, 'dir')
        addpath(miceSrcPath);
    end

    % Add MICE compiled MEX binary
    miceLibPath = fullfile(root, 'external', 'mice', 'lib');
    if exist(miceLibPath, 'dir')
        addpath(miceLibPath);
    end

    rehash toolboxcache

    disp('astroToolbox path added.');

    % Optional diagnostics
    if isempty(which('mice'))
        warning('MICE binary not found on path.');
    end
    if isempty(which('cspice_kclear'))
        warning('CSPICE wrappers not found on path.');
    end
end