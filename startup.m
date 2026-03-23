root = fileparts(mfilename('fullpath'));
addpath(root);

% Add MICE to path
micePath = fullfile(root, 'external', 'mice', 'src', 'mice');
if exist(micePath, 'dir')
    addpath(micePath);
end

disp('astroToolbox path added.');