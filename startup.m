function startup()
%STARTUP  Put the project on the MATLAB path.
%
%   Run this once per session, from anywhere:
%
%       run('/path/to/SoundFieldUQ/startup.m')
%
%   The +sfuq package folder is found automatically as long as the repository
%   ROOT is on the path -- MATLAB resolves package folders relative to path
%   entries, so the +sfuq directory itself must NOT be added.

    root = fileparts(mfilename('fullpath'));

    addpath(root);                              % makes +sfuq visible
    addpath(fullfile(root, 'experiments'));
    addpath(fullfile(root, 'tests'));

    for d = {'results', 'figures'}
        p = fullfile(root, d{1});
        if ~isfolder(p)
            mkdir(p);
        end
    end

    fprintf('SoundFieldUQ ready.  MATLAB %s\n', version('-release'));
    fprintf('  run_all    -- every experiment and figure\n');
    fprintf('  run_tests  -- the test suite\n');
end
