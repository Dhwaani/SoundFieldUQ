function results = run_tests()
%RUN_TESTS  Run the whole test suite.
%
%   results = RUN_TESTS() executes every test class in tests/ and prints a
%   summary. Run this before trusting any number the experiments produce, and
%   again before pushing.
%
%   The suite needs no data, no downloads and no toolboxes beyond base MATLAB.
%   Six test classes; the two that matter most are TCONFORMAL (the unweighted
%   guarantee and the Rayleigh factor) and TWEIGHTEDCONFORMAL (the new weighted
%   guarantee and the abstention behaviour). Runs in well under a minute.

    root = fileparts(mfilename('fullpath'));
    addpath(root);
    addpath(fullfile(root, 'tests'));

    import matlab.unittest.TestSuite
    import matlab.unittest.TestRunner
    import matlab.unittest.plugins.TestRunProgressPlugin

    suite  = TestSuite.fromFolder(fullfile(root, 'tests'));
    runner = TestRunner.withTextOutput('Verbosity', 1);
    runner.addPlugin(TestRunProgressPlugin.withVerbosity(2));

    results = runner.run(suite);

    fprintf('\n----------------------------------------------\n');
    fprintf('  %d passed, %d failed, %d incomplete  (%.1f s)\n', ...
            nnz([results.Passed]), nnz([results.Failed]), ...
            nnz([results.Incomplete]), sum([results.Duration]));
    fprintf('----------------------------------------------\n');

    if any([results.Failed])
        fprintf(2, '\nFailed tests:\n');
        failed = results([results.Failed]);
        for i = 1:numel(failed)
            fprintf(2, '  %s\n', failed(i).Name);
        end
    end
end
