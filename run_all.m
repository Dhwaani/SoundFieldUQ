function run_all(quick)
%RUN_ALL  Reproduce every experiment, table and figure in the paper.
%
%   RUN_ALL()      full run, roughly 15-30 minutes on a laptop
%   RUN_ALL(true)  quick run with fewer trials, roughly 3 minutes
%
%   One entry point on purpose: "these figures came from this code" should be
%   checkable by running one command, not by following instructions in the
%   right order.
%
%   Everything is synthetic. No dataset is downloaded, no hardware is used, no
%   measurement is required. Re-running on any machine with base MATLAB
%   reproduces the same numbers.
%
%   Experiments 01-05 establish the setting and audit the existing methods.
%   Experiments 06-08 are the contribution.

    if nargin < 1 || isempty(quick), quick = false; end

    root = fileparts(mfilename('fullpath'));
    run(fullfile(root, 'startup.m'));

    if quick
        n = struct('e1',3,'e2',4,'e3',4,'e4',4,'e5',3,'e6',4,'e7',3,'e8',4);
        fprintf('\n*** QUICK MODE: fewer trials, noisier numbers ***\n');
    else
        n = struct('e1',8,'e2',10,'e3',12,'e4',10,'e5',8,'e6',10,'e7',8,'e8',10);
    end

    t0 = tic;

    fprintf('\n--- Part I: the setting and the audit ---\n');
    exp01_accuracy(n.e1);
    exp02_marginal_calibration(n.e2);
    R3 = exp03_conditional_coverage(n.e3);
    exp04_aperture_shift(n.e4);
    exp05_frequency_sweep(n.e5);

    fprintf('\n--- Part II: geometry-weighted conformal (the contribution) ---\n');
    R6 = exp06_weighted_conformal(n.e6);
    R7 = exp07_support_failure(n.e7);
    R8 = exp08_placement(n.e8);

    fprintf('\n==============================================\n');
    fprintf('  Completed in %.1f minutes\n', toc(t0) / 60);
    fprintf('  Figures : %s\n', fullfile(root, 'figures'));
    fprintf('  Results : %s\n', fullfile(root, 'results'));
    fprintf('==============================================\n');

    fprintf('\nSUMMARY OF FINDINGS\n\n');

    fprintf('1. Split conformal is conditionally invalid (exp 03).\n');
    fprintf('   Coverage varies by %.3f across space; the weighted-normalized\n', ...
            max(R3.covCP) - min(R3.covCP));
    fprintf('   variant varies by only %.3f.\n\n', ...
            max(R3.covNCP) - min(R3.covNCP));

    fprintf('2. Under a DESIGNED density shift the exact geometric likelihood\n');
    fprintf('   ratio restores validity (exp 06).\n');
    fprintf('   Worst split-conformal coverage %.3f, worst weighted %.3f,\n', ...
            min(R6.covCP), min(R6.covWCP));
    fprintf('   against a nominal %.2f. Cost: ESS falls to %.2f of n.\n\n', ...
            1 - R6.alpha, min(R6.essFrac));

    fprintf('3. Beyond the calibration support NOTHING restores validity (exp 07).\n');
    fprintf('   Split conformal answers confidently and is wrong: %.3f coverage.\n', ...
            R7.covCP(end));
    fprintf('   Weighted conformal abstains on %.0f%% of points instead, and is\n', ...
            100 * R7.abstention(end));
    fprintf('   correct on the rest. Abstention beats a confident wrong answer.\n\n');

    fprintf('4. Placement follows (exp 08). All layouts are VALID; they differ\n');
    fprintf('   in cost. Compact vs matched: ESS %.2f vs %.2f, width %.5f vs %.5f.\n', ...
            R8.essFrac(1), R8.essFrac(3), R8.widthWCP(1), R8.widthWCP(3));
    fprintf('   Put the calibration microphones where the listeners are.\n');
end
