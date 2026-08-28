function c = colors()
%COLORS  Okabe-Ito palette, keyed by role rather than by index.
%
%   Colour-blind safe and legible in greyscale print. Referring to
%   c.conformal rather than to a row number means a reader of the figure code
%   can tell which series is which without cross-referencing.

    c.gp         = [0.00 0.45 0.70];   % blue    - Bayesian baseline
    c.conformal  = [0.90 0.62 0.00];   % orange  - split conformal
    c.normalized = [0.00 0.62 0.45];   % teal    - normalized conformal
    c.krr        = [0.80 0.47 0.65];   % pink    - kernel ridge point estimate
    c.ideal      = [0.30 0.30 0.30];   % grey    - the diagonal / target
    c.weighted   = [0.60 0.30 0.70];   % purple  - weighted conformal (new)
    c.warn       = [0.80 0.25 0.20];   % red     - violated guarantee
end
