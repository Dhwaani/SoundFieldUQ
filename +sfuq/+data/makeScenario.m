function S = makeScenario(varargin)
%MAKESCENARIO  Build a reproducible sound-field reconstruction scenario.
%
%   S = MAKESCENARIO() returns the default scenario used throughout the
%   experiments. Name-value pairs override any field.
%
%   Name-value options
%     'Room'        [6 5 3]        room dimensions, metres
%     'Source'      [1.2 3.4 1.5]  source position
%     'Beta'        0.7            wall reflection coefficient
%     'Order'       3              image-source order
%     'Centre'      [3 2.5 1.5]    centre of the measurement region
%     'HalfWidth'   0.5            half-width of the measurement aperture
%     'Frequency'   800            Hz
%     'SoundSpeed'  343            m/s
%     'NoiseRel'    0.02           measurement noise, relative to field RMS
%     'Seed'        0              base seed
%
%   The returned struct carries everything needed to regenerate the data, so
%   a results file plus its scenario struct is a complete description of an
%   experiment. Every experiment saves its scenario alongside its numbers.

    p = inputParser;
    p.addParameter('Room',       [6 5 3]);
    p.addParameter('Source',     [1.2 3.4 1.5]);
    p.addParameter('Beta',       0.7);
    p.addParameter('Order',      3);
    p.addParameter('Centre',     [3 2.5 1.5]);
    p.addParameter('HalfWidth',  0.5);
    p.addParameter('Frequency',  800);
    p.addParameter('SoundSpeed', 343);
    p.addParameter('NoiseRel',   0.02);
    p.addParameter('Seed',       0);
    p.parse(varargin{:});

    S = p.Results;
    S.k = 2 * pi * S.Frequency / S.SoundSpeed;
    [S.imgPos, S.imgAmp] = sfuq.data.imageSources(S.Source, S.Room, ...
                                                  S.Beta, S.Order);

    % Schroeder frequency and the half-wavelength spacing rule are printed by
    % the experiments as a sanity check: reconstructing above the spatial
    % Nyquist limit of the microphone layout is expected to fail, and knowing
    % where that limit falls stops a genuine failure being mistaken for a bug.
    S.wavelength = S.SoundSpeed / S.Frequency;
    S.nyquistSpacing = S.wavelength / 2;
end
