% Plot the synapse morphology - position of vesicles and channels.
% This will also estimate the channel location probability density
% as a heatmap - runs for several minutes.

% Auditory nerve fibers
fiber_properties = fiber_properties_0522();

fiber_properties = [ ...
    % fiber_properties.test_1
    fiber_properties.test_6_n4
];

%%

[args, opt, memopt, paropt] = common_opts(struct(), ...
    do_not_change_settings=true, ...
    silent=true);

antopt = antOpt( ...
    'script', 'ANT', ...
    'ant', 'v4', ...
    'fiber', {{'regular_HSR_v1'}});

% disp(antopt)

for i = 1:numel(fiber_properties)
    antopt.transductionopt = transductionOpt_v4_1(antopt.ant, fiber_properties{i});
    antopt.transductionopt.visualize_density()
end

%%

% Evan's Notes
% Gray Spheres: Vesicles not currently releasing neurotransmitters.
% Blue Spheres: Vesicles releasing neurotransmitters 

% Red things: Calcium channels