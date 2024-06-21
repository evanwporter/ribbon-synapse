%% Define Time & Voltage
t = linspace(0, .1, 1000); % Time from 0 to 1 second
dt = t(2) - t(1); % Time step
IHCVoltage = -70 * ones(size(t)); % Constant voltage or replace with actual data
n = 1;

V_steady_state = IHCVoltage(1,:)';


% Auditory nerve fibers
fiber_properties = fiber_properties_0522();

fiber_properties = [ ...
    % fiber_properties.test_1
    fiber_properties.test_6_n4
];

antopt = antOpt( ...
    'script', 'ANT', ...
    'ant', 'v4', ...
    'fiber', {{'regular_HSR_v1'}});


tropt = transductionOpt_v4_1(antopt.ant, fiber_properties{1});
x_pos = 0.5;


%% https://github.com/evanwporter/cochlea-nerve/blob/cc845a8870e4825796b05a13568a15e4361ce6cf/IHC/Transduction_v4.m#L64-L67
channels_open_ss_parameters_normal = tropt.generate_channels_open_ss_parameters('normal', x_pos);
channels_open_ss_parameters_burst = tropt.generate_channels_open_ss_parameters('burst', x_pos);

transmitter_release_parameters = tropt.generate_transmitter_release_parameters(x_pos);


%% https://github.com/evanwporter/cochlea-nerve/blob/cc845a8870e4825796b05a13568a15e4361ce6cf/IHC/Transduction_v4.m#L75-L81
C_Ca_background_base = 17.1e-6; % M
C_Ca_background_apex = 40.6e-6; % M

% interpolated
C_Ca_background = C_Ca_background_base + (C_Ca_background_apex - C_Ca_background_base) * x_pos;

C_Ca_background = C_Ca_background * tropt.C_Ca_background_factor;


%% https://github.com/evanwporter/cochlea-nerve/blob/cc845a8870e4825796b05a13568a15e4361ce6cf/IHC/Transduction_v4.m#L116-L121
num_CaV13 = tropt.num_CaV13;
        
tau_CaV13 = tropt.tau_CaV13.s;
tau_CaV13_blocked =  tropt.tau_CaV13_blocked.s;

num_vesicles = tropt.num_release_sites - tropt.num_inactive_release_sites;


% https://github.com/evanwporter/cochlea-nerve/blob/cc845a8870e4825796b05a13568a15e4361ce6cf/IHC/Transduction_v4.m#L125-L129
rs = RibbonSynapse_v4(...
    'num_release_sites', tropt.num_release_sites, ...
    'num_channels', tropt.num_CaV13, ...
    'distance_vesicle_membrane', tropt.ribbon_synapse_properties.distance_vesicle_membrane, ...
    'intervesicle_distance', tropt.ribbon_synapse_properties.intervesicle_distance, ...
    'channel_radius', tropt.ribbon_synapse_properties.channel_radius, ...
    'vesicle_radius', tropt.ribbon_synapse_properties.vesicle_radius, ...
    'channel_distribution_method', tropt.ribbon_synapse_properties.channel_distribution_method, ...
    'channel_distribution_parameters', tropt.ribbon_synapse_properties.channel_distribution_parameters, ...
    'plotflag', false ...
);


% https://github.com/evanwporter/cochlea-nerve/blob/cc845a8870e4825796b05a13568a15e4361ce6cf/IHC/Transduction_v4.m#L147-L148
channels = Channels(num_CaV13, tau_CaV13, tau_CaV13_blocked);



% https://github.com/evanwporter/cochlea-nerve/blob/cc845a8870e4825796b05a13568a15e4361ce6cf/IHC/Transduction_v4.m#L169-L175
V0t = channels_open_ss_parameters_normal{1};
S0t = channels_open_ss_parameters_normal{2};

channels.alpha = tropt.channels_markov_properties.alpha;
channels.kp0 = tropt.channels_markov_properties.kp0;
channels.V0t = V0t.V;
channels.S0t = S0t.V;



% https://github.com/evanwporter/cochlea-nerve/blob/cc845a8870e4825796b05a13568a15e4361ce6cf/IHC/Transduction_v4.m#L177-L178
vesicles = Vesicles(num_vesicles);


%% Size Info
% https://github.com/evanwporter/cochlea-nerve/blob/4e8b4f18f20782bfc39d88589db9f4e04dbcf507/IHC/Transduction_v4.m#L317-L328
size_info = create_size_info( ...
    "CaV13_channels_fraction", 1, ...
    "CaV13_channels_blocked", 1, ...
    "Ca_current", 1, ...
    "Ca_concentration", num_vesicles, ...
    "NT_free", num_vesicles, ...
    "NT_cleft", 1, ...
    "NT_reprocessing", 1, ...
    "proton_cleft", 1, ...
    "CaV13_num_inactivated", 1, ...
    "CaV13_num_normal", 1, ...
    "CaV13_num_burst", 1);


%% Initial conditions
m0 = sum([channels.state] == 'o') / channels.num;
C0 = repmat(C_Ca_background, 1, num_vesicles); % Row vector
q0 = ones(n, num_vesicles); % All NT in free pool
w0 = zeros(n, 1); % Empty
I0 = 0; % Initial current
c_proton0 = zeros(n, 1); % Proton cleft

y0 = [m0, sum([channels.state] == 'b') / channels.num, I0, C0, q0, 0, w0, c_proton0]';


% Solver options
solveropt = solverOpt('TimeStep', dt);
tspan = [t(1), t(end)];

% https://github.com/evanwporter/cochlea-nerve/blob/cc845a8870e4825796b05a13568a15e4361ce6cf/IHC/Transduction_v4.m#L108-L109
[numSteps, numSamples] = odeEuler_tspan(tspan, dt);


%% https://github.com/evanwporter/cochlea-nerve/blob/cc845a8870e4825796b05a13568a15e4361ce6cf/IHC/Transduction_v4.m#L194-L207
d = 3; %           ...  geometry factor
D = tropt.Ca_diffusion_coefficient;
rel_tol = tropt.Ca_conc_rel_tol;

for i = 1:num_CaV13
    
    rho = rs.rho(:,i) * 1e-9; % channel mouth--vesicle membrane
    psi = rs.psi(:,i) * 1e-9; % channel mouth--channel mouth

    psi_nernst = sqrt(psi.^2 + tropt.d_nernst.^2); % channel mouth--channel nernst point (above channel mouth)
    
    r = [psi_nernst; rho];

    ps(i) = Diffusion.PointSource(d, D, r, dt, 0:numSamples, 'rel_tol', rel_tol);
end


% https://github.com/evanwporter/cochlea-nerve/blob/cc845a8870e4825796b05a13568a15e4361ce6cf/IHC/Transduction_v4.m#L343-L349
if isscalar(tropt.x), tropt.x = repmat(tropt.x,[n,1]); end
if isscalar(tropt.y), tropt.y = repmat(tropt.y,[n,1]); end
if isscalar(tropt.r), tropt.r = repmat(tropt.r,[n,1]); end
if isscalar(tropt.l), tropt.l = repmat(tropt.l,[n,1]); end

uconv = @(var) Unit.batch_convert(simulation_units, var);

%% Differential
[y_out, t_out] = odeEuler(@TransductionRHS_v5, tspan, y0, solveropt, ...
    size_info, channels, vesicles, ps, tropt.y.Hz, tropt.l.Hz, tropt.x.Hz, tropt.r.Hz, ...
    tropt.G_Ca, C_Ca_background, ...
    channels_open_ss_parameters_normal, ...
    channels_open_ss_parameters_burst, ...
    transmitter_release_parameters, ...
    -70, dt, IHCVoltage);

% Voltage output
V_out = interp1(t, IHCVoltage, t_out); % Voltage interpolated at output times

% Dynamical output placeholder
y_dyn_out = []; % Modify this if additional dynamic simulations are required
