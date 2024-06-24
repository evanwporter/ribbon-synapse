% test.m

% https://github.com/evanwporter/cochlea-nerve/blob/4e8b4f18f20782bfc39d88589db9f4e04dbcf507/Wrapper/Options/transductionOpt_v4_1.m#L354-L355
num_release_sites = 14;
num_CaV13 = 84; %72;        


rs = RibbonSynapse_v4( ...
    'plotflag', false, ...
    'num_release_sites', num_release_sites, ...
    'num_channels',  num_CaV13, ...
    'channel_distribution_method', 'regular_band', ...
    'channel_distribution_parameters', struct('width', 100, 'length', 400), ...
    'plot_histograms', false ...
);

t = 0;      % Initial time
dt = .5;  % Time step


%%
% num_vesicles = rs_regular.num_release_sites;

% https://github.com/evanwporter/cochlea-nerve/blob/cc845a8870e4825796b05a13568a15e4361ce6cf/IHC/Transduction_v4.m#L121-L122
% num_vesicles = tropt.num_release_sites - tropt.num_inactive_release_sites;


% https://github.com/evanwporter/cochlea-nerve/blob/4e8b4f18f20782bfc39d88589db9f4e04dbcf507/Wrapper/Options/transductionOpt_v4_1.m#L335-L336
tau_CaV13 = Time(500, 'us'); % see Zampini 2010 Table 1
tau_CaV13_blocked = Time(1, 'ms');

% https://github.com/evanwporter/cochlea-nerve/blob/cc845a8870e4825796b05a13568a15e4361ce6cf/IHC/Transduction_v4.m#L147-L148
channels = Channels(num_CaV13, tau_CaV13.s, tau_CaV13_blocked.s)


%%
Vt = Voltage(-70, 'mV'); % Initial voltage
V_steady_state = Voltage(-70, 'mV').V; % Steady-state voltage



%%
% https://github.com/evanwporter/cochlea-nerve/blob/4e8b4f18f20782bfc39d88589db9f4e04dbcf507/IHC/Transduction_v4.m#L268-L284
% m0: Percent of open calcium channels
%   # of open states / by the total # of channels.
% C0: [Background calcium concentration] * num vesicles
% q0: # ready-to-release neurotransmitter in vesicles, set to ones
% w0: # neurotransmitters in the process of being reprocessed or recycled

% m0 = sum([channels.state] == 'o') / channels.num;

% C0 = C_Ca_background;
% C0 = repmat(C0, 1, num_vesicles); % row

% q0 = ones(n, num_vesicles);  % all NT in free pool

% w0 = zeros(n,1);           % empty



%%
% y0 = [q0, c_glut0, w0, c_prot0]';



%%
% https://github.com/evanwporter/cochlea-nerve/blob/4e8b4f18f20782bfc39d88589db9f4e04dbcf507/IHC/Transduction_v4.m#L514-L515
vesicles = Vesicles(num_vesicles);



%%
% https://github.com/evanwporter/cochlea-nerve/blob/4e8b4f18f20782bfc39d88589db9f4e04dbcf507/Wrapper/Options/transductionOpt_v4_1.m#L341-L350
transmitter_release_parameters = {500, 5.0, 100e-6};



%%
% https://github.com/evanwporter/cochlea-nerve/blob/4e8b4f18f20782bfc39d88589db9f4e04dbcf507/Wrapper/Options/transductionOpt_v4_1.m#L38-L50
% manufacture rate
y = Frequency(10, 'Hz'); % Sumner 2002

% loss rate
l = Frequency(1290, 'Hz'); % to roughly match the time-constant of H&H
% l = Frequency(2580, 'Hz'); % Sumner 2002

% reprocessing rate
x = Frequency(66.3, 'Hz'); % Sumner 2002

% re-uptake rate
r = Frequency(3290, 'Hz'); % to roughly match the time-constant of H&H
% r = Frequency(6580, 'Hz'); % Sumner 2002



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

% https://github.com/evanwporter/cochlea-nerve/blob/4e8b4f18f20782bfc39d88589db9f4e04dbcf507/IHC/Transduction_v4.m#L317-L328
size_info_2 = create_size_info( ...
    "NT_free", num_vesicles, ...
    "NT_cleft", 1, ...
    "NT_reprocessing", 1, ...
    "proton_cleft", 1 ...
);


%% 
C_vesicles = rand(num_vesicles, 1) * 1e-3; % Initial calcium concentration (in M)
q = ones(num_vesicles, 1); % Vesicle has transmitter (1) or not (0)
c = 0; % Initial cleft concentration
w = 0; % Initial reprocessing concentration 
c_proton = 0; % Initial proton concentration

% Need to convert x, y, l, r to Hz

% Initialize m (fraction of open Ca2+ channels)
% Let's assume that a certain percentage of channels are open initially.
initial_open_fraction = 0.2; % Initial fraction of open channels (20% as an example)
channels.state = repmat('c', channels.num, 1); % All channels start closed
num_open_channels = round(channels.num * initial_open_fraction);
channels.state(1:num_open_channels) = 'o'; % Open some channels
m = sum(channels.state == 'o') / channels.num; % Calculate initial m

% Initial CaV13_channels_blocked (Ca_blocked)
Ca_blocked = sum(channels.state == 'b') / channels.num;

% Initial Ca_current (I)
I = 0; % No initial current

z = [m, Ca_blocked, I, C_vesicles, q, c, w, c_proton]'; % Initial state vector



%%
% Calculate the neurotransmitter release rate
[dq, dc, dw, dc_proton, vesicles] = NTdynamicsRHS_v5_core( t, q, ...
    c, w, c_proton, ...
    vesicles, y.Hz, l.Hz, x.Hz, r.Hz, ...
    transmitter_release_parameters, ...
    dt, C_vesicles);


dz = TransductionRHS_v5(t, z, size_info, channels, vesicles, ps, rate_y, rate_l, rate_x, rate_r, G_Ca, C_Ca_background, ...
    rs_regular.channels_open_ss_parameters_normal, rs_regular.channels_open_ss_parameters_burst, ...
    transmitter_release_parameters, V_steady_state, dt, Vt);



% Output the rate of change of neurotransmitter release
fprintf('Rate of change of neurotransmitter release (dq):\n');
disp(dq);a