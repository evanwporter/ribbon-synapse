% Calc rate of change of z

function [ dz ] = Trans2( t, z, opts, V_steady_state, Vt)
arguments
    t
    z
    opts
    V_steady_state
    Vt
end

persistent last_t
persistent current_index

if isempty(last_t)
    last_t = 0;
    current_index = 1;
end 

% Compute timestep
dt = t - last_t;
if dt <= 0
    disp(t)
    dt = 1e-12;
end
last_t = t;


% opts.dt= dt;

Vt = Vt(:);

n = size(Vt,1);

m = decompose_z(z, 'CaV13_channels_fraction', opts.size_info);
Ca_blocked = decompose_z(z, 'CaV13_channels_blocked', opts.size_info);
C_vesicles = decompose_z(z, 'Ca_concentration', opts.size_info);
I = decompose_z(z, 'Ca_current', opts.size_info);
q = decompose_z(z, 'NT_free', opts.size_info);
c = decompose_z(z, 'NT_cleft', opts.size_info);
w = decompose_z(z, 'NT_reprocessing', opts.size_info);

c_proton = decompose_z(z, 'proton_cleft', opts.size_info);

m_old = m;
I_old = I;
Ca_blocked_old = Ca_blocked;
C_old = C_vesicles;

channels = opts.channels;
vesicles = opts.vesicles;
ps = opts.ps;


%%

[d_state_inactivated, d_state_normal, d_state_burst] = deal(zeros(1,n));

% number of Ca_V1.3 channels in the vicinity of the synapse
num_CaV13 = channels.num;


%%

p_block = CaV13.CaProtonBlock(c_proton);

assert(dt / channels.tau_blocked <= 1);


%% Calculate rates

S0t = channels.S0t;
V0t = channels.V0t;

alpha = channels.alpha;
beta = alpha + 1/S0t;

kp0 = channels.kp0;

km0 = kp0 * exp(V0t * (beta - alpha));


%% Compute conductance

kp = kp0 * exp(-alpha * Vt); % E: Opening rate
km = km0 * exp(-beta * Vt);  % E: Closing rate

Q = [1 - kp*dt, km *dt; 
     kp*dt, 1 - km*dt];


% EVAN
% Q prob of transitioning between the closed (c) and open (o) states during a time step 
%   1st row: prob of remaining closed (1 - k_p * dt) or transitioning to open (k_m \* dt).
%   2nd row: prob of transitioning to closed (k_p * dt) or remaining open (1 - k_m * dt).

cQ = cumsum(Q,1);

ST = ['c', 'o'];

%% Determine Whether Channel is Open/Closed/Blocked
%   Based on membrane potential
for ii = 1:channels.num

    if channels.state(ii) == 'i' % inactivated
        continue
    end

    % Unblock channel
    if channels.state(ii) == 'b' % blocked
        if t >= channels.tblocked(ii)
            channels.state(ii) = 'c';
        else
            continue
        end
    end

    % Block channel
    if rand(1) < p_block * dt / channels.tau_blocked

        % generate random close time from a gamma distribution
        % such that the mean (= a*b) is equal to the channels.tau_blocked
        %     a ... shape
        %     b ... scale

        a = 4;
        b = channels.tau_blocked / a;
        tau = b * gammaincinv(rand(1,1), a); % gamma random number without toolbox

        channels.state(ii) = 'b';
        channels.tblocked(ii) = t + tau;

    else
        r = rand(1);
    
        if channels.state(ii) == 'c' % closed
            st = 1;
        elseif channels.state(ii) == 'o' % open
            channels.topen(ii) = t; % last open time
            st = 2;
        else
            error('unknown state %s', channels.state(ii))
        end
        
        st_new = find(r < cQ(:,st), 1);
        try
            channels.state(ii) = ST(st_new);
        catch err
            throw(err)
        end

    end
end


m = sum(channels.state == 'o') / num_CaV13;
Ca_blocked = sum(channels.state == 'b') / num_CaV13;

%% Calculate concentration
% currently we calculate at two points:
%   1. close to the vesicle
%   2. close to the ion channel

all_channel_switch = true;
log_partial_concentrations = false;

if all_channel_switch
    ind_r_channel = 1:numel(ps);
else
    ind_r_channel = 1;
end

C_channels = zeros(numel(ps),1);
if current_index > 1 % skip first index to allow concentration to be set
    for i = 1:numel(ps)
        if all_channel_switch
            if log_partial_concentrations
                C = ps{i}.concentration(ind_r_channel, current_index - 1);
            else
                C = ps{i}.concentration(ind_r_channel);
            end
            
            % Debugging Output
            if any(isinf(C))
                warning('Infinite concentration detected at iteration %d for ps{%d}', current_index, i);
                disp('Concentration values:');
                disp(C);
            end
            
            C_channels(i) = sum(C);
            
            % Check for unrealistic values and cap them
            % if C_channels(i) > 1e3  % Example cap, adjust based on realistic limits
            %     warning('Capping concentration at iteration %d for ps{%d}, C = %d', it, i, C_channels(i));
            %     C_channels(i) = 1e3;
            % end
        else
            C_channels{i} = ps{i}.concentration(ind_r_channel, current_index - 1);
        end
    end
end

C_channels = C_channels + opts.C_Ca_background;

open_channels = channels.state == 'o';

% Evan -
% I_GHK: Used for calculating calcium current flowing into hair cell
%   Needed for determing local concentration of Calcium
%       Which in tern is needed for calculating neurotransmitter release rate

% I_ohm = calcium_current(Vt, G_Ca, open_channels, C_channels);
% I_GHK = calcium_current_GHK(-70, G_Ca, open_channels, C_channels);
I_GHK = calcium_current_GHK(Vt, opts.G_Ca, open_channels, C_channels);


I = I_GHK;

for i = 1:numel(ps)
    ps{i}.current(current_index) = - I(i); % minus because of convention
    ps{i}.dt = dt;
end

[C_vesicles, C_all] = calcium_concentration(vesicles, ps, current_index, all_channel_switch);


for i = 1:numel(ps)
    if log_partial_concentrations
        ps{i}.concentration(:, current_index) = C_all(:, i);
    else
        ps{i}.concentration = C_all(:, i);
    end
end

C_vesicles = C_vesicles + opts.C_Ca_background;


%% Transmitter Release and Recycling

[dq, dc, dw, dc_proton, vesicles] = NTdynamicsRHS_v5_core( t, ...
    q, c, w, c_proton, ...
    vesicles, opts.y.Hz, opts.l.Hz, opts.x.Hz, opts.r.Hz, ...
    opts.transmitter_release_parameters, ...
    dt, C_vesicles);


    
%% Build dz

dm = (m - m_old) / dt;

dCa_blocked = (Ca_blocked - Ca_blocked_old) / dt;

charge = 2;
amp_to_electron_per_second = 6.242e18;

I = sum(I) / amp_to_electron_per_second * charge;
dI = (I - I_old) / dt;

dC = (C_vesicles - C_old) / dt;

dz = [dm; dCa_blocked; dI; dC; dq; dc; dw; dc_proton];

d_state = [d_state_inactivated; d_state_normal; d_state_burst] / dt;

current_index = current_index + 1;

% dz = [dz; d_state];

if any(isnan(dz))
    error('NaN encountered')
end

end


% - Evan - 
% Vm: Membrane potential
% G: Channel conductance
% area: membrane area
% Used for calculating calcium current flowing into hair cell

function I = calcium_current_GHK(Vm, G, m, C)
arguments
    Vm (1,1) double
    G (1,1) double
    m (:,1) double  % num_channel x 1
    C (:,1) double  % num_channel x 1
end


% Ca2+ ion charge
charge = 2;

% extracellular concentration
C_extracellular = 1.3e-3; % M

T = 310; % [K] ~ 37C

area = 1;

me = 9.109383632e-31; % kg ... mass of electron
ce = -1.602e-19; % C ... charge of electron

fac = 6.0e-15 ; % m^3 / s^2 ... whatever this is
P = G .* m / ce^2 / charge^2 * me * fac;
P = P / area;

Cin = C*1000; % M to mol/m3
Cout = C_extracellular*1000; % M to mol/m3

Phi = GHKflux_fast(Vm, Cin, Cout, P, charge, T);

I = Phi * area;

amp_to_electron_per_second = 6.242e18;

I = I * amp_to_electron_per_second / charge; % ions/sec

% figure
% plot(t*1e3, chi.*j/1e3);
% ylabel('Ca2+ flux (ion/ms)')
% xlabel('time (ms)')

end

% vesicles: Struc w/ info about vesicles
% ps: Struct w/ channel properties and states
% it: Current iteration index
% all_channel_switch: Boolean to decide on full channel update
function [C_vesicle, C] = calcium_concentration(vesicles, ps, t, all_channel_switch)
    arguments
        vesicles
        ps
        t
        all_channel_switch
    end

    if all_channel_switch
        ind_r_channel = numel(ps);
    else
        ind_r_channel = 1;
    end
    
    N_A = 6.02214076e23;  % Avogadro constant, mol^-1
    
    num_vesicles = vesicles.num;
    num_channels = numel(ps);
    
    % Concentration at each vesicle
    C_vesicle = zeros(num_vesicles, 1);

    % Concentration Matrix
    %   C(:, i) -> Particular channel i (column)
    %   C(1:num_channels, :) -> Channel i to all other Channels concentration
    %   C(num_channels:end, :) -> Channel i to all other Vesiclse concentration
    C = zeros(ps{1}.nr, num_channels);   % Initialize concentration matrix
    
    % Calculate concentration from each channel
    for i = 1:num_channels
        C(:,i) = ps{i}.iterate(t); % (num_vesicles + num_channels) x 1 array
    end
    
    % Aggregate concentration for each vesicle
    for jj = 1:num_vesicles
        % Sum concentrations from channels close to this vesicle
        close_channels = vesicles.close_channels{jj};
        if isempty(close_channels)
            continue;
        end
        
        % Sum the concentrations from all relevant channels for this vesicle
        C_vesicle(jj) = sum(C(ind_r_channel + jj, close_channels));
    end
    
    % Normalize concentration
    C_vesicle = C_vesicle / N_A;  % Convert to mol/L (M)
    C = C / N_A;                  % Convert channel concentrations to mol/L (M)
end
    
    

function v = decompose_z(z, variable, size_info)

si = size_info.(variable);
v = z(si.start : si.end);

end
