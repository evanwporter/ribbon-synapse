% Calc rate of change of z

function [ dz ] = TransductionRHS_v6( t, z, opts, Vt, dt)
arguments
    t
    z
    opts
    Vt
    dt = 0
end

% if any(t)
%     disp(t)
% end

if opts.ode15s_
    if t == opts.tspan_array(end)
        dz = zeros(size(z));
        warning("Time did not increase.")
        return
    end
    opts.current_index = opts.current_index + 1;
    it = opts.current_index;
    dt_ = t - opts.tspan_array(end);
    opts.tspan_array(opts.current_index + 1) = t;
else
    dt_ = dt;
    % iteration index
    it = ceil(t/dt_);
end


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

% assert(dt_ / channels.tau_blocked <= 1);


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

Q = [1 - kp*dt_, km *dt_; 
     kp*dt_, 1 - km*dt_];


% EVAN
% Q prob of transitioning between the closed (c) and open (o) states during a time step 
%   1st row: prob of remaining closed (1 - k_p * dt_) or transitioning to open (k_m \* dt_).
%   2nd row: prob of transitioning to closed (k_p * dt_) or remaining open (1 - k_m * dt_).

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
    % Check if channel becomes unblocked
    if rand(1) < p_block * dt_ / channels.tau_blocked

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
if it > 1
    for i = 1:numel(ps)
        if all_channel_switch
            if log_partial_concentrations
                C = ps{i}.concentration(ind_r_channel, it - 1);
            else
                C = ps{i}.concentration(ind_r_channel);
            end
            
            if any(isinf(C))
                warning('Infinite concentration detected at iteration %d for ps{%d}', it, i);
                disp('Concentration values:');
                disp(C);
                disp('Time: '); disp(t)
            end
            
            C_channels(i) = sum(C);
            
            % Check for unrealistic values and cap them
            % if C_channels(i) > 1e3  % Example cap, adjust based on realistic limits
            %     warning('Capping concentration at iteration %d for ps{%d}, C = %d', it, i, C_channels(i));
            %     C_channels(i) = 1e3;
            % end
        else
            C_channels{i} = ps{i}.concentration(ind_r_channel, it - 1);
        end
    end
end

C_channels = C_channels + opts.C_Ca_background;

open_channels = channels.state == 'o';


%% Calculate current

% Evan -
% I_GHK: Used for calculating calcium current flowing into hair cell
%   Needed for determing local concentration of Calcium
%       Which in tern is needed for calculating neurotransmitter release rate

% I_ohm = calcium_current(Vt, G_Ca, open_channels, C_channels);
% I_GHK = calcium_current_GHK(-70, G_Ca, open_channels, C_channels);
I_GHK = calcium_current_GHK(Vt, opts.G_Ca, open_channels, C_channels);

I = I_GHK;


for i = 1:numel(ps)
    ps{i}.current(it) = - I(i); % minus because of convention
    ps{i}.lastopen = channels.topen(i);
end

[C_vesicles, C_all] = calcium_concentration(vesicles, ps, it, all_channel_switch, opts);

for i = 1:numel(ps)
    if log_partial_concentrations
        ps{i}.concentration(:, it) = C_all(:, i);
    else
        ps{i}.concentration = C_all(:, i);
    end
end

C_vesicles = C_vesicles + opts.C_Ca_background;



%% Transmitter Release and Recycling

% dq : RoC of Neurotransmitter contents in the vesicle; Num_Vesicle x 1
% dc : RoC of Neurotransmitter concentration in the synaptic cleft ??
% dw : Neurotransmitter Reprocessing rate
% dc_proton : RoC of proton concentration in the synaptic cleft
% vesicles : Updated vesicle state with NT release events logged (not used)

[dq, dc, dw, dc_proton, vesicles] = NTdynamicsRHS_v6_core( t, ...
    q, c, w, c_proton, ...
    vesicles, opts, ...
    dt_, C_vesicles);


    
%% Build dz

% dm : RoC of fraction of open channels (m)
% dCa_blocked : RoC of fraction of blocked channels (Ca_blocked)
% dI : RoC of calcium current (I)
% dC : RoC of calcium concentration at the vesicles (C_vesicles)

dm = (m - m_old) / dt_;

dCa_blocked = (Ca_blocked - Ca_blocked_old) / dt_;

charge = 2;
amp_to_electron_per_second = 6.242e18;

I = sum(I) / amp_to_electron_per_second * charge;
dI = (I - I_old) / dt_;

dC = (C_vesicles - C_old) / dt_;

dz = [dm; dCa_blocked; dI; dC; dq; dc; dw; dc_proton];

d_state = [d_state_inactivated; d_state_normal; d_state_burst] / dt_;

% dz = [dz; d_state];

if any(isnan(dz))
    error('NaN encountered')
end


% disp(['Time: ', num2str(t), ' dt: ', num2str(dt_)]);
% disp(['State z: ', num2str(z')]);

% disp(['Rate of change dz: ', num2str(dz')]);


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
    
    Cin = C * 1000; % M to mol/m3
    Cout = C_extracellular * 1000; % M to mol/m3
    
    
    % charge ... charge [K]
    % T ... temperature [K]
    %
    % Phi is the current density (flux) across the membrane carried by ion S, measured in amperes per square meter (A·m−2)
    % arguments
    %     Vm double % transmembrane potential in volts
    %     Cin double % intracellular concentration of ion S, measured in mol·m−3 or mmol·l−1
    %     Cout double % extracellular concentration of ion S, measured in mol·m−3 or mmol·l−1
    %     P double = 1; % permeability of the membrane for ion S measured in m·s−1
    % end
    
    F = 96485.3328959;   % [C/mol]
    R = 8.314459848;     % [J/mol/K]
    
    % Phi = P * charge^2 * (Vm*F^2 / (R*T)) * (Cin - Cout*exp(-charge*Vm*F/(R*T))) / (1 - exp(-charge*Vm*F/(R*T)));
    
    U = charge * Vm * F / R / T;
    emU = exp(-U);
    % epU = exp(U);
    
    Phi = P .* charge * F * U .* (Cin - Cout.*emU) ./ (1 - emU);
    
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
function [C_vesicle, C] = calcium_concentration(vesicles, ps, it, all_channel_switch, opts)
    arguments
        vesicles
        ps
        it
        all_channel_switch
        opts
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
        if opts.ode15s_
            C(:, i) = ps{i}.dt_iterate(opts.tspan_array) / 1e3;
        else
            if opts.e_iter
                C(:,i) = ps{i}.e_iterate(it) / 1e3;
            else
                C(:,i) = ps{i}.iterate(it) / 1e3; % (num_vesicles + num_channels) x 1 array
            end
        end
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