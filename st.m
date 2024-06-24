function st()
    
    % Voltage steps to simulate
    voltage_steps = linspace(-100, 100, 100); % Voltage range from -50 mV to 50 mV

    % Initialize results
    release_rates = zeros(length(voltage_steps), length(configurations));
    
    % Simulation parameters
    total_time = 1;  % Total simulation time
    tspan = [0 1];
    
    opts = SynapseOptions();  % Initialize options for this configuration
    
    for v = 1:length(voltage_steps)
        Vt = voltage_steps(v);
        disp(Vt)
        
        % Initialize synapse state (assuming initial conditions similar for all)
        initial_state = initialize_synapse_state(opts);

        options = odeset('RelTol',1e-5, 'AbsTol',1e-8, 'InitialStep',1e-12, 'MaxStep',1e-4);
        
        % Solve ODEs using ode15s
        [T, state] = ode15s(@(t, y) synapse_dynamics(t, y, Vt, opts), tspan, initial_state, options);
        
        % Extract final neurotransmitter release rate
        release_rates(v, c) = calculate_release_rate(state(end, :), opts);
    end

end

    disp("NO ERRORS")

    % Plot results
    figure;
    hold on;
    for c = 1:length(configurations)
        plot(voltage_steps, release_rates(:, c), 'DisplayName', configurations{c});
    end
    xlabel('Voltage Step (V)');
    ylabel('Neurotransmitter Release Rate');
    legend();
    hold off;
end

function dz = synapse_dynamics(t, z, Vt, opts)
    persistent last_t
    persistent current_index

    if isempty(last_t)
        last_t = -1e-12;  % Initialize on first call
        current_index = 1;
    end

    % Compute timestep
    dt = t - last_t;
    if dt <= 1e-12  % Prevent division by zero or too small step sizes
        dt = 1e-12;
    end
    last_t = t;

    opts.dt= dt;

    % Extract size info from opts
    size_info = opts.size_info;
    
    % Decompose the state vector z
    m = decompose_z(z, 'CaV13_channels_fraction', size_info);
    Ca_blocked = decompose_z(z, 'CaV13_channels_blocked', size_info);
    C_vesicles = decompose_z(z, 'Ca_concentration', size_info);
    I = decompose_z(z, 'Ca_current', size_info);
    q = decompose_z(z, 'NT_free', size_info);
    c = decompose_z(z, 'NT_cleft', size_info);
    w = decompose_z(z, 'NT_reprocessing', size_info);
    c_proton = decompose_z(z, 'proton_cleft', size_info);

    m_old = m;
    I_old = I;
    Ca_blocked_old = Ca_blocked;
    C_old = C_vesicles;

    % Calculate rates
    S0t = opts.channels.S0t;
    V0t = opts.channels.V0t;

    alpha = opts.channels.alpha;
    beta = alpha + 1/S0t;

    kp0 = opts.channels.kp0;
    km0 = kp0 * exp(V0t * (beta - alpha));

    kp = kp0 * exp(-alpha * Vt); % Opening rate
    km = km0 * exp(-beta * Vt);  % Closing rate

    % Transition matrix for channel states
    Q = [1 - kp*dt, km*dt; 
         kp*dt, 1 - km*dt];
    cQ = cumsum(Q, 1);
    ST = ['c', 'o'];

    p_block = CaV13.CaProtonBlock(c_proton);

    %% Update channel states
    for ii = 1:opts.num_CaV13
        % If inactive, skip
        if opts.channels.state(ii) == 'i'  % Inactivated
            continue;
        end

        % Check if blocked time has passed
        %   If yes then it becomes closed
        if opts.channels.state(ii) == 'b'  % Blocked
            if t >= opts.channels.tblocked(ii)
                opts.channels.state(ii) = 'c';
            else
                continue;
            end
        end

        % Checks whether channel becomes blocked due to proton binding
        % Probability of this event occurring is p_block
        %   timestep (dt) is compared against a random number. 
        % If the channel becomes blocked:
        %   A gamma random variable is generated using gammaincinv to 
        %       determine the blocking duration (tau).
        %   The channel's state is set to blocked ('b').
        %   The block end time (tblocked) is updated to the current 
        %       time plus the blocking duration (tau).
        if rand(1) < p_block * dt / opts.channels.tau_blocked
            a = 4;
            b = opts.channels.tau_blocked / a;
            tau = b * gammaincinv(rand(1,1), a);
            opts.channels.state(ii) = 'b';
            opts.channels.tblocked(ii) = t + tau;

        % Channel becomes transitions to closed (if open) or open (if
        %   closed)
        else
            r = rand(1);
        
            if opts.channels.state(ii) == 'c' % closed
                st = 1;
            elseif opts.channels.state(ii) == 'o' % open
                opts.channels.topen(ii) = t; % last open time
                st = 2;
            else
                error('unknown state %s', opts.channels.state(ii))
            end
            
            st_new = find(r < cQ(:,st), 1);
            opts.channels.state(ii) = ST(st_new);

        end
    end

    %% Calculate fractions
    % Fraction of calcium channels that are in the open state
    m = sum(opts.channels.state == 'o') / opts.num_CaV13;

    % Fraction of calcium channels that are in the close state
    Ca_blocked = sum(opts.channels.state == 'b') / opts.num_CaV13;



    % Calculate calcium concentration at vesicles
    [C_vesicles, ~] = calcium_concentration(opts.vesicles, opts.ps, t);
    C_vesicles = C_vesicles + opts.C_Ca_background;

    
    % Calculate calcium current for open channels
    open_channels = opts.channels.state == 'o';
    C_channels = calculate_concentration(opts, t);  % Update with correct concentrations for channels

    I_GHK = calcium_current_GHK(Vt, opts.G_Ca, open_channels, C_channels);

    open_channels = find(opts.channels.state == 'o');

    % Update PointSource with the calculated current
    for ch = 1:numel(open_channels)
        i = open_channels(ch);
        if i <= numel(opts.ps)
            opts.ps{i}.current(current_index) = - I_GHK(ch); % minus because of convention
            opts.ps{i}.lastopen = opts.channels.topen(i);
        end
    end

    % Calculate neurotransmitter release dynamics
    [dq, dc, dw, dc_proton, ~] = NTdynamicsRHS_v5_core(t, q, c, w, c_proton, ...
        opts.vesicles, opts.y.Hz, opts.l.Hz, opts.x.Hz, opts.r.Hz, ...
        opts.transmitter_release_parameters, dt, C_vesicles);

    
    %% Build dz
    dm = (m - m_old) / dt;
    dCa_blocked = (Ca_blocked - Ca_blocked_old) / dt;

    charge = 2;
    amp_to_electron_per_second = 6.242e18;
    I = sum(I) / amp_to_electron_per_second * charge;
    dI = (I - I_old) / dt;

    dC = (C_vesicles - C_old) / dt;

    dz = [dm; dCa_blocked; dI; dC; dq; dc; dw; dc_proton];

    % d_state = [d_state_inactivated; d_state_normal; d_state_burst] / dt;
    % dz = [dz; d_state];

    if any(isnan(dz))
        error('NaN encountered')
    end

    current_index = current_index + 1;

end

function release_rate = calculate_release_rate(state, opts)
    % Extract the neurotransmitter release rate from the state
    % Customize this based on your state variables and desired output
    release_rate = sum(state(opts.num_vesicles + 2:opts.num_vesicles + 1 + opts.num_vesicles));  % Sum of the neurotransmitter free pool
end

function initial_state = initialize_synapse_state(opts)    
    m_initial = 0;    
    Ca_blocked_initial = 0;    
    I_initial = 0;    
    C_vesicles_initial = opts.C_initial * ones(opts.num_vesicles, 1);
    q_initial = ones(opts.num_vesicles, 1);  % Example non-zero initial value
    c_initial = 0;    
    w_initial = 0;
    c_proton_initial = 0;

    initial_state = [m_initial; Ca_blocked_initial; I_initial; C_vesicles_initial; ...
                     q_initial; c_initial; w_initial; c_proton_initial];
end

function I = calcium_current_GHK(Vm, G, m, C)
    arguments
        Vm (1,1) double
        G (1,1) double
        m (:,1) double  % num_channel x 1
        C (:,1) double  % num_channel x 1
    end

    % Constants
    charge = 2;                  % Ca2+ ion charge
    C_extracellular = 1.3e-3;    % Extracellular concentration of calcium (in mol/L or M)
    T = 310;                     % Absolute temperature (in Kelvin), approximately 37°C
    area = 1;                    % Membrane area (in m²), assumed to be 1 for simplicity
    me = 9.109383632e-31;        % Mass of electron (in kg)
    ce = -1.602e-19;             % Charge of electron (in coulombs)
    F = 96485.3328959;           % Faraday constant (in C/mol)
    R = 8.314459848;             % Universal gas constant (in J/mol/K)

    % Permeability
    fac = 6.0e-15;                % Factor for permeability conversion
    P = G .* m / ce^2 / charge^2 * me * fac;
    P = P / area;                 % Permeability (P) in appropriate units

    % Concentration
    Cin = C * 1000;               % Intracellular concentration (from M to mol/m³)
    Cout = C_extracellular * 1000;% Extracellular concentration (from M to mol/m³)
    
    % GHK Flux Calculation
    U = charge * Vm * F / R / T;      % Dimensionless voltage term
    emU = exp(-U);                    % Exponential term for flux calculation
    Phi = P .* charge * F * U .* (Cin - Cout.*emU) ./ (1 - emU);

    % Current Calculation
    I = Phi * area;                   % Flux converted to current
    amp_to_electron_per_second = 6.242e18; % Conversion factor for ions/sec to current
    I = I * amp_to_electron_per_second / charge; % Final current in ions per second

end

function [C_vesicle, C] = calcium_concentration(vesicles, ps, t)
    N_A = 6.02214076e23;  % Avogadro's constant, mol^-1
    C = zeros(ps{1}.nr, vesicles.num);  % Initialize concentration matrix
    
    for i = 1:vesicles.num
        close_channels = vesicles.close_channels{i};
        for ch = close_channels
            C(:, i) = C(:, i) + ps{ch}.iterate(t);
        end
    end

    C = C / N_A;  % Convert to mol/L (M)
    C_vesicle = sum(C, 1)';  % Sum concentrations at each vesicle and transpose to column vector

    C = C / 1e3;
    C_vesicle = C_vesicle / 1e3;
end

function C_channels = calculate_concentration(opts, t)
    C_channels = zeros(numel(opts.ps), 1);
    for i = 1:numel(opts.ps)
        C_channels(i) = sum(opts.ps{i}.iterate(t));
    end
end

function v = decompose_z(z, variable, size_info)
    si = size_info.(variable);
    v = z(si.start : si.end);
end
