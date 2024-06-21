function st()
    % Define parameters for each configuration
    configurations = {'r1'}; % Add more configurations as needed
    
    % Voltage steps to simulate
    voltage_steps = linspace(-0.1, 0.1, 100);  % Example range from -100mV to 100mV

    % Initialize results
    release_rates = zeros(length(voltage_steps), length(configurations));
    
    % Simulation parameters
    total_time = 0.01;  % Total simulation time
    tspan = [0 total_time];
    
    for c = 1:length(configurations)
        config = configurations{c};
        opts = SynapseOptions(config);  % Initialize options for this configuration
        
        for v = 1:length(voltage_steps)
            Vt = voltage_steps(v);
            
            % Initialize synapse state (assuming initial conditions similar for all)
            initial_state = initialize_synapse_state(opts);
            
            % Solve ODEs using ode15s
            [~, state] = ode15s(@(t, y) synapse_dynamics(t, y, Vt, opts), tspan, initial_state);
            
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

function dy = synapse_dynamics(t, y, Vt, opts)
    % Extract state variables
    C_vesicles = y(1:10);    % Calcium concentration at vesicles
    c_proton = y(11);        % Proton concentration
    q = y(12:21);            % Neurotransmitter free pool
    c = y(22);               % Neurotransmitter cleft concentration
    w = y(23);               % Neurotransmitter reprocessing

    % Calculate calcium current
    I_GHK = calcium_current_GHK(Vt, opts.G_Ca, C_vesicles, opts.C_Ca_background);
    
    % Calculate calcium concentration
    [C_vesicles, ~] = calcium_concentration(opts.vesicles, opts.ps, t, opts.all_channel_switch);
    
    % Calculate neurotransmitter release dynamics
    [dq, dc, dw, ~, ~] = NTdynamicsRHS_v5_core(t, q, c, w, c_proton, ...
        opts.vesicles, opts.rate_y, opts.rate_l, opts.rate_x, opts.rate_r, ...
        opts.transmitter_release_parameters, opts.dt, C_vesicles);
    
    % Build the derivative of the state vector
    dy = [C_vesicles; c_proton; dq; dc; dw];  % Adjust based on state variables
end

function release_rate = calculate_release_rate(state, opts)
    % Extract the neurotransmitter release rate from the state
    % Customize this based on your state variables and desired output
    release_rate = sum(state(12:21));  % Sum of the neurotransmitter free pool
end

function initial_state = initialize_synapse_state(opts)
    % Define the initial state of the synapse
    % Customize based on the synapse configuration
    initial_state = [opts.C_initial * ones(10, 1); 0; ones(10, 1); 0; 0];  % Example initial state
end

function I = calcium_current_GHK(Vm, G, C_in, C_background)
    % Calculates the calcium current using the Goldman-Hodgkin-Katz equation
    % Vm: Membrane potential (V)
    % G: Conductance (S)
    % C_in: Intracellular calcium concentration (M)
    % C_background: Extracellular calcium concentration (M)
    
    % Constants
    charge = 2;  % Ca2+ ion charge
    R = 8.314;  % Universal gas constant (J/(mol·K))
    T = 310;  % Absolute temperature (K)
    F = 96485;  % Faraday constant (C/mol)

    % Calculate the reversal potential
    E_Ca = (R * T) / (charge * F) * log(C_background / C_in);
    
    % Calculate the current using the GHK equation
    I = G * (Vm - E_Ca);
end

function [C_vesicle, C] = calcium_concentration(vesicles, ps, t, all_channel_switch)
    % Calculates the calcium concentration at vesicles and channels
    % vesicles: Vesicle structure
    % ps: Channel properties structure
    % t: Time (s)
    % all_channel_switch: Boolean to decide on full channel update

    N_A = 6.02214076e23;  % Avogadro's constant (mol^-1)
    C = zeros(ps.num, vesicles.num);  % Initialize concentration matrix
    
    for i = 1:vesicles.num
        % Sum concentrations of close channels
        close_channels = vesicles.close_channels(i);
        C(:,i) = sum([ps.concentration(:, close_channels)], 2);
    end

    % Normalize concentration
    C = C / N_A;  % Convert to mol/L (M)
    C_vesicle = sum(C, 2);  % Sum concentrations at each vesicle
end
