function SynapseDynamics()

    % Voltage step to simulate
    Vt = -.01; % Voltage in mV
    voltage_steps = linspace(-.07, -.005, 20);
    
    opts = SynapseOptions();

    % To reproduce results
    rng(opts.seed, "twister") 
    
    tspan = opts.tspan_array;

    % Initialize results
    release_rates = zeros(1, length(voltage_steps));

    % rs_s = ["r1"; "r2"; "i1"; "i2"; "i3"; "i4"; "i5"];
    % 
    % for rs = 1:length(rs_s)
    %     opts.rs = getSynapse(rs_s(rs));
    %     disp(rs_s(rs))
    for v = 1:length(voltage_steps)
        Vt = voltage_steps(v);
        initial_state = initialize_synapse_state(opts);
        solveropt = solverOpt('TimeStep', opts.dt);
        V_steady_state = -70; % Example value in mV
        [t_out, y_out] = odeEuler(@TransductionRHS_v6, tspan, initial_state, solveropt, ...
                                opts, opts.dt, Vt);
    
        release_rates(v) = calc_q_released(t_out, y_out, opts);
    end
    % 
    % end

    disp(release_rates)

    initial_state = initialize_synapse_state(opts);
    solveropt = solverOpt('TimeStep', opts.dt);

    % Simulation using odeEuler with TransductionRHS_v5
    [t_out, y_out] = odeEuler(@TransductionRHS_v6, tspan, initial_state, solveropt, ...
                              opts, opts.dt, Vt);

    calc_q_released(t_out, y_out, opts);

    % For debugging: Display output (y_out, t_out)
    disp("Simulation complete.");
    plot(t_out, decompose_z(y_out, 'NT_free', opts.size_info));
    xlabel('Time (s)');
    ylabel('State Variables');
    title('Synapse Dynamics Simulation');
end

% Helper function to initialize the state vector
function initial_state = initialize_synapse_state(opts)
    m_initial = 0;    
    Ca_blocked_initial = 0;    
    I_initial = 0;    
    C_vesicles_initial = opts.C_initial * ones(opts.num_vesicles, 1);
    q_initial = ones(opts.num_vesicles, 1);
    c_initial = 0;    
    w_initial = 0;
    c_proton_initial = 0;

    initial_state = [m_initial; Ca_blocked_initial; I_initial; C_vesicles_initial; ...
                     q_initial; c_initial; w_initial; c_proton_initial];
end


function total_release = calc_q_released(t, z_array, opts)
    q_start = opts.size_info.NT_free.start;
    q_end = opts.size_info.NT_free.end;

    q_values = z_array(:, q_start:q_end);

    % Calculate total release by summing the decreases in q
    total_release = 0;
    for i = 2:length(t)
        try
            released_this_step = q_values(i - 1, :) - q_values(i, :);
        catch e
            throw(e)
        end
        total_release = total_release + sum(released_this_step(released_this_step > 0));
    end

    disp("Total neurotransmitter released: " + total_release);
end

function v = decompose_z(z, variable, size_info)

si = size_info.(variable);
v = z(:, si.start : si.end);

end