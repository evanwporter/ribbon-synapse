function SynapseDynamics()

    % Voltage step to simulate
    % Vt = -.0495; % Voltage in mV
    voltage_steps = linspace(-.07, -.005, 20);
    voltage_steps = voltage_steps(12);
    
    opts = SynapseOptions();

    % To reproduce results
    rng(opts.seed, "twister") 
    
    if opts.ode15s_
        tspan = opts.tspan;
    else
        tspan = opts.tspan_array;
        solveropt = solverOpt('TimeStep', opts.dt);
    end

    % Initialize results
    release_rates = zeros(1, length(voltage_steps));
    % results = zeros(1001, 34, length(voltage_steps));

    opts.initial_state = initialize_synapse_state(opts);
    

    % rs_s = ["r1"; "r2"; "i1"; "i2"; "i3"; "i4"; "i5"];
    % 
    % for rs = 1:length(rs_s)
    %     opts.rs = getSynapse(rs_s(rs));
    %     disp(rs_s(rs))
    for v = 1:length(voltage_steps)
        Vt = voltage_steps(v);

        if opts.ode15s_
            warning("Ode15s is  experimental at this moment.")
            options = odeset('RelTol',1e-3, 'AbsTol',1e-6, 'InitialStep', 1e-4);
            [t_out, y_out] = ode15s(@(t, z) TransductionRHS_v6(t,z,opts,Vt), opts.tspan, opts.initial_state, options);
        else
            [t_out, y_out] = odeEuler(@TransductionRHS_v6, tspan, opts.initial_state, solveropt, ...
                                    opts, Vt, opts.dt);
        end
        
        v_s = num2str(Vt);
        v_s = v_s(4:end);
        v_s(end+1:end+10) = "_mv_e_1e_3";

        save(v_s, "y_out");
        % results(:,:,v) = y_out;
        release_rates(v) = calc_q_released(t_out, y_out, opts);
    end
    % 
    % end

    % Plot release rates against voltage steps
    figure;
    plot(voltage_steps, release_rates, '-o');
    xlabel('Voltage (V)');
    ylabel('Neurotransmitter Release Rate');
    title('Neurotransmitter Release Rate vs Voltage');
    grid on;



    disp(release_rates)

    % initial_state = initialize_synapse_state(opts);
    % solveropt = solverOpt('TimeStep', opts.dt);

    % % Simulation using odeEuler with TransductionRHS_v5
    % [t_out, y_out] = odeEuler(@TransductionRHS_v6, tspan, initial_state, solveropt, ...
    %                           opts, opts.dt, Vt);

    % calc_q_released(t_out, y_out, opts);

    % % For debugging: Display output (y_out, t_out)
    % disp("Simulation complete.");
    % plot(t_out, decompose_z(y_out, 'NT_free', opts.size_info));
    % xlabel('Time (s)');
    % ylabel('State Variables');
    % title('Synapse Dynamics Simulation');
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