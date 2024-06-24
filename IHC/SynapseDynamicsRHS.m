function SynapseDynamicsRHS()

tspan = [0 100];

opts = SynapseOptions();
initial_state = initialize_synapse_state(opts);

Vt = 10;

dynamics(tspan(0), initial_state, Vt, opts)

disp(dz)
    
end

function dz = dynamics(t, z, Vt, opts)
    %% Time Step & Index
    persistent last_t
    persistent current_index

    if isempty(last_t)
        last_t = 0;  % Initialize on first call
        current_index = 1;
    end

    % Compute timestep
    dt = t - last_t;
    last_t = t;

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

    % Update channel states
    for ii = 1:opts.num_CaV13
        if opts.channels.state(ii) == 'i'  % Inactivated
            continue;
        end

        if opts.channels.state(ii) == 'b'  % Blocked
            if t >= opts.channels.tblocked(ii)
                opts.channels.state(ii) = 'c';
            else
                continue;
            end
        end

        if rand(1) < p_block * dt / opts.channels.tau_blocked
            a = 4;
            b = opts.channels.tau_blocked / a;
            tau = b * gammaincinv(rand(1,1), a);  % Generate gamma random number
            opts.channels.state(ii) = 'b';
            opts.channels.tblocked(ii) = t + tau;
        else
            r = rand(1);
            st = find(opts.channels.state(ii) == ST);
            st_new = find(r < cQ(:,st), 1);
            opts.channels.state(ii) = ST(st_new);
            if st_new == 2  % Open state
                opts.channels.topen(ii) = t;
            end
        end
    end

    % Calculate open channel fraction
    m = sum(opts.channels.state == 'o') / opts.num_CaV13;
    Ca_blocked = sum(opts.channels.state == 'b') / opts.num_CaV13;

    % Calculate calcium concentration at vesicles
    [C_vesicles, ~] = calcium_concentration(opts.vesicles, opts.ps, t);
    C_vesicles = C_vesicles + opts.C_Ca_background;

    % Calculate calcium current for open channels
    open_channels = find(opts.channels.state == 'o');
    C_channels = calculate_concentration(opts, t);  % Update with correct concentrations for channels
    if ~isempty(open_channels)
        I_GHK = calcium_current_GHK(Vt, opts.G_Ca, ones(size(open_channels)), C_channels(open_channels));
    else
        I_GHK = zeros(size(C_channels));
    end

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

    % Build dz
    dm = (m - m_old) / dt;
    dCa_blocked = (Ca_blocked - Ca_blocked_old) / dt;
    dI = (sum(I_GHK) - I_old) / dt;
    dC = (C_vesicles - C_old) / dt;

    dz = [dm; dCa_blocked; dI; dC; dq; dc; dw; dc_proton];
    current_index = current_index + 1;
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