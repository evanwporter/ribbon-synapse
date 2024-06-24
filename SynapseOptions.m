classdef SynapseOptions
    properties (SetAccess = immutable)
        % https://github.com/evanwporter/cochlea-nerve/blob/4e8b4f18f20782bfc39d88589db9f4e04dbcf507/Wrapper/Options/transductionOpt_v4_1.m#L354-L355
        num_release_sites (1,1) {mustBePositive, mustBeInteger} = 14;
        num_CaV13 (1,1) {mustBePositive, mustBeInteger} = 84; %72

        % https://github.com/evanwporter/cochlea-nerve/blob/50d39b91828e149530a871f8f2cffa432c6c53f0/Wrapper/Options/transductionOpt_v4_1.m#L15-L18
        Ca_diffusion_coefficient (1,1) double {mustBePositive} = 5.2e-10;

        % https://github.com/evanwporter/cochlea-nerve/blob/50d39b91828e149530a871f8f2cffa432c6c53f0/Wrapper/Options/transductionOpt_v4_1.m#L8-L11
        Ca_conc_rel_tol (1,1) double {mustBeNonnegative} = 10e-2;

        % https://github.com/evanwporter/cochlea-nerve/blob/50d39b91828e149530a871f8f2cffa432c6c53f0/Wrapper/Options/transductionOpt_v4_1.m#L339
        d_nernst (1,1) double {mustBeNonnegative} = 20e-9;

        num_inactive_release_sites (1,1) {mustBeNonnegative, mustBeInteger} = 0;

        num_vesicles (1,1) {mustBeNonnegative, mustBeInteger};

        tau_CaV13  % see Zampini 2010 Table 1
        tau_CaV13_blocked

        channels Channels;
        vesicles Vesicles;

        y Frequency;
        l Frequency;
        x Frequency;
        r Frequency;

        transmitter_release_parameters

        % https://github.com/evanwporter/cochlea-nerve/blob/50d39b91828e149530a871f8f2cffa432c6c53f0/Wrapper/Options/transductionOpt_v4_1.m#L80-L82
        channels_markov_properties struct = struct('alpha', -133, 'kp0', 5e5);

        channels_open_ss_parameters_normal cell;
        channels_open_ss_parameters_burst cell;

        ps

        C_initial (1,1) double {mustBePositive} = 1e-7; % Initial calcium concentration in M

        %% https://github.com/evanwporter/cochlea-nerve/blob/cc845a8870e4825796b05a13568a15e4361ce6cf/IHC/Transduction_v4.m#L75-L81
        C_Ca_background_base (1,1) double {mustBePositive} = 17.1e-6; % M
        C_Ca_background_apex (1,1) double {mustBePositive} = 40.6e-6; % M

        %% https://github.com/evanwporter/cochlea-nerve/blob/50d39b91828e149530a871f8f2cffa432c6c53f0/Wrapper/Options/transductionOpt_v4_1.m#L27-L28
        C_Ca_background_factor (1,1) double {mustBePositive} = 1;

        C_Ca_background (1,1) double;

        % Position property, must be between 0 and 1 (I think)
        x_pos (1,1) double {mustBeGreaterThanOrEqual(x_pos, 0), mustBeLessThanOrEqual(x_pos, 1)} = 0.5;

        size_info struct;

        simulation_units struct = struct( ...
            'Voltage', 'V', ...
            'Frequency', 'Hz', ...
            'Concentration', 'M', ...
            'Time', 's'...
        );

        % Calcium Conductance
        % https://github.com/evanwporter/cochlea-nerve/blob/50d39b91828e149530a871f8f2cffa432c6c53f0/Wrapper/Options/transductionOpt_v4_1.m#L333
        G_Ca (1,1) double {mustBePositive} = 15e-12;

        dt (1,1) double {mustBePositive} = 1e-4;
        tspan double {mustBeReal, mustBeFinite, mustBeNonnegative};

    end

    properties
        rs RibbonSynapse_v4;
    end
    
    methods
        function obj = SynapseOptions()

            uconv = @(var) Unit.batch_convert(obj.simulation_units, var);

            obj.num_vesicles = obj.num_release_sites - obj.num_inactive_release_sites;

            
            %% Set Time

            obj.tspan = [0 obj.dt * 10];

            [numSteps, numSamples] = odeEuler_tspan(obj.tspan, obj.dt);

            %% Channels

            % https://github.com/evanwporter/cochlea-nerve/blob/4e8b4f18f20782bfc39d88589db9f4e04dbcf507/Wrapper/Options/transductionOpt_v4_1.m#L335-L336
            obj.tau_CaV13 = Time(500, 'us');        % Calcium channel time constant
            obj.tau_CaV13_blocked = Time(1, 'ms');  % Blocked state time constant

            % https://github.com/evanwporter/cochlea-nerve/blob/cc845a8870e4825796b05a13568a15e4361ce6cf/IHC/Transduction_v4.m#L147-L148
            obj.channels = Channels(obj.num_CaV13, obj.tau_CaV13.s, obj.tau_CaV13_blocked.s);

            obj.x_pos = .5;

            obj.transmitter_release_parameters = uconv(obj.generate_transmitter_release_parameters(obj.x_pos));
            
            % https://github.com/evanwporter/cochlea-nerve/blob/50d39b91828e149530a871f8f2cffa432c6c53f0/IHC/Transduction_v4.m#L64-L65
            obj.channels_open_ss_parameters_normal = uconv(obj.generate_channels_open_ss_parameters('normal', obj.x_pos));
            obj.channels_open_ss_parameters_burst = uconv(obj.generate_channels_open_ss_parameters('burst', obj.x_pos));
            
            V0t = obj.channels_open_ss_parameters_normal{1};
            S0t = obj.channels_open_ss_parameters_normal{2};
    
            obj.channels.alpha = obj.channels_markov_properties.alpha;
            obj.channels.kp0 = obj.channels_markov_properties.kp0;
            obj.channels.V0t = V0t;
            obj.channels.S0t = S0t;

            % https://github.com/evanwporter/cochlea-nerve/blob/4e8b4f18f20782bfc39d88589db9f4e04dbcf507/IHC/Transduction_v4.m#L514-L515
            obj.vesicles = Vesicles(obj.num_vesicles);

            % https://github.com/evanwporter/cochlea-nerve/blob/50d39b91828e149530a871f8f2cffa432c6c53f0/IHC/Transduction_v4.m#L213-L228
            for i = 1:obj.num_vesicles
                obj.vesicles.close_channels{i} = 1:obj.num_CaV13;
            end

            %% https://github.com/evanwporter/cochlea-nerve/blob/4e8b4f18f20782bfc39d88589db9f4e04dbcf507/Wrapper/Options/transductionOpt_v4_1.m#L38-L50
            obj.y = Frequency(10, 'Hz'); % Sumner 2002
            obj.l = Frequency(1290, 'Hz'); % to roughly match the time-constant of H&H
            obj.x = Frequency(66.3, 'Hz'); % Sumner 2002
            obj.r = Frequency(3290, 'Hz'); % to roughly match the time-constant of H&H

            % Initialize the RibbonSynapse
            obj.rs = RibbonSynapse_v4('num_channels', obj.num_CaV13, 'num_release_sites', obj.num_release_sites);

            d = 3; %           ...  geometry factor

            %% Create PointSources
            for i = 1:obj.num_CaV13
                rho = obj.rs.rho(:,i) * 1e-9; % channel mouth--vesicle membrane distance
                psi = obj.rs.psi(:,i) * 1e-9; % channel mouth--channel mouth distance
                psi_nernst = sqrt(psi.^2 + obj.d_nernst.^2); % channel mouth--channel nernst point (above channel mouth)
                r = [psi_nernst; rho];
                obj.ps{i} = Diffusion.PointSource(d, obj.Ca_diffusion_coefficient, r, obj.dt, 0:numSamples, 'rel_tol', obj.Ca_conc_rel_tol);
            end

            %% https://github.com/evanwporter/cochlea-nerve/blob/cc845a8870e4825796b05a13568a15e4361ce6cf/IHC/Transduction_v4.m#L75-L81
            obj.C_Ca_background = obj.C_Ca_background_base + (obj.C_Ca_background_apex - obj.C_Ca_background_base) * obj.x_pos;
            obj.C_Ca_background = obj.C_Ca_background * obj.C_Ca_background_factor;

            % https://github.com/evanwporter/cochlea-nerve/blob/4e8b4f18f20782bfc39d88589db9f4e04dbcf507/IHC/Transduction_v4.m#L317-L328
            obj.size_info = create_size_info( ...
                "CaV13_channels_fraction", 1, ...         % m
                "CaV13_channels_blocked", 1, ...          % Ca_blocked
                "Ca_current", 1, ...                      % I
                "Ca_concentration", obj.num_vesicles, ... % C_vesicles
                "NT_free", obj.num_vesicles, ...          % q
                "NT_cleft", 1, ...                        % c
                "NT_reprocessing", 1, ...                 % w
                "proton_cleft", 1, ...                    % c_proton
                "CaV13_num_inactivated", 1, ...           % CaV13_num_inactivated
                "CaV13_num_normal", 1, ...                % CaV13_num_normal
                "CaV13_num_burst", 1);                    % CaV13_num_burst

        end

        % https://github.com/evanwporter/cochlea-nerve/blob/50d39b91828e149530a871f8f2cffa432c6c53f0/Wrapper/Options/transductionOpt_v4_1.m#L130-L140
        function params = generate_transmitter_release_parameters(obj, x_pos)
            
            % https://github.com/evanwporter/cochlea-nerve/blob/50d39b91828e149530a871f8f2cffa432c6c53f0/Wrapper/Options/transductionOpt_v4_1.m#L341-L350
            transmitter_release_parameters = struct( ...
                'apex', {{ ...
                    Frequency(500, 'Hz'), ...
                    5.0, ...      % Hill's coefficient
                    100e-6}}, ... % Half activation [Ca2+ M]
                'base', {{ ...
                    Frequency(500, 'Hz'), ...
                    5.0, ...
                    100e-6}} ...
            );

            par_base = transmitter_release_parameters.base;
            par_apex = transmitter_release_parameters.apex;

            params = cell(numel(par_base), 1);
            for i = 1:numel(par_base)
                params{i} = par_base{i} + (par_apex{i} - par_base{i}) * x_pos;
            end

        end

        % https://github.com/evanwporter/cochlea-nerve/blob/50d39b91828e149530a871f8f2cffa432c6c53f0/Wrapper/Options/transductionOpt_v4_1.m#L142-L153
        function params = generate_channels_open_ss_parameters(obj, mode, x_pos)
            
            % https://github.com/evanwporter/cochlea-nerve/blob/50d39b91828e149530a871f8f2cffa432c6c53f0/Wrapper/Options/transductionOpt_v4_1.m#L53-L77
            channels_open_ss_parameters = struct( ...
                'normal', struct( ...
                    'apex', {{ ...
                        Voltage(-23.4, 'mV'), ...   % [V] - IHC voltage - half-activation boltzmann function parameter
                        Voltage(8, 'mV'), ... % [V] - IHC voltage - slope boltzmann function parameter
                        0.0, 1, ...
                        }}, ...
                    'base', {{ ...
                        Voltage(-25.9, 'mV'), ...   % [V] - IHC voltage - half-activation boltzmann function parameter
                        Voltage(7.4, 'mV'), ...  % [V] - IHC voltage - slope boltzmann function parameter
                        0.0, 1, ...
                        }} ...
                ), ...
                'burst', struct( ...
                    'apex', {{ ...
                        Voltage(0, 'mV'), ...
                        Voltage(1, 'mV'), ... % [V] - IHC voltage - slope boltzmann function parameter
                        0, 1, }}, ...
                    'base', {{ ...
                        Voltage(0, 'mV'), ...
                        Voltage(1, 'mV'), ...  % [V] - IHC voltage - slope boltzmann function parameter
                        0, 1, }} ...
                ) ...
            );
            
            par_base = channels_open_ss_parameters.(mode).base;
            par_apex = channels_open_ss_parameters.(mode).apex;
            
            % linear model
            params = cell(numel(par_base), 1);
            for i = 1:numel(par_base)
                params{i} = par_base{i} + (par_apex{i} - par_base{i}) * x_pos;
            end
        end
        
        function obj = set.dt(obj, val)
            for i = 1:obj.num_CaV13
                obj.ps{i}.dt = val;
            end
            obj.dt = val;
        end

        
    end
end

