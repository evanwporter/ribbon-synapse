classdef SynapseOptions
    properties
        G_Ca                % Conductance in Siemens
        C_Ca_background     % Background calcium concentration in M
        vesicles            % Vesicle structure
        ps                  % Channel properties structure
        rate_y              % Manufacturing rate
        rate_l              % Loss rate
        rate_x              % Reprocessing rate
        rate_r              % Re-uptake rate
        transmitter_release_parameters % Parameters for transmitter release
        all_channel_switch  % Boolean for channel update
        dt                  % Time step
        C_initial           % Initial calcium concentration in M
    end
    
    methods
        function obj = SynapseOptions(config)
            % Initialize properties based on configuration
            switch config
                case 'r1'
                    obj.G_Ca = 1e-12;  % Example
                    obj.C_Ca_background = 1e-6;  % Example
                    obj.vesicles = struct('num', 10, 'close_channels', {1:10});  % Example
                    obj.ps = struct('num', 10, 'concentration', {1:10});  % Example
                    obj.rate_y = 0.5;  % Example
                    obj.rate_l = 0.3;  % Example
                    obj.rate_x = 0.2;  % Example
                    obj.rate_r = 0.1;  % Example
                    obj.transmitter_release_parameters = {1, 5, 0.5};  % Example
                    obj.all_channel_switch = true;  % Example
                    obj.dt = 1e-4;  % Time step
                    obj.C_initial = 1e-7;  % Initial calcium concentration in M
                % Add other cases for different configurations here
                otherwise
                    error('Unknown configuration: %s', config);
            end
        end
    end
end
