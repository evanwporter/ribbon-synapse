d = 3;  % Dimensionality
D = 1e-6;  % Diffusion coefficient (m^2/s)
r = [1e-6, 2e-6, 5e-6];  % Distances from the point source (meters)
dt = 10;  % Time step (seconds)
total_time = 1000;  % Total simulation time (seconds)
it = 0:dt:total_time;  % Time steps

frequency = 1;  % Frequency in Hz
amplitude = 1e-9;  % 

% Initialize Channels object
num_channels = 2;  % Total number of channels
tau_CaV13 = 1e-3;  % Channel time constant (seconds)
tau_CaV13_blocked = 1e-2;  % Blocked state time constant (seconds)
channels = Channels(num_channels, tau_CaV13, tau_CaV13_blocked);

% Ensure two channels are always open
channels.state(:) = 'o';
channels.topen(:) = 0;  % Initially open at time 0

% Current Input
current_input = amplitude * sin(2 * pi * frequency * it);


ps = cell(num_channels, 1);
for i = 1:num_channels
    ps{i} = PointSource2(d, D, r, dt, it(1), 'rel_tol', 1e-3);
end

concentrations = zeros(length(r), length(it), num_channels);
for t = 1:length(it)
    for ch = 1:num_channels
        if channels.state(ch) == 'o'
            ps{ch}.current(t) = current_input(t);  % Set current at time step t if channel is open
            channels.topen(ch) = it(t);  % Update last open time
            ps{ch}.lastopen = channels.topen(ch);
        end
        concentrations(:, t, ch) = ps{ch}.iterate(it(t));  % Compute concentration
    end
end

% Combine concentrations from all channels
total_concentrations = sum(concentrations, 3);

figure;
hold on;
for i = 1:length(r)
    plot(it, total_concentrations(i, :), 'DisplayName', sprintf('New - Distance = %.1f um', r(i) * 1e6));
end
xlabel('Time (s)');
ylabel('Concentration (mol/m^3)');
legend;
title('New: Calcium Concentration Over Time at Different Distances');
hold off;
