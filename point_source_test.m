method = "evans"

d = 3;  % Dimensionality
D = 5.2e-10;  % Diffusion coefficient (m^2/s)
r = [2.5e-07 4e-7];%, 5e-6];  % Distances from the point source (meters)
dt = 1e-6;  % Time step (seconds)
total_time = 10e-3;  % Total simulation time (seconds)
it = 0:dt:total_time;  % Time steps

frequency = 1;  % Frequency in Hz
amplitude = 40e-12;

% Channels object
num_channels = 1; 
tau_CaV13 = 1e-3;
tau_CaV13_blocked = 1e-2;
channels = Channels(num_channels, tau_CaV13, tau_CaV13_blocked);

% Ensure two channels are always open
channels.state(:) = 'o';
channels.topen(:) = 0;  % Initially open at time 0

% Current input
current_input = amplitude * abs(sin(2 * pi * frequency * it)); % sinusoidal current
% current_input = amplitude * ones(size(it)); % constant current
% current_input = amplitude * (1 + it / max(it)); % simple increasing current
% current_input = amplitude * (1 - it / max(it));   % decreasing
% current_input = opts.ps{1}.current;


% Initialize PointSource objects for each channel
ps = cell(num_channels, 1);
for i = 1:num_channels
    ps{i} = PointSource(d, D, r, dt, it);%, 'rel_tol', 1e-3);
    % ps{i}.N = 100;
end


concentrations = zeros(length(r), length(it), num_channels);
for t = 2:length(it)
    for ch = 1:num_channels
        % Needs to be done for iterate method to work
        if channels.state(ch) == 'o'
            ps{ch}.current(t) = current_input(t);
            channels.topen(ch) = it(t);
            ps{ch}.lastopen = channels.topen(ch);
        end
        % if t == 3
        %     channels.state(1) = 'c';
        % end

        if method == "evans"
            concentrations(:, t, ch) = ps{ch}.e_iterate(t);  % Compute concentration
        else
            concentrations(:, t, ch) = ps{ch}.iterate(t);% / 1000;
        end
    end
end

% Combine concentrations from all channels
% There's only one so this doesn't do anything
total_concentrations = sum(concentrations, 3);

figure;
hold on;
for i = 1:length(r)
    plot(it, total_concentrations(i, :), 'DisplayName', sprintf('Distance = %.1f um', r(i) * 1e6));
end
xlabel('Time (s)');
ylabel('Concentration (mol/m^3)');
legend;
title('Ondras Method, N=100');
hold off;


% Calculate at a single time point
% t = 101;
% ps{1}.current(t) = current_input(t);
% channels.topen(ch) = it(t);
% ps{ch}.lastopen = channels.topen(ch);
% C = ps{ch}.iterate(t);
% disp(["Concnetration of PointSource class ", num2str(C(1))]);