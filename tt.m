% % Parameters
D = 1e-6;
dt = 10; % Time step (seconds)
total_time = 1000;  % Total simulation time (seconds)
r = 2e-6;
t = flip(0:dt:total_time);  % Time steps
amplitude = 1; % Assuming some amplitude
frequency = 0.1; % Assuming some frequency
I = @(t) amplitude * sin(2 * pi * frequency * t);

% Preallocate concentration array
C = zeros(length(t), length(r));

% Calculate concentration for each combination of t and r
for i = 1:length(t)
    for j = 1:length(r)
        % Sum contributions from all previous time steps
        for k = 1:i
            if t(i) ~= t(k)  % Avoid division by zero
                contribution = (I(t(k)) / (4 * pi * D * r(j))) * erfc(r(j) / sqrt(4 * D * (t(k) - t(i) + dt)));
                C(i, j) = C(i, j) + contribution;
            end
        end
    end
end

% Plot the results
figure;
hold on;
plot(t, C(:,1));
% plot(t, C(:,2));
xlabel('Time t');
ylabel('Concentration C');
zlabel('Concentration C');
title('Concentration of Diffusing Substance from a Point Source');
hold off;
