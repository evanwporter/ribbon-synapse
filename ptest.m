% Parameters
D = 1e-6;
dt = 1e-4; % seconds
total_time = 1;  % seconds
r = [1e-6 2e6];
t = 0:dt:total_time;  

amplitude = 1;
frequency = 0.1;
I = amplitude * abs(sin(2 * pi * frequency * t));
% I = amplitude * ones(size(t)); % constant
% I = amplitude * (1 + t / max(t)); % increasing
% I = amplitude * (1 - t / max(t));   % decreasing


% Preallocate concentration array
C = zeros(length(t), length(r));

% GreensFunction = @(t) ( sqrt(pi) * erfc(obj.r ./ sqrt(4*obj.D*t)) ) ./ ( 4 * obj.D * pi^(obj.d/2) * r );


% for it = 1:length(t)
%     for j = 1:length(r)
%         C(it, j) = 2 * (I(it) / (4 * pi * D * r(j))) * erfc(r(j) / sqrt(4 * D * t(it)));
%     end
% end


% % Exponential decay function
% tau = 100;  % Decay constant (seconds)

% % Compute concentration with history of current values
% for it = 1:length(t)
%     for j = 1:length(r)
%         sum_term = 0;
%         for k = 1:it
%             time_diff = t(it) - t(k);
%             if time_diff > 0
%                 sum_term = sum_term + I(k) * exp(-time_diff / tau) * erfc(r(j) / sqrt(4 * D * time_diff));
%             end
%         end
%         C(it, j) = 2 * sum_term / (4 * pi * D * r(j));
%     end
% end

% Calculate concentration for each combination of t and r
for i = 1:length(t)
    conc(C, I, r, t, D, dt, i-1);
end

% for it = 1:length(t)
%     % Loop over each distance
%     for i = 1:length(r)
%         % Sum contributions from all previous time steps
%         for n = 1:it
%             % Time elapsed since time step n
%             elapsed_t = (it - n + 1) * dt;
%             % Contribution from time step n
%             C(it, i) = C(it, i) + 2 * (I(n) / (4 * pi * D * r(i))) * erfc(r(i) / sqrt(4 * D * elapsed_t));
%         end
%     end
% end

% Call the C function for each time step
% for i = 2:length(t)
%      conc(r, dt, D, C(i,:), I, length(r), i, t(1:i));
% end

% Plot the results
figure;
hold on;
plot(t, C(:,1));
plot(t, C(:,2));
xlabel('Time t');
ylabel('Concentration C');
zlabel('Concentration C');
title('Concentration of Diffusing Substance from a Point Source');
hold off;



% % Parameters
% D = 1e-6;
% dt = 10; % Time step (seconds)
% total_time = 1000;  % Total simulation time (seconds)
% r = 2e-6;
% t = 0:dt:total_time;  % Time steps
% amplitude = 1; % Assuming some amplitude
% frequency = 0.1; % Assuming some frequency
% I = @(t) amplitude * sin(2 * pi * frequency * t);
% % Preallocate concentration array
% C = zeros(length(t), length(r));
% % Precompute weights
% weights = @(delta_t, r, D) exp(-r.^2 / (4 * D * delta_t)) ./ (4 * pi * D * delta_t).^(3/2);
% % Calculate concentration for each combination of t and r
% for i = 1:length(t)
%     for j = 1:length(r)
%         % Integrate weighted contributions from all past time steps
%         for k = 1:i
%             tau = t(k);
%             delta_t = t(i) - tau;
%             if delta_t > 0
%                 weight = weights(delta_t, r(j), D);
%                 C(i, j) = C(i, j) + I(tau) * weight * dt;
%             end
%         end
%     end
% end
% % Plot the results
% figure;
% hold on;
% plot(t, C(:,1));
% % plot(t, C(:,2));
% xlabel('Time t');
% ylabel('Concentration C');
% zlabel('Concentration C');
% title('Concentration of Diffusing Substance from a Point Source');
% hold off;


% % Parameters
% D = 1e-6;
% dt = 10; % Time step (seconds)
% total_time = 1000;  % Total simulation time (seconds)
% r = 2e-6;
% t = 0:dt:total_time;  % Time steps
% amplitude = 1; % Assuming some amplitude
% frequency = 0.1; % Assuming some frequency
% I = @(t) amplitude * sin(2 * pi * frequency * t);
% 
% % Preallocate concentration array
% C = zeros(length(t), length(r));
% 
% % Calculate concentration for each combination of t and r
% for i = 1:length(t)
%     for j = 1:length(r)
%         % Integrate contributions from all past time steps
%         for k = 1:i
%             tau = t(k);
%             delta_t = t(i) - tau;
%             if delta_t > 0
%                 weight = exp(-r(j)^2 / (4 * D * delta_t)) / (4 * pi * D * delta_t)^(3/2);
%                 C(i, j) = C(i, j) + I(tau) * weight * dt;
%             end
%         end
%     end
% end
% 
% % Plot the results
% figure;
% hold on;
% plot(t, C(:,1));
% % plot(t, C(:,2));
% xlabel('Time t');
% ylabel('Concentration C');
% zlabel('Concentration C');
% title('Concentration of Diffusing Substance from a Point Source');
% hold off;