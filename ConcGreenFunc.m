%% Parameters
t_total = 1000; % total time to compute the concentration
dt = 10; % time step size in seconds
D = 5.2e-10; % Diffusion coefficient (m^2/s)
r = [1e-6, 5e-6]; % Distance from the source
nr = 2;

frequency = 1; % Frequency in Hz
amplitude = 40e-12; % Amplitude of the source term

I = @(t) amplitude * abs(sin(2 .* pi .* frequency .* t));

it = 0:dt:t_total;
nt = length(it);

%% Precompute Green's Function
% G_array has 3 dim
% dim 1 : distance (nr)
% dim 2 : current time dim (nt)
% dim 3 : prev time dim (t')
%   forall <= dim 2 == 0

% G_array = zeros(length(r), nt, nt);
% for i = 1:length(r)
%     for t = 1:nt
%         for t_prime = 1:t-1
%             G_array(i, t, t_prime) = greens_function(it(t), it(t_prime), r(i), D);
%         end
%     end
% end

G_array = PrecomputeGreen(r, it, D);

%% Compute Concentration
C = zeros(length(r), nt);

% for i = 1:nt
%     t = it(i);
%     C(:,i) = approximate_solution(t, r, I, dt, G_array, it);
% end

I_array = I(it);
for i = 1:nt
    t = it(i);
    C(:,i) = CalcG(i, nr, I_array, dt, G_array, nt);
end

%% Display Results
figure;
hold on;
for i = 1:length(r)
    plot(it, C(i, :), 'DisplayName', sprintf('Distance = %.1f um', r(i) * 1e6));
end
xlabel('Time (s)');
ylabel('Concentration (mol/m^3)');
legend;
title('Calcium Concentration Over Time at Different Distances');
hold off;


function C = approximate_solution(t, r, I, dt, G_array, it)
    % t : time at which to evaluate the concentration
    % r : array of distances from the source
    % D : diffusion coefficient
    % I : function handle for the source term I(t')
    % dt : time step for the approximation
    % G_array : precomputed Green's function array
    % it : time array
    
    % Number of steps
    ns = floor(t / dt);

    nr = length(r);

    C = zeros(1, nr);

    for i = 1:nr % distance
        for j = 1:ns % time
            t_prime = it(j);
            G = G_array(i, ns+1, j);
            C(i) = C(i) + I(t_prime) * G * dt;
        end
    end
end
