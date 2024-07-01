%% Parameters
t_total = 1000; % total time to compute the concentration
dt = 10; % time step size in seconds
D = 5.2e-10; % Diffusion coefficient (m^2/s)
r = [1e-6, 5e-6]; % Distance from the source

frequency = 1; % Frequency in Hz
amplitude = 1e-9; % Amplitude of the source term

I = @(t) amplitude * abs(sin(2 * pi * frequency * t));

it = 0:dt:t_total;
nt = length(it);

%% Compute
C = zeros(length(r), nt);

for i = 1:nt
    t = it(i);
    C(:,i) = approximate_solution(t, r, D, I, dt);
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

% C1 = analytical_solution(t_total, r, D, I);
% C2 = concG(t_total, r, D, I, dt);
% C3 = approximate_solution(t_total, r, D, I, dt);
% 
% disp(['Integral Concentration at time t = ' num2str(t) ' and distance r = ' num2str(r) ': ' num2str(C1)]);
% disp(['C Approximation Concentration at time t = ' num2str(t) ' and distance r = ' num2str(r) ' and time step of 10 : ' num2str(C2)]);
% disp(['Approximation Concentration at time t = ' num2str(t) ' and distance r = ' num2str(r) ' and time step of 10 : ' num2str(C3)]);


%% Functions
function C = analytical_solution(t, r, D, I)
    t_prime_min = 0;
    t_prime_max = t;

    integrand = @(t_prime) (I(t_prime) ./ (4 * pi * D * (t - t_prime)).^(3/2)) ...
                            .* exp(-r^2 ./ (4 * D * (t - t_prime)));

    C = integral(@(t_prime) integrand(t_prime), t_prime_min, t_prime_max, 'RelTol', 1e-8, 'AbsTol', 1e-12, 'ArrayValued', true);
end

function G = greens_function(t, t_prime, r, D)
    if t == t_prime
        G = 0; % Green's function is zero when t equals t_prime
    else
        G = (4 * pi * D * (t - t_prime))^(-3/2) * exp(-r^2 / (4 * D * (t - t_prime)));
    end
end

function C = approximate_solution(t, r, D, I, dt)
    % t : time at which to evaluate the concentration
    % r : array of distances from the source
    % D : diffusion coefficient
    % I : function handle for the source term I(t')
    % dt : time step for the approximation
    
    % Number steps 
    ns = floor(t / dt);

    nr = length(r);

    C = zeros(1, nr);

    for i = 1:nr
        for j = 0:ns - 1
            t_prime = j * dt;
            G = greens_function(t, t_prime, r(i), D);
            C(i) = C(i) + I(t_prime) * G * dt;
        end
    end
end
