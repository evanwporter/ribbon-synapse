function [time_vector, channel_state_matrix] = GenerateChannelUpdateMatrix(opts, num_time_points, dt, Vt, to_csv)
arguments
    opts
    num_time_points = 10000;
    dt = 1e-4;
    Vt = -0.04;
    to_csv = true;
end

% initial_channel_states = randi([0, 1], [1, 72]);
initial_channel_states = zeros(1, 72);

channel_state_matrix = zeros(num_time_points, 72);
channel_state_matrix(1, :) = initial_channel_states;

% V_half = -0.03;	
% k = 0.006;
% P_top = 0.35; 	
% P_bottom = 0;	
% boltzmann = @(Voltage) P_bottom + ((P_top - P_bottom) / (1+exp((V_half-Voltage)/k)));

for t = 2:num_time_points
    % if rand(1) < boltzmann(Vt)
    %     channel_state_matrix(t, :) = 1;
    % else
    %     channel_state_matrix(t, :) = 0;
    % end
    prev_channel_states = channel_state_matrix(t-1, :);
    channel_state_matrix(t, :) = ChannelStateUpdate(opts, dt, prev_channel_states, Vt);
end

% disp(channel_state_matrix);

time_vector = (0:num_time_points-1)' * dt;

% channel_state_matrix = [ones(1,72); channel_state_matrix];
% time_vector = [-1; time_vector];
 
% assignin("base", "channel_state_matrix", channel_state_matrix);
% assignin("base", "time_vector", time_vector);

% Write the matrix to a CSV file
if to_csv
    output_matrix = [time_vector, channel_state_matrix];
    header = ['time', arrayfun(@(x) sprintf('channel%d', x), 1:72, 'UniformOutput', false)];
    header_str = strjoin(header, ',');
    csv_filename = 'channel_state_matrix.csv';
    fid = fopen(csv_filename, 'w');
    fprintf(fid, '%s', header_str);
    fclose(fid);
    writematrix(output_matrix, csv_filename, 'WriteMode', 'append');
end

% end

function channel_states = ChannelStateUpdate(opts, dt, prev_channel_states, Vt)
arguments
    opts
    dt = 1e-4;
    prev_channel_states = randi([0, 1], [1, 72]);
    Vt = -0.04
end

channel_states = zeros(size(prev_channel_states));

%% Calculate rates

S0t = opts.channels.S0t;
V0t = opts.channels.V0t;

alpha = opts.channels.alpha;
beta = alpha + 1/S0t;

kp0 = opts.channels.kp0;

km0 = kp0 * exp(V0t * (beta - alpha));


%% Compute conductance

kp = kp0 * exp(-alpha * Vt); % E: Opening rate
km = km0 * exp(-beta * Vt);  % E: Closing rate

Q = [1 - kp*dt, km *dt; 
     kp*dt, 1 - km*dt];

cQ = cumsum(Q,1);

prev_channel_states = prev_channel_states + 1;

for i = 1:length(prev_channel_states)

    % TODO: Implement blocked channels
    r = rand(1);

    state = prev_channel_states(i);
    
    
    new_state = find(r < cQ(:,state), 1);

    channel_states(i) = new_state - 1;

end

% channel_states = randi([0, 1], [1, 72]);

% if sum(channel_states) ~= 0
%     disp(sum(channel_states))
% end


% kp = 1.2685e6;
% km = 20806;
% Vm = 0.007;
% dt = .1;
% prev_channel_states = randi([0, 1], [1, 72]);