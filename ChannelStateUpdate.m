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