function channel_state = ChannelBlock(current_state, channel_num, t, varargin)
persistent q w c_proton;
persistent channels;
persistent last_t;

transmitter_release_parameters = {500; 5; 1e-4};
l = 1290; % Hz
tau_blocked = 1e-3;
dt = 1e-4;

nVesicles = 16;
nChannels = 72;

if isempty(q)
    q = ones(nVesicles, 1);
    w = 0;
    c_proton = 0;

    last_t = 0;

    % These parameters are from SynapseOptions.channels
    channels = Channels(72, 5e-4, tau_blocked);
end
    
if last_t ~= t 
    last_t = t;
    
    C_vesicle = cell2mat(varargin);
    
    k = RibbonSynapse_v4.TransmitterRelease(C_vesicle, transmitter_release_parameters{:});
    
    NTTransport = q & (rand(nVesicles,1) < k*dt);
    q = q - NTTransport;
    
    NTReprocessing = ~q & (rand(nVesicles,1) < floor(w) * opts.x.Hz * dt);
    q = q - NTReprocessing;
    
    NTManufacture = ~q & (rand(nVesicles,1) < opts.y.Hz*dt);
    q = q + NTManufacture;
    
    dc_proton = NTTransport * nVesicles / dt - 0.9 * l .* c_proton;
    c_proton = c_proton + dc_proton * dt;
end

p_block = CaV13.CaProtonBlock(c_proton);

channel_state = 0;

if current_state == 0 % inactivated
    return
end

% Unblock channel
if channels.state(channel_num) == 'b' % blocked
    if t >= channels.tblocked(channel_num)
        channels.state(channel_num) = 'c';
    else
        return
    end
end

% Block channel
% Check if channel becomes unblocked
if rand(1, 1) < p_block * dt_ / tau_blocked

    % generate random close time from a gamma distribution
    % such that the mean (= a*b) is equal to the opts.channels.tau_blocked
    %     a ... shape
    %     b ... scale

    a = 4;
    b = tau_blocked / a;
    tau = b * gammaincinv(rand(1,1), a); % gamma random number without toolbox

    channels.state(channel_num) = 'b';
    channels.tblocked(channel_num) = t + tau;
else
    channel_state = current_state;
end
    
end

