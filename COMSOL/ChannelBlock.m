function channel_state = ChannelBlock(current_state, varargin)
persistent q w c_proton;
persistent channels;

transmitter_release_parameters = {500; 5; 1e-4};
l = 1290; % Hz
tau_blocked = 1e-3;

nVesicles = 16;
nChannels = 72;

if isempty(q)
    q = ones(nVesicles, 1);
    w = 0;
    c_proton = 0;

    % These parameters are from SynapseOptions.channels
    channels = Channels(72, 5e-4, tau_blocked);
    channels.alpha = -133;
    channels.kp0 = 500000;
    channels.V0t = -0.0239;
    channels.S0t = 0.0079;
end

dt = 1e-4;

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

p_block = CaV13.CaProtonBlock(c_proton);

channel_state = 0;

for ii = 1:nChannels

    if channels.state(ii) == 'i' % inactivated
        continue
    end

    % Unblock channel
    if channels.state(ii) == 'b' % blocked
        if t >= channels.tblocked(ii)
            channels.state(ii) = 'c';
        else
            continue
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

        channels.state(ii) = 'b';
        channels.tblocked(ii) = t + tau;
    else
        channel_state = current_state;
    end
    
end

