function channel_state = ChannelBlock(Cin, current_state)
persistent c_proton
persistent channels tblocked;

transmitter_release_parameters = {500; 5; 1e-4};
l = 1290; % Hz
tau_blocked = 1e-3;

nVesicles = 16;
nChannels = 72;


dt = 1e-4;

k = RibbonSynapse_v4.TransmitterRelease(Cin, transmitter_release_parameters{:});
rho = k*dt;

NTTransport = rand(1,1) < rho;

dc_proton = NTTransport * nVesicles / dt - 0.9 * l .* c_proton;

c_proton = c_proton + dc_proton * dt;

p_block = CaV13.CaProtonBlock(c_proton);

channel_state = 0;

for ii = 1:nChannels

    if channels(ii) == 'i' % inactivated
        continue
    end

    % Unblock channel
    if channels(ii) == 'b' % blocked
        if t >= tblocked(ii)
            channels(ii) = 'c';
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

        channels(ii) = 'b';
        tblocked(ii) = t + tau;
    else
        channel_state = current_state;
    end
    
end

