function [ dq, dc, dw, dc_proton, Nqk] = NTdynamics( q, c, w, M, opts, dt, C_vesicles)
arguments
    q % vesicule has transmitter (1) or not (0); 1xM array
    c 
    w
    M % num_vesicles
    opts
    dt % time step
    C_vesicles % calcium concentration of vesicle
end

C_vesicles = C_vesicles(:);

%% Transmitter Release and Recycling

% k1 ... single spot release rate
% k1 = kmax .* Hill_Langmuir_A(C_vesicles, n_Hill, KA_Hill);
k1 = RibbonSynapse_v4.TransmitterRelease(C_vesicles, opts.transmitter_release_parameters{:});

q_old = q;

% Nqk = number released from vesicle (random) <- dq : RoC of release
% Nwx = number reprocessed <- dc RoC of repocess
% NMqt = number created <- dw RoC of creation

% dc_proton <- RoC of proton concentration

% k1*dt is the amount of neurotransmitters release in the timeframe
Nqk = NTTransport(q, k1*dt, M); % number released from vesicle (random)
q = q - Nqk;

% reprocessed
Nwx = NTReprocessing(~q, floor(w) * opts.x.Hz * dt, M); % number reprocessed
q = q + Nwx;

% manufactured
NMqy = NTManufacture(~q, opts.y.Hz*dt, M);
q = q + NMqy;

dq = (q - q_old) / dt;
dc = sum(Nqk) / dt - opts.l.Hz .* c - opts.r.Hz .* c;
dw = opts.r.Hz .* c - sum(Nwx) / dt;

dc_proton = sum(Nqk)/dt - 0.9 * opts.l.Hz .* c_proton;

% https://github.com/evanwporter/cochlea-nerve/blob/4e8b4f18f20782bfc39d88589db9f4e04dbcf507/IHC/Transduction_v4.m#L147-L148

%% Stochastic NT Transport

% n: # vesicles

function N = NTReprocessing(q, rho, n)
    n_empty = sum(q);
    try
        if n_empty > 0 && rho > 0
            N = q & (rand(n,1) < (rho / n_empty));
        else
            N = zeros(size(q));
        end
    catch exception
        disp("n_empty");disp(n_empty)
        disp("rho");disp(rho)
        disp("q");disp(q)
        disp("n");disp(n)
        throw(exception)
    end
end

function N = NTManufacture(q, rho, n)
    n_empty = sum(q);
    if n_empty > 0
        try     
            N = q & (rand(n,1) < (rho / n_empty));
        catch exception
            disp("n_empty");disp(n_empty)
            disp("rho");disp(rho)
            disp("q");disp(q)
            disp("n");disp(n)
            throw(exception)
        end
    else
        N = zeros(size(q));
    end
end

function N = NTTransport(q, rho, n)
    n_occupied = sum(q); % number occupied
    % disp("n_occupied"); disp(n_occupied)
    if n_occupied > 0
        try
            N = q & (rand(n,1) < rho);
        catch exception
            disp(q);
            disp(rho);
            disp(n)
            throw(exception)
        end
    else
        N = zeros(size(q));
    end
end


end
