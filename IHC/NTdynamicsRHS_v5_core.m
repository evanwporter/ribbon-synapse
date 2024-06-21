function [ dq, dc, dw, dc_proton, vesicles] = NTdynamicsRHS_v5_core( t, ...
    q, c, w, c_proton, ...
    vesicles, y, l, x, r, ...
    transmitter_release_parameters, ...
    dt, C_vesicles)
arguments
    t
    q % vesicule has transmitter (1) or not (0); 1xM array
    c 
    w
    c_proton
    vesicles
    y
    l
    x
    r
    transmitter_release_parameters
    dt % time step
    C_vesicles % calcium concentration of vesicle
end
%TRANSDUCTIONRHS

C_vesicles = C_vesicles(:);

M = vesicles.num;

% disp("W"); disp(w);
% disp("c"); disp(c)

%% Transmitter Release and Recycling

kmax = transmitter_release_parameters{1};
% n_Hill = transmitter_release_parameters{2};
% KA_Hill = transmitter_release_parameters{3};

% assert(kmax*dt <= 1)

% k1 ... single spot release rate
% k1 = kmax .* Hill_Langmuir_A(C_vesicles, n_Hill, KA_Hill);
k1 = RibbonSynapse_v4.TransmitterRelease(C_vesicles, transmitter_release_parameters{:});

q_old = q;

% Nqk = number released from vesicle (random) <- dq : RoC of release
% Nwx = number reprocessed <- dc RoC of repocess
% NMqt = number created <- dw RoC of creation

% dc_proton <- RoC of proton concentration

% k1*dt is the amount of neurotransmitters release in the timeframe
Nqk = NTTransport(q, k1*dt, M); % number released from vesicle (random)
q = q - Nqk;

% disp("k1*dt"); disp(k1*dt)
% disp("q");disp(q)

% reprocessed
Nwx = NTReprocessing(~q, floor(w)*x*dt, M); % number reprocessed
q = q + Nwx;

% disp("# reproc"); disp(Nwx);

% manufactured
NMqy = NTManufacture(~q, y*dt, M);
q = q + NMqy;

% disp(q); disp(q_old); disp(dt);

% dq = (Nwx + NMqy - Nqk)/dt;
dq = (q - q_old) / dt;
dc = sum(Nqk) / dt - l .* c - r .* c;
dw = r .* c - sum(Nwx) / dt;

dc_proton = sum(Nqk)/dt - 0.9*l.*c_proton;

for i = 1:sum(Nqk)
    vesicles.logReleaseEvent(t);
end

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
        N = q & (rand(n,1) < rho);
    else
        N = zeros(size(q));
    end
end


end
