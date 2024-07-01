function V = u_ondra(obj, J, E, n, N, nt, nr)
    % Evan's MATLAB version of the u_ondra.c
    % For debugging purposes

    %  n : current time index
    %  N : max history length
    % nt : number time steps in history
    % nr : number of distances

    % Effective history length
    % Use the number of timesteps (nt) if N is larger this number
    % Ensures we don't go beyond max history
    mm = min(n, N);

    V = zeros(1, nr);

    % Most recent J, is multiplied by earliest diffusion values
    % ie: J7 * E4, J8 * E3, J9 * E2, J10 * E1
    for i = 1:nr        
        V(i) = dot(E(i, width(E) - mm + 1:end), J(1:mm));
    end

    
end