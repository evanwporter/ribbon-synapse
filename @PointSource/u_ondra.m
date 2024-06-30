function V = u_ondra(obj, J, E, n, nn, N, nt, nr)
    % Evan's MATLAB version of the u_ondra.c
    % For debugging purposes

    %  n : current time index
    % nn : total number of time steps
    %  N : max history length
    % nt : number time steps in history
    % nr : number of distances

    % Effective history length
    % Use the number of timesteps (nt) if N is larger this number
    % Ensures we don't go beyond max history
    mm = min(n, N);

    % Offset index for history length.
    m = N - mm;

    % Offset index for current index
    k = nn - mm;

    % k + j + 1 == it


    for i = 1:nr
        V(i) = 0.0;
        
        for j = 0:(mm - 1)
            ind = i + (m + j) * nr; 
            tmp = V(i) + J(k + j + 1) * E(ind);  
            if isinf(tmp)
                disp("err");
            end
            V(i) = tmp;
        end
    end
end