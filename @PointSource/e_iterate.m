function C = e_iterate(obj, it) % t, r, D, I, dt)
    % t : time at which to evaluate the concentration
    % r : array of distances from the source
    % D : diffusion coefficient
    % I : function handle for the source term I(t')
    % dt : time step for the approximation
    %
    % C = zeros(1, obj.nr);
    % 
    % t = it * obj.dt;
    % 
    % GF = @(t, t_prime, r, D) (4 * pi * D * (t - t_prime))^(-3/2) * exp(-r^2 / (4 * D * (t - t_prime)));
    % 
    % for i = 1:obj.nr
    %     for j = 0:it - 1
    %         t_prime = j * obj.dt;
    %         G = GF(t, t_prime, obj.r(i), obj.D);
    %         C(i) = C(i) + obj.current(j + 1) * G * obj.dt;
    %     end
    % end

    % C = CalcConc(obj.r, obj.nr, obj.D, obj.dt, obj.current, it);
    C = CalcG(it, obj.nr, obj.current, obj.dt, obj.G_array, obj.nt);
end