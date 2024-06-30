function c = iterate(obj, it)
    arguments
        obj
        it (1,1) double    
    end

    c = zeros(obj.nr, 1);
    current_t = it * obj.dt;

    for i = numel(it)
        last_it_open = obj.lastopen / obj.dt;
        if last_it_open < it(i) - obj.N
            continue
        end
        n = it(i);
    
        cc = zeros(obj.nr, 1);
        u_ondra_mex(obj.current, obj.e_pre_rev, n, n, obj.N, obj.nt, cc);
        % cc = obj.u_ondra(obj.current, obj.e_pre_rev, n, n, obj.N, obj.nt, obj.nr);
        c(:, i) = cc;
    
        % Error prevention
        if any(isinf(c(:, 1)) | isnan(c(:, 1)))
            disp(obj.u_ondra(obj.current, obj.e_pre_rev, n, n, obj.N, obj.nt, obj.nr));
            fprintf('Iteration: %d, Concentration: %f\n', i, c(:, i));
            fprintf('Current: %f, e_pre_rev: %f\n', obj.current, obj.e_pre_rev);
            error('Divergence detected at iteration %d\n', i);
        end
    end
    
    
    
        % c(:, i) = obj.u(obj.current,n,n);
        % c(:, i) = obj.u_mex(obj.current,n,n);
    
        % cc(:, i) = obj.u(obj.current,n,n);
        % if c(:,i) ~= 0
        %     assert(all(abs((c(:,i)-cc(:,i))./c(:,i)) < 1e-14))
        % end
    
    % From the paper
    % For d = 1 and d = 3, Eq. (1) assumes that the ions are diffusing
    % into open space on both sides of the membrane. However, in the 
    % ion channel situation, we are only interested in the ions 
    % diffusing into the compartment on one side of the membrane. 
    % Therefore, in these two geometries, one must use twice the 
    % channel flux in the equations.
    
    if obj.d == 1 || obj.d == 3                
        % since u is linear with respect to j, we can multiply c
        c = c * 2;
    end
    
end