classdef PointSource < handle
    %POINTSOURCE

    % Evan's Notes
    % d : dimensionality
    % D : diffusion coefficient
    % r : vector of distances from the point source to the point of interest

    
    properties
        d
        D
        r
        
        nr
        nt
        
        dt
        
        N (1,1) double = +inf; % max history
        
        u1_const
        
        rel_tol
    end
    properties % (Transient) % transient property acces is slow (R2022a)
        e_pre
        e_pre_rev

        lastopen
        
        current
        concentration
    end
    
    methods
        function obj = PointSource(d, D, r, dt, it, args)
            arguments
                d  (1,1) double {mustBePositive, mustBeInteger, mustBeLessThanOrEqual(d, 3)}
                D  (1,1) double {mustBePositive}
                r  (:,1) double {mustBePositive}
                dt (1,1) double {mustBePositive}
                it (1,:) double {mustBeNonnegative, mustBeInteger}

                args.max_history (1,1) double {mustBePositive} = +inf
                args.rel_tol (1,1) double {mustBeNonnegative} = 0
                args.log_concentration (1,1) logical = false
            end
            %POINTSOURCE
            
            obj.d = d;
            obj.D = D;
            obj.r = r;
            
            obj.nr = numel(r);
            obj.nt = numel(it);

            obj.lastopen = -inf;
            
            obj.dt = dt;
            
            obj.rel_tol = args.rel_tol;

            if args.rel_tol > 0 && args.max_history < +inf
                error('only one argument can be set to non-default value')
            end
            
            if args.rel_tol > 0 || args.max_history < +inf
                approximation = true;
            else
                approximation = false;
            end
            
            %%
            
            obj.u1_const = pi^(-d/2)/2/D * r.^(2-d);
            
            % obj.e_pre = obj.e(it);
            obj.e_pre_rev = obj.e(flip(it)); % reversed order
            
            
            %%
            
            if args.rel_tol > 0
                N = obj.find_max_history(args.rel_tol);
                % fprintf('History of %d steps required for %g rel error.\n', N, args.rel_tol)
            elseif isfinite(args.max_history)
                N = args.max_history;
            else
                N = obj.nt;
            end
            
            mustBeInteger(N)
            obj.N = N;
            
            tdim = 2;
            
            debug = false;
            
            if approximation == true
                approx_error_est = sum(obj.e_pre_rev(:,1:obj.nt - N), tdim) ./ sum(obj.e_pre_rev, tdim);

                if debug == true
                    dbprintf('Approximation error estimate\n');
                    for i = 1:obj.nr
                        dbprintf('%8.2f%%   ... for r = %g\n', approx_error_est(i) * 100, r(i));
                    end
                end
                
                assert(all(approx_error_est < obj.rel_tol), sprintf('Error estimate > %g%%', obj.rel_tol * 100));

                itt = flip(it(1:N));
                obj.e_pre_rev = obj.e(itt); % reversed order

            end
            
            obj.current = zeros(1, obj.nt);

            if args.log_concentration
                obj.concentration = zeros(obj.nr, obj.nt);
            end
        end
        % -----------------------------------------------------------------
        function N = find_max_history(obj, tol)
            
            tdim = 2;
            
            err_est = cumsum(obj.e_pre_rev, tdim);
            rel_err_est = err_est ./ err_est(:,end);
            
            ind = zeros(obj.nr, 1);
            for i = 1:obj.nr
                ind(i) = find(rel_err_est(i,:) > tol, 1) - 1;
            end
            assert(all(ind > 0));
            N = obj.nt - min(ind);
            
        end
        % -----------------------------------------------------------------
        
        function c = iterate(obj, it)
            arguments
                obj
                it (1,1) double    
            end
        
            c = zeros(obj.nr, 1);
            current_t = it * obj.dt;
            
            % for i = 1:obj.nr
            %     current_r = obj.r(i);
            %     c(i) = (1 / (4 * pi * obj.D * current_t)^(obj.d / 2)) * exp(-current_r^2 / (4 * obj.D * current_t));
            % end

        
            for i = numel(it)
                last_it_open = obj.lastopen / obj.dt;
                if last_it_open < it(i) - obj.N
                    continue
                end
                n = it(i);
            
                cc = zeros(obj.nr, 1);
                u_ondra_mex(obj.current, obj.e_pre_rev, n, n, obj.N, obj.nt, cc);
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

        function V = u_ondra(obj, J, E, n, nn, N, nt, nr)
            % Evan's MATLAB version of the u_ondra.c
            % For debugging purposes

            mm = min(n, N);
            m = N - mm;
            k = nn - mm;
        
            V = zeros(nr, 1);

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
        
        function val = u(obj, j, n, nn)
            arguments
                obj
                j double
                n double
                nn double
            end

            % -- Original (without approx.) --
            % for m = 0:n-1
            %     val = val + chi(n-m) * j(n-m) * e(m,r,dt,d,D);
            % end

            % -- Precalc v1 --
            % val = 0;
            % for m = 0:min(n-1,N-1)
            %     val = val + j(n-m) * e_pre(m+1);
            % end
            
            % exact for N := n but very comp. demanding. N = 100 usually good approx.
            
            mm = min(n, obj.N);
            
            % m = obj.nt - mm + 1; % v1
            m = obj.N - mm + 1; % v2
            
            % nn = numel(j);
            
            tdim = 2;

            % je = j(nn:-1:(nn-mm+1)) .* obj.e_pre(:,1:mm);
            % je = flip(je, tdim);
            % val = sum(je, tdim);
            
            % we sum in the reverse order because of finite computer
            % precision; the array e_prev_rev is increasing.
            
            % je = j((nn-mm+1):nn) .* obj.e_pre_rev(:,m:obj.nt); % v1
            je = j((nn-mm+1):nn) .* obj.e_pre_rev(:,m:obj.N); % v2

            val = sum(je, tdim);
            
            % val2 = zeros(size(obj.e_pre_rev, 1), 1);
            % for k = 1:size(obj.e_pre_rev, 1)
            %     val2(k) = obj.e_pre_rev(k,m:obj.nt)*j((nn-mm+1):nn)';
            % end
            % assert(norm(val-val2)/(1 + norm(val)) < 100*eps);

        end

        function val = u_mex(obj, j, n, nn)
            arguments
                obj
                j double
                n double
                nn double
            end

            val = u_mex(j, n, nn, obj.N, obj.nt, obj.e_pre_rev);
            
        end
        % -----------------------------------------------------------------

        function val = e(obj, m)
            % 
            val = obj.u1((m+1) * obj.dt) - obj.u1(m * obj.dt);

        end

        % -----------------------------------------------------------------
        
        function val = u1(obj, t)
            % Equation 8 of the paper

            % val = pi^(-d/2)/2/D * r^(2-d) * obj.U(r./sqrt(4*D*t), d);
            val = obj.u1_const .* obj.U(obj.r ./ sqrt(4*obj.D*t));

        end
        
        % -----------------------------------------------------------------

        function val = U(obj, alpha)
            % Equation 9 of the paper

            if obj.d == 1
                % form entire membrane
                val = exp(1).^(-alpha.^2) ./ alpha - sqrt(pi) * erfc(alpha);
            elseif obj.d == 2
                % space between cells
                error('Check gamma definition');
                val = gamma(0, alpha.^2) / 2;
            elseif obj.d == 3
                % radially from a channel
                val = sqrt(pi) * erfc(alpha) / 2;    
            else
                error('Value d = %d not valid', obj.d)
            end
        end
        
        % -----------------------------------------------------------------
    end
end

