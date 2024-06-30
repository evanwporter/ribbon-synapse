classdef PointSource2 < handle
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
        tspan
        dt
    end
    
    methods
        function obj = PointSource2(d, D, r, dt, tspan, args)
            arguments
                d  (1,1) double {mustBePositive, mustBeInteger, mustBeLessThanOrEqual(d, 3)}
                D  (1,1) double {mustBePositive}
                r  (:,1) double {mustBePositive}
                dt (1,1) double {mustBePositive}
                tspan (1,:) double {mustBeNonnegative, mustBeNonempty}

                args.max_history (1,1) double {mustBePositive} = +inf
                args.rel_tol (1,1) double {mustBeNonnegative} = 1e-3
                args.log_concentration (1,1) logical = false
            end
            
            obj.d = d;
            obj.D = D;
            obj.r = r;
            obj.dt = dt;
            
            obj.nr = numel(r);
            obj.tspan = tspan(1); % Initialize with the first time point

            obj.lastopen = -inf;
            
            obj.rel_tol = args.rel_tol;

            if args.rel_tol > 0 && args.max_history < +inf
                error('only one argument can be set to non-default value')
            end
            
            if args.rel_tol > 0 || args.max_history < +inf
                approximation = true;
            else
                approximation = false;
            end
            
            obj.u1_const = pi^(-d/2)/2/D * r.^(2-d);
            
            % obj.e_pre = obj.e(1); 
            obj.e_pre_rev = obj.e(1); % Initialize with the first element
            
            if args.rel_tol > 0
                N = obj.find_max_history(args.rel_tol);
            elseif isfinite(args.max_history)
                N = args.max_history;
            else
                N = 1; % Start with one element
            end
            
            mustBeInteger(N)
            obj.N = N;
            
            if approximation == true
                approx_error_est = sum(obj.e_pre_rev(:,1:N), 2) ./ sum(obj.e_pre_rev, 2);

                assert(all(approx_error_est < obj.rel_tol), sprintf('Error estimate > %g%%', obj.rel_tol * 100));

                obj.e_pre_rev = obj.e(1:N); % Initialize with the first N elements
            end
            
            obj.current = zeros(1, 1); % Initialize with one element

            if args.log_concentration
                obj.concentration = zeros(obj.nr, 1); % Initialize with one element
            end
        end
        % -----------------------------------------------------------------
        function N = find_max_history(obj, tol)
            arguments
                obj
                tol (1,1) double {mustBeNonnegative}
            end

            tdim = 2;
            
            err_est = cumsum(obj.e_pre_rev, tdim);
            rel_err_est = err_est ./ err_est(:,end);
            
            ind = zeros(obj.nr, 1);
            for i = 1:obj.nr
                % Find first index where relative error estimate exceeds the tolerance
                %   Need to subtract 1 to get index
                ind_val = find(rel_err_est(i,:) > tol, 1); 
                if isempty(ind_val) || ind_val <= 2
                    ind(i) = numel(obj.tspan); % If not found, use the maximum history
                else
                    ind(i) = ind_val - 1;
                end
            end
            
            % Make sure we have an indice
            assert(all(ind > 0), 'Relative error estimate failed to exceed tolerance in all cases.');
            N = numel(obj.tspan) - min(ind);
        end
        % -----------------------------------------------------------------
        
        function c = iterate(obj, t)
            arguments
                obj
                t (1,1) double    
            end
        
            c = zeros(obj.nr, 1);
            
            current_t = t;
            dt = current_t - obj.tspan(end);

            obj.tspan(end + 1) = current_t; % Append new time to tspan

            % Update e_pre_rev dynamically
            obj.e_pre_rev = [obj.e(t) obj.e_pre_rev]; 
            
            if obj.rel_tol > 0
                N = obj.find_max_history(obj.rel_tol);
            else
                N = obj.N;
            end

            N = 1;
            
            % last_it_open = obj.lastopen / dt;
            % if last_it_open < numel(obj.tspan) - N
            %     return
            % end
            n = numel(obj.tspan);
        
            % u_ondra_mex(obj.current, obj.e_pre_rev, n, n, N, numel(obj.tspan), c);
            cc = obj.u_ondra(obj.current, obj.e_pre_rev, n, n, N, numel(obj.tspan), obj.nr);

            c = cc            
            % Error prevention
            if any(isinf(c(:, 1)) | isnan(c(:, 1)))
                % disp(obj.u_ondra(obj.current, obj.e_pre_rev, n, n, N, numel(obj.tspan), obj.nr));
                fprintf('Iteration: %d, Concentration: %f\n', numel(obj.tspan), c(:, end));
                fprintf('Current: %f, e_pre_rev: %f\n', obj.current, obj.e_pre_rev);
                error('Divergence detected at iteration %d\n', i);
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

        function V = U_ONDRA(J, E, n, nn, N, nt, nr)
            % Function to compute the weighted sum of J and E and update V.
            %
            % Inputs:
            %   J  - Array of current values
            %   E  - Array of precomputed values
            %   n  - Current time index
            %   nn - Total number of time steps
            %   N  - Maximum history length
            %   nt - Number of time steps in the history
            %   nr - Number of spatial points (distances r)
            %
            % Output:
            %   V  - Updated array based on weighted sum
        
            mm = min(n, N);
            m = N - mm; % Adjusted start index for history length
            k = nn - mm;
        
            V = zeros(nr, 1);
        
            for i = 1:nr
                V(i) = 0.0;
                for j = 0:(mm - 1)
                    ind = i + (m + j) * nr;
                    if (k + j + 1) > 0 && (k + j + 1) <= numel(J)
                        V(i) = V(i) + J(k + j + 1) * E(ind);
                    else
                        warning('Index out of bounds: k + j + 1 = %d', k + j + 1);
                    end
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

            val = u_mex(j, n, nn, obj.N, numel(obj.tspan), obj.e_pre_rev);
            
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

