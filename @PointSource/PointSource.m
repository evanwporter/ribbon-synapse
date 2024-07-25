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

        method
        
        N (1,1) double = +inf; % max history
        
        u1_const
        
        rel_tol

        test_var
        time_array % evan
        G_array

    end
    properties % (Transient) % transient property acces is slow (R2022a)
        e_pre
        e_pre_rev

        lastopen
        
        current (1, :) double;
        concentration

        E
    end
    
    methods
        function obj = PointSource(d, D, r, dt, it, args)
            arguments
                d  (1,1) double {mustBePositive, mustBeInteger, mustBeLessThanOrEqual(d, 3)}
                D  (1,1) double {mustBePositive}
                r  (:,1) double {mustBePositive}

                dt (1,1) double {mustBeNonnegative} = 0
                it (1,:) double {mustBeNonnegative} = 0

                args.method string = "ondrej"; % Other methods are evan and ode15s

                args.max_history (1,1) double {mustBePositive} = +inf
                args.rel_tol (1,1) double {mustBeNonnegative} = 0
                args.log_concentration (1,1) logical = false

                % args.tspan (1, 2) double{mustBeNonnegative} = 0
            end
            %POINTSOURCE
            
            obj.d = d;
            obj.D = D;
            obj.r = r;
            
            obj.nr = numel(r);

            obj.method = args.method;

            if not(obj.method == "ode15s")
                    
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
                
                if args.log_concentration
                    obj.concentration = zeros(obj.nr, obj.nt);
                end

                obj.current = zeros(1, obj.nt);

            end

            obj.G_array = PrecomputeGreen(reshape(obj.r,1,[]), it, D);

            % obj.test_var = (2 / (4 * pi * obj.D .* obj.r)) * erfc(obj.r ./ sqrt(4 * obj.D .* it));
        end
        % -----------------------------------------------------------------
        function N = find_max_history(obj, tol)
            
            tdim = 2;
            
            err_est = cumsum(obj.e_pre_rev, tdim);

            % Gives the fraction of the total cumulative effect up to each time step
            rel_err_est = err_est ./ err_est(:,end);
            
            ind = zeros(obj.nr, 1);
            for i = 1:obj.nr
                ind(i) = find(rel_err_est(i,:) > tol, 1) - 1;
            end
            assert(all(ind > 0));
            N = obj.nt - min(ind);
            
        end

        c = e_iterate(obj, it)        

        function C = dt_iterate(obj, time_array)
            arguments
                obj
                time_array (1, :) double {mustBeNonnegative}
            end
            % t : time at which to evaluate the concentration
            % r : array of distances from the source
            % D : diffusion coefficient
            % I : function handle for the source term I(t')
            % dt : time step for the approximation

            C = zeros(1, obj.nr);
            
            t = time_array(end);
            
            GF = @(t, t_prime, r) (4 * pi * obj.D * (t - t_prime))^(-3/2) * exp(-r^2 / (4 * obj.D * (t - t_prime)));
            
            for i = 1:obj.nr
                for j = 2:(length(time_array) - 1)
                    t_prime = time_array(j);
                    dt = t_prime - time_array(j - 1);
                    G = GF(t, t_prime, obj.r(i));
                    C(i) = C(i) + obj.current(j) * G * dt;
                end
            end
        end

        function c = simple_iterate(obj, it)
            % Evan's Simple version of iterate
            % Same as PointSource.iterate when PointSource.N == 1

            assert(obj.d == 3, "Requires dimensionality of 3")
            c = zeros(obj.nr, 1);

            current_t = it * obj.dt;

            for i = 1:length(obj.r)
                c(i) = 2 * (obj.current(it) / (4 * pi * obj.D * obj.r(i))) * erfc(obj.r(i) / sqrt(4 * obj.D * current_t));
                % c(i) = obj.test_var(it) .* obj.current(it);
            end

            % Precompute the following:
            % 2 / (4 * pi * obj.D .* obj.r) * erfc(obj.r ./ sqrt(4 * obj.D .* it))
        end

        c = iterate(obj, it)

        V = u_ondra(obj, J, E, n, nn, N, nt, nr)

        % -----------------------------------------------------------------

        function val = e(obj, m)
            val = obj.u2((m+1) * obj.dt) - obj.u2(m * obj.dt);
        end

        function u2 = u2(obj, t)
            assert(obj.d == 3, "Dimensionality must be 3 to use u2.")
            u2 = (obj.r.^(2-obj.d) / (4*obj.D*pi^(obj.d/2))) .* ... 
                 (sqrt(pi) * erfc(obj.r ./ sqrt(4*obj.D*t)));
        end
    end
end

