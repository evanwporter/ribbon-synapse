classdef PointSource < handle
    %POINTSOURCE

    % Evan's Notes
    % d : dimensionality
    % D : diffusion coefficient
    % r : vector of distances from the point source to the point of interest 
    %   (where concentration needs to be calculated)

    % nr : # distances
    % nt : # time intervals

    
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

        time_array (1,:) double {mustBeNonnegative, mustBeInteger} = [0]
    end
    properties % (Transient) % transient property acces is slow (R2022a)
        e_pre
        e_pre_rev

        lastopen
        
        current
        concentration
    end
    
    methods
        function obj = PointSource(d, D, r)
            arguments
                d  (1,1) double {mustBePositive, mustBeInteger, mustBeLessThanOrEqual(d, 3)}
                D  (1,1) double {mustBePositive}
                r  (:,1) double {mustBePositive}
            end
            %POINTSOURCE
            
            obj.d = d;
            obj.D = D;
            obj.r = r;
            
            obj.nr = numel(r);
            obj.nt = 0;
                                    
            %%
            
            obj.u1_const = pi^(-d/2)/2/D * r.^(2-d);
                        
            obj.current = [];
        end
        % -----------------------------------------------------------------
        
        function c = iterate(obj, t)
            arguments
                obj
                t (1, 1) double       % current time
            end

            c = zeros(obj.nr, 1);

            obj.time_array = [obj.time_array, t];
            obj.nt = obj.nt + 1;

            e_pre_rev = obj.e(flip(obj.time_array)); % reversed order

            cc = zeros(obj.nr,1);
            u_ondra_mex(obj.current, e_pre_rev, t, t, obj.nt, obj.nt, cc);

            c(:, 1) = cc;

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
        
        % -----------------------------------------------------------------
        
        function val = u_orig(r, chi, j, n, N, dt, d, D)
            arguments
                r (1,1) double      % distance from the channel
                chi 
                j double
                n (1,1) double      % number of time intervals
                N (1,1) double      % max number of time intervals to use in the computation
                dt (1,1) double     % time step
                d (1,1) double      % channel space setup parameter
                D (1,1) double      % diffusion coefficient
            end
            %
            % exact for N := n but very comp. demanding. N = 100 usually good approx.
            %

            % -- Original (without approx.) --
            for m = 0:n-1 % min(n-1,N-1)
                val = val + chi(n-m) * j(n-m) * e(m,r,dt,d,D);
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

