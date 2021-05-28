classdef spacecraft
    %spacecraft generalized spacecraft class
    %   Ideally can be used for launch vehicles and satellites
    
    properties
        p       % constants for dynamics, etc
        rv      % INITIAL position and velocity of vehicle
        oe      % INITIAL orbital elements of vehicle
        state_vector % udpated in integration
        Ic      % Inertia matrix of sysem, may change with stage seps.
        n = 13; % state vector dimension
    end
    
    methods
        function self = spacecraft(state, type, params)
            %spacecraft Constructor for this class
            %  state must be 1) [r;v], 2) an oe struct, 3) a full 13x1
            %  state vector [m r v mrp omega]
            
            if ~exist('params', 'var')
                self.p = struct();
            else
                self.p = params;
            end
            
            % Check for required value
            if ~isfield(self.p, 'mu')
                self.p.mu = 3.986004415e+05; %km3s?2
            end
            if ~isfield(self.p, 'mu_m')
                self.p.mu_m = 3.986004418E14; %m3s?2
            end
            if ~isfield(self.p, 'j2')
                self.p.j2 = 0.0010826269;
            end
            if ~isfield(self.p, 'Re')
                self.p.Re = 6371; % km
            end
            if ~isfield(self.p, 'control_flag')
                self.p.control_flag = 0;
            end
            if ~isfield(self.p, 'perturb_flag')
                self.p.perturb_flag = 0;
            end
            
            
            % compute oe2rv and state vector as needed
            if strcmp(type, 'oe')
                p_out = schaub_elements(state, self.p,0,2,1);
                self.rv = [p_out.r; p_out.v];
                self.state_vector = zeros(13,1);
                self.state_vector(2:7) = self.rv;
                % computes more interesting values and replaces...
                p_out = schaub_elements(self.rv,self.p,0,1,2);
                self.oe = state;
                self.oe.h = p_out.h;
                self.oe.E = p_out.E;
                self.oe.nu= p_out.nu;
                self.oe.n = p_out.n;
            elseif strcmp(type, 'rv')
                self.oe = schaub_elements(state,self.p,0,1,2);
                self.rv = state;
                self.state_vector = zeros(13,1);
                self.state_vector(2:7) = self.rv;
            elseif strcmp(type, 'statevector')
                self.rv = state(2:7);
                self.oe = schaub_elements(self.rv,self.p,0,1,2);
                self.state_vector = state;
            end
            
        end
        
        %% helper functions, etc
        function update_inertia(self, Ic_new)
            % other operations here.
            self.Ic = Ic_new*(self.p.mu/self.p.mu);
        end
        
        % just incase the internal state needs to be updated
        function update_state(self, state, type)
            if strcmp(type, 'oe')
                self.oe = state;
            elseif strcmp(type, 'rv')
                self.rv = state;
            end
        end
        
        % J2 Effect from oblate spheroid
        function a_J2 = J2(self, x_N)
            r1 = norm(x_N);
            j2_const = (-3/2)* ...
                ((self.p.mu*self.p.j2*(self.p.Re)^2)/(r1^4));
            x = x_N(1,:);
            y = x_N(2,:);
            z = x_N(3,:);
            a_J2 = j2_const* ...
               [(1-(5*((z^2)/(r1^2))))*(x/r1); ...
                (1-(5*((z^2)/(r1^2))))*(y/r1); ...
                (3-(5*((z^2)/(r1^2))))*(z/r1)];
            if r1-self.p.Re < 0
                a_J2 = [0 0 0].';
            end
        end
        
        function a_exo = atmo_drag_a(self, x_N, v_N, sigma_BN, m)
            % constants
            omega_E = [0 0 7.29211505392569e-05].';
                % equatorial rotation of earth
            H = 7.250;                      % km
            rho_naught = 1.225 * (1000^3);  % kg/km^3
            A = 1/(1000^2);                 % km^2
            
            % make this dependant upon the MRP at some point 
            alpha = pi/2;
            r = norm(x_N);
            vatm = v_N - skew((self.p.Re/r)*omega_E)*x_N;
            % taking care if singular cases
            if r-self.p.Re >= 0
                rho = rho_naught*exp(-(r-self.p.Re)/H);
            else
                rho = 0;
            end
            Cd = 2*(sin(alpha)^3);
            v = norm(vatm);
            v_hat = vatm./v;
            q = (rho*v*v)/2;
            a_exo = -(q*Cd*A*v_hat)/m;
        end
        
        function tau_exo = atmo_drag_tau(self)
            tau_exo = [0 0 0].'.*self.p.mu;
        end
        
        % Gravity Gradient Torque, in body frame
        function tau_exo = grav_grad_tau(self, x_N, sigma_BN)
            BN = MRP2C(sigma_BN);
            x_B = BN*(x_N.*1000); % must convert to meters
            rc = norm(x_B);
            tau_exo = (((3*self.p.mu_m)/(rc^5))* ...
                (skew(x_B)*(self.Ic*x_B)));
        end
        
        % find the perturbing accelerations and torques for post procesing
        function p_exo = perturbing(self, time, state)
            p_exo.tau = zeros(3,length(time));
            p_exo.a   = zeros(3, length(time));
            if self.p.perturb_flag
                % Takes a time history of rototranslational state and time and
                % converts to a time history of expected perturbations
                % instantaneously at these times
                
                for i = 1: length(time)
                    state_in = state(:,i);
                    m =         state_in(1);
                    x_N =       state_in(2:4);
                    v_N =       state_in(5:7);
                    sigma_BN =  state_in(8:10);
                    p_exo.tau(:,i) = grav_grad_tau(self, x_N, sigma_BN);
                    p_exo.a(:,i) = J2(self, x_N) ...
                        + atmo_drag_a(self, x_N, v_N, sigma_BN, m);
                end
            else
                return
            end
        end
        
        % Gravity Gradient Torque in body frame, accelerations are in the
        % the inertial frame
        function u = compute_control(self, t, state_in)
            m =         state_in(1);
            x_N =       state_in(2:4);
            v_N =       state_in(5:7);
            sigma_BN =  state_in(8:10);
            omega_BN =  state_in(11:13);
            
            % gains for detumble law
            P = max(diag(self.Ic).*(2/120));
            K = (P^2)./self.Ic(2,2);
            
            % Zanetti guidance...
            
            % body frame torques and inertial frame accelerations
            % ACS solution is naive, just quick and dirty
            u.tau = -K*sigma_BN - P*omega_BN;
            u.a = NH*HR*[ur ut uz].';
            u.mdot = 0;
        end
        
        %% Dynamics function: can activate controls or peturbations
        function f_x = dynamics(self, t, state_in)
            % distribute the states for use! these are row vectors.
            m =         state_in(1);
            x_N =       state_in(2:4);
            v_N =       state_in(5:7);
            sigma_BN =  state_in(8:10);
            omega_BN =  state_in(11:13);
            
            % Should I make the control compute an external function?
            if self.p.control_flag
                u = compute_control_formation(self, t, state_in);
            else
                u.tau = [0 0 0].';
                u.a = [0 0 0].';
                u.mdot = 0;
            end
            
            % Check for pertubation flag, and compute those
            if self.p.perturb_flag
                % Compute exogenous torques and accelerations
                % still need to add attitude torques from drag effects
                tau_exo = grav_grad_tau(self, x_N, sigma_BN);
                a_exo = J2(self, x_N) ...
                    + atmo_drag_a(self, x_N, v_N, sigma_BN, m);
            else
                tau_exo = [0 0 0].';
                a_exo = [0 0 0].';
            end

            % Mass depletion dynamics (kg/s)
            m_dot = -u.mdot;
            % Inertial velocity (km/s)
            x_dot = v_N;
            % Inertial acceleration with input acceleration and perturbs
            v_dot = -(self.p.mu/(norm(x_N)^3)).*x_N ...
                + a_exo ...
                + u.a;
            
            % MRP integration (rad/s for angular rates)
            sigma_dot = (0.25.*((1 -(sigma_BN.'*sigma_BN))*eye(3) ...
                + 2*skew(sigma_BN) ...
                +(2*sigma_BN*(sigma_BN.'))))*omega_BN;
            
            % Angular Rate Integration with input torques
            omega_dot = self.Ic\((-skew(omega_BN)*(self.Ic*omega_BN)) ...
                + tau_exo ...
                + u.tau);
            
            % output the derivative
            f_x = [m_dot x_dot.' v_dot.' sigma_dot.' omega_dot.'].';
            
        end
        
        
        %% ODE45 function call
        function [time, state_out] = integrate_dynamics(self, time, X_0, odesettings)
            if ~exist('odesettings', 'var')
                odesettings = odeset('AbsTol', 1E-12, 'RelTol', 1E-12);
            end
            f_dot = @(t_in, state_in) self.dynamics(t_in, state_in);
            tic
            [time, state_out] = ode45(f_dot, time, X_0, odesettings);
            toc 
        end
        
        %% Rk4 integrator call to the 6DoF dynamics with control
        function [time, state_out] = integrate_dynamics_rk4(self, time, X_0)
            f_dot = @(t_in, state_in) self.dynamics(t_in, state_in);
            dt = time(2) - time(1);
            npoints = length(time);
            state_out = zeros(self.n, npoints);
            state_out(:,1) = X_0;
            tic
            for i = 1:npoints-1
                k_1 = f_dot(time(i), state_out(:,i));
                k_2 = f_dot(time(i)+0.5*dt, state_out(:,i)+0.5*dt*k_1);
                k_3 = f_dot((time(i)+0.5*dt), (state_out(:,i)+0.5*dt*k_2));
                k_4 = f_dot((time(i)+dt), (state_out(:,i)+k_3*dt));
                state_out(:,i+1) = state_out(:,i) + (1/6)*(k_1+(2*k_2)+(2*k_3)+k_4)*dt;
                s = norm(state_out(8:10, i+1));
                if s > 1
                    state_out(8:10,i+1) = -(state_out(8:10,i+1) ./(s^2));
                end
            end
            toc
        end
    end
end

