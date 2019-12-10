function f_x = formation_RHS(t_in, state_in, p)
    num_sats = length(p.formation);
    state_n = length(state_in)/num_sats;
    f_x = [];
    % state vector is [m, r_vec, v_vec, sigma_vec, omega_vec]
    %                   1, 2:4,     5:7,    8:10,   11:14

    for k = 1:num_sats
        % distribute the states for use! these are row vectors.
        ind_val = (k-1)*state_n;
        m =         state_in(ind_val + 1);
        x_N =       state_in(ind_val+2 : ind_val+4);
        v_N =       state_in(ind_val + 5: ind_val + 7);
        sigma_BN =  state_in(ind_val + 8: ind_val + 10);
        omega_BN =  state_in(ind_val + 11:ind_val + 13);

        if p.formation(k).p.control_flag
            u = compute_control(t_in, state_in, p, k);
        else
            u.tau = [0 0 0].';
            u.a = [0 0 0].';
            u.mdot = 0;
        end

        % Check for pertubation flag, and compute those
        if p.formation(k).p.perturb_flag
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
        v_dot = -(p.mu/(norm(x_N)^3)).*x_N ...
            + a_exo ...
            + u.a;

        % MRP integration (rad/s for angular rates)
        sigma_dot = (0.25.*((1 -(sigma_BN.'*sigma_BN))*eye(3) ...
            + 2*skew(sigma_BN) ...
            +(2*sigma_BN*(sigma_BN.'))))*omega_BN;

        % Angular Rate Integration with input torques
        omega_dot = p.formation(k).Ic\((-skew(omega_BN)*(p.formation(k).Ic*omega_BN)) ...
            + tau_exo ...
            + u.tau);

        % output the derivative
        f_x = [f_x; [m_dot x_dot.' v_dot.' sigma_dot.' omega_dot.'].'];
    end
end























