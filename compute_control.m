function u = compute_control(t, state_in, p, k)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
    
    % pull out state for the satellite we are controlling
    num_sats = length(p.formation);
    state_n = length(state_in)/num_sats;
    ind_val = (k-1)*state_n;
    m =         state_in(ind_val + 1);
    x_N =       state_in(ind_val + 2 : ind_val+4);
    v_N =       state_in(ind_val + 5: ind_val + 7);
    sigma_BN =  state_in(ind_val + 8: ind_val + 10);
    omega_BN =  state_in(ind_val + 11:ind_val + 13);
    
    Kp = 0.001;
    Kd = 0.003;
    Kz = 0.003;
    
    % zanetti guidance model
    if p.control_mode == 1
        theta = p.glideslope;
        n = p.formation(1).oe.n;
        % find the hill frame vectors
        x_c = state_in(2:4); v_c = state_in(5:7);
        HN = Cart2Hill([x_c.' v_c.']);
        rho_H = HN*(x_N-x_c);
        f_c = (norm(p.formation(1).oe.h))/(norm(x_c)^2);
        rhod_H= HN*(v_N-v_c) - skew([0 0 f_c].')*rho_H; % via transport
        % find the DCM from hill to glideslope frame
        RH = [sin(theta) -cos(theta) 0; cos(theta) sin(theta) 0; 0 0 1];
        rtz = RH*rho_H;
        rtzd= RH*rhod_H;
        % compute control and STM
        r = rtz(1);
        tee = rtz(2);
        r_dot = rtzd(1);
        t_dot = rtzd(2);
        z_dot = rtzd(3);
        
        if t < p.T_f_zanetti
            A = [0 1 0 0; ...
                3*n*n*(sin(theta)^2) 0 0 -1; ...
                -9*(n^4)*(cos(theta)^2)*(sin(theta)^2) 6*(n^3)*cos(theta)*sin(theta) 0 ...
                -3*(n^2)*(sin(theta)^2); ...
                6*(n^3)*sin(theta)*cos(theta) -4*(n^2) -1 0];
            u_r_star = [0 0 0 -1]*expm(A*t)*p.z_0;
            u_t_star = 2*n*r_dot - 3*(n^2)*r*sin(theta)*cos(theta);
        else
            u_r_star = 0; u_t_star = 0;
        end
        u_r = u_r_star - 2*n*t_dot - 3*(n^2)*sin(theta)*cos(theta)*tee;
        u_t = u_t_star - 3*(n^2)*(cos(theta)^2)*tee - Kp*tee - Kd*t_dot;
        u_z = -Kz*z_dot;
        u_H = RH.'*[u_r u_t u_z].';
        u_N = HN.'*u_H;
    end

    
    % gains for detumble law
%     P = max(diag(p.formation(2).Ic).*(1/100000));
%     K = ((P^2)./p.formation(2).Ic(2,2))/100000;

    % body frame torques and inertial frame accelerations
    % ACS solution is naive, just quick and dirty
%     u.tau = -K*sigma_BN - P*omega_BN;
    u.tau = [0 0 0].';
    u.a = u_N;
    u.mdot = 0; 
end

