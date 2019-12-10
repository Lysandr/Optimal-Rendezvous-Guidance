function [time, state_out] = general_rk4(f_dot, time, X_0, p)
% General RK4 (DOES NOT HAVE THE MRP CONSTRAINTS !)
%     f_dot = @(t_in, state_in) self.nlrel(t_in, state_in);
    dt = time(2) - time(1);
    npoints = length(time);
    state_out = zeros(length(X_0), npoints);
    state_out(:,1) = X_0;
    tic
    for i = 1:npoints-1
        k_1 = f_dot(time(i), state_out(:,i),p);
        k_2 = f_dot(time(i)+0.5*dt, state_out(:,i) + 0.5*dt*k_1,p);
        k_3 = f_dot((time(i)+0.5*dt), (state_out(:,i)+0.5*dt*k_2),p);
        k_4 = f_dot((time(i)+dt), (state_out(:,i)+k_3*dt),p);
        state_out(:,i+1) = state_out(:,i) + (1/6)*(k_1+(2*k_2)+(2*k_3)+k_4)*dt;
        
        % keep MRPs within norm
        s = norm(state_out(8:10, i+1));
        if s > 1
            state_out(8:10,i+1) = -(state_out(8:10,i+1) ./(s^2));
        end
        s = norm(state_out(13+8:13+10, i+1));
        if s > 1
            state_out(13+8:13+10,i+1) = -(state_out(13+8:13+10,i+1) ./(s^2));
        end
    end
    toc
end