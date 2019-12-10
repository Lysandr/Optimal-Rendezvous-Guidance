function plot_two_sats(time,chief, deputy)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

    % Pull out the states
    m_c =         chief(:, 1);
    x_N_c =       chief(:, 2:4);
    v_N_c =       chief(:, 5:7);
    sigma_BN_c =  chief(:, 8:10);
    omega_BN_c =  chief(:, 11:13);
    m_c =         deputy(:, 1);
    x_N_d =       deputy(:, 2:4);
    v_N_d =       deputy(:, 5:7);
    sigma_BN_d =  deputy(:, 8:10);
    omega_BN_d =  deputy(:, 11:13);
    
    %% plot the angular momentum along the entire orbit (debug)
    h_norm_d = zeros(1, length(time));
    h_norm_c = zeros(1, length(time));
    h1 = norm(skew(x_N_c(1,:).')*v_N_c(1,:).');
    h2 = norm(skew(x_N_d(1,:).')*v_N_d(1,:).');
    for i = 1:length(time)
        h_norm_c(i) = norm(skew(x_N_c(i,:).')*v_N_c(i,:).') - h1;
        h_norm_d(i) = norm(skew(x_N_d(i,:).')*v_N_d(i,:).') - h2;
    end
    figure;
    plot(time, h_norm_c); hold on;
    plot(time, h_norm_d);
    title('Normed ang mom difference over time')
    
    
    %% plot the orbital energy over time (debug)
    energy_orb_c = zeros(1, length(time));
    energy_orb_d = zeros(1, length(time));
    for i = 1:length(time)
        energy_orb_c(i) = (norm(v_N_c(i,:))^2)/2 - (398600/norm(x_N_c(i,:)));
        energy_orb_d(i) = (norm(v_N_d(i,:))^2)/2 - (398600/norm(x_N_d(i,:)));
    end
    figure;
    plot(time, energy_orb_c); hold on;
    plot(time, energy_orb_d);
    title('Orbital energy over time')
   

    %% Plot the orbit history onto the earth sphere
    figure;
    plot3(x_N_c(:,1), x_N_c(:,2), x_N_c(:,3)); hold on;
    plot3(x_N_c(1,1), x_N_c(1,2), x_N_c(1,3), 'go');
    plot3(x_N_d(:,1), x_N_d(:,2), x_N_d(:,3)); hold on;
    plot3(x_N_d(1,1), x_N_d(1,2), x_N_d(1,3), 'bo');
    earth_sphere();
    title('Orbit over time')
    set(gca,'Color','k')
    grid on;
    
    %% Plot the norm of the radial vector over time
    figure;
    plot(time, vecnorm(x_N_c.')); hold on;
    plot(time, vecnorm(x_N_d.'));
    title('Orbit radius over time')
    set(gca,'Color','k')
    grid on;
    
    %% Plot velocity norm history over time
    figure;
    plot(time, vecnorm(v_N_c.')); hold on;
    plot(time, vecnorm(v_N_d.'));
    title('Velocity Magnitude over time')
    set(gca,'Color','k')
    grid on;


end

