function plot_all_the_things(time,state)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

    % Pull out the states
    m =         state(:, 1);
    x_N =       state(:, 2:4);
    v_N =       state(:, 5:7);
    sigma_BN =  state(:, 8:10);
    omega_BN =  state(:, 11:13);
    
    %% plot the angular momentum along the entire orbit (debug)
%     h_norm = zeros(1, length(time));
%     h1 = norm(skew(x_N(1,:).')*v_N(1,:).');
%     for i = 1:length(time)
%         h_norm(i) = norm(skew(x_N(i,:).')*v_N(i,:).') - h1;
%     end
%     figure;
%     plot(time, h_norm);
%     title('Normed ang mom difference over time')
    
    
    %% plot the orbital energy over time (debug)
%     energy_orb = zeros(1, length(time));
%     for i = 1:length(time)
%         energy_orb(i) = (norm(v_N(i,:))^2)/2 - (398600/norm(x_N(i,:)));
%     end
%     figure;
%     plot(time, energy_orb);
%     title('Orbital energy over time')
%    

    %% Plot the orbit history onto the earth sphere
    figure;
    plot3(x_N(:,1), x_N(:,2), x_N(:,3), 'r'); hold on;
    plot3(x_N(1,1), x_N(1,2), x_N(1,3), 'go');
    earth_sphere(); hold off;
    title('Orbit over time')
    set(gca,'Color','k')
    grid on;
    
    %% Plot the norm of the radial vector over time
%     figure;
%     plot(time, vecnorm(x_N.'));
%     title('Orbit radius over time')
%     set(gca,'Color','k')
%     grid on;
    
    %% Plot velocity norm history over time
%     figure;
%     plot(time, vecnorm(v_N.'));
%     title('Velocity Magnitude over time')
%     set(gca,'Color','k')
%     grid on;

    %% Plot the body frame angular rate over time
    figure;
    plot(time, omega_BN(:,1)); hold on;
    plot(time, omega_BN(:,2)); hold on;
    plot(time, omega_BN(:,3)); hold on;
    plot(time, vecnorm(omega_BN.'));
    legend('omega1', 'omega2', 'omega3', 'omeganorm')
    title('Orbit \omega_{BN} over time')
    grid on;

    %% Plot the attitude history
    figure;
    plot(time, sigma_BN(:,1)); hold on;
    plot(time, sigma_BN(:,2)); hold on;
    plot(time, sigma_BN(:,3)); hold on;
    title('Orbit \sigma_{BN} over time')
    grid on;

    %% Plot inertial frame omega over time
    % [BN]T * B_omega_BN = N_omega_BN
    N_omega_BN = zeros(3,length(time));
    for i = 1:length(time)
        N_omega_BN(:,i) = MRP2C(sigma_BN(i,:).').'*(omega_BN(i,:).');
    end
    figure;
    plot(time, (N_omega_BN)); grid on; hold on;
    plot(time, vecnorm(N_omega_BN))
    legend('omega1', 'omega2', 'omega3', 'omeganorm')
    title('Inertial frame \omega_{BN} over time')
    
    %% plot the boresight (+x) vector over time
    b = zeros(3,length(time));
    for i = 1:length(time)
        NB = MRP2C(sigma_BN(i,:).').';
        HN = Cart2Hill([x_N(i,:) v_N(i,:)]);
        b(:,i) = HN*NB*([0 0 1].');
    end
    figure;
    plot3(b(1,:),b(2,:),b(3,:)); hold on;
    line_to_ic = b(:,1)*(0:0.01:1);
    plot3(line_to_ic(1,:),line_to_ic(2,:),line_to_ic(3,:),'r');
    plot3(b(1,1),b(2,1),b(3,1),'ro');
    plot3(0:.01:1,zeros(1,101),zeros(1,101),'Color',[0.8500 0.3250 0.0980]);
    plot3(zeros(1,101),0:.01:1,zeros(1,101),'Color',[0.8500 0.3250 0.0980]);
    plot3(zeros(1,101),zeros(1,101),0:.01:1,'Color',[0.8500 0.3250 0.0980]);
    grid on; axis equal;
%     fv = stlread('femur.stl');
%     patch(fv,'FaceColor',       [0.8 0.8 1.0], ...
%          'EdgeColor',       'none',        ...
%          'FaceLighting',    'gouraud',     ...
%          'AmbientStrength', 0.15);
% 
%     % Add a camera light, and tone down the specular highlighting
%     camlight('headlight');
%     material('dull');
    % Fix the axes scaling, and set a nice view angle
    axis('image');
    view([-71 35]);
    xlabel('Or (radial)');
    ylabel('Otheta (local tangent)');
    zlabel('Oh (radial normal)');
    title('Hill frame boresight vector over time')
    
    %% plot THE NOODLE
%     figure;
%     for i=1:length(time)
%         att_vect = MRP2C(sigma_BN(i,:).')*[0;0;1];
%         plot3(att_vect(1), att_vect(2), att_vect(3), ...
%             'or','MarkerSize',5,'MarkerFaceColor','r');
%         axis([-1.1 1.1 -1.1 1.1 -1.1 1.1]);
%         hold on; grid on;
%         pause(0.0000000001)
%     end
end

