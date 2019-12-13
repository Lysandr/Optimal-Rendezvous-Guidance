%% padraig lysandrou
clc; close all; clear all;

mu =  3.986004415e+05;      p.mu = mu;

% orbital elements of chief		orbital elements of deputy
oe_c.a 		= 6371 + 400;		oe_d.a 		= oe_c.a - 0.2;
oe_c.e 		= 0;				oe_d.e 		= oe_c.e;
oe_c.inc 	= deg2rad(51.6);	oe_d.inc 	= oe_c.inc;
oe_c.Omega 	= deg2rad(0);		oe_d.Omega 	= oe_c.Omega;
oe_c.omega 	= deg2rad(0);		oe_d.omega 	= oe_c.omega;
oe_c.Mo 	= deg2rad(0);		oe_d.Mo 	= oe_c.Mo + deg2rad(.0001);

% Plan out the timing
orbits = 1/5;         T = 2*pi*sqrt((oe_c.a^3)/mu);
T_tot = 600;

% Guidance law timing
T_f_zanetti = T_tot;  p.T_f_zanetti = T_f_zanetti;

% Sim timing
dt = 0.01;
time = linspace(0,T_tot,T_tot/dt);
control_time = time(1:end-1);
npoints = length(time);

% get other parameters
p_out_c = schaub_elements(oe_c, p,0,2,1);
p_out_d = schaub_elements(oe_d, p,0,2,1);
rc = p_out_c.r;         vc = p_out_c.v;
rd = p_out_d.r;         vd = p_out_d.v;

% get the initial hill frame coordinates
HN = Cart2Hill([rc.' vc.']);
H_omega_ON = [0; 0; p_out_c.n];
rho_N = (rd-rc);
rho_Hill = HN*rho_N;
rhod_Hill = HN*(vd-vc) - skew(H_omega_ON)*rho_Hill;

% the things that matter [or otheta oz]
% rho_Hill = [-.2 0.01 0].';
% rhod_Hill= [0 0 0].';
theta = deg2rad(-90);

% from this, get the transverse and
p.theta = theta;
RH = [sin(theta) -cos(theta) 0; cos(theta) sin(theta) 0; 0 0 1];
rtz = RH*rho_Hill;
rtzd= RH*rhod_Hill;

% find the costate initial condition
n = p_out_c.n;          p.n = n;
A = [0 1 0 0; ...
    3*(n^2)*(sin(theta)^2) 0 0 -1; ...
    -9*(n^4)*(cos(theta)^2)*(sin(theta)^2) 6*(n^3)*cos(theta)*sin(theta) 0 ...
    -3*(n^2)*(sin(theta)^2); ...
    6*(n^3)*sin(theta)*cos(theta) -4*(n^2) -1 0];
p.A = A;
Phi_tf_0 = expm(A*T_f_zanetti);
Phi_rl = Phi_tf_0(1:2,3:4);
Phi_rr = Phi_tf_0(1:2,1:2);
costate_0 = Phi_rl\([0 0].' - Phi_rr*([rtz(1) rtzd(1)].'));

% setup the initial conditions and other bluushuit
X_0 = [rtz; rtzd; costate_0];
state_out = zeros(length(X_0), npoints);
control_hist = zeros(3, npoints-1);
control_hist_H = zeros(3, npoints-1);
state_out(:,1) = X_0;
f_dot = @(t_in, state_in, p) zanetti_RHS(t_in, state_in, p);

% gains for things
Kp = 0.0005;
Kd = 0.05;
Kz = 0.003;

% calculate the dt STM
Phi_dt = expm(A*dt);
eig(Phi_dt)

tic
% integration loop
for i = 1:npoints-1
    rtz = state_out(1:3,i);
    rtzd= state_out(4:6,i);
    costate = state_out(7:8,i);
    
    % compute the control
    x_t = [rtz(1) rtzd(1) costate.'].';
    ur_star = [0 0 0 -1]*Phi_dt*x_t ...
        - 2*n*rtzd(2) - 3*(n^2)*sin(theta)*cos(theta)*rtz(2);
    ut_star = 2*n*rtzd(1) - 3*(n^2)*rtz(1)*sin(theta)*cos(theta) ...
        - 3*(n^2)*(cos(theta)^2)*rtz(2) - Kp*rtz(2) - Kd*rtzd(2);
    uz = -Kz*rtzd(3);
    control_hist(:,i) = [ur_star ut_star uz].';
    p.u = control_hist(:,i);
    control_hist_H(:,i) = RH.'*p.u;
%     p.u = [0 0 0].';
    
    % integrate the dynamics
    k_1 = f_dot(time(i), state_out(:,i),p);
    k_2 = f_dot(time(i)+0.5*dt, state_out(:,i) + 0.5*dt*k_1,p);
    k_3 = f_dot((time(i)+0.5*dt), (state_out(:,i)+0.5*dt*k_2),p);
    k_4 = f_dot((time(i)+dt), (state_out(:,i)+k_3*dt),p);
    state_out(:,i+1) = state_out(:,i) + (1/6)*(k_1+(2*k_2)+(2*k_3)+k_4)*dt;
end
toc


state_plot = zeros(length(X_0), npoints);
input_hill = zeros(3, npoints -1);
for i = 1:npoints
    state_plot(1:3,i) = RH.'*state_out(1:3,i);
    state_plot(4:6,i) = RH.'*state_out(4:6,i);
end

for i = 1:npoints-1
   input_hill(:,i) =   RH.'*control_hist(:,i); 
end


%% plotting
close all
width=4; height=4;

figure('Units','inches','Position',[0 0 width height],'PaperPositionMode','auto');
plot(-state_plot(2,:), state_plot(1,:));  hold on;
plot(-state_plot(2,1), state_plot(1,1),'go');
plot(-state_plot(2,end), state_plot(1,end),'ro');
axis equal; grid on;
xlabel('o_\theta alongtrack km');
ylabel('o_r radial km')
legend('traj','IC','FC','location','best')
title('Hill Frame Trajectory of Deputy over time')
set(gca,...
'Units','normalized',...
'FontUnits','points',...    
'FontWeight','normal',...
'FontSize',9,...
'FontName','Times')
print -depsc2 traj.eps

% plot the same thing with quivers
figure('Units','inches','Position',[0 0 4 6],'PaperPositionMode','auto');
plot(-state_plot(2,:), state_plot(1,:));  hold on;
plot(-state_plot(2,1), state_plot(1,1),'go');
plot(-state_plot(2,end), state_plot(1,end),'ro');
axis equal; grid on;
% sfact = 300;
% pfact = 1/40;
q = quiver(-state_plot(2,1:1000:end-1), state_plot(1,1:1000:end-1), ...
    -input_hill(2,1:1000:end), input_hill(1,1:1000:end));
q.Color = 'red';
q.LineWidth = 0.5;
q.MaxHeadSize = 0.01;
% q.AutoScale = 'off';
% hold on
xlabel('o_\theta alongtrack km');
ylabel('o_r radial km')
legend('traj','IC','FC','Thrust Vector','location','best')
title('Hill Frame Trajectory of Deputy over time')
set(gca,...
'Units','normalized',...
'FontUnits','points',...    
'FontWeight','normal',...
'FontSize',9,...
'FontName','Times')
print -depsc2 traj_quiv.eps

figure('Units','inches','Position',[0 0 width height],'PaperPositionMode','auto');
plot(time, state_plot(4:6,:)); grid on;
title('Hill frame velocities over time')
xlabel('Time, s')
ylabel('Velocity km/s')
legend('xdot','ydot','zdot','location','best')
set(gca,...
'Units','normalized',...
'FontUnits','points',...    
'FontWeight','normal',...
'FontSize',9,...
'FontName','Times')
print -depsc2 hillvels.eps

figure('Units','inches','Position',[0 0 width height],'PaperPositionMode','auto');
plot(control_time, control_hist); grid on;
title('Control History vs Time')
xlabel('Time, s')
ylabel('Acceleration km/s2')
legend('ur','ut','uz')
set(gca,...
'Units','normalized',...
'FontUnits','points',...    
'FontWeight','normal',...
'FontSize',9,...
'FontName','Times')
print -depsc2 controls.eps

figure('Units','inches','Position',[0 0 4 2],'PaperPositionMode','auto');
plot(control_time, control_hist_H); grid on;
title('Optimal Control vs Time')
xlabel('Time, s')
ylabel('Acceleration km/s2')
legend('ux (Hill)','uy','uz','location','best')
set(gca,...
'Units','normalized',...
'FontUnits','points',...    
'FontWeight','normal',...
'FontSize',9,...
'FontName','Times')
print -depsc2 controlsH.eps

figure('Units','inches','Position',[0 0 width height],'PaperPositionMode','auto');
plot(time, state_out(4,:)); hold on
plot(time, state_out(5,:));
title('RADIAL TRANSVERSAL frame velocity ')
legend('r','t'); hold off; grid on;
set(gca,...
'Units','normalized',...
'FontUnits','points',...    
'FontWeight','normal',...
'FontSize',9,...
'FontName','Times')
print -depsc2 transvels.eps





dv_over_time = zeros(1, npoints-1);
dv_over_time(1) = dt*norm(control_hist_H(:,1));
for i = 1:npoints-2
   dv_over_time(i+1) = (dt*norm(control_hist_H(:,i+1))) + dv_over_time(i); 
end

figure('Units','inches','Position',[0 0 4 2],'PaperPositionMode','auto');
plot(control_time, dv_over_time); grid on;
title('Optimal Control dV over Time')
xlabel('Time, s')
ylabel('dV km/s')
set(gca,...
'Units','normalized',...
'FontUnits','points',...    
'FontWeight','normal',...
'FontSize',9,...
'FontName','Times')
print -depsc2 dVOptimal.eps








%%
f_dot = @(t_in, state_in, p) two_body(t_in, state_in, p);

init_conds = [rc; rd; vc; vd];
state_cart = zeros(12, npoints);
state_cart(:,1) = init_conds;
rdd_hist = zeros(3, npoints-1);
rho_H_hist= zeros(3, npoints-1);
u_d_hist = zeros(3,npoints-1);


% Kp = 5e-4;
% Kd = 1e-2;
K1 = .0005*eye(3);
K2 = .05*eye(3);
tic
% integration loop
for i = 1:npoints-1
    rc = state_cart(1:3,i);
    rd = state_cart(4:6,i);
    vc = state_cart(7:9,i);
    vd = state_cart(10:12,i);
    
    % compute the control
    rdd = rc;
    vdd = vc;
    a_dd= -(p.mu/(norm(rdd)^3))*rdd;
    delta_r = rd - rdd;
    delta_r_dot = vd - vdd;
    p.u_d = (p.mu/norm(rd)^3)*rd + a_dd - K1*delta_r - K2*delta_r_dot;
    
    u_d_hist(:,i) = p.u_d;
    
    % integrate the dynamics
    k_1 = f_dot(time(i), state_cart(:,i),p);
    k_2 = f_dot(time(i)+0.5*dt, state_cart(:,i) + 0.5*dt*k_1,p);
    k_3 = f_dot((time(i)+0.5*dt), (state_cart(:,i)+0.5*dt*k_2),p);
    k_4 = f_dot((time(i)+dt), (state_cart(:,i)+k_3*dt),p);
    state_cart(:,i+1) = state_cart(:,i) + (1/6)*(k_1+(2*k_2)+(2*k_3)+k_4)*dt;
end
toc

rho_cart = zeros(3, npoints);
hframvels = rho_cart;
input_cart = zeros(3, npoints -1);
dv_over_time = zeros(1,npoints-1);
dv_over_time(1) = dt*norm(u_d_hist(:,1));
for i = 1:npoints
    rc_loop = state_cart(1:3,i);
    vc_loop = state_cart(7:9,i);
    HN = Cart2Hill([rc_loop.' vc_loop.']);
    rho_cart(:,i) = HN*(state_cart(4:6,i) - rc_loop);
    hframvels(:,i) = HN*(state_cart(10:12,i)-vc_loop) ...
        - skew(H_omega_ON)*rho_cart(:,i);
end

for i = 1:npoints-2
   dv_over_time(i+1) = (dt*norm(u_d_hist(:,i+1))) + dv_over_time(i); 
end

for i = 1:npoints-1
    rc_loop = state_cart(1:3,i);
    vc_loop = state_cart(7:9,i);
    HN = Cart2Hill([rc_loop.' vc_loop.']);
    input_cart(:,i) = HN*u_d_hist(:,i); 
end

figure('Units','inches','Position',[0 0 4 2],'PaperPositionMode','auto');
plot(control_time, input_cart); grid on;
title('Lyapunov Control vs Time')
xlabel('Time, s')
ylabel('Acceleration km/s2')
legend('ux (Hill)','uy','uz','location','best')
set(gca,...
'Units','normalized',...
'FontUnits','points',...    
'FontWeight','normal',...
'FontSize',9,...
'FontName','Times')
print -depsc2 controlsCart.eps


figure('Units','inches','Position',[0 0 width height],'PaperPositionMode','auto');
plot(-rho_cart(2,:), rho_cart(1,:));  hold on;
plot(-rho_cart(2,1), rho_cart(1,1),'go');
plot(-rho_cart(2,end), rho_cart(1,end),'ro');
axis equal; grid on;
xlabel('o_\theta alongtrack km');
ylabel('o_r radial km')
legend('traj','IC','FC','location','best')
title('Hill Frame Trajectory of Deputy over time')
set(gca,...
'Units','normalized',...
'FontUnits','points',...    
'FontWeight','normal',...
'FontSize',9,...
'FontName','Times')
print -depsc2 trajLyap.eps


figure('Units','inches','Position',[0 0 width height],'PaperPositionMode','auto');
plot(time, hframvels); grid on;
title('Hill frame velocities over time')
xlabel('Time, s')
ylabel('Velocity km/s')
legend('xdot','ydot','zdot','location','best')
set(gca,...
'Units','normalized',...
'FontUnits','points',...    
'FontWeight','normal',...
'FontSize',9,...
'FontName','Times')
print -depsc2 velsLyap.eps


figure('Units','inches','Position',[0 0 4 2],'PaperPositionMode','auto');
plot(control_time, dv_over_time); grid on;
title('Lyapunov dV over Time')
xlabel('Time, s')
ylabel('dV km/s')
set(gca,...
'Units','normalized',...
'FontUnits','points',...    
'FontWeight','normal',...
'FontSize',9,...
'FontName','Times')
print -depsc2 dVCart.eps
















