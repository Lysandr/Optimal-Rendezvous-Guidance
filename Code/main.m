%% Padraig Lysandrou ASEN6014 Project 1
clc; close all; clear all;

%% Begin the beguine
mu =  3.986004415e+05;

% orbital elements of chief
oe_c.a 		= 6371 + 400;
oe_c.e 		= 0;
oe_c.inc 	= deg2rad(51.6);
oe_c.Omega 	= deg2rad(0);
oe_c.omega 	= deg2rad(0);
oe_c.Mo 	= deg2rad(0);

% orbital elements of deputy
oe_d.a 		= oe_c.a + 0.001;
oe_d.e 		= oe_c.e;
oe_d.inc 	= oe_c.inc;
oe_d.Omega 	= oe_c.Omega;
oe_d.omega 	= oe_c.omega;
oe_d.Mo 	= oe_c.Mo;

% Plan out the timing
orbits = 1;
T = 2*pi*sqrt((oe_c.a^3)/mu);
T_f_zanetti = T;
p.T_f_zanetti = T_f_zanetti;
T_tot = orbits*T;
dt = 1;
t_0 = 0;
time = t_0:dt:T_tot;
npoints = length(time);

% Insantiate a spacecraft object, set some parameters
chief = spacecraft(oe_c, 'oe');
deputy= spacecraft(oe_d, 'oe');

% Turn on our pertubations
chief.p.perturb_flag    = 0;
deputy.p.perturb_flag   = 0;
chief.p.control_flag    = 0;
deputy.p.control_flag   = 1;

% determine other initial conditions
sigma_0_c = C2MRP(angle2dcm(0,-pi/1.9,0)*Cart2Hill(chief.rv'));
sigma_0_d = C2MRP(angle2dcm(0,-pi/1.9,0)*Cart2Hill(deputy.rv'));
omega_0 = [2*pi/T+0.0005 -0.0001 0.00].';
m_0 = 100;
r= 0.5; h= 200;
Ic = diag([(1/12)*(m_0*(3*r*r + h*h)) ...
    (1/12)*(m_0*(3*r*r + h*h)) (m_0*r*r)/2]);
chief.Ic = Ic;
deputy.Ic= Ic;

% initial condition vectors for both spacecraft
state_initial_c = [m_0 chief.rv.' sigma_0_c.' omega_0.'].';
chief.state_vector = state_initial_c;
state_initial_d = [m_0 deputy.rv.' sigma_0_d.' omega_0.'].';
deputy.state_vector = state_initial_d;

% store these mofos
p.formation(1) = chief;
p.formation(2) = deputy;
p.mu = mu;
p.control_mode = 1; % zanetti control mode
p.dt = dt;
p.glideslope = deg2rad(180);


theta = p.glideslope;
n = p.formation(1).oe.n;
% find the hill frame vectors
x_c = chief.rv(1:3); v_c = chief.rv(4:6);
x_d = deputy.rv(1:3); v_d = deputy.rv(4:6);
HN = Cart2Hill(chief.rv.');
rho_H = HN*(x_d-x_c);
f_c = (norm(p.formation(1).oe.h))/(norm(x_c)^2);
rhod_H= HN*(v_d-v_c) - skew([0 0 f_c].')*rho_H; % via transport
% find the DCM from hill to glideslope frame
RH = [sin(theta) -cos(theta) 0; cos(theta) sin(theta) 0; 0 0 1];
rtz = RH*rho_H; rtzd= RH*rhod_H;
% compute control and STM
r0 = rtz(1);
v0 = rtzd(1);

A = [0 1 0 0; ...
    3*n*n*(sin(theta)^2) 0 0 -1; ...
    -9*(n^4)*(cos(theta)^2)*(sin(theta)^2) 6*(n^3)*cos(theta)*sin(theta) 0 ...
    -3*(n^2)*(sin(theta)^2); ...
    6*(n^3)*sin(theta)*cos(theta) -4*(n^2) -1 0];
Phi_tf_0 = expm(A*T_f_zanetti);
Phi_rl = Phi_tf_0(1:2,3:4);
Phi_rr = Phi_tf_0(1:2,1:2);
costate_0 = Phi_rl\([0.001 0].' - Phi_rr*[r0 v0].');
p.z_0 = [r0 v0 costate_0.'].';


%% INTEGRATE THE 6DOF MODEL
%
f_x = @(t_in, state_in, p) formation_RHS(t_in, state_in, p);
X_0 = [state_initial_c.' state_initial_d.'].';
[time_rk, state_rk] = general_rk4(f_x, time, X_0, p);


% recall:
% state vector is [m, r_vec, v_vec, sigma_vec, omega_vec]
%                   1, 2:4,     5:7,    8:10,   11:14

figure;
subplot(3,1,1);
plot(time_rk, state_rk(2,:) - state_rk(13+2,:))
title('CARTESIAN POSITION ERROR')
subplot(3,1,2);
plot(time_rk, state_rk(3,:) - state_rk(13+3,:))
subplot(3,1,3);
plot(time_rk, state_rk(4,:) - state_rk(13+4,:))

figure;
subplot(3,1,1);
plot(time_rk, state_rk(5,:) - state_rk(13+5,:))
title('CARTESIAN VELOCITY ERROR')
subplot(3,1,2);
plot(time_rk, state_rk(6,:) - state_rk(13+6,:))
subplot(3,1,3);
plot(time_rk, state_rk(7,:) - state_rk(13+7,:))


figure;
subplot(3,1,1);
plot(time_rk, state_rk(13+2,:))
title('DEPUTY CARTESIAN position')
subplot(3,1,2);
plot(time_rk, state_rk(13+3,:))
subplot(3,1,3);
plot(time_rk, state_rk(13+4,:))


figure;
subplot(3,1,1);
plot(time_rk, state_rk(2,:))
title('CHIEF CARTESIAN position')
subplot(3,1,2);
plot(time_rk, state_rk(3,:))
subplot(3,1,3);
plot(time_rk, state_rk(4,:))


figure;
subplot(3,1,1);
plot(time_rk, state_rk(13+5,:))
title('DEPUTY CARTESIAN VELOCITY')
subplot(3,1,2);
plot(time_rk, state_rk(13+6,:))
subplot(3,1,3);
plot(time_rk, state_rk(13+7,:))


figure;
subplot(3,1,1);
plot(time_rk, state_rk(5,:))
title('CHIEF CARTESIAN VELOCITY')
subplot(3,1,2);
plot(time_rk, state_rk(6,:))
subplot(3,1,3);
plot(time_rk, state_rk(7,:))




% figure;
% subplot(3,1,1);
% plot(time_rk, state_rk(13 + 11,:))
% title('DEPUTY ANGULAR RATE')
% subplot(3,1,2);
% plot(time_rk, state_rk(13+12,:))
% subplot(3,1,3);
% plot(time_rk, state_rk(13 + 13,:))
























