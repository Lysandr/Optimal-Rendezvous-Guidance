%% lysandrou prj2 symsolve
clc; clear all; close all


syms x y z xd yd zd ux uy uz n theta real
syms r t rd td ur ut real

xdd = 2*n*yd + 3*n*n*x + ux;
ydd = -2*n*xd + uy;
zdd = -n*n*z+ uz;

R = [sin(theta) -cos(theta); cos(theta) sin(theta)];
Rinv = inv(R).';
Rinv = subs(Rinv, cos(theta)^2 + sin(theta)^2, 1)

% substitution shit
states = [x;y];
states_d = [xd;yd];
statesg = [r;t];
statesg_d = [rd;td];

rad_transv = Rinv*statesg;
rad_transv_d = Rinv*statesg_d;
rad_transv_dd = R*[xdd;ydd];

rdd = rad_transv_dd(1);
tdd = rad_transv_dd(2);

rdd = subs(rdd, states, rad_transv);
rdd = subs(rdd, states_d, rad_transv_d);
rdd = simplify(rdd);
latex(rdd)

tdd = subs(tdd, states, rad_transv);
tdd = subs(tdd, states_d, rad_transv_d);
tdd = simplify(tdd);
latex(tdd)

% do the inputs now
inputs = [ux; uy];
inputsg = [ur; ut];
rt_inputs = Rinv*inputsg;

rdd = subs(rdd, inputs, rt_inputs);
rdd = simplify(rdd);
latex(rdd)

tdd = subs(tdd, inputs, rt_inputs);
tdd = simplify(tdd);
latex(tdd)



thingforhamiltonion = (2*n*rd - 3*n*n*r*sin(theta)*cos(theta))^2
latex(diff(thingforhamiltonion,r)/2)



























