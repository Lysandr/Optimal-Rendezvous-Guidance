function [r,v] = oe2rv(mu,oe)
%orbital elements to rv cartesian vectors, needs h,e,omega,Omega,i,f0
% Vallado page 118
    h = norm(oe.h);
    e = norm(oe.e);
    omega = oe.omega;
    Omega = oe.Omega;
    inc = oe.inc;
    nu = oe.nu;

    rp = (h^2/mu)*(1/(1 + e*cos(nu)))*(cos(nu)*[1;0;0] + sin(nu)*[0;1;0]);
    vp = (mu/h)*(-sin(nu)*[1;0;0]+(e + cos(nu))*[0;1;0]);

    % R3_omega, DCM rotation about z, omega amount
    % R3_Omega, DCM rotation about the z, Omega amount
    % R1_inc,  DCM rotation about the x, inc amount
    % i wonder if anyone has written a quaternion form for this
    R3_omega = [ cos(omega)  sin(omega)  0 
                -sin(omega)  cos(omega)  0
                0       0     1];
    R3_Omega = [ cos(Omega)  sin(Omega)  0
                -sin(Omega)  cos(Omega)  0
                 0        0     1];
    R1_inc = [1       0          0
                0   cos(inc)  sin(inc)
                0  -sin(inc)  cos(inc)];
    % total rotations SO3 group, can be multiplied any way
    perif_to_ECI = (R3_omega*R1_inc*R3_Omega).';
    r = perif_to_ECI*rp;
    v = perif_to_ECI*vp;
end

