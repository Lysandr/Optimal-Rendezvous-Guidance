
function p_out = schaub_elements(state, p,delta_t,inputFlag,outputFlag)
%schaub_elements performs operations: rv2oe and oe2rv for keplerian
%elements and other trashy routines
% state must be a [r v] for [1, 2] flags, and oe struct for [2,1] flags
% todo: [1,1] flags being integate rv from t1 to t2, same for [2,2]

if inputFlag == 1 && outputFlag == 2
    % compute the RV2OE transform
    % first extract the important parameters
    mu = p.mu;
    r = state(1:3);
    v = state(4:6);
    
    rhat = r./norm(r);
    h = cross(r,v);
    hhat = h./norm(h);
    e = (1/mu).*((v.'*v - (mu./norm(r))).*r - (r.'*v).*v);
    ehat = e./norm(e);
    ehat_p = cross([0 0 1],ehat).';
    p = norm(h)^2/mu;
    energy = (0.5*norm(v)^2) - (mu/norm(r));
    inc = acos(h(3)/norm(h));
    zhat = [0 0 1].';
    nhat_O = (cross(zhat,hhat))/(norm(cross(zhat,hhat)));
    Omega = atan2(dot([0 1 0].',nhat_O),dot([1 0 0].',nhat_O));
    if inc == 0
        Omega = 0;
    end
    nhat_Op = cross(hhat,nhat_O);
    omega = atan2(dot(e,nhat_Op),dot(e,nhat_O));
    a = p/(1-norm(e)^2);
    % calculate the initial true anomaly
    f0 = atan2(dot(r,ehat_p),dot(r,ehat));
    T = 2*pi*sqrt(abs(a^3)/mu);
    n = (2*pi)/T;
    
    % package up interesting parameters
    oe.h = h;
    oe.a = a;
    oe.e = norm(e);
    oe.inc = inc;
    oe.omega = omega;
    oe.Omega = Omega;
    oe.p = p;
    oe.f0 = f0;
    oe.T = T;
    oe.n = n;
    oe.energy = energy;
    % calculate the actual true anomaly    
%     f0 = atan((cross(ehat, rhat).'*hhat)/(ehat.'*rhat));
    if norm(e) <= 1
       E = 2*atan(sqrt((1-norm(e))/(1+norm(e)))*tan(f0/2));
       M = E - (norm(e)*sin(E));
    elseif norm(e) > 1
       H = 2*atanh(sqrt((norm(e)-1)/(1+norm(e)))*tan(f0/2));
       M = (norm(e)*sinh(H) - H);
    end
    if M < 0
        M = M + 2*pi;
    elseif M > 2*pi
        M = mod(M, 2*pi);
    end
    oe.Mo = M - (sqrt(mu/abs(a^3))*delta_t);
    [oe.nu, Ecc_anom] = findE(delta_t,norm(e), sqrt(mu/abs(a^3)), oe.Mo);
    oe.E = Ecc_anom;
    p_out = oe;
    
    
elseif inputFlag == 2 && outputFlag == 1
    % compute the OE2RV transform
    mu = p.mu;
    oe = state;
    a = oe.a;
    e = norm(oe.e);
    h = sqrt(mu*a*(1-(e^2)));
    omega = oe.omega;
    Omega = oe.Omega;
    inc = oe.inc;
    n = sqrt(mu/abs(a^3));
    [nu, Ecc_anom] = findE(delta_t, e, n, oe.Mo);
    
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

    % total rotations SO3 group
    perif_to_ECI = (R3_omega*R1_inc*R3_Omega).';
    r = perif_to_ECI*rp;
    v = perif_to_ECI*vp;
    
    % package up and send out
    p_out.E = Ecc_anom;
    p_out.r = r;
    p_out.v = v;
    p_out.nu = nu;
    p_out.h = h;
    p_out.n = n;
    p_out.a = a;
    p_out.e = e;
    p_out.inc =  inc;
    p_out.omega = omega;
    p_out.Omega = Omega;
    p_out.Mo = oe.Mo;
else
    
end


end

