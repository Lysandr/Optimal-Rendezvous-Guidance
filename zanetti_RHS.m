function f_x = zanetti_RHS(t_in, state_in, p)
    rtz = state_in(1:3);
    rtzd= state_in(4:6);
    costate = state_in(7:8);
    
    ur = p.u(1);
    ut = p.u(2);
    uz = p.u(3);
    
    rdd = 2*p.n*rtzd(2) + 3*(p.n^2)*rtz(1)*(sin(p.theta)^2) ...
        + 3*(p.n^2)*rtz(2)*sin(p.theta)*cos(p.theta) + ur;
    tdd = -2*p.n*rtzd(1)+ 3*(p.n^2)*rtz(2)*(cos(p.theta)^2) ...
        + 3*(p.n^2)*rtz(1)*sin(p.theta)*cos(p.theta) + ut;
    zdd = -(p.n^2)*rtz(3) + uz;
    
    xd = p.A*([rtz(1) rtzd(1) costate.'].');
    lambda_d = xd(3:4);
    
    f_x = [rtzd.' rdd tdd zdd lambda_d.'].';
end























