function dxdt = cw_hill_ode(~, X, p)
%cw_hill_ode Linear Clohessy Wiltshire equations in the hill frame
% CW equations assume circular orbit for the chief
    % pull out constants
    oe = p.oe_c;
    n = sqrt(p.mu/(oe.a^3));

    % pull out states  
    x 	= X(1);
    y 	= X(2);
    z 	= X(3);
    xd	= X(4);
    yd	= X(5);
    zd	= X(6);
    

    % ODEs to solve (all linear...)
    xdd = (2*n*yd) + (3*n*n*x);
    ydd = (-2*n*xd);
    zdd = (-n*n*z);

    dxdt =[xd; yd; zd; xdd; ydd; zdd];
end
