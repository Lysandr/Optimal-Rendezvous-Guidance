function dxdt = nl_rel_ode(~, X, p)
%nl_rel_ode List of the noninear relative motion differential equations
%   State vector to be integrated to understand the relative motion of two
%   satellite (one chief and one deputy)

    % pull out constants
    mu = p.mu;
    oe = p.oe_c; % orbital elements of the chief
    rectum = oe.a*(1 - (oe.e)^2);
    h = sqrt(rectum * mu);

    % pull out states  
    x = X(1);
    y = X(2);
    z = X(3);
    xd= X(4);
    yd= X(5);
    zd= X(6);
    rc= X(7);
    rcd=X(8);
    
    % intermediate variables
    fd = h/(rc^2);
    rd = norm([(rc+x), y, z]);

    % ODEs to solve
    xdd = 2*fd*(yd - (y*(rcd/rc))) + (x*fd*fd) + (mu/(rc*rc)) + ...
        - ((mu/(rd^3))*(rc + x));
    ydd = -2*fd*(xd - (x*(rcd/rc))) + (y*fd*fd) - ((mu/(rd^3))*y);
    zdd = -((mu/(rd^3))*z);
    rcdd = (rc*fd*fd)*(1-(rc/rectum));

    dxdt =[xd; yd; zd; xdd; ydd; zdd; rcd; rcdd];
end
