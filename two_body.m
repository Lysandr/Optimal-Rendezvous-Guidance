function dxdt = two_body(~, x, p)
% twobody function for two satellites
	r_v_1 = x(1:3);
    r_v_2 = x(4:6);
	v_v_1 = x(7:9);
    v_v_2 = x(10:12);
	r1 = norm(r_v_1);
    r2 = norm(r_v_2);

    % reveal the deputy control
	a_v_1 = -(p.mu/r1^3)*r_v_1;
	a_v_2 = -(p.mu/r2^3)*r_v_2 + p.u_d;
    
	dxdt =[v_v_1; v_v_2; a_v_1; a_v_2];
end
