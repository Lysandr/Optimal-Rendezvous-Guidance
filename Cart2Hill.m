function [HN] = Cart2Hill(rv)
%Cart2Hill converts a set of R and V vectors into a hill frame DCM
%   hill frame order defined as [or otheta oh]

    % decompose
    rc = rv(1:3).';
    vc = rv(4:6).';
    
    hvec = skew(rc)*vc;
    hhat = hvec./norm(hvec);
    o_r = rc./norm(rc);
    o_h = hhat;
    o_theta = skew(o_h)*o_r;
    HN = [o_r.'; o_theta.'; o_h.'];
    
    
    
%     hvec = skew(rc)*vc;
%     hhat = hvec./norm(hvec);
%     o_r = rc./norm(rc);
%     o_h = hhat;
%     o_theta = skew(o_h)*o_r;
%     HN = [o_r.'; o_theta.'; o_h.'];
%     rho_N = (rd-rc);
%     rho_Hill = HN*rho_N;
end

