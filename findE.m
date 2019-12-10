function [nu,En_plus_1] = findE(t,e,n,Mo)
    if nargin == 4
        M = Mo + n*t;
    else
        M = n*t;
    end
    if e <= 1 % elliptical case
        En_plus_1 = M;
    % this was some weird sh Vallado recommends. I don't.
    %     if ( M >-pi && M<0) || M>pi
    %         En_plus_1 = M-e;     
    %     else
    %         En_plus_1 = M+e;
    %     end 
        En = 1E6; % make arbitrarily large
        while (abs(En_plus_1 - En) > 1e-12)
            En = En_plus_1;
            En_plus_1 = En - ((-M + En - e*sin(En))/(1-(e*cos(En))));
        end
        nu = 2*(atan2(real(sqrt(1+e)*(sin(En_plus_1/2))), ...
            real(sqrt(1-e)*(cos(En_plus_1/2)))));
    else
        Hn_plus_1 = M;
        Hn = 1E6; % make arbitrarily large
        while (abs(Hn_plus_1 - Hn) > 1e-12)
            Hn = Hn_plus_1;
            Hn_plus_1 = Hn - ((e*sinh(Hn)- Hn - M)/((e*cosh(Hn)-1)));
        end
        nu = 2*(atan2(real(sqrt(1+e)*(sinh(Hn_plus_1/2))), ...
            real(sqrt(e-1)*(cosh(Hn_plus_1/2)))));
        En_plus_1 = Hn_plus_1;
    end
end



