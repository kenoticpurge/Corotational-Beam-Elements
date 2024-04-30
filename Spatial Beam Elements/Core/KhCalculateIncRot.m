function Kh = KhCalculateIncRot(th1, th2, fg)

% Equation 77 in Battini's paper

for i = 1:2

    if i == 1
        th = th1;
    elseif i == 2
        th = th2;
    end

    mm = (i - 1)*6 + 4;
    v = fg(mm:mm+2);

    T = sqrt(th(1)^2 + th(2)^2 + th(3)^2);
    
        if T == 0
            Kh = zeros(3,3); %still not sure about this
        else
            e = th/T;
            vsk = skewsymmatrix(v);
        
            k1 = (sin(T)/T - (sin(T/2)/(T/2))^2) * cross(e,v) * e';
            
            k2 = 0.5 * ((sin(T/2)/(T/2))^2) * vsk;
        
            k3 = ((cos(T) - sin(T)/T)*(1/T)) * (v*e' - (e'*v)*(e*e'));
        
            k4 = ((1 - sin(T)/T)*(1/T)) * (e*v' - 2*(e'*v)*(e*e') + (e'*v)*eye(3,3));
        
            Kh = -k1 + k2 + k3 + k4;
        end

    if i == 1
        Kh1 = Kh;
    elseif i == 2
        Kh2 = Kh;
    end

end

Kh = [zeros(3) zeros(3) zeros(3) zeros(3);
      zeros(3)   Kh1    zeros(3) zeros(3);
      zeros(3) zeros(3) zeros(3) zeros(3);
      zeros(3) zeros(3) zeros(3)   Kh2];

end