function [kh, kh1, kh2, eta] = khCalculateRegRot(dl,fl)

% Equation 45 in the Battini paper

for i = 1:2

t = dl(3*i-1:3*i+1); % define for the node
a = norm(t);

    if a == 0 % this is something I'm still unsure on
        kh = zeros(3,3);
        eta = 0;
    else
        e = (2*sin(a) - a*(1 + cos(a)))/(2*a^2*sin(a));
        
        if i == 1
          eta = e;
        end
        
        mu = (a*(a + sin(a)) - 8*sin(a/2)^2)/(4*a^4*sin(a/2)^2);
        
        v = fl(3*i-1:3*i+1);
        
        vsk = skewsymmatrix(v);
        tsk = skewsymmatrix(t);
        
        Tsi = TsInverse(t);
        
        kh = (e * (t*v' - 2*v*t' + dot(t',v)*eye(3)) + mu*tsk^2*(v*t') - 0.5*vsk) * Tsi;
    end

    if i == 1
        kh1 = kh;
    elseif i == 2
        kh2 = kh;
    end

end

kh = [0 zeros(1,3) zeros(1,3);
      zeros(3,1) kh1 zeros(3,3);
      zeros(3,1) zeros(3,3) kh2];

end