function Ts = TsCalculate(phi)

% Function to calculate Ts based on the rotation vector

% Equation 12 in the first paper by Batinni

phisk = skewsymmatrix(phi);

p = sqrt(phi(1)^2 + phi(2)^2 + phi(3)^2);

if p == 0

    Ts = eye(3);

else

e = phi/p;

Ts = (sin(p)/p)*eye(3,3) + (1 - (sin(p)/p))*(e*e') + 0.5*((sin(p/2)/(p/2))^2)*phisk;

end

end
