function Ths = EndRotation(t,Ts,ws)

% Function which calculates the end rotation at time t

% Equation 2.1 in the mentioned paper

if t <= Ts

    Ths = (ws/Ts) * (0.5*t^2 + (Ts/(2*pi))^2 * (cos((2*pi*t)/Ts) - 1));

elseif t > Ts

    Ths = ws * (t - Ts/2);

end

end