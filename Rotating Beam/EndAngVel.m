function ts_t = EndAngVel(t, ws, Ts)

% Function which calculates the end angular velocity at time t

if t <= Ts

    ts_t = (ws/Ts) * (t - (Ts/2*pi)*sin(2*pi*t/Ts));

elseif t > Ts
    
    ts_t = ws;

end