function ts_tt = EndAngAcc(t, ws, Ts)

% Function which calculates the end angular acceleration at time t

if t <= Ts
    
    ts_tt = (ws/Ts) * (1 - cos(2*pi*t/Ts));
        
elseif t > Ts
    
    ts_tt = 0;

end

end