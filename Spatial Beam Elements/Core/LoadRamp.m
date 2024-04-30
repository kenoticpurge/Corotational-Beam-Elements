function P = LoadRamp(t,T,P0)

% Function to calculate the load at each pseudo time step

% The load is treated as a ramp which is what ABAQUS uses for its
% formulation

if t < T

    P = (t/T) * P0;

else
    
    P = P0;

end