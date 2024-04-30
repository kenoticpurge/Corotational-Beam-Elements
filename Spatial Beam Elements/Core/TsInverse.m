function Ts_in = TsInverse(theta)

% Function to calculate the inverse of Ts for a given rotation vector
% Inverse of Ts is a matrix. Some calculations require the transpose of the
% inverse so keep that in mind

% Eqn 14 in the first paper by Batinni

% theta must be of size 3 x 1
assert(size(theta,1) == 3,"Input Vector is wrong size");

t = (theta'*theta)^0.5;


if t == 0

    Ts_in = eye(3);

else

    tsk = [0 -theta(3) theta(2);
           theta(3) 0 -theta(1);
           -theta(2) theta(1) 0];
    
    
    T1 = ((t/2)/tan(t/2)) * eye(3);
    
    T2 = (1 - ((t/2)/tan(t/2))) * (theta*theta')/t^2;
    
    T3 = 0.5 * tsk;
    
    Ts_in = T1 + T2 - T3;

end

end