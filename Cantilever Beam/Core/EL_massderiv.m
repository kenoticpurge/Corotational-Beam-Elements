function [M_beta, M_t1bar, M_t2bar, Mlb, Mlt1, Mlt2] = EL_massderiv(params,ML,T,Lo)

% This function calculates the derivatives of the mass matrix

% These derivatives are the ones seen in Eqns 54 to 56

% Since the mass matrix is a function of the state, derivatives are
% required for some of the calculations

rho = params.rho;
A = params.A;
L0 = Lo;
    
% Expression for M_beta - Eqn 52 and 54

I = reshape([0 1 0 0 0 0; ...
            -1 0 0 0 0 0; ...
             0 0 0 0 0 0; ...
             0 0 0 0 1 0; ...
             0 0 0 -1 0 0; ...
             0 0 0 0 0 0], [6,6]); 

Mlb = (I'*ML + ML*I);

M_beta = T' * Mlb * T;

% Expression for M_t1bar - Eqn 55

c1 = (rho * A * L0)/60;

Mlt1 = c1 * reshape([0 3 0 0 -3 0; ...
                     3 0 0 2 0 0; ...
                     0 0 0 0 0 0; ...
                     0 2 0 0 -2 0; ...
                     -3 0 0 -2 0 0; ...
                     0 0 0 0 0 0],[6,6]);

M_t1bar = T' * Mlt1 * T;

% Expression for M_t2bar - Eqn 56

c2 = (rho * A * L0)/60;

Mlt2 = c2 * reshape([0 -2 0 0 2 0; ...
                     -2 0 0 -3 0 0; ...
                     0 0 0 0 0 0; ...
                     0 -3 0 0 3 0; ...
                     2 0 0 3 0 0; ...
                     0 0 0 0 0 0],[6,6]);

M_t2bar = T' * Mlt2 * T;

end

 
