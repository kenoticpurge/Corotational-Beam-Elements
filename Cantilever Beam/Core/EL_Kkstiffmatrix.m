function Kk = EL_Kkstiffmatrix(params,ML,qd,qdd,T,B,r,z,Ln,Lo)

% This function calculates the Cetrifugal Stiffness Matrix using Eqn 64 to
% Eqn 70

assert(size(qd,1) == 6, 'qd needs to be 6x1 in size');
assert(size(qdd,1) == 6, 'qdd needs to be 6x1 in size');

% Mass Matrix Derivatives
[M_beta, M_t1bar, M_t2bar, Mlb, Mlt1, Mlt2] = ...
    EL_massderiv(params,ML,T,Lo);

I = reshape([0 1 0 0 0 0; ...
            -1 0 0 0 0 0; ...
             0 0 0 0 0 0; ...
             0 0 0 0 1 0; ...
             0 0 0 -1 0 0; ...
             0 0 0 0 0 0], [6,6]);

M_b_b = T'*(I'*Mlb + Mlb*I)*T;

M_b_t1 = T'*(I'*Mlt1 + Mlt1*I)*T;

M_b_t2 = T'*(I'*Mlt2 + Mlt2*I)*T;

M_t1_b = M_b_t1;

M_t2_b = M_b_t2;

b2 = B(2,:)';
b3 = B(3,:)';

% Centrifugal Stiffness Matrix Componenets
K1 = M_beta*qdd*(z'/Ln) + M_t1bar*qdd*b2' + M_t2bar*qdd*b3';

K2 = ((z'/Ln) * qd) * (M_b_b*qd*(z'/Ln) + M_b_t1*qd*b2' + M_b_t2*qd*b3')...
    + (b2'*qd)*M_t1_b*qd*(z'/Ln) + (b3'*qd)*M_t2_b*qd*(z'/Ln) ...
    - (M_beta - M_t1bar - M_t2bar)*qd*qd'*((r*z' + r*z')/Ln^2);

K3t = (qd'*M_b_b*qd)*((z*z')/Ln^2) + (qd'*M_b_t1*qd)*(z/Ln)*b2' ...
    + (qd'*M_b_t2*qd)*(z/Ln)*b3' + (qd'*M_b_t1*qd)*b2*(z'/Ln) ...
    + (qd'*M_b_t2*qd)*b3*(z'/Ln) - qd'*(M_beta - M_t1bar - M_t2bar)*qd*((r*z' + r*z')/Ln^2);

K3 = 0.5 * K3t;

% Element Centrifugal Stiffness Matrix
Kk = K1 + K2 - K3;

end


