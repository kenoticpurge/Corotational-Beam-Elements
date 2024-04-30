function fK = EL_InertiaForceVector(params,ML,T,B,z,Ln,Lo,M,qd,qdd)

% The function calculates the Inertia Force Vector of each element

assert(size(qd,1) == 6, 'qd needs to be 6x1 in size');
assert(size(qdd,1) == 6, 'qdd needs to be 6x1 in size');

% We evaluate the inertia force vector for each element then assemble it
% for the overall structure in matrix assembly

% This force vector is required to evaluate the dynamic force balance
% during the NR iterations

% Mass matrix derivatives which are required for fK
[M_beta, M_t1bar, M_t2bar, ~, ~, ~] = EL_massderiv(params,ML,T,Lo);

fk1 = M * qdd;

b2 = B(2,:)';
b3 = B(3,:)';

fk2 = (M_beta*((z'/Ln)*qd) + M_t1bar*(b2'*qd) + M_t2bar*(b3'*qd)) * qd;

fk3 = ((0.5*qd'*M_beta*qd)*(z/Ln)) + ((0.5*qd'*M_t1bar*qd)*b2) + ((0.5*qd'*M_t2bar*qd)*b3);

fK = fk1 + fk2 - fk3;

end