function CK = EL_Ckmatrix(params,qd,T,B,z,Ln,Lo,ML)

% This function calculates the Gyroscopic Matrix Ck

[M_beta, M_t1bar, M_t2bar,~,~,~] = EL_massderiv(params,ML,T,Lo);

b2 = B(2,:)';
b3 = B(3,:)';

% Again, I don't think the transposes are right in the paper

Mdot = M_beta*((z'/Ln)*qd) + M_t1bar*(b2'*qd) + M_t2bar*(b3'*qd);

C1 = M_beta*(qd*(z'/Ln)) + M_t1bar*(qd*b2') + M_t2bar*(qd*b3');

CK = Mdot + C1 - C1';

end