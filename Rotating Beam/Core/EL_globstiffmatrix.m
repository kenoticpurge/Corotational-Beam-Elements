function [fl, Kg] = EL_globstiffmatrix(params,B,r,z,Ln,Lo,qbar)

% This function calculates the Local Force Vector and the Global Stiffness
% matrix for the element

% The functions for the Local Force Vector and Local Stiffness Matrix were
% derived in Maple

ubar = qbar(1);
tbar = qbar(2:3);

Omega = params.Omega;
EA = params.EA; EI = params.EI;

% Local Force Vector Function
fl = EL_localforcevector(tbar(1),tbar(2),Lo,ubar,EA,EI,Omega);

% Local Tangent Stiffness Matrix
kl = EL_stiffmatrix(tbar(1),tbar(2),Lo,ubar,EA,EI,Omega);

% Components of the Local Force Vector
N = fl(1);
M1 = fl(2);
M2 = fl(3);

% Calculate the Global Stiffness Matrix
Kg = B'*kl*B + ((z*z')/Ln)*N + (1/Ln^2)*(r*z' + z*r')*(M1 + M2);

end
