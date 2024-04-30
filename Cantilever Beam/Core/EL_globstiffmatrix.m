function [fl, Kg] = EL_globstiffmatrix(params,B,r,z,Ln,Lo,qbar)

% This function calculates the Local Force Vector and the Global Stiffness
% matrix for the element

% This function corresponds to Eqn 14 and Eqn 16

% While the Maple code is provided in the paper (see the appendix), there is an error.
% The error has been adjusted for in this code

% Error in provided Maple Code : "t + dw" should be "dw - t"

% Establish Element Parameters
E = params.E;
A = params.A;
I = params.I;
Omega = params.Omega;

% Local Elastic Displacements
ubar = qbar(1);
tbar = qbar(2:3);

EA = E*A; EI = E*I;

% Local Force Vector Function - Eqn 14
fl = EL_localforcevector(tbar(1),tbar(2),Lo,ubar,EA,EI,Omega);

% Local Tangent Stiffness Matrix - Kl in Eqn 16
kl = EL_stiffmatrix(tbar(1),tbar(2),Lo,ubar,EA,EI,Omega);

% Components of the Local Force Vector
N = fl(1);
M1 = fl(2);
M2 = fl(3);

% Calculate the Global Element Stiffness Matrix - Eqn 16 
Kg = B'*kl*B + ((z*z')/Ln)*N + (1/Ln^2)*(r*z' + z*r')*(M1 + M2);

end
