function [ML] = EL_massmatrix(params,Lo,tbar)

% This function calculates the mass matrix of each element in its
% corotational frame

% Note this formulation means that the mass matrix is dependent on the
% configuration of the beam element as opposed to being constant

rho = params.rho;
A = params.A;
I = params.I;
L0 = Lo;

t1bar = tbar(1);
t2bar = tbar(2);

m1 = 21 * t1bar - 14 * t2bar;
m2 = 14 * t1bar - 21 * t2bar;

t1 = 140;
t2 = 70;
t3 = 156;
t4 = 22 * L0;
t5 = 54;
t6 = 13 * L0;
t7 = 4 * L0^2;
t8 = 3 * L0^2;

ML1 = reshape([t1,m1,0,t2,-m1,0,...
              m1,t3,t4,m2,t5,-t6,...
              0,t4,t7,0,t6,-t8,...
              t2,m2,0,t1,-m2,0,...
              -m1,t5,t6,-m2,t3,-t4,...
              0,-t6,-t8,0,-t4,t7],[6,6]);

r1 = 36;
r2 = 3 * L0;
r3 = 4 * L0^2;
r4 = L0^2;

ML2 = reshape([0,0,0,0,0,0,...
               0,r1,r2,0,-r1,r2,...
               0,r2,r3,0,-r2,-r4,...
               0,0,0,0,0,0,...
               0,-r1,-r2,0,r1,-r2,...
               0,r2,-r4,0,-r2,r3],[6,6]);

% Element Mass Matrix in the Local Corotational Coordinate System
ML = ((rho * A * L0)/420) * ML1 + ((rho * I)/(30 * L0)) * ML2;

end