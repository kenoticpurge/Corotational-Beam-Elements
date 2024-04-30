function kl = EL_stiffmatrix(tl1,tl2,L,u,EA,EI,Omega)

% Function to calculate the local stiffness matrix for each element

% This function corresponds to the Kl term in Eqn 16

% This function was generated using Maple - a software for symbollic
% programming

% While the Maple code is provided in the paper (see the appendix), there is an error.
% The error has been adjusted for in this code

% Error in provided Maple Code : "t + dw" should be "dw - t"

t1 = tl2 - tl1;
t2 = 0.1e1 / 0.90e2;
t3 = 0.1e1 / 0.360e3;
t4 = t1 * Omega * (-Omega - 0.1e1 / 0.6e1);
t5 = 0.1e1 + 0.12e2 * Omega;
t6 = 0.1e1 / L;
t5 = 0.1e1 / t5 ^ 2;
t7 = t5 ^ 2;
t8 = 0.12e2 * EA * (t2 * tl1 - t3 * tl2 + t4) * t5;
t2 = (-0.12e2 * t2 * tl2 + 0.12e2 * t3 * tl1 + 0.12e2 * t4) * EA * t5;
t1 = t1 ^ 2;
t3 = EA * L;
t4 = (6220800 * EI) + t3 * (0.64800e5 * t1 * L + (518400 * u));
t5 = tl1 ^ 2;
t9 = 3240;
t10 = 23760 * u;
t11 = 648000 * EI;
t12 = 240;
t13 = 1560 * u;
t14 = 46800 * EI;
t15 = (Omega * t4 + (3628800 * EI) + t3 * (0.21600e5 * t1 * L + (172800 * u))) * Omega;
t16 = 8;
t17 = 40 * u;
t18 = 1200 * EI;
t19 = tl2 ^ 2 + t5;
t20 = 0.1e1 / 0.150e3;
t21 = 0.1e1 / 0.300e3;
t1 = t20 * ((300 * EI) + (((-Omega * t4 / 0.2e1 - (259200 * EI) + t3 * (-0.10800e5 * t1 * L - (86400 * u))) * Omega + (64800 * EI) + t3 * (-(9720 * u) + (0.2520e4 * tl1 * tl2 - 0.1350e4 * t19) * L)) * Omega + (9000 * EI) + t3 * (-(420 * u) + (0.120e3 * tl1 * tl2 - 0.75e2 * t19) * L)) * Omega - EA * (-0.3e1 * tl1 * tl2 + t19) * L ^ 2 - 0.5e1 * t3 * u) * t6 * t7;
t4 = 0.2520e4;
t19 = 0.120e3;

% Local Element Stiffness Matrix - Kl in Eqn 16
kl = [EA * t6 t8 -t2; t8 t21 * (((t3 * (t9 * ((0.7e1 / 0.9e1 * tl2 - 0.5e1 / 0.3e1 * tl1) * tl2 + t5) * L + t10) + t11 + t15) * Omega + t3 * (t12 * ((tl2 / 0.2e1 - 0.5e1 / 0.4e1 * tl1) * tl2 + t5) * L + t13) + t14) * Omega + t3 * (t16 * ((0.3e1 / 0.8e1 * tl2 - tl1 / 0.2e1) * tl2 + t5) * L + t17) + t18) * t6 * t7 t1; -t2 t1 t21 * (((t3 * (t4 * ((0.9e1 / 0.7e1 * tl2 - 0.15e2 / 0.7e1 * tl1) * tl2 + t5) * L + t10) + t11 + t15) * Omega + t3 * (t19 * ((-0.5e1 / 0.2e1 * tl1 + 0.2e1 * tl2) * tl2 + t5) * L + t13) + t14) * Omega + t3 * ((0.3e1 * (0.8e1 / 0.3e1 * tl2 - 0.4e1 / 0.3e1 * tl1) * tl2 + 0.3e1 * t5) * L + t17) + t18) * t6 * t7;];


end