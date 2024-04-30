function fl = EL_localforcevector(tl1,tl2,L,u,EA,EI,Omega)

% Function to calculate the local force vector for each element

% This function corresponds to the fl term in Eqn 14

% This function was generated using Maple - a software for symbollic
% programming

% While the Maple code is provided in the paper (see the appendix), there is an error.
% The error has been adjusted for in this code

% Error in provided Maple Code : "t + dw" should be "dw - t"

t1 = (tl2 - tl1);
t2 = tl1 ^ 2;
t3 = tl2 ^ 2 + t2;
t4 = t1 ^ 2;
t5 = Omega + 0.1e1 / 0.12e2;
t6 = 24 * u;
t7 = 0.12e2 * t5;
t8 = (t4 * L);
t9 = (EA * L);
t10 = t9 * t1 * (518400 * u + 21600 * t8);
t11 = tl2 / 0.2e1;
t12 = -t11 * tl1 + t3;
t13 = 3600;
t14 = Omega ^ 2;
t8 = (0.64800e5 * t1 * (288 * EI + t9 * (t6 + t8)) * t14 ^ 2);
t14 = -0.30e2;
t15 = 0.1e1 / L;
t16 = 0.1e1 / 0.900e3;
t7 = 0.1e1 / t7 ^ 2;
t17 = t7 ^ 2;
t11 = t16 * (-t8 + t9 * (tl1 - tl2 / 0.4e1) * (0.8e1 * t12 * L + (120 * u)) + t13 * (t11 + tl1) * EI + (((-t10 + 0.10886400e8 * EI * (tl1 - tl2 / 0.7e1)) * Omega + t9 * ((0.71280e5 * tl1 - 0.58320e5 * tl2) * u - 0.3240e4 * t1 * ((0.5e1 / 0.6e1 * tl2 - 0.3e1 / 0.2e1 * tl1) * tl2 + t2) * L) + 0.1944000e7 * EI * (tl1 + tl2 / 0.5e1)) * Omega + (0.140400e6 * tl1 + 0.54000e5 * tl2) * EI + t9 * (-0.240e3 * t1 * ((0.5e1 / 0.8e1 * tl2 - 0.7e1 / 0.8e1 * tl1) * tl2 + t2) * L + (0.4680e4 * tl1 - 0.2520e4 * tl2) * u)) * Omega) * t15 * t17;

fl = [0.6e1 * EA * ((t3 / 0.90e2 - tl1 * tl2 / 0.180e3 + t4 * Omega * (Omega + 0.1e1 / 0.6e1)) * L + t6 * t5 ^ 2) * t7 * t15 t11 t16 * (t8 + t9 * (tl1 - 0.4e1 * tl2) * (-0.2e1 * t12 * L + t14 * u) + (((-0.1555200e7 * EI * (tl1 - 0.7e1 * tl2) + t10) * Omega + t9 * (0.2700e4 * t1 * ((0.6e1 / 0.5e1 * tl2 - 0.9e1 / 0.5e1 * tl1) * tl2 + t2) * L + (-0.58320e5 * tl1 + 0.71280e5 * tl2) * u) + 0.388800e6 * EI * (tl1 + 0.5e1 * tl2)) * Omega + t9 * (0.150e3 * t1 * ((0.8e1 / 0.5e1 * tl2 - 0.7e1 / 0.5e1 * tl1) * tl2 + t2) * L + (-0.2520e4 * tl1 + 0.4680e4 * tl2) * u) + 0.54000e5 * EI * (tl1 + 0.13e2 / 0.5e1 * tl2)) * Omega + 0.1800e4 * EI * (tl1 + 0.2e1 * tl2)) * t15 * t17]';

end