function kl = EL_LocalStiffMatrix(dl,Lo,params)

% Function to obtain the local element stiffness matrix obtained using
% Maple

E = params.E;
A = params.A;
G = params.G;
J = params.J;
Io = params.Io; Iyy = params.Iyy; Izz = params.Izz; Irr = params.Irr;
L = Lo;

u = dl(1);
th11 = dl(2); th21 = dl(3); th31 = dl(4);
th12 = dl(5); th22 = dl(6); th32 = dl(7);

t1 = 1 / L;
t2 = (t1 * (th11 - th12));
t3 = t1 ^ 2;
t4 = t3 * (th21 + th22);
t5 = t1 * (4 * th21 + 2 * th22);
t6 = -24;
t7 = 6 * t5 * t1;
t8 = 0.18e2 / 0.5e1 * t4;
t9 = 2 * th21;
t10 = (((((t1 * (12 * th21 + 6 * th22) + 8 * t5) * t1) / 0.3e1 + t8 + t1 * (t4 * t6 - t7) * L / 0.4e1) * L - (t1 * (16 * th21 + 4 * th22)) / 0.2e1) * L + t9) * L;
t11 = t3 * (th31 + th32);
t12 = -t1 * (4 * th31 + 2 * th32);
t13 = 6 * t12 * t1;
t14 = 0.18e2 / 0.5e1 * t11;
t15 = 2 * th31;
t6 = ((((((t1 * (12 * th31 + 6 * th32) - 8 * t12) * t1) / 0.3e1 + t14 + t1 * (t11 * t6 + t13) * L / 0.4e1) * L - (t1 * (16 * th31 + 4 * th32)) / 0.2e1) * L + t15) * L);
t16 = 6 * t1;
t17 = L ^ 2;
t7 = ((t17 * ((t1 * (t16 * th21 + 4 * t5)) / 0.3e1 + t8 - t1 * (0.12e2 * t4 + t7) * L / 0.4e1) - t9) * L);
t8 = (t17 * ((t1 * (t16 * th31 - 4 * t12)) / 0.3e1 + t14 + t1 * (-0.12e2 * t11 + t13) * L / 0.4e1) - t15) * L;
t9 = E * t1;
t13 = t9 * A;
t14 = t13 / 0.2e1;
t15 = (t14 * t10);
t18 = t14 * t6;
t19 = t14 * t7;
t14 = t14 * t8;
t9 = t9 * Io * t2;
t20 = 1 / A;
t2 = t2 ^ 2;
t21 = Io * t20 * t2;
t22 = 0.3e1 / 0.2e1;
t23 = 0.9e1 / 0.5e1;
t4 = t1 * (u + L * t17 * (((t11 * t12 - t4 * t5) * t22 + (t11 ^ 2 + t4 ^ 2) * t23 * L) * L + 0.2e1 * t11 * th31 + (t12 ^ 2) / 0.3e1 + 0.2e1 * t4 * th21 + (t5 ^ 2) / 0.3e1) / 0.2e1 + ((t12 * th31 - t5 * th21) * L + (th21 ^ 2) + (th31 ^ 2) + t21) * L / 0.2e1);
t5 = -3;
t2 = (0.2e1 * t3 * Io * (t4 + t21) + t5 * (Io ^ 2 * t20 - Irr) * t2 * t3);
t5 = G * J * t1;
t11 = E / 0.2e1;
t12 = (t11 * t2 * L + t5);
t2 = (-t11 * t2 * L - t5);
t5 = (t9 / 0.2e1);
t17 = (t5 * t6);
t20 = t5 * t7;
t21 = (t5 * t8);
t5 = (t5 * t10);
t22 = 0.4e1 / 0.15e2 * A * t4;
t4 = -A * t4 / 0.15e2;
t23 = E * ((t3 * ((16 * Izz) + A * t10 * t7 / 0.2e1) / 0.2e1 + t4 / 0.2e1) * L - (t16 * Izz));
t24 = t13 / 0.4e1;
t25 = t24 * t6;
t26 = t25 * t10;
t24 = t24 * t8;
t27 = t24 * t10;
t4 = E * ((t3 * ((16 * Iyy) + A * t6 * t8 / 0.2e1) / 0.2e1 + t4 / 0.2e1) * L - (t16 * Iyy));
t16 = (t25 * t7);
t24 = t24 * t7;

kl = [t13 t9 t15 t18 -t9 t19 t14; t9 t12 t5 t17 t2 t20 t21; t15 t5 E * ((t3 * ((32 * Izz) + A * t10 ^ 2 / 0.2e1) / 0.2e1 + t22 / 0.2e1) * L - (12 * Izz * t1)) t26 -t5 t23 t27; t18 t17 t26 E * ((t3 * ((32 * Iyy) + (A * t6 ^ 2) / 0.2e1) / 0.2e1 + t22 / 0.2e1) * L - (12 * Iyy * t1)) -t17 t16 t4; -t9 t2 -t5 -t17 t12 -t20 -t21; t19 t20 t23 t16 -t20 t11 * (t3 * ((8 * Izz) + (A * t7 ^ 2) / 0.2e1) + t22) * L t24; t14 t21 t27 t4 -t21 t24 t11 * (t3 * ((8 * Iyy) + A * t8 ^ 2 / 0.2e1) + t22) * L;];

end
