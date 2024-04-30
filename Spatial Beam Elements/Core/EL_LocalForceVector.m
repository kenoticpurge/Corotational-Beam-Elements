function fl = EL_LocalForceVector(dl,Lo,params)

% Function to obtain the local element force vector obtained using
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
t2 = (t1 ^ 2);
t3 = (t2 * (th31 + th32));
t2 = (t2 * (th21 + th22));
t4 = -t1 * (4 * th31 + 2 * th32);
t5 = t1 * (4 * th21 + 2 * th22);
t6 = 6;
t7 = t1 * (th11 - th12);
t8 = 1 / A;
t9 = t7 ^ 2;
t10 = 0.1e1 / 0.3e1;
t11 = 0.3e1 / 0.2e1;
t12 = 0.9e1 / 0.5e1;
t13 = L ^ 2;
t14 = 0.1e1 / 0.2e1;
t11 = t1 * (t14 * ((((-t2 * t5 + t3 * t4) * t11 + (t2 ^ 2 + t3 ^ 2) * t12 * L) * L + (t6 * (t2 * th21 + t3 * th31) + t4 ^ 2 + t5 ^ 2) * t10) * L * t13 + ((t4 * th31 - t5 * th21) * L + (th21 ^ 2) + (th31 ^ 2) + (Io * t9 * t8)) * L) + u);
t8 = (t7 * t1 * (-(t9 * (Io ^ 2 * t8 - Irr)) + 0.2e1 * t11 * Io));
t7 = G * J * t7;
t9 = (t14 * E);
t12 = -t5;
t15 = -48;
t16 = 0.12e2 * t12 * t1;
t17 = Izz * t1;
t18 = t6 * t5 * t1;
t19 = -24 * t2;
t20 = 0.18e2 / 0.5e1 * t2;
t21 = 2 * th21;
t11 = A * t11;
t22 = 12 * t2;
t23 = t22 * Izz;
t24 = L * E;
t25 = -t4;
t26 = 12 * t25 * t1;
t27 = Iyy * t1;
t28 = t6 * t4 * t1;
t29 = -24 * t3;
t30 = 0.18e2 / 0.5e1 * t3;
t31 = 2 * th31;
t32 = 12 * t3;
t33 = t32 * Iyy;
t34 = t6 * t1;
t32 = (t24 * (t14 * t1 * ((4 * Iyy * t25) + t11 * (t13 * (t10 * t1 * (t34 * th31 - 4 * t4) + t30 + t1 * (-t32 + t28) * L / 0.4e1) - t31) * L) + t33 + t27 * (-t26 + t29) * L / 0.4e1));

fl = [t11 * E t9 * t8 * L + t7 t24 * (t14 * t1 * (-0.8e1 * Izz * t12 + t11 * (((t10 * (t1 * (t6 * th22 + 12 * th21) + 8 * t5) * t1 + t20 + t1 * (t19 - t18) * L / 0.4e1) * L + t14 * t1 * (-16 * th21 - 4 * th22)) * L + t21) * L) + t23 + t17 * ((t15 * t2) + t16) * L / 0.4e1) t24 * (t14 * t1 * ((8 * Iyy * t25) + t11 * (((t10 * (t1 * (t6 * th32 + 12 * th31) - 8 * t4) * t1 + t30 + t1 * (t28 + t29) * L / 0.4e1) * L + t14 * t1 * (-16 * th31 - 4 * th32)) * L + t31) * L) + t33 + t27 * (t15 * t3 - t26) * L / 0.4e1) -t9 * t8 * L - t7 t24 * (t14 * t1 * (-0.4e1 * Izz * t12 + t11 * (t13 * (t10 * t1 * (t34 * th21 + 4 * t5) + t20 - t1 * (t22 + t18) * L / 0.4e1) - t21) * L) + t23 + t17 * (t19 + t16) * L / 0.4e1) t32]';

end

