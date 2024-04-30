function [dl,Ln,Lo,GT,E,P,Rr,R1g,R2g,r] = EL_getKinematics3D_IncRot(u0,u,th,Rgn,R0,params)

% In this function we get all the necessary kinematic data of the 3D beam
% element including the Local Nodal Displacement Vector

% The details of the expressions can be found in Section 3 of the paper
% "Co-rotational beam elements with warping eﬀects
% in instability problems"

x1g = u0(1:3); u1g = u(1:3);
x2g = u0(4:6); u2g = u(4:6);

Lo = params.Le; % Initial Element Length
Ln = norm(x2g + u2g - x1g - u1g); % Current Element Length

% Calculate the Global Rotation matrix for each node at the next increment
% based on the previous rotation matrix and the global incremental rotation
% vector
th1 = th(1:3); th2 = th(4:6);
R1gn = Rgn{1}; R2gn = Rgn{2};
% Update nodal rotation matrix based on Equation 81
R1g = expm(skewsymmatrix(th1))*R1gn;
R2g = expm(skewsymmatrix(th2))*R2gn;

% Calculate p - we set R0 before hand to describe the orientation of the
% local frame at the initial configuration
p1  = R1g * R0 * [0 1 0]';
p2 = R2g * R0 * [0 1 0]';

p = 0.5 * (p1 + p2);

% Calculate the components of Rr

r1 = (x2g + u2g - x1g - u1g)/Ln;
r3 = cross(r1,p)/norm(cross(r1,p));
r2 = cross(r3,r1);

r = [-r1' zeros(1,3) r1' zeros(1,3)];

Rr = [r1 r2 r3]; % Rigid Rotation matrix

% Calculate bar matrices i.e. the local rotations of the nodes - Eqn 18
R1bar = Rr' * R1g * R0;
R2bar = Rr' * R2g * R0;

% Now compute the components of the local displacement vector
ubar = Ln - Lo;
t1bar = logm(R1bar);
t1bar = [t1bar(3,2); t1bar(1,3); t1bar(2,1)];
t2bar = logm(R2bar);
t2bar = [t2bar(3,2); t2bar(1,3); t2bar(2,1)];

% Local Displacement vector components are needed for the local stiffness
% and force vector calculations
dl = [ubar t1bar' t2bar']';
assert(size(dl,1) == 7,'dl is not a 7 x 1 vector');

% Compute other important matrices etc.

E = [Rr zeros(3) zeros(3) zeros(3);
     zeros(3) Rr zeros(3) zeros(3);
     zeros(3) zeros(3) Rr zeros(3);
     zeros(3) zeros(3) zeros(3) Rr];

% Vectors and components to obtain GT and P
gt = Rr' * p; %pj as components
gt1 = Rr' * p1; %pij as components
gt2 = Rr' * p2;

e = gt(1)/gt(2); 
e11 = gt1(1)/gt(2); e12 = gt1(2)/gt(2);
e21 = gt2(1)/gt(2); e22 = gt2(2)/gt(2);

GT = [0 0 e/Ln e12/2 -e11/2 0 0 0 -e/Ln e22/2 -e21/2 0;
      0  0 1/Ln   0      0   0 0 0 -1/Ln   0      0   0;
      0 -1/Ln 0 0 0 0 0 1/Ln 0 0 0 0];

P = ([zeros(3) eye(3) zeros(3) zeros(3); ...
      zeros(3) zeros(3) zeros(3) eye(3);] - [GT; GT]);

end
