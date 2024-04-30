function [T,B,r,z,Ln,Lo,qbar] = EL_getKinematics(q0,q)

% This function establishes the important kinematic parameters for the beam
% element

% The equations in this function follow the equations in Section 2 of the
% paper
% The paper omits some details. The extra details can be found in the
% thesis by Thanh-Nam Le. % Note that there are errors in paper which have
% been corrected for in this code.

assert(size(q0,1) == 6, 'q0 needs to be 6x1 in size');
assert(size(q,1) == 6, 'q needs to be 6x1 in size');

% q0 - Initial Position of Nodes 
% q - Elastic Displacement Vector [u1 w1 th1 u2 w2 th2]'

% Initial Nodal Position
x1 = q0(1); z1 = q0(2);
x2 = q0(4); z2 = q0(5);

% Initial Element Length - Eqn 6
Lo = sqrt((x2 - x1)^2 + (z2 - z1)^2);
% Initial Angle beta0
co = (x2 - x1)/Lo;
so = (z2 - z1)/Lo;

% Global Elastic Displacements
u1 = q(1); u2 = q(4);
w1 = q(2); w2 = q(5);

% Current Length of the Element - Eqn 7
Ln = sqrt((x2 + u2 - x1 - u1)^2 + (z2 + w2 - z1 - w1)^2); 

% Current Angle beta in the global configuration - Eqns 8 and 9
c = (x2 + u2 - x1 - u1)/Ln;
s = (z2 + w2 - z1 - w1)/Ln;

% Planar Rotation Matrix
A = [c s 0;
     -s c 0;
     0 0 1];

% Element Transformation Matrices and Vectors
T = [A zeros(3,3);
     zeros(3,3) A];

r = [-c -s 0 c s 0]'; % Eqn 17

z = [s -c 0 -s c 0]'; % Eqn 18

B = [-c -s 0 c s 0;
     -s/Ln c/Ln 1 s/Ln -c/Ln 0;
     -s/Ln c/Ln 0 s/Ln -c/Ln 1]; % Eqn 12

% Rigid Rotation Alpha - see thesis for details
sa = s*co - c*so;
ca = c*co + s*so;

% Local Rotation Calculations i.e. t1bar, t2bar - see thesis for details
t1 = q(3); t2 = q(6);

st1bar = sin(t1)*ca - cos(t1)*sa;
ct1bar = cos(t1)*ca + sin(t1)*sa;

st2bar = sin(t2)*ca - cos(t2)*sa;
ct2bar = cos(t2)*ca + sin(t2)*sa;

% Conditonal based derivation of tbar. The if statements are needed due to
% inverse trigonometry being used - see thesis for details
% t1bar
if st1bar >= 0 && ct1bar >= 0
    t1bar = asin(st1bar);
elseif st1bar < 0 && ct1bar >= 0
    t1bar = asin(st1bar);
elseif st1bar >= 0 && ct1bar < 0
    t1bar = acos(ct1bar);
elseif st1bar < 0 && ct1bar < 0
    t1bar = -acos(ct1bar);
end

% t2bar
if st2bar >= 0 && ct2bar >= 0
    t2bar = asin(st2bar);
elseif st2bar < 0 && ct2bar >= 0
    t2bar = asin(st2bar);
elseif st2bar >= 0 && ct2bar < 0
    t2bar = acos(ct2bar);
elseif st2bar < 0 && ct2bar < 0
    t2bar = -acos(ct2bar);
end

ubar = Ln - Lo; % Eqn 3

% Current Configuration of the element in the corotational frame - Eqn 2
qbar = [ubar t1bar t2bar]';

end







