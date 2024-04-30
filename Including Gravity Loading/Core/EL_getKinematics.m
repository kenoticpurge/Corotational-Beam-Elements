function [T,B,r,z,Ln,Lo,qbar,c] = EL_getKinematics(q0,q)

% This function establishes the important kinematic parameters of the beam
% element

assert(size(q0,1) == 6, 'q0 needs to be 6x1 in size');
assert(size(q,1) == 6, 'q needs to be 6x1 in size');

%q0 - vector of intial positon [x1 z1 bet0 x2 z2 bet0] 
%q - vector of global displacements [u1 w1 t1 u2 w2 t2]

% want to keep the notation as close a possible to the paper

% Initial Position
x1 = q0(1); z1 = q0(2);
x2 = q0(4); z2 = q0(5);

% Initial Element Length
Lo = sqrt((x2 - x1)^2 + (z2 - z1)^2);
% Initial Angle beta0
co = (x2 - x1)/Lo;
so = (z2 - z1)/Lo;

% Global Displacements
u1 = q(1); u2 = q(4);
w1 = q(2); w2 = q(5);

% Current Length of the Element
Ln = sqrt((x2 + u2 - x1 - u1)^2 + (z2 + w2 - z1 - w1)^2);
% Current Angle beta (in the global configuration)
c = (x2 + u2 - x1 - u1)/Ln;

s = (z2 + w2 - z1 - w1)/Ln;

A = [c s 0;
     -s c 0;
     0 0 1];

T = [A zeros(3,3);
     zeros(3,3) A];

r = [-c -s 0 c s 0]';

z = [s -c 0 -s c 0]';

B = [-c -s 0 c s 0;
     -s/Ln c/Ln 1 s/Ln -c/Ln 0;
     -s/Ln c/Ln 0 s/Ln -c/Ln 1];

% Rigid Rotation Alpha
sa = s*co - c*so;
ca = c*co + s*so;

% Local Rotation Calculations i.e. t1bar, t2bar
t1 = q(3); t2 = q(6);

st1bar = sin(t1)*ca - cos(t1)*sa;
ct1bar = cos(t1)*ca + sin(t1)*sa;

st2bar = sin(t2)*ca - cos(t2)*sa;
ct2bar = cos(t2)*ca + sin(t2)*sa;

% Conditonal based derivation of tbar
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

%t2bar
if st2bar >= 0 && ct2bar >= 0
    t2bar = asin(st2bar);
elseif st2bar < 0 && ct2bar >= 0
    t2bar = asin(st2bar);
elseif st2bar >= 0 && ct2bar < 0
    t2bar = acos(ct2bar);
elseif st2bar < 0 && ct2bar < 0
    t2bar = -acos(ct2bar);
end

ubar = Ln - Lo;

% Current Configuration of the element in the corotational frame
qbar = [ubar t1bar t2bar]';

end







