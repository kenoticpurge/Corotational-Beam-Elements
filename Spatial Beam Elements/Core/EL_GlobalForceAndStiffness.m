function [fg, Kg] = EL_GlobalForceAndStiffness(Rr,P,E,GT,Ln,Lo,dl,params)

% Unsure:
% eta - is alpha calculated for each node? I think it's for node 1 based on
% a vector


% Gotta be careful with this since there are changes of variables in the
% formulation

% We follow the method presented in the Battini paper

% the formulations for the force and stiffness matrices seem fine

t1bar = dl(2:4); %Assuming these are available
t2bar = dl(5:7);

Ts1 = TsInverse(t1bar);
Ts2 = TsInverse(t2bar);

B1 = [1 zeros(1,3) zeros(1,3);
      zeros(3,1) Ts1 zeros(3,3);
      zeros(3,1) zeros(3,3) Ts2];

fl = EL_LocalForceVector(dl,Lo,params);

% fa in the Battini paper - Eqn 40
fa = B1' * fl;
m = fa(2:7);

% Calculate the local stiffness matrix
kl = EL_LocalStiffMatrix(dl,Lo,params);

[kh,~,~,eta] = khCalculateRegRot(dl,fl);

ka = B1'*kl*B1 + kh; % Equation 43

e1 = Rr(:,1);
r = [-e1' zeros(1,3) e1' zeros(1,3)];

Ba = [r; P*E'];

% Element internal force vector in the global coordinates
fg = Ba'*fa;

% Element stiffness matrix in the global coordinates

D3 = (1/Ln) * (eye(3) - (e1*e1'));

D = [D3 zeros(3) -D3 zeros(3);
     zeros(3) zeros(3) zeros(3) zeros(3);
     -D3 zeros(3) D3 zeros(3);
     zeros(3) zeros(3) zeros(3) zeros(3)];

a = [                  0    ;
     eta*(m(1) + m(4))/Ln - (m(2) + m(5))/Ln;
               (m(3) + m(6))/Ln];

Q = P'* m;

Q1 = Q(1:3); Q2 = Q(4:6); Q3 = Q(7:9); Q4 = Q(10:12);

Q = [skewsymmatrix(Q1);
     skewsymmatrix(Q2);
     skewsymmatrix(Q3);
     skewsymmatrix(Q4);];

Km = D*fa(1) - E*Q*GT*E' + E*GT'*a*r;

Kg = Ba'*ka*Ba + Km; % Element stiffness matrix in the global coordinates

end






