function [M, Ck, Kk, Kg, Fg, FK] = matrix_assembly(q0,q,qd,qdd,params)

% This function assembles the matrices and vectors for the overall
% beam structure using standard Finite Element Assembly procedure

% M - Mass Matrix
% Ck, Kk - Gyroscopic Matrices
% Kg - Global Stiffness Matrix
% Fg - Global Elastic Force Vector
% FK - Global Inertia Force Vector

% q0 - Initial Position of Nodes
% q - Elastic Displacement Vector
% qd - Elastic Velocity Vector
% qdd - Elastic Acceleration Vector
% params - Beam parameters

% Matrix reset to stop redundant addtion
Kk = zeros(params.n);
M = zeros(params.n);
Ck = zeros(params.n);
Kg = zeros(params.n);
Fg = zeros(params.n,1); 
FK = zeros(params.n,1); 

% Loop over all elements
for ie = 1:params.ne
    % Find range of nodal coordinates for this element
    istart = 3 * ie - 2;
    iend = 3 * ie + 3;
    % Isolating important variables for each element
    q0e = q0(istart:iend); qe = q(istart:iend);
    qde = qd(istart:iend); qdde = qdd(istart:iend);

    % Obtaining the important kinematic variables based on the elements
    % configuration - Section 2
    [T,B,r,z,Ln,Lo,qbar] = EL_getKinematics(q0e,qe);
    
    tbar = qbar(2:3);

    % Obtain the Element Mass matrix in the Corotational Frame - Eqn 45
    ML = EL_massmatrix(params,Lo,tbar);
    
    % Transform Element Mass Matrix into the Global Frame - Eqn 43
    Me = T' * ML * T;

    % Element Inertia Force Vector - Eqn 51
    fKe = EL_InertiaForceVector(params,ML,T,B,z,Ln,Lo,Me,qde,qdde);

    % Obtain the Global Stiffness Matrix and the Local Elastic Force Vector for the element - Eqn 14 and 16
    [fle, Kge] = EL_globstiffmatrix(params,B,r,z,Ln,Lo,qbar);
    % Element Global Elastic Force Vector - Eqn 14
    fge = B'*fle; 

    % Obtain the inertia associated Ck matrix - Gyroscopic Matrix - for the
    % element - Eqn 62
    Cke = EL_Ckmatrix(params,qde,T,B,z,Ln,Lo,ML);

    % Obtain the inertia associated Kk matrix - Centrifugal Matrix - for
    % the element - Eqn 67
    Kke = EL_Kkstiffmatrix(params,ML,qde,qdde,T,B,r,z,Ln,Lo);

    % Matrix partitioning
    Kg00 = Kge(1:3,1:3); Kg01 = Kge(1:3,4:6); % Kge = [Kg00 Kg01;  6 x 6
    Kg10 = Kge(4:6,1:3); Kg11 = Kge(4:6,4:6); %       [Kg10 Kg11]

    Kk00 = Kke(1:3,1:3); Kk01 = Kke(1:3,4:6); % Kk = [Kk00 Kk01;  6 x 6
    Kk10 = Kke(4:6,1:3); Kk11 = Kke(4:6,4:6); %      [Kk10 Kk11]

    M00 = Me(1:3,1:3); M01 = Me(1:3,4:6); % Me = [M00 M01; 6 x 6
    M10 = Me(4:6,1:3); M11 = Me(4:6,4:6); %      [M10 M11]

    Ck00 = Cke(1:3,1:3); Ck01 = Cke(1:3,4:6); % Ck = [Ck00 Ck01; 6 x 6
    Ck10 = Cke(4:6,1:3); Ck11 = Cke(4:6,4:6); %      [Ck10 Ck11]

    % Assemble matrices of the overall structure using Standard Finite
    % Element Assembly Procedure
    ii = (ie-1)*3 + 1;

    % Mass Matrix Assembly
    M(ii:ii+2,ii+3:ii+5) = M01; M(ii+3:ii+5,ii:ii+2) = M10;
    M(ii:ii+2,ii:ii+2) = M(ii:ii+2,ii:ii+2) + M00;
    M(ii+3:ii+5,ii+3:ii+5) = M(ii+3:ii+5,ii+3:ii+5) + M11;

    % Global Stiffness Matrix Assembly
    Kg(ii:ii+2,ii+3:ii+5) = Kg01; Kg(ii+3:ii+5,ii:ii+2) = Kg10;
    Kg(ii:ii+2,ii:ii+2) = Kg(ii:ii+2,ii:ii+2) + Kg00;
    Kg(ii+3:ii+5,ii+3:ii+5) = Kg(ii+3:ii+5,ii+3:ii+5) + Kg11;

    % Inertia Centrifugal Matrix Assembly
    Kk(ii:ii+2,ii+3:ii+5) = Kk01; Kk(ii+3:ii+5,ii:ii+2) = Kk10;
    Kk(ii:ii+2,ii:ii+2) = Kk(ii:ii+2,ii:ii+2) + Kk00;
    Kk(ii+3:ii+5,ii+3:ii+5) = Kk(ii+3:ii+5,ii+3:ii+5) + Kk11;
    % 
    % Inertia Gyroscopic Matrix Assembly
    Ck(ii:ii+2,ii+3:ii+5) = Ck01; Ck(ii+3:ii+5,ii:ii+2) = Ck10;
    Ck(ii:ii+2,ii:ii+2) = Ck(ii:ii+2,ii:ii+2) + Ck00;
    Ck(ii+3:ii+5,ii+3:ii+5) = Ck(ii+3:ii+5,ii+3:ii+5) + Ck11;

    % Global Elastic Force Vector
    Fg(istart:iend) = Fg(istart:iend) + fge;
    
    % Global Inertia Force Vector
    FK(istart:iend) = FK(istart:iend) + fKe;

end

end








