function [Kg, Fg, Rg] = MatrixAssembly(params,u0,R0,u,th,Rg,Rgn)

Kg = zeros(params.n); % Reset matrices and vectors to prevent unneccessary addition
Fg = zeros(params.n,1);
    
  % Loop over all elements
  for ie = 1:params.ne
        istart = 3*ie - 2; iend = 3*ie + 3;
        
        % Extract important terms for each element
        u0e = u0(istart:iend); 
        ue = u(istart:iend); the = th(istart:iend);
        % Nodal Rotation matrix for each element
        Rgne = Rgn(ie:ie+1);
        
        % Obtain important Kinematic terms - Section 3 from paper
        [dl,Ln,Lo,GT,E,P,Rr,R1g,R2g,~] = EL_getKinematics3D_IncRot(u0e,ue,the,Rgne,R0,params);

        % Storing the updated global nodal rotation matrices for each
        % element
        Rg(ie:ie+1) = {R1g ; R2g};

        % Calculate the Global Force Vector and Global Stiffness Matrix for
        % each element
        [fge_loc, kge_loc] = EL_GlobalForceAndStiffness(Rr,P,E,GT,Ln,Lo,dl,params);

        % Obtain theta for the 2 nodes of each element
        th1 = the(1:3);
        th2 = the(4:6);
        Ts1 = TsCalculate(th1);
        Ts2 = TsCalculate(th2);
        
        % Change of variables transformation - Eqn 80 and Eqn 77
        H = [eye(3) zeros(3) zeros(3) zeros(3);
             zeros(3) Ts1 zeros(3) zeros(3);
             zeros(3) zeros(3) eye(3) zeros(3);
             zeros(3) zeros(3) zeros(3) Ts2];

        K_h = KhCalculateIncRot(th1,th2,fge_loc);

        % Element force and stiffness matrix after change of variables -
        % Eqn 80
        fg_glob = H'*fge_loc;
        kg_glob = H'*kge_loc*H + K_h;

        % Matrix Assembly using Standard Finite Element Procedure
        Kg00 = kg_glob(1:6,1:6) ; Kg01 = kg_glob(1:6,7:12);
        Kg10 = kg_glob(7:12,1:6) ; Kg11 = kg_glob(7:12,7:12);

        ii = (ie - 1)*6 + 1;
    
        % Global Stiffness Matrix for whole structure
        Kg(ii:ii+5,ii+6:ii+11) = Kg01; Kg(ii+6:ii+11,ii:ii+5) = Kg10;
        Kg(ii:ii+5,ii:ii+5) = Kg(ii:ii+5,ii:ii+5) + Kg00;
        Kg(ii+6:ii+11,ii+6:ii+11) = Kg(ii+6:ii+11,ii+6:ii+11) + Kg11;
        
        % Global Internal Force Vector for whole structure
        Fg(6*ie-5:6*ie+6) = Fg(6*ie-5:6*ie+6) + fg_glob;

  end