% This code simulates the dynamics of a Cantilever Beam under large
% deformations.

% The beam is modelled using "Corotational Beam Elements"

% The formulation of these elements can be found in the paper and thesis: 
    % Efﬁcient formulation for dynamics of corotational 2D beams
    % Nonlinear dynamics of flexible structures using corotational beam elements

% This code mirrors the example in the paper

% To solve the Equations of Motion, the Hilber Hughes Taylor (HHT) Alpha Method
% is used 
% Details for the HHT implementation can be found in the following
% materials:
    % "Non-linear Finite Element Analysis of Solids and Structures Volume
    % 2: Advanced Topics"
    % "Dynamics of 3D beam elements in a corotational context: A comparitve
    % study of established and new formulations"

% This code is not sophisticated or advanced. Use this code to help
% understand the formulations seen in the papers and books mentioned.

clear; clc; close all;
addpath Core\
tic

%% Beam Properties
% Geometric properites - rectangular cross section
params.b = 0.3;
params.h = 0.25;
params.A = params.b * params.h;
params.I = (params.b * params.h^3)/12;
params.Lt = 10; % Total length
% Material Properties
params.rho = 7850;
params.E = 210e9;
params.nu = 0.3;
params.G = params.E/(2*(1+params.nu));
params.Ks = 5/6;
% Element Properties
params.ne = 5;
params.n = 3 * (params.ne + 1);
params.Le = params.Lt/params.ne; % Define params.Le as the element length

params.Omega = (params.E*params.I)/(params.G*params.A*params.Ks*(params.Le^2));
%% Loading Details
Po = 6e6; % Magnitude
om = 50; % Frequency 

%% q0 here is the initial position of the element nodes etc.
q0 = zeros(params.n,1);
for ie = 1:params.ne+1
    q0(3*ie - 2) = (ie - 1) * params.Le;
end
q0(1:3) = zeros(3,1); % Cantilever Boundary Conditions

%% Initialize global matrices
p = zeros(params.n,1); % External Load Vector P
% Elastic Displacements and Derivatives assumed to be zeros at initial
% state
q = zeros(params.n,1);
qd = zeros(params.n,1);
qdd = zeros(params.n,1);
Dq = zeros(params.n,1);
deltq = zeros(params.n,1);

%% Solver Parameters
% Newton Raphson Parameters
m = 0; maxit = 200;
er = 1e-5; e = 1; 

% HHT Alpha Integration Method
alph = -0.05;
gam = 0.5 * (1 - 2*alph); beta = 0.25 * (1 - alph)^2;
dt = 1e-4;
t = 0; tend = 1;

% Data storage parameters
hout = 1e-3; % Saving data every 1e-3 seconds
nout = floor(tend/hout) + 1;
data.t = zeros(1, nout); 
data.q = zeros(params.n, nout);
data.t(1) = 0;
data.q(:,1) = q;
iout = 1;

%% Start Solver

% The matrices M, Kg, Ck, Kk are needed for the Newton Raphson (NR) method. The appear after
% linearasing the nonlinear dynamic EOM. Look at the Wikipedia page on
% "Linearization".

while t < tend
    m = 0; %  Reset iteration count
    e = 1; % Reset e

    % Storing the values from the previous timestep i.e. N
    q_prev = q; qd_prev = qd; qdd_prev = qdd;
    p_p = zeros(params.n-3,1);
    p_p(end-1) = Po * sin(om * t);

    % Calculating matrices for the previous timestep - N
    [M_p, Ck_p, Kk_p, Kg_p, Fg_p, Fk_p] = matrix_assembly(q0,q_prev,qd_prev,qdd_prev,params);
    M_p = M_p(4:params.n,4:params.n);
    Ck_p = Ck_p(4:params.n,4:params.n);
    Kk_p = Kk_p(4:params.n,4:params.n);
    Fk_p = Fk_p(4:params.n,1);
    Kg_p = Kg_p(4:params.n,4:params.n);
    Fg_p = Fg_p(4:params.n,1);
    
    % Calculating the external force for the next timestep - N + 1
    p = zeros(params.n-3,1); p(end-1) = Po * sin(om * (t + dt));

    % Calculating the predictor Dq
    Df = (1 + alph)*p - Fg_p - Fk_p - alph*p_p ... 
        + M_p/(beta*dt^2) * (dt*qd_prev(4:params.n) + dt^2*0.5*qdd_prev(4:params.n)) ...
         + Ck_p * ((gam/beta)*qd_prev(4:params.n) - (dt*((2*beta - gam)/(2*beta))*qdd_prev(4:params.n)));

    Kbar_p = (1 + alph)*Kg_p + Kk_p + gam/(beta*dt) * Ck_p...
         + 1/(beta*dt^2) * M_p;

    Dq = Kbar_p\Df; % Predictor improves NR convergence speed

    q = q_prev + [zeros(3,1); Dq]; % Initial estimate for displacement

    Dqdd = 1/(beta*dt^2) * ([zeros(3,1); Dq] - dt*qd_prev - dt^2*0.5*qdd_prev);
    Dqd = gam/(beta*dt) * [zeros(3,1); Dq] - (gam/beta)*qd_prev + dt*((2*beta - gam)/(2*beta))*qdd_prev;

    % Initial Estimate for Velocity and Acceleration at time instance N + 1
    qd = qd_prev + Dqd;
    qdd = qdd_prev + Dqdd;

    % Force term based on previous timestep - needed for HHT equilibrium
    % equation
    F_prev = (Fg_p - p_p); 

    while (m < maxit & e > er) % Newton Raphson Loop for Corrector term
        m = m + 1;
        % Assembling matrices for the overall beam structure
        [M, Ck, Kk, Kg, Fg, Fk] = matrix_assembly(q0,q,qd,qdd,params);
        % Adjust matrices for Cantilever Boundary conditions
        Mred = M(4:params.n,4:params.n);
        Ckred = Ck(4:params.n,4:params.n);
        Kkred = Kk(4:params.n,4:params.n);
        Kgred = Kg(4:params.n,4:params.n);
        Fgred = Fg(4:params.n,1);
        Fkred = Fk(4:params.n,1);
   
        % HHT Nonlinear Dynamic Equilibrium Equation
        phi = (1 + alph)*(Fgred - p) - alph*F_prev + Fkred;
    
        if m == 1
            phif = phi;
        end
        % Error term for convergence check
        e = norm(phi)/norm(phif);
        
        % Tangent Stiffness Matrix
        Ktbar = (1 + alph)*Kgred + Kkred + gam/(beta*dt) * Ckred + 1/(beta*dt^2) * Mred;
        
        % Calculation of Corrector term
        deltq = -(Ktbar\phi);
    
        % Displacments etc. for N+1 at the next iteration
        q = q + [zeros(3,1); deltq];
        qd = qd + gam/(beta*dt) * [zeros(3,1); deltq];
        qdd = qdd + 1/(beta*dt^2) * [zeros(3,1); deltq];
        % Values at time instance N+1 are constantly updated until
        % convergence obtained

    end

    t = t + dt; % Move to next timestep N+1

    if m == maxit
        % Convergence checker
        fprintf("No Convergence at t = %f \n",t);
    end
    
    if t >= iout * hout
        % Save Output
        data.t(iout+1) = t;
        data.q(:,iout+1) = q;
        iout = iout + 1;
    end

end

U = data.q(end-2,:); % Horizontal Tip Displacement
V = data.q(end-1,:); % Vertical Tip Displacement

figure; % Tip Displacement Plot
tiledlayout(1,2);
nexttile
plot(data.t,U,'LineWidth',2,'Color','k');
title('Horizontal Tip Displacement'); grid on;
xlabel('Time [s]'); ylabel('U [m]'); xlim([0,tend]);
nexttile
plot(data.t,V,'LineWidth',2,'Color','k'); 
title('Vertical Tip Displacement'); grid on;
xlabel('Time [s]'); ylabel('V [m]'); xlim([0,tend]);

toc

% load("AbaqusT.mat")
% 
% figure;
% plot(AbaqusU2T(:,1),AbaqusU2T(:,2),"LineWidth",2,"Color","k","LineStyle","--");
% hold on; grid on;
% plot(data.t,V,"LineWidth",2,"Color","b");
% xlabel("Time [s]"); ylabel("Vertical Displacement [m]"); xlim([0 0.5]);
% legend("Abaqus - 50 Elements","Corotational - 5 Elements");

