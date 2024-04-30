% We can take the Corotational Beam Formulation seen in the Cantilever Case
% and apply it to other scenarios

% Make sure you understand the Cantilever Example before looking at this
% one

% In this code, we consider the deformation of very flexible pendulum under
% its own weight

% The parameters for this example are found in the paper:
    % Development of simple models for the elastic forces 
    % in the absolute nodal co-ordinate formulation

% The effect of gravity is considered using the method in:
  % "Consistent co-rotational framework for Euler-Bernoulli 
  % and Timoshenko beam-column elements under distributed member loads"

clear; clc; close all;
tic
addpath Core\

%% Beam Properties
% Geometric properites
params.A = 0.0018;
params.I = 1.215e-8;
params.Lt = 1.2; % Total length
% Material Properties
params.rho = 5540;
params.E = 0.7e6;
params.nu = 0.3;
params.G = params.E/(2*(1+params.nu));
params.Ks = 5/6;
% Element Properties
params.ne = 10;
params.n = 3 * (params.ne + 1);
params.Le = params.Lt/params.ne; % Define params.Le as the element length

params.Omega = (params.E*params.I)/(params.G*params.A*params.Ks*(params.Le^2));

% Gravity Terms
params.g = -9.81;
params.p = params.g * params.rho * params.A;

%% q0 here is the initial position of the element nodes etc.
q0 = zeros(params.n,1);
for ie = 1:params.ne+1
    q0(3*ie - 2) = (ie - 1) * params.Le;
end
q0(1:3) = zeros(3,1); % Boundary conditions

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
er = 1e-6; e = 1; 

% HHT Alpha Integration Method
alph = 0.0;
gam = 0.5 * (1 - 2*alph); beta = 0.25 * (1 - alph)^2;
dt = 1e-4;
t = 0; tend = 1.1;

% Data storage parameters
hout = 1e-3; % Saving data every 1e-3 seconds
nout = floor(tend/hout) + 1;
data.t = zeros(1, nout); 
data.q = zeros(params.n, nout);
data.t(1) = 0;
data.q(:,1) = q;
iout = 1;

%% Start Solver

while t < tend
    m = 0; %  Reset iteration count
    e = 1; % Reset e

    % Storing the values from the previous timestep i.e. N
    q_prev = q; qd_prev = qd; qdd_prev = qdd;
    p_p = zeros(params.n-2,1);
    
    % Calculating matrices for the previous timestep - N
    [M_p, Ck_p, Kk_p, Kg_p, Fg_p, Fk_p] = matrix_assembly(q0,q_prev,qd_prev,qdd_prev,params);
    
    M_p = M_p(3:params.n,3:params.n);
    Ck_p = Ck_p(3:params.n,3:params.n);
    Kk_p = Kk_p(3:params.n,3:params.n);
    Fk_p = Fk_p(3:params.n,1);
    Kg_p = Kg_p(3:params.n,3:params.n);
    Fg_p = Fg_p(3:params.n,1);
    
    % Calculating the external force for the next timestep - N + 1
    p = zeros(params.n-2,1);

    % Calculating the predictor Dq
    Df = (1 + alph)*p - Fg_p - Fk_p - alph*p_p ... 
        + M_p/(beta*dt^2) * (dt*qd_prev(3:params.n) + dt^2*0.5*qdd_prev(3:params.n)) ...
         + Ck_p * ((gam/beta)*qd_prev(3:params.n) - (dt*((2*beta - gam)/(2*beta))*qdd_prev(3:params.n)));

    Kbar_p = (1 + alph)*Kg_p + Kk_p + gam/(beta*dt) * Ck_p...
         + 1/(beta*dt^2) * M_p;

    Dq = Kbar_p\Df; % Predictor improves NR convergence speed

    q = q_prev + [zeros(2,1); Dq]; % Initial estimate for displacement

    Dqdd = 1/(beta*dt^2) * ([zeros(2,1); Dq] - dt*qd_prev - dt^2*0.5*qdd_prev);
    Dqd = gam/(beta*dt) * [zeros(2,1); Dq] - (gam/beta)*qd_prev + dt*((2*beta - gam)/(2*beta))*qdd_prev;

    % Initial Estimate for Velocity and Acceleration at time instance N + 1
    qd = qd_prev + Dqd;
    qdd = qdd_prev + Dqdd;

    % Force term based on previous timestep - needed for HHT equilibrium
    % equation
    F_prev = (Fg_p - p_p); 

    while (m < maxit & e > er) % Newton Raphson Loop for Corrector term
        m = m + 1;
        % Assembling matrices for the overall structure
        [M, Ck, Kk, Kg, Fg, Fk] = matrix_assembly(q0,q,qd,qdd,params);
        % Adjust matrices for Pendulu, Boundary Conditions
        Mred = M(3:params.n,3:params.n);
        Ckred = Ck(3:params.n,3:params.n);
        Kkred = Kk(3:params.n,3:params.n);
        Kgred = Kg(3:params.n,3:params.n);
        Fgred = Fg(3:params.n,1);
        Fkred = Fk(3:params.n,1);
   
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
        q = q + [zeros(2,1); deltq];
        qd = qd + gam/(beta*dt) * [zeros(2,1); deltq];
        qdd = qdd + 1/(beta*dt^2) * [zeros(2,1); deltq];
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

%%
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

