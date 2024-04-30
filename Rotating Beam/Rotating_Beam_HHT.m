% The rotating beam is a classic example in validating large deformation
% beam formulations

% This code applies the Corotational Beam Formulation from the Cantilever
% case to the Rotating Beam Example

% Make sure you are familiar with the formulation in the Cantilever case so
% you can understand the nuances and differences in this Rotating Beam
% case

% The parameters for the rotating beam example can be found in the paper:
    % On the Dynamics of Flexible Beams Under Large Overall Motions-The Plane Case: Part II

% Keep in mind that this code can take a while to run and might have some
% initial convergence problems. The results obtained are validated.

clear; clc; close all;
tic
addpath Core\
%% Beam Properties
% Geometric + Material Properties
params.EA = 2.8e7;
params.EI = 1.4e4;
params.Arho = 1.2;
params.Irho = 6e-4;
params.Lt = 10;
params.nu = 0.3;

params.GA = params.EA/(2*(1+params.nu));
params.Ks = 5/6;
% Element Properties
params.ne = 10;
params.n = 3 * (params.ne + 1);
params.Le = params.Lt/params.ne; % Define params.Le as the element length

params.Omega = params.EI/(params.GA*params.Ks*(params.Le^2));

Ts = 15; ws = 6; % Rotation Quantities

%% q0 here is the initial position of the element nodes etc.
q0 = zeros(params.n,1);
for ie = 1:params.ne+1
    q0(3*ie - 2) = (ie - 1) * params.Le;
end
q0(1:3) = zeros(3,1); % Boundary conditions

%% Initialize global matrices
p = zeros(params.n,1); % External Load Vector P
% Displacement and Derivatives assumed to be zeros
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
alph = -0.01;
gam = 0.5 * (1 - 2*alph); beta = 0.25 * (1 - alph)^2;
dt = 1e-3;
t = 0; tend = 20;

% Data storage parameters
hout = 1e-3; % Saving data every 1e-3 seconds
nout = floor(tend/hout) + 1;
data.t = zeros(1, nout); 
data.q = zeros(params.n, nout);
data.qd = zeros(params.n, nout);
data.qdd = zeros(params.n, nout);
data.t(1) = 0;
data.q(:,1) = q;
data.qd(:,1) = qd;
data.qdd(:,1) = qdd;
iout = 1;

%% Start Solver
while t < tend
    m = 0; %  Reset iteration count
    e = 1; % Reset e

    % Storing the values from the previous timestep i.e. N
    q_prev = q; qd_prev = qd; qdd_prev = qdd; p_p = zeros(params.n-3,1);
    
    [~, ~, ~, ~, Fg_p, ~] = matrix_assembly(q0,q_prev,qd_prev,qdd_prev,params);
    Fg_p = Fg_p(4:params.n,1);

    % Initial Estimate for the velocity and acceleration - predictor seen
    % in Simo and Vu-Quoc Paper
    qdd = -(qd_prev/(dt*beta) + ((0.5 - beta)/beta)*qdd_prev);
    qd = qd_prev + dt*((1 - gam)*qdd_prev + gam*qdd);
        
    % Imposing Rotation Quantities at the next timestep
    q(3) = EndRotation(t+dt,Ts,ws);
    qd(3) = EndAngVel(t+dt,ws,Ts);
    qdd(3) = EndAngAcc(t+dt,ws,Ts);

    p = zeros(params.n-3,1);
          
    % Force term based on previous timestep - needed for HHT equilibrium
    % equation
    F_prev = (Fg_p - p_p); 

    while (m < maxit & e > er) % Newton Raphson Loop for Corrector term
        m = m + 1;
        % Assembling matrices for the overall structure
        [M, Ck, Kk, Kg, Fg, Fk] = matrix_assembly(q0,q,qd,qdd,params);
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
        data.qd(:,iout+1) = qd;
        data.qdd(:,iout+1) = qdd;
        iout = iout + 1;
    end
end

U = data.q(end-2,:); % Horizontal Tip Displacement
V = data.q(end-1,:); % Vertical Tip Displacement

toc

% Calculate Relative Tip Displacements
r = params.Lt;

Xr = zeros(size(data.t));
Yr = zeros(size(data.t));
th = zeros(size(data.t));

for i = 1:size(data.t,2)
    
    ts = data.t(i);
    tht = EndRotation(ts,Ts,ws);
    Xr(i) = r * cos(tht);
    Yr(i) = r * sin(tht);
    th(i) = tht;

end

Ur = r - Xr ;
Vr = Yr;

Vt = V - Vr;
Ut = Ur + U;
Urel = zeros(size(data.t));
Vrel = zeros(size(data.t));

for i = 1:size(data.t,2)
    Vrel(i) = Vt(i) * cos(th(i)) - Ut(i) * sin(th(i));
    Urel(i) = Ut(i) * cos(th(i)) + Vt(i) * sin(th(i));
end
%%
% Relative Tip Displacements
figure(2);
tiledlayout(1,2);
nexttile;
yline(0,'k--');
hold on;
plot(data.t(1:end-1),Urel(1:end-1),'Linewidth',2,'Color','k'); grid on; xlim([0,tend])
xlabel('Time [s]'); ylabel('Tip Displacement U [m]');
hold off;
nexttile;
yline(0,'k--'); hold on;
plot(data.t(1:end-1),Vrel(1:end-1),'LineWidth',2,'Color','k'); grid on; xlim([0,tend])
xlabel('Time [s]'); ylabel('Tip Displacement V [m]'); hold off;