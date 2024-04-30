% We can take the Corotational Beam Formulation seen in the Cantilever Case
% and apply it to other scenarios

% Make sure you understand the Cantilever Examples first, before looking at
% this code

% The effect of gravity is considered using the method in:
  % "Consistent co-rotational framework for Euler-Bernoulli 
  % and Timoshenko beam-column elements under distributed member loads"

clear; clc; close all;
addpath Core\
tic
%% Beam Properties
% Geometric properites
params.b = 0.5;
params.h = 0.25;
params.A = params.b * params.h;
params.I = (params.b * params.h^3)/12;
params.Lt = 10; %total length
% Material Properties
params.rho = 7850;
params.E = 210e9;
params.nu = 0.3;
params.G = params.E/(2*(1+params.nu));
params.Ks = 5/6;
% Element Properties
params.ne = 50;
params.n = 3 * (params.ne + 1);
params.Le = params.Lt/params.ne; % Define params.Le as the element length

params.Omega = (params.E*params.I)/(params.G*params.A*params.Ks*(params.Le^2));

% Gravity Terms
params.g = -9.81;
params.p = params.g * params.rho * params.A;

%% Loading - Normalized Loading
Lf = 1; %Loading Factor
P = (Lf * params.E * params.I)/(params.Lt^2);
p = zeros(params.n,1); %External Force Vector
p(end-1) = -P; 

%% q0 here is the initial position of the element nodes etc.
% Need for the length calculations
q0 = zeros(params.n,1);
for ie = 1:params.ne+1
    q0(3*ie - 2) = (ie - 1) * params.Le;
end
q0(1:3) = zeros(3,1); % Boundary conditions

%% Initialize global matrices
% Elastic Displacements and Derivatives assumed to be zeros at initial
% state
q = zeros(params.n,1);
qd = zeros(params.n,1);
qdd = zeros(params.n,1);
Dq = zeros(params.n,1);
deltq = zeros(params.n,1);

%% Solver Parameters
% Newton Raphson Parameters
m = 0; maxit = 300;
er = 1e-6; e = 1; 

%% Start Solver

while (m < maxit & e > er) % Newton Raphson Loop for Corrector term
        m = m + 1;
        % Assembling matrices for the overall structure
        [~, ~, ~, Kg, Fg, ~] = matrix_assembly(q0,q,qd,qdd,params);
        % Adjust for Cantilever Boundary Conditions
        Kgred = Kg(4:params.n,4:params.n);
        Fgred = Fg(4:params.n,1);
        pred = p(4:params.n,1);

        % Nonlinear Static Equilibrium Equation
        phi = (Fgred - pred);
    
        if m == 1
            phif = phi;
        end
        % Error term for convergence check
        e = norm(phi)/norm(phif);
        
        % Calculation of Corrector term
        deltq = -(Kgred\phi);
    
        % Displacments etc. for N+1 at the next iteration
        q = q + [zeros(3,1); deltq];
      
end

if m == maxit
        % Convergence checker
        fprintf("No Convergence for the current load case \n");
end

% Horizontal and Vertical Tip Displacements
U = q(end-2); V = q(end-1);



