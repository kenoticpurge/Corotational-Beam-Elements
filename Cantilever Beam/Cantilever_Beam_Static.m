% This code simulates the statics of a Cantilever Beam under large
% deformations.

% The beam is modelled using "Corotational Beam Elements"

% The formulation of these elements can be found in the paper: 
    % Efﬁcient formulation for dynamics of corotational 2D beams

% This code adapts the example in the paper for static solving.

% The static equilbrium equations are solved using the Newton Raphson
% method

% This code is not sophisticated or advanced. Use this code to help
% understand the formulations seen in the papers and books mentioned.

clear; clc; close all;
addpath Core\
addpath Data\

%% Beam Properties
tic
% Geometric properites
params.b = 0.25;
params.h = 0.3;
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
params.ne = 50; % Number of Elements
params.n = 3 * (params.ne + 1);
params.Le = params.Lt/params.ne; % Define params.Le as the element length

params.Omega = (params.E*params.I)/(params.G*params.A*params.Ks*(params.Le^2));

%% Loading - Normalized Loading
Lf = 1; % Loading Factor
P = (Lf * params.E * params.I)/(params.Lt^2);
p = zeros(params.n,1);
p(end-1) = P; % Load in the Positive Y direction

%% q0 here is the initial position of the element nodes etc.
% Need for the length calculations
q0 = zeros(params.n,1);
for ie = 1:params.ne+1
    q0(3*ie - 2) = (ie - 1) * params.Le;
end
q0(1:3) = zeros(3,1); % Boundary conditions

%% Initialize global matrices
% Elastic Displacement and Derivatives assumed to be zeros
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
        % Adjust Matrices and Vector for Cantilever Boundary Conditions
        Kgred = Kg(4:params.n,4:params.n);
        Fgred = Fg(4:params.n,1);
        pred = p(4:params.n,1);

        % Nonlinear Equilibrium Equation
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



