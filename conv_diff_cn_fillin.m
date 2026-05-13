% Course: CFD Lab
% TU-Muenchen, SS 2020
%
% This is the Matlab script for the unsteady 1D convection diffusion equation
%
% You must fill in the missing parts by yourself!
% Missing parts are marked by ???
%
% Tianshi Sun, tianshi.sun@tum.de
% 
%
%--------------------------------------------------------------------------
%
% time integration of
%
%   d phi            d phi             d^2 phi
%  ------- =  -U0 * ------- + Gamma * ---------
%    dt               dx                dx^2 
%
% 0 <= x <= 2pi
%
% periodic boundary condition
% phi(0) = phi(2pi)
%
% initial condition
% t = 0  ==>  phi = sin(x)
%
% Central difference scheme (CDS) for spatial discretization
% Crank Nicolson scheme for time advancement
%
%--------------------------------------------------------------------------
%
% Clear all variables and plots.
format long;
clear;
hold off;

% Set convection velocity
U0 = 1.0;

% Set diffusion coefficient
Gamma = 1.0;

% Discrete spacing in space
xend   = 2.0 * pi;
points = 40; 
dx     = xend / ( points - 1 );
% Grid with x locations:
x = 0.0 : dx : xend;

% Discrete spacing in time
% tstep = number of discrete timesteps
tsteps = 1000;
dt     = 0.1;
tend   = dt * tsteps;

% Initialise coefficient matrix A, constant vector b
% and solution vector phi
A   = zeros(points,points);
b   = zeros(points,1);
phi = zeros(points,1);
phi_a = zeros(points,1);

% Initialise the solution (initial condition)
% Loop over grid points in space:
for j = 1 : points
   phi(j) = sin(x(j));
end

% Check initial field:
plot(x, phi, 'r');
hold on;
pause(3);

% Compute coefficients of matrix A
a_w = - U0 * dt / (4.0 * dx) ...
      - Gamma * dt / (2.0 * dx^2);
a_p = 1.0 ...
      + Gamma * dt / (dx^2);
a_e =   U0 * dt / (4.0 * dx) ...
      - Gamma * dt / (2.0 * dx^2);

% Crank Nicolson:
%----------------
%
% Loop over timesteps:
for i = 1 : tsteps

  A(:,:) = 0.0;
  b(:)   = 0.0;

  % Periodic boundary conditions at x=0:
  A(1,points-1) = a_w;
  A(1,1)        = a_p;
  A(1,2)        = a_e;
  b(1) = ...
      ( U0 * dt / (4.0 * dx) + Gamma * dt / (2.0 * dx^2) ) * phi(points-1) ...
    + ( 1.0 - Gamma * dt / (dx^2) ) * phi(1) ...
    + ( - U0 * dt / (4.0 * dx) + Gamma * dt / (2.0 * dx^2) ) * phi(2);

  % Loop over grid points in space:
  for j = 2 : points - 1

    A(j,j-1) = a_w;
    A(j,j)   = a_p;
    A(j,j+1) = a_e;
    b(j) = ...
        ( U0 * dt / (4.0 * dx) + Gamma * dt / (2.0 * dx^2) ) * phi(j-1) ...
      + ( 1.0 - Gamma * dt / (dx^2) ) * phi(j) ...
      + ( - U0 * dt / (4.0 * dx) + Gamma * dt / (2.0 * dx^2) ) * phi(j+1);

  end

  % Periodic boundary conditions at x=2*pi:
  A(points,points-1) = a_w;
  A(points,points)   = a_p;
  A(points,2)        = a_e;
  b(points) = ...
      ( U0 * dt / (4.0 * dx) + Gamma * dt / (2.0 * dx^2) ) * phi(points-1) ...
    + ( 1.0 - Gamma * dt / (dx^2) ) * phi(points) ...
    + ( - U0 * dt / (4.0 * dx) + Gamma * dt / (2.0 * dx^2) ) * phi(2);

  % Solve the linear system of equations
  phi = A\b;

  % Analytical solution
  t = i * dt;
  for j = 1 : points
    phi_a(j) = exp(-Gamma * t) * sin(x(j) - U0 * t);
    % hint: t(i) = i * dt
  end

  % Plot transported wave for each timestep
  plot(x, phi, 'r', x, phi_a, 'g');
  hold off;
  pause(0.003);

end
