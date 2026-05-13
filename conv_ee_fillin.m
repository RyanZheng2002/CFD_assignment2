% Course: CFD Lab
% TU Muenchen, Summer term 2020
%
% This is the Matlab script for the unsteady 1D convection equation
%
% You must fill in the missing parts by yourself!
% Missing parts are marked by ???
%
% Author: Tianshi Sun, tianshi.sun@tum.de
% 
%
%--------------------------------------------------------------------------
%
% time integration of
%
%   d phi            d phi
%  ------- =  -U0 * -------
%    dt               dx
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
% Explicit Euler scheme for time advancement
%
%--------------------------------------------------------------------------
%
% Clear all variables and plots.
format long;
clear;
hold off;

% Set convection velocity
U0 = 1.0;

% Discrete spacing in space
xend   = 2.0 * pi;
points = 80; 
dx     = xend / ( points - 1 );
% Grid with x locations:
x = 0.0 : dx : xend;

% Discrete spacing in time
% tstep = number of discrete timesteps
tsteps = 1000;
dt     = 0.1;
tend   = dt * tsteps;

% Initialise the solution (initial condition)
phi    = zeros(points,1);
phinew = zeros(points,1);
phi_a  = zeros(points,1);

% Loop over grid points in space:
for j = 1 : points
   phi(j) = sin(x(j));
end

% Check initial field:
plot(x, phi, 'r');
hold on;
pause(3);

% Explicit Euler:
%----------------
%
% phinew is phi at new time level
% phinew must be written back to phi for new timestep
%
% Loop over timesteps:
for i = 1 : tsteps

  % Periodic boundary conditions at x=0:
  phinew(1) = phi(1) ...
            - U0 * dt / (2.0 * dx) * (phi(2) - phi(points-1));

  % Loop over grid points in space:
  for j = 2 : points - 1

    phinew(j) = phi(j) ...
              - U0 * dt / (2.0 * dx) * (phi(j+1) - phi(j-1));

  end

  % Periodic boundary conditions at x=2*pi:
  phinew(points) = phi(points) ...
                 - U0 * dt / (2.0 * dx) * (phi(2) - phi(points-1));

  % Write new field back to old field:
  phi = phinew;

  % Analytical solution
  t = i * dt;
  for j = 1 : points
    phi_a(j) = sin(x(j) - U0 * t);
    % hint: t(i) = i * dt
  end

  % Plot transported wave for each timestep
  plot(x, phi, 'r', x, phi_a, 'g');
  hold off;
  pause(0.003);

end
