%{
    A pseudo-spectral solver for Dysthe-Lo-Mei equation.
    See Equations (2.6) - (2.9) from Lo & Mei (1985, JFM)
    paper for more details.

    Here we test the code by simulating the NLS ground state evolution
    and deformation under higher order Dysthe terms. This solution is
    exact if the parameter e = 0.

    Copyright (C) 2020 Denys DUTYKH

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

%}

%%% -------------------------------------------------- %%%
%%% The main file for pseudo-spectral Dysthe solver    %%%
%%% -------------------------------------------------- %%%
%%% Author: Denys Dutykh, CNRS -- LAMA, USMB           %%%
%%% E-mail: Denys.Dutykh@univ-smb.fr                   %%%
%%% Web:    http://www.denys-dutykh.com/               %%%
%%% Blog:   http://dutykh.github.io/                   %%%
%%% GitHub: https://github.com/dutykh/                 %%%
%%% -------------------------------------------------- %%%
%%% Distributed under GNU General Public License       %%%
%%% -------------------------------------------------- %%%

clear
close all
format longE
warning('on', 'verbose');

%%% Other source files we use:
addpath('sources/');

%%% Global variables:
global e j g k op K2i

%%% Model parameters:
e = 0.05;    % epsilon (steepness), if e = 0 : NLS equation
g = 1.0;    % gamma
h = 10;     % unperturbed water depth (if > 20: deep water)

%%% Numerical parameters:
N  = 2048;  % number of Fourier modes in the discrete solution
l  = 40.0;  % half-length of the domain in x-space
dx = 2*l/N; % distance between two x-points
x  = (1-N/2:N/2)'*dx; % x-space discretization
k  = [0:N/2-1 0 1-N/2:-1]'*pi/l; % vector of wavenumbers

% Composition of d/dx (i*k) and Hilbert transform (in finite depth, becomes Hilbert if h > 20)
op = -0.5*k.*tanh(2.0*k*h);

% antialising treatment
j  = (N/4+2:N/4*3);  % the frequencies we sacrify
k(j) = 0;

% Define some useful operators:
K2i = -1i*g*g*k.^2;  % symbol of the linear operator

TimeStepLoad;   % we load in memory the coefficients of the embedded RK scheme

fftw('planner', 'exhaustive');

%%% Initial condition specification:
om = 1.0;       % breathing frequency
a  = sqrt(2.0*om);
u0 = a*sech(sqrt(om)/g*x);
v  = fft(u0);

tol = 1e-13;            % tolerance for the time stepping
er1 = tol; er2 = tol;   % local errors for adaptive time stepping
rho1 = 1.0; rho2 = 1.0; % time step multipliers

omega = inline ('1+atan(r-1)', 'r'); % time step limiter function
% (kappa = 1)

%%% Time stepping parameters:
t = 0.0;        % the discrete time variable
Tf = 30.0;      % final simulation time
tloc = 0.0;     % local time between two plots
tdraw = 0.1;    % we plot results each 'tdraw' seconds

dt2 = 1.0e-3;   % initial guess of the time step

figure;
set(gcf, 'pos', [744 600 731 422]);

while (t < Tf)
   er1 = 2*tol;
   while (er1 > tol)
       rho1 = rho2; dt = dt2;
       
       k1 = RHS(v, 0);
       k2 = RHS(v + a2*dt*k1, a2*dt);
       k3 = RHS(v + dt*(b31*k1 + b32*k2), a3*dt);
       k4 = RHS(v + dt*(b41*k1 + b43*k3), a4*dt);
       k5 = RHS(v + dt*(b51*k1 + b53*k3 + b54*k4), a5*dt);
       k6 = RHS(v + dt*(b61*k1 + b64*k4 + b65*k5), a6*dt);
       k7 = RHS(v + dt*(b71*k1 + b74*k4 + b75*k5 + b76*k6), a7*dt);
       k8 = RHS(v + dt*(b81*k1 + b86*k6 + b87*k7), a8*dt);
       k9 = RHS(v + dt*(b91*k1 + b96*k6 + b97*k7 + b98*k8), a9*dt);
       k10 = RHS(v + dt*(b10_1*k1 + b10_6*k6 + b10_7*k7 + b10_8*k8 + b10_9*k9), a10*dt);
       k11 = RHS(v + dt*(b11_1*k1 + b11_7*k7 + b11_8*k8 + b11_9*k9 + b11_10 * k10), a11*dt);
       k12 = RHS(v + dt*(b12_1*k1 + b12_6*k6 + b12_7*k7 + b12_8*k8 + b12_9*k9 + b12_10 * k10 + b12_11 * k11), a12*dt);
       k13 = RHS(v + dt*(b13_1*k1 + b13_6*k6 + b13_7*k7 + b13_8*k8 + b13_9*k9 + b13_10*k10 + b13_11*k11 + b13_12*k12), a13*dt);
       k14 = RHS(v + dt*(b14_1*k1 + b14_6*k6 + b14_7*k7 + b14_8*k8 + b14_9*k9 + b14_10*k10 + b14_11*k11 + b14_12*k12 + b14_13*k13), dt);
       k15 = RHS(v + dt*(b15_1*k1 + b15_6*k6 + b15_7*k7 + b15_8*k8 + b15_9*k9 + b15_10*k10 + b15_11*k11 + b15_12*k12 + b15_13*k13), a15*dt);
       k16 = RHS(v + dt*(b16_1*k1 + b16_6*k6 + b16_7*k7 + b16_8*k8 + b16_9*k9 + b16_10*k10 + b16_11*k11 + b16_12*k12 + b16_13*k13 + b16_15*k15), dt);
       
       er2 = dt*norm(e1*k1 + e8*k8 + e9*k9 + e10*k10 + e11*k11 +...
           e12*k12 + e13*k13 + e14*k14 + e15*k15 + e16*k16, inf);
       rho2 = (tol/er2)^C36*(tol/er1)^C36*(rho1^(-0.25));
       dt2 = max(omega(rho2)*dt, 1e-6); er1 = er2;
  end
  w = v + dt*(c1*k1 + c8*k8 + c9*k9 + c10*k10 + c11*k11 +...
      c12*k12 + c13*k13 + c14*k14);
  v = turn(w, -dt);
  
  if (tloc + dt > tdraw)
    dt = tdraw - tloc;
  end % if ()

  if (t + dt > Tf)
    dt = Tf - t;
  end % if ()

  tloc = tloc + dt;
  t = t + dt;

  % we plot numerical results
  if ((tloc == tdraw) || (t == Tf))

      tloc = 0;
      
      u_hat = v;
      u = ifft(u_hat);
      
      subplot (2,1,1)
      plot (x, abs(u), 'b-', 'LineWidth', 2.0), grid on
      axis([-l l -0.1 1.1*a])
      xlabel('$x$', 'interpreter', 'LaTeX');
      ylabel ('$|A(x,t)|$', 'interpreter', 'LaTeX');
      title (['Absolute value of A(x,t) at t = ',num2str(t,'%5.2f'),...
          '; dt = ',num2str(dt,'%8.7f')]);

      subplot (2,1,2)
      loglog(1:N/2, abs(v(1:N/2,1)/N), 'k.-'), grid on
      axis([1 N/2 1e-20 1e1])
      title ('Fourier half-spectrum of the solution')
      xlabel('Fourier mode number, $n$', 'interpreter', 'LaTeX');
      ylabel('$|\hat{A}(k,t)|$', 'interpreter', 'LaTeX');
      set(gcf, 'color', 'w');
      drawnow
  end
end