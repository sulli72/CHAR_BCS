function [UNP1] = SSPRK3(U,dx,dt,wvsp)
% 3rd Order TVD Runge-Kutta Explicit Time Marching
% Follows scheme developed and proved by Gottlieb and Shu:
% "TOTAL VARIATION DIMINISHING RUNGE-KUTTA SCHEMES - 1998
% (Journal: Mathematics of Computation)

% Inputs:
% U == solution at current timestep
% dx == grid spacing
% dt == time step
% wvsp == wavespeed for linear wave eqn

% Outputs:
% UNP1 == solution after RK3 algorithm marches forward in time

%%% General Conservative formulation:
%  UNP1(j) = UN(j) - dt/dx ( F(i+1/2) - F(i-1/2) )

F = PRIM2FLUX(U,wvsp);                 % <-- create flux functions from solution variables
for i=1:size(F,2)                      % <-- loop over each variable
  DFDX(:,i) = COMPACTDIFF4(F(:,i),dx); % <-- get d/dx f(u) 
end
DFDX = CHARBC(U,DFDX,dx,wvsp);         % <-- Modify fluxes at boundaries
U1 = U - dt*DFDX;                      % <-- March to next RK step


% Second RK3 step
F = PRIM2FLUX(U1,wvsp);
for i=1:size(F,2)
  DFDX(:,i) = COMPACTDIFF4(F(:,i),dx);
end
DFDX = CHARBC(U1,DFDX,dx,wvsp);
U2 = 3/4*U + 1/4*(U1 - dt*DFDX);

% Third RK3 step
F = PRIM2FLUX(U2,wvsp);
for i=1:size(F,2)
  DFDX(:,i) = COMPACTDIFF4(F(:,i),dx);
end
DFDX = CHARBC(U2,DFDX,dx,wvsp);
UNP1 = 1/3*(U) + 2/3*(U2 - dt*DFDX);

end