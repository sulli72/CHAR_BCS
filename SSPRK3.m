function [UNP1] = SSPRK3(U,dx,dt,wvsp,ifilt,alphafilt,tv,n)
% 3rd Order TVD Runge-Kutta Explicit Time Marching
% Follows scheme developed and proved by Gottlieb and Shu:
% "TOTAL VARIATION DIMINISHING RUNGE-KUTTA SCHEMES - 1998
% (Journal: Mathematics of Computation)

% Inputs:
% Q == solution at current timestep
% dx == grid spacing
% dt == time step
% wvsp == wavespeed for linear wave eqn

% Outputs:
% QNP1 == solution after RK3 algorithm marches forward in time

%%% General Conservative formulation:
%  UNP1(j) = UN(j) - dt/dx ( F(i+1/2) - F(i-1/2) )

F = PRIM2FLUX(U,wvsp);
for i=1:size(F,2)
  DFDX(:,i) = COMPACTDIFF4(F(:,i),dx);
end

% TREAT DFDX LIKE RHS
DFDX = CHARBC(U,DFDX,dx,wvsp,tv,n);

U1 = U - dt*DFDX;
% u1 = ui - dt*U(:,2);


F = PRIM2FLUX(U1,wvsp);
for i=1:size(F,2)
  DFDX(:,i) = COMPACTDIFF4(F(:,i),dx);
end
% TREAT DFDX LIKE RHS
DFDX = CHARBC(U1,DFDX,dx,wvsp,tv,n);

U2 = 3/4*U + 1/4*(U1 - dt*DFDX);

% u2 = 3/4*ui + 1/4*(u1 - dt*U2(:,2));


F = PRIM2FLUX(U2,wvsp);
for i=1:size(F,2)
  DFDX(:,i) = COMPACTDIFF4(F(:,i),dx);
end
% TREAT DFDX LIKE RHS
DFDX = CHARBC(U2,DFDX,dx,wvsp,tv,n);

UNP1 = 1/3*(U) + 2/3*(U2 - dt*DFDX);

% uout = 1/3*(ui) + 2/3*(u2 - dt*UNP1(:,2));

%%% Filter once after SSP RK3 is done
if ifilt==1
  for i=1:2
    % UNP1(:,i) = PADEFILT6(UNP1(:,i),alphafilt);
    UNP1(:,i) = PADEFILT2(UNP1(:,i),alphafilt);
  end
  % UNP1 = PADEFILT4(UNP1,alphafilt);
end


end