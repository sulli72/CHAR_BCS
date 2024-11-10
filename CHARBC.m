function DFDXMOD = CHARBC(U,DFDX,dx,c,tv,n)

% Unpack variables
w=U(:,1);
v=U(:,2);
il=length(w);
DFDXMOD=DFDX;

% i=1 point -- use proper upwinding for derivatives
% wx = (-3*w(1)+4*w(2)-1*w(3) )/(2*dx);
% vx = (-3*v(1)+4*v(2)-1*v(3) )/(2*dx);
% 
% wx = (2*w(1)-9*w(2)+18*w(3)-11*w(4))/(6*dx);
% vx = (2*w(1)-9*w(2)+18*w(3)-11*w(4))/(6*dx);

wx = -(-25*w(1)+48*w(2)-36*w(3)+16*w(4)-3*w(5))/(12*dx);
vx = -(-25*v(1)+48*v(2)-36*v(3)+16*v(4)-3*v(5))/(12*dx);

L1 = -c^2*wx+c*vx; % <-- +c wave
L2 =  c^2*wx+c*vx; % <-- -c wave

%%% at i=1, +c coming into domain
L1 = -L2;  % <-- Neumann BC   (value based on IC u(x,0)
L1 = L2;   % <-- Dirichlet BC (value based on IC u(x,0)
% L1 = 0;      % <-- Non-reflective BC


% Map to fluxes for each variable
flx1 = L1+L2;
flx2 = c*(-L1+L2);


DFDXMOD(1,:) = [flx1,flx2]; % <-- replace RHS at i=1




% i=il point -- use proper upwinding for derivatives
wx = (w(il-2)-4*w(il-1)+3*w(il) )/(2*dx);
vx = (v(il-2)-4*v(il-1)+3*v(il) )/(2*dx);

wx = (-2*w(il-3)+9*w(il-2)-18*w(il-1)+11*w(il))/(6*dx);
vx = (-2*w(il-3)+9*w(il-2)-18*w(il-1)+11*w(il))/(6*dx);

wx = (-25*w(il)+48*w(il-1)-36*w(il-2)+16*w(il-3)-3*w(il-4))/(12*dx);
vx = (-25*v(il)+48*v(il-1)-36*v(il-2)+16*v(il-3)-3*v(il-4))/(12*dx);

L1 = -c^2*wx+c*vx; % <-- +c wave
L2 =  c^2*wx+c*vx; % <-- -c wave

%%% at i=il, -c coming into domain
L2 = -L1;  % <-- Neumann BC   (value based on IC u(x,0)
L2 = L1;   % <-- Dirichlet BC (value based on IC u(x,0)
L2 = 0;      % <-- Non-reflective BC
% L2 = 0.5*(1 - cos(10*pi*tv(n)));

% Map to fluxes for each variable
flx1 = L1+L2;
flx2 = c*(-L1+L2);


DFDXMOD(il,:) = [flx1,flx2]; % <-- replace RHS at i=il


db=1;

end