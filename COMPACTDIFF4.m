function PHIDERIV = COMPACTDIFF4(PHI,dx)
% 4th Order Compact FD scheme
il=length(PHI);

% Scheme coeffs
ALPHASCHEME=1/4;
a=3/2; %14/9;
b=0; %1/9;

% Tri-diag matrix vectors
ap1=zeros(il-1,1); % super diag
ac0=ones(il,1);   % main diag -- hard fixed for dirichlet BCs
am1=zeros(il-1,1); % sub diag
kv=zeros(il,1);    % RHS vector

% Build tri-diag vectors
for ii=2:il-1
  ap1(ii)=ALPHASCHEME;
  ac0(ii)=1;
  am1(ii-1)=ALPHASCHEME;
  kv(ii)=a*(PHI(ii+1)-PHI(ii-1))/(2*dx);
end

% C4 at first boundary point
ii=1;
ap1(1)=3;
a1=-17/6;
b1=3/2;
c1=3/2;
d1=-1/6;
kv(ii) = ( a1*PHI(ii) + b1*PHI(ii+1) + c1*PHI(ii+2) + d1*PHI(ii+3)   )/(1*dx);

% C4 at last boundary point
ii=il;
am1(ii-1)=3;
a1=17/6;
b1=-3/2;
c1=-3/2;
d1=1/6;
kv(ii) = ( a1*PHI(ii) + b1*PHI(ii-1) + c1*PHI(ii-2) + d1*PHI(ii-3)   )/(1*dx);

% Compute derivatives using Thomas algorithm
PHIDERIV = THOMAS(ac0,ap1,am1,kv);
end