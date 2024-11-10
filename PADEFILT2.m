function PHIFIL = PADEFILT2(PHI,ALPHAFIL)
% 2nd Order Pade Filter for Stabilizing Compact FD scheme
il=length(PHI);

% Filter coefficients
af0=1/2 + ALPHAFIL;
af1=1/2 + ALPHAFIL;

afv=[af0 af1];

% Tri-diag matrix vectors
ap1=zeros(il-1,1); % super diag
ac0=ones(il,1);   % main diag -- hard fix for Dirichlet BCs
am1=zeros(il-1,1); % sub diag
kv=zeros(il,1);    % RHS vector

for ii=2:il-1

  ap1(ii)=ALPHAFIL;
  ac0(ii)=1;
  am1(ii-1)=ALPHAFIL;
  kv(ii)=0;

  for aa=1:length(afv)
    amap=aa-1;
    kv(ii) = kv(ii) + afv(aa)/2*(PHI(ii-amap) + PHI(ii+amap));
  end


end

PHIFIL = THOMAS(ac0,ap1,am1,kv);


end