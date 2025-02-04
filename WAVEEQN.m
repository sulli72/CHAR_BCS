%%% Header
clearvars; close all; clc;
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
wavepos=[50 50 800 700];
nfg=0;
%%% End of Header

%% SECOND ORDER WAVE EQUATION
% u_tt - c^2 u_xx = 0
%
% rewrite as sytem of first order eqns
% [w v]'_t + f( [w v]')_x = 0
% 
% w = u_x , v = u_t

%%% ----------- NOTE ----------- %%%
%%% Boundary conditions are set  %%%
%%% by altering L1,L2 at end pts %%%
%%%       in CHARBC routine      %%%   
%%% ---------------------------- %%%

%% Mesh and timestep set up 
% Mesh
il=401;
x0=0; xf=1;
norm=2*pi/(xf-x0); % <-- helps normalize wavenumbers
dx=(xf-x0)/(il-1);
xv=x0:dx:xf;

% Wave Speed
c=1;

% Timestep
t0=0; tf=2;
CFL=.7;
dt=CFL*dx/c;
tv=t0:dt:tf;
nt=length(tv);

%% Initial Spatial Condition

% icpick=1; % <-- Single harmonic standing wave
% icpick=2; % <-- Smaller single gaussian pulse
% icpick=3; % <-- Single tall gaussian pulse
icpick=4; % <-- undisturb, u0(x) = 0
u0=ICMAKER(icpick,xv);

% u_t IC
dudt0 = zeros(il,1); % <-- no initial impulse

%% FORCING FOR UNDISTURBED IC
if icpick==4

  % Spatial distribution of forcing function
  f0x = cos(5*norm*xv).^2;
  xst = interp1(xv,xv,.45,'nearest');
  xfn = interp1(xv,xv,.55,'nearest');
  ist=find(xv==xst);
  ifn=find(xv==xfn);
  f0x(1:ist-1)=0;
  f0x(ifn+1:end)=0;
  f0x=0.5*f0x;
  tdyn=zeros(1,nt);
  % tdyn=sin(omega*2*pi*tv); % monochromatic forcing

  % Apply harmonic forcing over finite time interval
  omegav=[10];
  tst = interp1(tv,tv,0,'nearest');
  tfn = interp1(tv,tv,.3,'nearest');
  nst=find(tv==tst);
  nfn=find(tv==tfn);
  for w=1:length(omegav)
    omega=omegav(w);
    tdyn(nst:nfn)=tdyn(nst:nfn)+sin(omega*2*pi*tv(nst:nfn));
  end
  fxt=f0x'*tdyn; % outer product to get forcing function over space and time
else
  fxt=zeros(il,nt); % <-- no forcing for other IC cases 1,2,3
end

%% Get variables w0,v0 from IC of wave equation
w0=COMPACTDIFF4(u0,dx);
v0=dudt0;


%% Allocate variables and assign IC
w=zeros(il,nt); v=zeros(il,nt); u=zeros(il,nt);
w(:,1)=w0; v(:,1)=v0; u(:,1)=u0;

%% Numerically solve equation system
for n=1:nt-1

  % Instantaneous fields
  wmy=w(:,n);  vmy=v(:,n);  

  % Add forcing to time derivative variable 
  vmy = vmy + fxt(:,n);

  % Create solution vector
  Q = [wmy,vmy];

  % Integrate using RK3 algorithm
  [QNP1] = SSPRK3(Q,dx,dt,c);

  % Fill solution at next timestep
  w(:,n+1)=QNP1(:,1);  v(:,n+1)=QNP1(:,2);
end

% Integrate v(x,t) = u_t to get actual solution
for i=1:il
  vi=v(i,:);
  u(i,:)=u0(i) + cumtrapz(tv,vi);
end

%% CHARACTERISTIC PLOT
[X,T] = meshgrid(xv,tv);
nfg=nfg+1;
figure(nfg)
imagesc(xv,tv,u')
colormap jet;
set(gca,'YDir','normal')
title('$u(x,t)$')
xlabel('$x$')
ylabel('$t$')
fontsize(gcf,16,'points')
ax=gca;
borderpos = tightPosition(ax);
annotation("rectangle",borderpos,Color="black",LineWidth=1.5)
%% ANIMATE FIELD
nfg=nfg+1;
f=figure(nfg);
for n=1:nt
  umy=u(:,n);
  

  plot(xv,umy,'k-','LineWidth',2)
  xlabel('$x$')
  ylabel('$u(x,t)$')
  ylim([-1.2 1.2]);
  fontsize(gcf,16,'points')
  xlim([x0 xf]);
  grid off
  ax=gca;
  borderpos = tightPosition(ax);
  annotation("rectangle",borderpos,Color="black",LineWidth=1.5)
  drawnow

end
