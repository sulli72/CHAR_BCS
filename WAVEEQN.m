clearvars; close all; clc;
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
addpath 'C:\Users\jacks\Documents\MATLAB\Add-Ons\Collections\linspecer';
addpath 'C:\Users\jacks\Documents\MATLAB\Add-Ons\Collections\cmap\cmap-master\cmap-master\';
m1pos=[50 50 700 700];
m2pos=[-2500 -900 1800 1200];
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

% FORCING
% famp=.01;
% omega=10;
% f0x=zeros(il,1);
% f0x(ceil(il/2))=famp; % add midpoint forcing
% % f0x(1)=famp; % add endpoint forcing
% tdyn=zeros(1,nt);
% % tdyn=sin(omega*2*pi*tv); % monochromatic forcing
% 
% 
% omegav=[10];
% for nf=1:length(omegav)
%   omega=omegav(nf);
%   tdyn=tdyn+sin(omega*2*pi*tv);
% end
% fxt=f0x*tdyn; % outer product to get forcing function over space and time


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
tf=2; t0=0;
CFL=.7;
dt=CFL*dx/c;
tv=t0:dt:tf;
nt=length(tv);

%% Initial Spatial Condition
icpick=4;
u0=ICMAKER(icpick,xv);

% u_t IC
dudt0 = zeros(il,1);
% dudt0(:)=xv;

% IC CHECK
% figure()
% plot(xv,u0,'k-','LineWidth',2)
% xlim([x0 xf]);
% ylim([-.2 1.2]);


% FORCING
famp=1;
omega=2;
f0x=zeros(il,1);
f0x(ceil(il/2))=famp; % add midpoint forcing
% f0x(1)=famp; % add endpoint forcing

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


omegav=[10];
tst = interp1(tv,tv,0,'nearest');
tfn = interp1(tv,tv,.3,'nearest');
nst=find(tv==tst);
nfn=find(tv==tfn);
for w=1:length(omegav)
  omega=omegav(w);
  tdyn(nst:nfn)=tdyn(nst:nfn)+sin(omega*2*pi*tv(nst:nfn));
  % tdyn(nst:nfn)=tdyn(nst:nfn)+2*omega*pi*cos(omega*2*pi*tv(nst:nfn)); % <-- v(t) forcing
end
fxt=f0x'*tdyn; % outer product to get forcing function over space and time


%% Initial Temporal Condition
ut0(1:il)=0; % <-- zero impulse IC (non-moving medium?)

%% Get variables w,v from IC
w0=COMPACTDIFF4(u0,dx);
v0=dudt0;

% figure()
% plot(xv,w0,'k-','LineWidth',2)
% xlim([x0 xf]);
% ylim([-.2 1.2]);

%% BUILD SOLUTION
w=zeros(il,nt);
v=zeros(il,nt);
u=zeros(il,nt);
w(:,1)=w0;
v(:,1)=v0;
u(:,1)=u0;
for n=1:nt-1
  wmy=w(:,n); 
  vmy=v(:,n);
  umy=u(:,n);

 % Add forcing to prim variable 
 vmy = vmy + fxt(:,n);

  Q = [wmy,vmy];

  ifilt=0;
  alphafilt=.3;

  [QNP1] = SSPRK3(Q,dx,dt,c,ifilt,alphafilt,tv,n);



  w(:,n+1)=QNP1(:,1);
  v(:,n+1)=QNP1(:,2);


end

% Integrate v(x,t) = u_t to get actual solution
for i=1:il
  vi=v(i,:);
  u(i,:)=u0(i) + cumtrapz(tv,vi);
end

% CHARACTERISTIC PLOT
[X,T] = meshgrid(xv,tv);
nfg=nfg+1;
figure(nfg)
imagesc(tv,xv,u')
colormap turbo;
set(gca,'YDir','normal')




%% ANIMATE FIELD
nfg=nfg+1;
f=figure(nfg);
nfnw=ceil(nt/2);
nfnw=nt;
for n=1:nfnw
  umy=u(:,n);
  

  plot(xv,umy,'k-','LineWidth',2)
  xlabel('$x$')
  ylabel('$u(x,t)$')
  ylim([-1 1]);
  fontsize(gcf,16,'points')
  % FORMATFIG(f,wavepos,0,0,0)


%   yp=.7;
%   ylim([-yp yp]);
  xlim([x0 xf]);
  grid off
  ax=gca;
borderpos = tightPosition(ax);
annotation("rectangle",borderpos,Color="black",LineWidth=1.5)

  drawnow
  frame = getframe(f);
  im{n} = frame2im(frame);

  % fontsize(gcf,scale=1.3)

  % exportgraphics(gcf,'test.gif','Append',true);
  fprintf(1,'DID FRAME WRITE %i/%i\n',n,nfnw)
end


%% GIF WRITING
filename = "dirichlet_non_ref_combo.gif"; % Specify the output file name
dtime = 0.025; % <-- delay time in seconds

for idx = 1:nfnw
    [A,map] = rgb2ind(im{idx},256);
    if idx == 1
        imwrite(A,map,filename,"gif",LoopCount=Inf, ...
                DelayTime=dtime)
    else
        imwrite(A,map,filename,"gif",WriteMode="append", ...
                DelayTime=dtime)
    end
      fprintf(1,'DID GIF FRAME %i/%i\n',idx,nfnw)
end



%%

% L2 NORM
% E0=u0*u0';
% for n=1:nt
%   Et(n)=u(:,n)'*u(:,n);
% end
% Et=Et/E0;
% 
% figure(6)
% plot(tv,Et,'k-','LineWidth',2)
% % ylim([-1 1]);
% %   yp=.7;
% %   ylim([-yp yp]);
% xlim([t0 tf]);
% xlabel('$t$')
% ylabel('$E(t)$')
% grid off


% Basis vectors for Ek spectrum
% kmax=10;
% for k=1:kmax
%   PHI(:,k) = sin(norm*k*xv);
% end
% 
% for n=1:nt
%   umy=u(:,n);
%   for k=1:kmax
%     EK(k,n) = umy'*PHI(:,k);
%   end
% 
% end

% kv=1:1:kmax;
% EK0=EK(:,1);
% figure(7)
% plot(kv,EK(:,1:10:100)./EK0,'o','LineWidth',2)
% % ylim([-1 1]);
% %   yp=.7;
% %   ylim([-yp yp]);
% % xlim([t0 tf]);
% xlabel('$\hat{k}$')
% ylabel('$E(k)$')
% grid off