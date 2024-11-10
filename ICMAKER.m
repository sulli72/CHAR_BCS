function [u0] = ICMAKER(pick,xv)

x0=min(xv);
xf=max(xv);
il=length(xv);
norm=2*pi/(xf-x0);

if pick==1
% SIN WAVE IC
alpha=1; %1/2;
u0 = sin(alpha*norm*xv);


elseif pick==2
% ISOLATED BUMP IC
u0 = cos(5*norm*xv).^2;
xst = interp1(xv,xv,.45,'nearest');
xfn = interp1(xv,xv,.55,'nearest');
ist=find(xv==xst);
ifn=find(xv==xfn);
u0(1:ist-1)=0;
u0(ifn+1:end)=0;
u0 = 0.5*u0;
uadd = 0; % linspace(0,0.5,il);
u0 = u0+uadd;
% u0(:)=0; % for forcing case


elseif pick==3
% GAUSSIAN IC
sig = .01; % std dev
mu = (x0+xf)/2; % mean (in x)
u0 = 1/(sig*sqrt(2*pi)) * exp(-0.5*( (xv-mu)/sig ).^2);
u0 = u0/max(u0);



elseif pick==4
% NO INITIAL PROFILE
u0(1:il)=0; % for forcing case


end