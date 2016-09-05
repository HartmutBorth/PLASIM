%-------------------------------------%
% Figure of the Laplacian Dissipation %
%-------------------------------------%
nx = 64;
dt = 1;

%--- polynomical of small-scale dissipation
ksig   =  [32];
psig   =  [1];
rtsig  =  [0.1];  % 1/tsig in units of dt 

%--- polynomical of large-scale dissipation
klam   =  [2];
plam   =  [-1];
rtlam  =  [0.001];  % 1/tlam in units of dt

%--- wave number grid
kmax   = fix(nx/3);

%============== internal part ============

%--- spectral grid
k =  1:0.1:kmax;

sig = rtsig*(1/ksig)^(2*psig);
lam = rtlam*(1/klam)^(2*plam);

r   = k;
r2  = k.^2;

D_sig   = -sig*r2.^psig;
D_lam   = -lam*r2.^plam;

D       = D_sig + D_lam;

ff = figure
ax = axes
set(ax,'FontSize',[16]);
set(ax,'FontWeight','bold');
set(ax,'XLim',[1,kmax]);
hold
grid on
pp = plot(r,exp(dt*(D_sig + D_lam)),'k');
set(pp,'LineWidth',[2]);
