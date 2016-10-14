%-------------------------------------%
% Figure of the Laplacian Dissipation %
%-------------------------------------%
nx = 64;
dt = 1;


%--- polynomical of small-scale dissipation
rsig    =  [64];
nsig2    =  [2];
nsig3    =  [3];
tsigm1  =  [0.1]; % tsig^(-1) in units of dt 

%--- polynomical of large-scale dissipation
rlam    =  [1];
nlam    =  [1];
tlamm1  =  [0.01]; % tlam^{-1} in units of dt

%--- wave number grid
kx   = fix(nx/3);

%============== internal part ============

%--- spectral grid
k   =  0:0.1:kx;
rk  =  k;

D_nsig2 = (rk./rsig).^(2*(nsig2-1));
D_nsig3 = (rk./rsig).^(2*(nsig3-1));

D_lamn = (rk./rlam).^(-2*nlam);


expDnsig2 = exp(-D_nsig2*(dt*tsigm1));
expDnsig3 = exp(-D_nsig3*(dt*tsigm1));
expDlamn = exp(-D_lamn*(dt*tlamm1));

ff = figure
ax = axes
set(ax,'FontSize',[16]);
set(ax,'FontWeight','bold');
set(ax,'XLim',[0,kx]);
set(ax,'YLim',[0.9,1]);
hold
grid on
pp1 = plot(rk,expDnsig2,'k');
set(pp1,'LineWidth',[2]);
pp2 = plot(rk,expDlamn,'k--');
set(pp2,'LineWidth',[2]);

pp3 = plot(rk,expDnsig3,'r');
set(pp1,'LineWidth',[2]);


