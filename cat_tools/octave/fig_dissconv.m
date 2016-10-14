%------------------------------------%
% Dissipation and Spectral Filtering %
%------------------------------------%

%--- grid
ngx = 64;    % number of grid points in x-direction
ngy = 64;    % number of grid points in y-direction

%ngx = 128;    % number of grid points in x-direction
%ngy = 128;    % number of grid points in y-direction

%--- control flags
tit_f = 0; % 0/1 print title no/yes
pri_f = 1; % 0/1 save figures
pdf_f = 1; % 0/1 convert figures to pdf

%--- gaussian vorticity field
%qtp = 'gauss';
qtp = 'rank';
%qtp = 'jet';

%dtp = 'trunc';
dtp = 'lap';

%------------%
% parameters %
%------------%


% gaussian vortices 'gauss'
%--------------------------


%--- positions of vortex centers (grid points)
%x0 = [8 16 24 32 40 48 56]; 
%y0 = [8 16 24 32 40 48 56];


%x0 = [32 64 96];
%y0 = [32 64 96];

%x0 = [64];
%y0 = [64];

x0 = [32];
y0 = [32];


%--- half-width in grid points
hw  = 4;

%--- vortex radius in grid points
vrad = 8;

%--- amplitude of vortices
am  = 1;


% jet 'jet' 
%------------------
hwjet    = 20;  % half width of jet
wsheet   = 10;  % width of vortex sheet


%--- dissipation by circular truncation
rkmax = 18;

%--- dissipation by Laplacian, e.g. gaussian filter
rkhalf = 18;      % half-width of filter
rts    = 1;       % dt/tscale 


%===============%
% internal part %
%===============%

%--- create physical grid
[X,Y] = meshgrid(1:ngx,1:ngy);


%--- truncation mask of CAT 
nkx = fix(ngx/3);
nky = fix(ngy/3);


%--- create spectral grid
kx  = [0:1:0.5*ngx -0.5*ngx+1:1:-1];
ky  = [0:1:0.5*ngy -0.5*ngy+1:1:-1];

[KX,KY] = meshgrid(kx,ky);
rk      = sqrt(KX.^2 + KY.^2);


%--- create vorticity field
switch qtp
  case 'gauss'
    fx    = zeros(size(X));
    fxtmp = zeros(size(X));
    nx0 = length(x0);
    for kk = 1:nx0
      x0tmp = x0(kk);
      y0tmp = y0(kk);
      nfac  = 1/(2*pi*hw^2); 
%     fxtmp = am*nfac*exp(-((X-x0tmp).^2 + (Y-y0tmp).^2)./(2*hw^2));
      fxtmp = am*exp(-((X-x0tmp).^2 + (Y-y0tmp).^2)./(2*hw^2));
      fx    = fx + fxtmp;
    end
  case 'rank'
    fx    = zeros(size(X));
    fxtmp = zeros(size(X));
    nx0 = length(x0);
    for kk = 1:nx0
      x0tmp = x0(kk);
      y0tmp = y0(kk);
      r2tmp = (X-x0tmp).^2 + (Y-y0tmp).^2;
      fxtmp(find(r2tmp <= vrad.^2)) = am;
      fx    = fx + fxtmp;
    end
  case 'jet'
    fx    = zeros(size(X));
    fxtmp = zeros(size(X));
    cs = 0.5*ngx-hwjet;
    cn = 0.5*ngx+1+hwjet;
 
    is = cs-wsheet:1:cs; % index south sheet
    in = cn:1:cn+wsheet; % index north sheet
  
    fx(is,:)  =  1;
    fx(in,:)  = -1;
endswitch


%--- create dissipation operator
switch dtp
  case 'trunc'
    diss_op = zeros(size(rk));
    diss_op(find(rk <= rkmax)) = 1;
  case 'lap'
    diss_op = exp(-(rk./rkhalf).^2 * rts);
endswitch


%--- Vorticity in spectral space and truncation
Ffx = fft2(fx);

%--- Truncation of CAT 
indx = nkx+2:1:ngx-nkx;
indy = nky+2:1:ngy-nky;
Ffx(indy,indx) = 0;

RFfx = real(Ffx);
IFfx = imag(Ffx);
AFfx = abs(Ffx);


%--- Vorticity in physical space after spectral filtering
IFFfx = ifft2(Ffx);


%--- Dissipation with spectral filtering
DissFfx = Ffx.*diss_op;
Dissfx  = ifft2(DissFfx);

%---------%
% Figures %
%---------%

%--- vorticity in phyiscal space
ff = figure
ax = axes
set(ax,'YLim',[1 ngy]);
set(ax,'XLim',[1 ngx]);
set(ax,'FontSize',[16]);
set(ax,'FontWeight','bold');

hold
if tit_f == 1
  title(['Vorticity in Physical Space']);
end
pcolor(fx);

cb=colorbar;
set(cb,'FontSize',[12]);
set(cb,'FontWeight',['bold']);

fname = [num2str(ngx) qtp dtp 'phys'];
if pri_f == 1
  print('-depsc', fname);
end
if pdf_f == 1
  system(['convert' ' ' fname '.eps' ' ' fname '.pdf']);
end

%--- vorticity in spectral space
ff = figure
ax = axes
set(ax,'YLim',[-0.5*ngy,0.5*ngy-1]);
set(ax,'XLim',[-0.5*ngx,0.5*ngx-1]);
set(ax,'FontSize',[16]);
set(ax,'FontWeight','bold');

hold
if tit_f == 1  
  title(['Vorticity in Spectral Space']);
end
ktmp = -0.5*ngx:1:0.5*ngx-1;
ATmp = fftshift(AFfx);
pp = pcolor(ktmp,ktmp,ATmp);

cb=colorbar;
set(cb,'FontSize',[12]);
set(cb,'FontWeight',['bold']);

fname = [num2str(ngx) qtp dtp 'spec'];
if pri_f == 1
  print('-depsc', fname);
end
if pdf_f == 1
  system(['convert' ' ' fname '.eps' ' ' fname '.pdf']);
end


%--- vorticity in physical space after spectral filtering
ff = figure
ax = axes
set(ax,'YLim',[1,ngy]);
set(ax,'XLim',[1,ngx]);
set(ax,'FontSize',[16]);
set(ax,'FontWeight','bold');

hold
if tit_f == 1  
  title(['Vorticity in Physical Space after spectral Filtering']);
end
pp = pcolor(IFFfx);

cb=colorbar;
set(cb,'FontSize',[12]);
set(cb,'FontWeight',['bold']);

fname = [num2str(ngx) qtp dtp 'filt' qtp dtp];
if pri_f == 1
  print('-depsc', fname);
end
if pdf_f == 1
  system(['convert' ' ' fname '.eps' ' ' fname '.pdf']);
end

%--- error ifft2(fft2(fx))  
ff = figure
ax = axes
set(ax,'YLim',[1,ngy]);
set(ax,'XLim',[1,ngx]);
set(ax,'FontSize',[16]);
set(ax,'FontWeight','bold');

hold
if tit_f == 1  
  title(['Difference (fx - IF(F(fx))']);
end
pp = pcolor(fx-IFFfx);

cb=colorbar;
set(cb,'FontSize',[12])
set(cb,'FontWeight',['bold'])

fname = [num2str(ngx) qtp dtp 'errfilt'];
if pri_f == 1
  print('-depsc', fname);
end
if pdf_f == 1
  system(['convert' ' ' fname '.eps' ' ' fname '.pdf']);
end

%--- Vorticity after application of Dissipation (spectral filtering)
ff = figure
ax = axes
set(ax,'YLim',[1 ngy]);
set(ax,'XLim',[1 ngx]);
set(ax,'FontSize',[16]);
set(ax,'FontWeight','bold');

hold
if tit_f == 1  
  title(['Vorticity in Physical Space after Dissipation']);
end
pcolor(Dissfx);

cb=colorbar;
set(cb,'FontSize',[12]);
set(cb,'FontWeight',['bold']);

fname = [num2str(ngx) qtp dtp 'diss'];
if pri_f == 1
  print('-depsc', fname);
end
if pdf_f == 1
  system(['convert' ' ' fname '.eps' ' ' fname '.pdf']);
end

%--- Difference
ff = figure
ax = axes
set(ax,'YLim',[1 ngy]);
set(ax,'XLim',[1 ngx]);
set(ax,'FontSize',[16]);
set(ax,'FontWeight','bold');

hold
if tit_f == 1  
  title(['Difference (fx-Dissfx)']);
end
Dtmp = fx-Dissfx;
pcolor(Dtmp);

cb=colorbar;
set(cb,'FontSize',[12]);
set(cb,'FontWeight',['bold']);

fname = [num2str(ngx) qtp dtp 'errdiss'];
if pri_f == 1
  print('-depsc', fname);
end
if pdf_f == 1
  system(['convert' ' ' fname '.eps' ' ' fname '.pdf']);
end
