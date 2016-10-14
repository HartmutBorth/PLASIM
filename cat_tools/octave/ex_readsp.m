%--------------------------------------%
% Example for reading a spectral field %  
%--------------------------------------%

%--- plotting flag
plot_flag = "conf";     %  plot options:
                        % --------------
                        %
                        %   "checker"  : checkerboard plot 
                        %   "conf"     : filled contours

%--- define resolution
nx     = 64;
ny     = 64;

%--- define i/o
infile = ["cat_sp.nc"];
invar  = ["var138"];

%--- read amplitudes and phases
[qamp,qarg] = f_readsp(infile,invar);

%--- define centered wave numbers
k_x    =  0:1:fix(nx/3);
k_y    = -fix(ny/3):1:fix(ny/3);

%--- sum up spectral field of the time coordinate
sumqamp = sum(qamp,3);

%--- rotate field
sumqamp = sumqamp';


figure
set(gca,'DataAspectRatio',[1 1 1]);
set(gca,'FontSize',[15]);
set(gca,'FontWeight','bold');
set(gca,'XLim',[min(k_x),max(k_x)]);
set(gca,'YLim',[min(k_y),max(k_y)]);

hold
title(["Cumulated Vorticity in Centered Spectral Space"]);
xlabel(["k_x"]);
ylabel(["k_y"]);

switch plot_flag
   case("checker")
      pcolor(k_x,k_y,sumqamp);
   case("conf")
      contourf(k_x,k_y,sumqamp);
      grid on;
   otherwise
      disp([plot_flag " is unknown option"]);
end

hold
