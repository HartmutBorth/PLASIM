function  [out_amp,out_arg] = f_readsp(fname,vname)
% F_READSP reads the cat output in spectral format and transforms it
% to a format centred in spectral space. For octave users f_readsp
% only works if package <netcdf> is installed and loaded. The variables
% to be read are assumed to be in the spectral CAT output format.
% For more details see below. 
%
% syntax
%  [out_amp,out_arg] = f_readsp(fname,vname)
%
% description
%   input arguments
%     fname    : string, specifying the netcdf input file.
%     vname    : string, specifying the name of the variable to
%                be loaded.
%   output arguments
%     out_amp : array with the amplitude of the complex spectral array.
%     out_arg : array with the phase of the complex spectral array as
%               fractions of pi.
%
% The input and output formats are specified as follows:
%
% Dimensions:
% [nfx+1,nfy+1,nt].
%
%  with:
% -------
%  nfx = 2*nkx+1 with nkx = ngx/3 (integer division)
%  nfy = 2*nky   with nky = ngy/3 (integer division)
%  nt  = number of time-steps
% 
%  in the first dimension 2*(nkx+1) alternating real and 
%  imaginary part of the positive wave-numbers k_x in 
%  x-direction including the zero wave number are stored. 
%  In the second dimension the wave numbers k_y in y-direction
%  are stored.  
%
%  The first input spectral coordinate is given by
%
%             | 
%  Re(0)      |
%  Im(0)      |
%  Re(1)      |
%  Im(1)      |
%             | 
%    .      
%    .       k_x 
%    .        
%             | 
%  Re(ngx/3)  |
%  Im(ngx/3)  |
%             | 
%             V
%
%  and the second spectral coordinate by
%
%
%  -------------------- k_y ------------------>
%
%  0,1,2, ... , ngy/3,-ngy/3,-ngy/3+1, ... ,-1,  
%
%  where ngx/3 and ngy/3 are the integer parts remaining after division
%
%  The third dimension is the time coordinate, which remains the same
%  for input and output.
%
%
%  The first output spectral coordinate is given by
%         
%    0   
%    1  
%    2      |
%    .      |
%    .       
%    .     k_x
%    .       
%    .      |
%    .      V 
%    .          
%  ngx/3-1     
%  ngx/3          
%
%  The second output spectral coordinate is given by
% 
%
%               -------- k_y ------->
%
%  -ngy/3,-ngy/3+1, ... ,-1,0,1,2, ... , ngy/3
%
%  On the output grid the amplitude and the phase of the complex
%  fields are given, i.e. we can write
%
%        z =  Re(F) + i Im(F) = amplitude(F) exp(i phase(F))
%
%  The spectral fields F(k_x,k_y) are characterized by the 
%  periodicity condition
%
%  F(k_x+N,k_y+M) = F(k_x,k_y) with the periods N,M
%
%  and the symmetry properties 
%
%  F(k_x,k_y) = F^*(-k_x,-k_y)  with F^* the conjugate complex 
%  value of F
%
%  This symmetry property allows spectral fields to be represented 
%  only on a part of the wave number space (The rest can be 
%  complemented by the symmetry properties)
%
%  It follows that F(0,0) (the average of F) is real and
%  that |F(0,k_y)| = |F(0,-k_y)|. Moreover we have that 
%  phase(F(0,k_y)) = - phase(F(0,-k_y)).
%
%  To complement the amplitudes, just take a point mirror image
%  at the origin (k_x,k_y) = (0,0). For the phases one has in 
%  addition to change the sign.
%
%  To guarantee that the original field f(x,y,t) is real one has to
%  satisfy the symmetry properties on the line (k_x = 0).
%
%--------------------------------------------------------------------

%--- read variable from netcdf-file
[var] = ncread(fname,vname);

%--- determine size of input field
nx = size(var,1);
ny = size(var,2);
nt = size(var,3);

ireal = 1:2:nx;
iim   = 2:2:nx;

varreal = var(ireal,:,:);
varim   = var(iim,:,:);

%--- determine amplitude
out_amp = sqrt(varreal.^2 + varim.^2);

%--- determine argument (phase)
for kk = 1:nt
   rtmp            = varreal(:,:,kk);
   itmp            = varim(:,:,kk);
   out_arg(:,:,kk) = atan2(itmp,rtmp)/pi;
end

%--- centralize wave-numbers in y-direction
nysplit = 0.5*(ny-1) + 1;

out_tmp1 = out_amp(:,1:nysplit,:);
out_tmp2 = out_amp(:,nysplit+1:ny,:);
out_amp  = [out_tmp2 out_tmp1];

out_tmp1 = out_arg(:,1:nysplit,:);
out_tmp2 = out_arg(:,nysplit+1:ny,:);
out_arg  = [out_tmp2 out_tmp1];
