	function zz=destagger(zz)

%	function zo=destagger(z)
%       Destaggers an Arakawa E field for plotting
%	(JvH 12/06)

	[nx,ny,nt]=size(zz);

	for it=1:nt
          for j=1:ny
              zz(:,j,it)=circshift(zz(:,j,it),-floor((j-1)/2)); 
          end
        end

