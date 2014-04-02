module splitmod
implicit none
integer :: nlon,nlat,ntru
integer :: icemon = 0
integer :: sstmon = 0
character(len=1024) :: ifile,odir
integer :: head(8)
real, parameter :: GA = 9.81 ! g
real, parameter :: TMELT = 273.16 
real, allocatable :: f(:,:)
character(len=1), allocatable :: dls(:,:)
character(len=1), allocatable :: dglac(:,:)
character(len=1), allocatable :: dice(:,:)
real,  allocatable :: doro(:,:)

end module splitmod

program surfsplit
use splitmod
implicit none

interface
   integer function iargc()
   end function iargc
end interface

integer :: nargc,iostat

nargc = iargc()

if (nargc /= 2) then
   print '("Usage surfsplit <surface_file> <output_directory>")'
   stop
endif

call getarg(1,ifile)
call getarg(2,odir)

open(10,file=ifile)

print '(46("*"))'
print '("* Surface file splitter for Planet Simulator *")'
print '(46("*"))'

read (10,'(8I10)',iostat=iostat) head(:)
if (iostat /= 0) then
   print '("Could not read from file <",a)',ifile
   stop
endif

nlon = head(5)
nlat = head(6)
ntru = (nlon-1) / 3

allocate(f(nlon,nlat))
allocate(dls(nlon,nlat))
allocate(dglac(nlon,nlat))
allocate(dice(nlon,nlat))
allocate(doro(nlon,nlat))

do while (iostat == 0)
   print '(6I10)',head(1:6)
   read (10,'(8E12.6)') f(:,:)
   if (head(1) == 129 ) call subdoro
   if (head(1) == 169 ) call subdsst
   if (head(1) == 172 ) call subdls
   if (head(1) == 173 ) call subdz0
   if (head(1) == 1730) call subdz0noveg
   if (head(1) == 174 ) call subdalb
   if (head(1) == 209 ) call subdsoil
   if (head(1) == 210 ) call subdice
   if (head(1) == 232 ) call subdglac
   read (10,'(8I10)',iostat=iostat) head(:)
enddo ! while (iostat == 0)

stop
end

subroutine subdls
use splitmod
implicit none
integer :: jlon,jlat
character(len=80) yform

character (len=1024) :: ofile
write (ofile,'(a,"/t",i3.3,"_landsea.txt")') trim(odir),ntru

open(11,file=ofile)
dls(:,:) = ' '
where (f(:,:) > 0.5) dls(:,:) = '*'
write(11,'("[PLASIM SURFACE ARRAY]")')
write(11,'("CODE    = 172")')
write(11,'("ARRAY   = dls  ")')
write(11,'("COMMENT = Land Sea Mask <*> = Land")')
write(11,'("NLON    = ",I4)') nlon
write(11,'("NLAT    = ",I4)') nlat
write(11,'("FORMAT  = (",I3,"A1)")') nlon
write(11,'("[BEGIN ARRAY]")')
write(yform,'("(",i,"a1)")') nlon
write(11,yform) dls(:,:)
write(11,'("[END ARRAY]")')
close(11)
return
end

subroutine subdglac
use splitmod
implicit none
integer :: jlon,jlat
character(len=80) yform

character (len=1024) :: ofile
write (ofile,'(a,"/t",i3.3,"_mask.txt")') trim(odir),ntru

open(11,file=ofile)
dglac(:,:) = dls(:,:)
where (f(:,:) > 0.5) dglac(:,:) = '='
write(11,'("[PLASIM SURFACE ARRAY]")')
write(11,'("CODE    = 232")')
write(11,'("ARRAY   = dglac")')
write(11,'("COMMENT = Glacier Mask <=>")')
write(11,'("NLON    = ",I4)') nlon
write(11,'("NLAT    = ",I4)') nlat
write(11,'("FORMAT  = (",I3,"A1)")') nlon
write(11,'("[BEGIN ARRAY]")')
write (yform,'("(",i,"a1)")') nlon
write (11,yform) dglac(:,:)
write(11,'("[END ARRAY]")')
close(11)
return
end

subroutine subdice
use splitmod
implicit none
integer :: jlon,jlat
character(len=80) yform
character (len=1024) :: ofile

icemon = icemon + 1
if (icemon == 1) then
   write (ofile,'(a,"/t",i3.3,"_ice.txt")') trim(odir),ntru
   open(11,file=ofile)
endif

dice(:,:) = dglac(:,:)
where (f(:,:) > 0.9) dice(:,:) = '.'
write (yform,'("(",i,"a1)")') nlon
write(11,'("[PLASIM SURFACE ARRAY]")')
write(11,'("CODE    = 210")')
write(11,'("ARRAY   = xclicec")')
write(11,'("COMMENT = Sea Ice Cover <.> And Glacier Mask <=>")')
write(11,'("NLON    = ",I4)') nlon
write(11,'("NLAT    = ",I4)') nlat
write(11,'("FORMAT  = (",I3,"A1)")') nlon
write(11,'("MONTH   = ",I2)') icemon
write(11,'("[BEGIN ARRAY]")')
write(11,yform) dice(:,:)
write(11,'("[END ARRAY]")')
return
end

subroutine subdoro
use splitmod
implicit none
integer :: jlon,jlat
character (len=1024) :: ofile
integer :: itopo(nlon,nlat)

write (ofile,'(a,"/t",i3.3,"_topo.txt")') trim(odir),ntru
open(11,file=ofile)

doro(:,:) = f(:,:)
itopo(:,:) = nint(f(:,:) / GA)

write(11,'("[PLASIM SURFACE ARRAY]")')
write(11,'("CODE    = 129")')
write(11,'("ARRAY   = doro")')
write(11,'("UNIT    = gpm")')
write(11,'("COMMENT = Topography = Orography / g")')
write(11,'("NLON    = ",I4)') nlon
write(11,'("NLAT    = ",I4)') nlat
write(11,'("FORMAT  = (16I5)")')
write(11,'("[BEGIN ARRAY]")')
write(11,'(16I5)') itopo(:,:)
write(11,'("[END ARRAY]")')
close(11)
return
end

subroutine subdz0
use splitmod
implicit none
integer :: jlon,jlat
character (len=1024) :: ofile
integer :: itopo(nlon,nlat)

write (ofile,'(a,"/t",i3.3,"_z0.txt")') trim(odir),ntru
open(11,file=ofile)


write(11,'("[PLASIM SURFACE ARRAY]")')
write(11,'("CODE    = 173")')
write(11,'("ARRAY   = dz0clim")')
write(11,'("COMMENT = Climatological roughness length")')
write(11,'("NLON    = ",I4)') nlon
write(11,'("NLAT    = ",I4)') nlat
write(11,'("FORMAT  = (8F10.6)")')
write(11,'("[BEGIN ARRAY]")')
write(11,'(8F10.6)') f(:,:)
write(11,'("[END ARRAY]")')
close(11)
return
end

subroutine subdz0noveg
use splitmod
implicit none
integer :: jlon,jlat
character (len=1024) :: ofile
integer :: itopo(nlon,nlat)

write (ofile,'(a,"/t",i3.3,"_z0noveg.txt")') trim(odir),ntru
open(11,file=ofile)


write(11,'("[PLASIM SURFACE ARRAY]")')
write(11,'("CODE    = 1730")')
write(11,'("ARRAY   = dz0climo")')
write(11,'("COMMENT = Clim. roughness length without vegetation")')
write(11,'("NLON    = ",I4)') nlon
write(11,'("NLAT    = ",I4)') nlat
write(11,'("FORMAT  = (8F10.6)")')
write(11,'("[BEGIN ARRAY]")')
write(11,'(8F10.6)') f(:,:)
write(11,'("[END ARRAY]")')
close(11)
return
end

subroutine subdsoil
use splitmod
implicit none
integer :: jlon,jlat
character (len=1024) :: ofile
integer :: itopo(nlon,nlat)

write (ofile,'(a,"/t",i3.3,"_soil.txt")') trim(odir),ntru
open(11,file=ofile)

where (f(:,:) > 100-0) f(:,:) = f(:,:) - TMELT

write(11,'("[PLASIM SURFACE ARRAY]")')
write(11,'("CODE    = 209")')
write(11,'("ARRAY   = dtclsoil")')
write(11,'("UNIT    = Celsius")')
write(11,'("COMMENT = Climatological Soil Temperature")')
write(11,'("NLON    = ",I4)') nlon
write(11,'("NLAT    = ",I4)') nlat
write(11,'("FORMAT  = (8F7.2)")')
write(11,'("[BEGIN ARRAY]")')
write(11,'(8F7.2)') f(:,:)
write(11,'("[END ARRAY]")')
close(11)
return
end

subroutine subdalb
use splitmod
implicit none
integer :: jlon,jlat
character (len=1024) :: ofile
integer :: itopo(nlon,nlat)

write (ofile,'(a,"/t",i3.3,"_albedo.txt")') trim(odir),ntru
open(11,file=ofile)


write(11,'("[PLASIM SURFACE ARRAY]")')
write(11,'("CODE    = 174")')
write(11,'("ARRAY   = dalbclim")')
write(11,'("COMMENT = Climatological Albedo")')
write(11,'("NLON    = ",I4)') nlon
write(11,'("NLAT    = ",I4)') nlat
write(11,'("FORMAT  = (8F8.4)")')
write(11,'("[BEGIN ARRAY]")')
write(11,'(8F8.4)') f(:,:)
write(11,'("[END ARRAY]")')
close(11)
return
end

subroutine subdsst
use splitmod
implicit none
integer :: jlon,jlat
character(len=80) yform
character (len=1024) :: ofile

sstmon = sstmon + 1
if (sstmon == 1) then
   write (ofile,'(a,"/t",i3.3,"_sst.txt")') trim(odir),ntru
   open(11,file=ofile)
endif

where (f(:,:) > 100-0) f(:,:) = f(:,:) - TMELT

write (yform,'("(",i,"a1)")') nlon
write(11,'("[PLASIM SURFACE ARRAY]")')
write(11,'("CODE    = 169")')
write(11,'("ARRAY   = xclsst")')
write(11,'("UNIT    = Celsius")')
write(11,'("COMMENT = SST and Surface Temperature for Land")')
write(11,'("NLON    = ",I4)') nlon
write(11,'("NLAT    = ",I4)') nlat
write(11,'("MONTH   = ",I4)') sstmon
write(11,'("FORMAT  = (8F7.2)")')
write(11,'("[BEGIN ARRAY]")')
write(11,'(8F7.2)') f(:,:)
write(11,'("[END ARRAY]")')
return
end
