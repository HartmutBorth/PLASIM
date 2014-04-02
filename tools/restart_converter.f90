module rcmod
implicit none 

integer :: verb =    1
integer :: nlat =   32
integer :: nlon =   64
integer :: nlev =   10
integer :: nrsp =  506
integer :: nhor = 2048
integer :: nlsoil =  5
integer :: nlev_oce = 1

integer :: nstep
integer :: naccuout

character (len=16) :: yname

real, allocatable :: sphor(:)
real, allocatable :: splev(:,:)
real, allocatable :: gphor(:)
real, allocatable :: gplev(:,:)
real, allocatable :: gpm14(:,:)
real, allocatable :: gpm12(:,:)
end module rcmod

program rescon
use rcmod

allocate(sphor(nrsp))
allocate(splev(nrsp,nlev))
allocate(gphor(nhor))
allocate(gplev(nhor,nlev))
allocate(gpm14(nhor,14))
allocate(gpm12(nhor,12))

open(10,file='puma_restart',form='unformatted')
open(20,file='plasim_restart',form='unformatted')

if (verb == 1) then
   print '(77("*"))'
   print '("* Planet Simulator Restart File Converter 1.0                               *")'
   print '(77("*"))'
   print '("* Part I      * puma_restart                                                *")'
   print '(77("*"))'
endif

! scalars

call scai

if (verb == 1) print '(77("*"))'

! spectral arrays

call sp3d('sz              ')
call sp3d('sd              ')
call sp3d('st              ')
call sp3d('sq              ')
call sp2d('sp              ')

call sp3d('szm             ')
call sp3d('sdm             ')
call sp3d('stm             ')
call sp3d('sqm             ')
call sp2d('spm             ')

call sp2d('so              ')
call sp3d('sr              ')

! gridpoint arrays

call gp2d('dls             ')
call gp2d('drhs            ')
call gp2d('dalb            ')
call gp2d('dz0             ')
call gp2d('dicec           ')
call gp2d('diced           ')
call gp2d('dwatc           ')
call gp2d('drunoff         ')
call gp2d('dveg            ')
call gp2d('dforest         ')

call gp3d('dcc             ')
call gp3d('dql             ')
call gp3d('dqsat           ')

call gp2d('dt              ')
call gp2d('dq              ')
call gp2d('dust3           ')

! accumulated diagnostics

call gp2d('aprl            ')
call gp2d('aprc            ')
call gp2d('aprs            ')
call gp2d('aevap           ')
call gp2d('ashfl           ')
call gp2d('alhfl           ')
call gp2d('aroff           ')
call gp2d('asmelt          ')
call gp2d('asndch          ')
call gp2d('acc             ')
call gp2d('assol           ')
call gp2d('asthr           ')
call gp2d('atsol           ')
call gp2d('atthr           ')
call gp2d('ataux           ')
call gp2d('atauy           ')
call gp2d('atsolu          ')
call gp2d('assolu          ')
call gp2d('asthru          ')
call gp2d('aqvi            ')
call gp2d('atsa            ')
call gp2d('ats0            ')

close (10)

open(10,file='land_restart',form='unformatted')

read (10) nstep_land

if (nstep /= nstep_land) then
   print *,'nstep in puma_restart = ',nstep
   print *,'ntspe in land_restart = ',nstep_land
   print *,'*** error stop ***'
   stop
endif


if (verb == 1) then
   print '(77("*"))'
   print '("* Part II     * land_restart                                                *")'
   print '(77("*"))'
endif

! landmod arrays

call gp2d('dtsl            ')
call gp2d('dtsm            ')
call gp2d('dqs             ')
call gp2d('driver          ')
call gp2d('duroff          ')
call gp2d('dvroff          ')
call gp2d('darea           ')
call gp2d('dwmax           ')
call gp12('dtcl            ')
call gp12('dwcl            ')
call gp2d('dsnowt          ')
call gp2d('dsnowz          ')
call gp3l('dsoilt          ')
call gp2d('dcveg           ')
call gp2d('dcsoil          ')
call gp2d('dglac           ')
call gp2d('dz0clim         ')
call gp2d('dz0climo        ')
call gp2d('dalbclim        ')

! simba arrays

call gp2d('agpp            ')
call gp2d('anpp            ')
call gp2d('agppl           ')
call gp2d('agppw           ')
call gp2d('anogrow         ')
call gp2d('aresh           ')
call gp2d('alitter         ')

close (10)

open(10,file='sea_restart',form='unformatted')

read (10) nstep_sea,naccua

if (nstep /= nstep_sea) then
   print *,'nstep in puma_restart = ',nstep
   print *,'ntspe in sea_restart  = ',nstep_sea
   print *,'*** error stop ***'
   stop
endif

yname = 'naccua'
write (20) yname
write (20) naccua

if (verb == 1) then
   print '(77("*"))'
   print '("* Part III    * sea_restart                                                 *")'
   print '(77("*"))'
   print 1000,'naccua',naccua
endif

! seamod arrays

call gp2d('dts             ')
call gp2d('cheata          ')
call gp2d('cpmea           ')
call gp2d('cprsa           ')
call gp2d('croffa          ')
call gp2d('ctauxa          ')
call gp2d('ctauya          ')
call gp2d('cust3a          ')
call gp2d('cshfla          ')
call gp2d('cshdta          ')
call gp2d('clhfla          ')
call gp2d('clhdta          ')
call gp2d('cswfla          ')
call gp2d('clwfla          ')

close (10)

open(10,file='ice_restart',form='unformatted')

read (10) nstep_ice,naccuice,naccuo

! Ignore nstep from icemod and use nstep from plasim

!if (nstep /= nstep_ice) then
!   print *,'nstep in puma_restart = ',nstep
!   print *,'nstep in ice_restart  = ',nstep_ice
!   print *,'*** error stop ***'
!   stop
!endif

yname = 'naccuice'
write (20) yname
write (20) naccuice
yname = 'naccuo'
write (20) yname
write (20) naccuo

if (verb == 1) then
   print '(77("*"))'
   print '("* Part IV     * ice_restart                                                 *")'
   print '(77("*"))'
   print 1000,'naccuice',naccuice
   print 1000,'naccuo',naccuo
endif

! icemod arrays

call gp2d('xls             ')
call gp2d('xts             ')
call gp2d('xicec           ')
call gp2d('xiced           ')
call gp2d('xsnow           ')

call gp14('xclicec         ')
call gp14('xcliced         ')
call gp14('xclsst          ')

call gp2d('cheat           ')
call gp2d('cfresh          ')
call gp2d('ctaux           ')
call gp2d('ctauy           ')
call gp2d('cust3           ')
call gp2d('xflxicea        ')
call gp2d('xheata          ')
call gp2d('xofluxa         ')
call gp2d('xqmelta         ')
call gp2d('xcfluxa         ')
call gp2d('xsmelta         ')
call gp2d('ximelta         ')
call gp2d('xtsfluxa        ')
call gp2d('xhnewa          ')

close (10)

open(10,file='ocean_restart',form='unformatted')

read (10) nstep_oce,naccuoce

yname = 'naccuoce'
write (20) yname
write (20) naccuice
yname = 'nlev_oce'
write (20) yname
write (20) nlev_oce

if (verb == 1) then
   print '(77("*"))'
   print '("* Part V      * ocean_restart                                               *")'
   print '(77("*"))'
   print 1000,'naccuoce',naccuoce
endif

! oceanmod arrays

call gp2d('yls             ')
call gp3o('ysst            ')
call gp2d('yiflux          ')
call gp14('yclsst          ')
call gp2d('yheata          ')
call gp2d('yfssta          ')
call gp2d('yifluxa         ')
call gp2d('ydssta          ')

close (20)

if (verb == 1) print '(77("*"))'

stop
1000 format("* Integer     * ",a8," * ",I12,36X," *")
end program rescon


subroutine scai
use rcmod

character (len=16) :: yn

read (10) nstep,naccuout
yn = 'nstep   '
write (20) yn
write (20) nstep
yn = 'naccuout'
write (20) yn
write (20) naccuout
yn = 'nlat    '
write (20) yn
write (20) nlat
yn = 'nlon    '
write (20) yn
write (20) nlon
yn = 'nlev    '
write (20) yn
write (20) nlev
yn = 'nrsp    '
write (20) yn
write (20) nrsp
yn = 'nlsoil  '
write (20) yn
write (20) nlsoil

if (verb == 1) then
   print 1000,'nstep   ',nstep
   print 1000,'naccuout',naccuout
   print 1000,'nlat    ',nlat 
   print 1000,'nlon    ',nlon    
   print 1000,'nlev    ',nlev 
   print 1000,'nrsp    ',nrsp    
   print 1000,'nlsoil  ',nlsoil  
endif
return
1000 format("* Integer     * ",a8," * ",I12,36X," *")
end subroutine scai


subroutine sp3d(yn)
use rcmod

character (len=16) :: yn

read (10) splev
write (20) yn
write (20) splev

if (verb == 1) print 1000,yn,splev(1,1:4)
return
1000 format("* Spectral 3D * ",a8," * ",4e12.3," *")
end subroutine sp3d


subroutine sp2d(yn)
use rcmod

character (len=16) :: yn

read (10) sphor
write (20) yn
write (20) sphor

if (verb == 1) print 1000,yn,sphor(1:4)
return
1000 format("* Spectral 2D * ",a8," * ",4e12.3," *")
end subroutine sp2d


subroutine gp3d(yn)
use rcmod

character (len=16) :: yn

read (10) gplev
write (20) yn
write (20) gplev

if (verb == 1) print 1000,yn,gplev(1,1:4)
return
1000 format("* Grid     3D * ",a8," * ",4e12.3," *")
end subroutine gp3d


subroutine gp2d(yn)
use rcmod

character (len=16) :: yn

read (10) gphor
write (20) yn
write (20) gphor

if (verb == 1) print 1000,yn,gphor(1:4)
return
1000 format("* Grid     2D * ",a8," * ",4e12.3," *")
end subroutine gp2d


subroutine gp3l(yn)
use rcmod

character (len=16) :: yn

read (10) gplev(:,1:nlsoil)
write (20) yn
write (20) gplev(:,1:nlsoil)

if (verb == 1) print 1000,yn,gplev(1,1:4)
return
1000 format("* Grid     3D * ",a8," * ",4e12.3," *")
end subroutine gp3l


subroutine gp14(yn)
use rcmod

character (len=16) :: yn

read (10) gpm14(:,2:13)
gpm14(:, 1) = gpm14(:,13)
gpm14(:,14) = gpm14(:, 2) 
write (20) yn
write (20) gpm14

if (verb == 1) print 1000,yn,gpm14(1,1:4)
return
1000 format("* Grid annual * ",a8," * ",4e12.3," *")
end subroutine gp14


subroutine gp12(yn)
use rcmod

character (len=16) :: yn

read (10) gpm12
write (20) yn
write (20) gpm12

if (verb == 1) print 1000,yn,gpm12(1,1:4)
return
1000 format("* Grid 12-mon * ",a8," * ",4e12.3," *")
end subroutine gp12


subroutine gp3o(yn)
use rcmod

character (len=16) :: yn

read (10) gplev(:,1:nlev_oce)
write (20) yn
write (20) gplev(:,1:nlev_oce)

if (verb == 1) print 1000,yn,gplev(1:4,1)
return
1000 format("* Grid     3D * ",a8," * ",4e12.3," *")
end subroutine gp3o

