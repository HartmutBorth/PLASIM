module simmod
use catmod

!--- flags and switches
logical :: lsimnl = .false.    ! .true. if sim_namelist is present in run 
                               ! directory

integer :: ios_simnl = 1       ! 0 if <sim_namelist> is readable

!--- i/o units
integer, parameter :: nusimnl    = 110  ! sim namelist
integer, parameter :: nusimgp    = 120  ! grid point fields
integer, parameter :: nusimsp    = 130  ! spectral fields

!--- i/o file names
character (256) :: sim_namelist = "sim_namelist"

!--- predefined experimets
character (256) :: sim = "iv00"  ! type of predefined simulation
                                 ! options are:

                                 !----------------------------------!
                                 ! Initial Value Problems (decaying !
                                 ! flows)                           !
                                 !----------------------------------!
           
                                 ! "iv00"   top hat jet

                                 !-----------------------! 
                                 ! Forced decaying Flows !
                                 !-----------------------! 
          
                                 ! "fd00"   top hat wind forcing 



!--- parameters of iv01 (initial top hat jet)
integer :: iv00w1    = 8      ! half width of jet center (in grid points)
integer :: iv00w2    = 4      ! width of vortex sheet (in grid points)
integer :: iv00scl   = 1      ! horizontal scale of jet
real(8) :: iv00qmax  = 1.0    ! amplitude of vortex sheets

!--- parameters of fd01 (forced decaying top hat jet)
integer :: fd00w1    = 8      ! half width of jet center (in grid points)
integer :: fd00w2    = 4      ! width of vortex sheet (in grid points)
integer :: fd00scl   = 1      ! horizontal scale of jet
real(8) :: fd00qmax  = 1.0d-4 ! amplitude of vortex sheets

end module simmod


! ***********************
! * SUBROUTINE SIMSTART *
! ***********************
subroutine simstart
use simmod
implicit none

call sim_readnl
call sim_cases

return
end subroutine simstart


! *************************
! * SUBROUTINE SIM_READNL *
! *************************
subroutine sim_readnl
use simmod
implicit none

!--- define sim_namelist
namelist /sim_nl/ sim        ,                                 &
                  iv00qmax   ,iv00w1     ,iv00w2     ,iv00scl, & 
                  fd00qmax   ,fd00w1     ,fd00w2     ,fd00scl   

!--- check if sim_namelist is present
inquire(file=sim_namelist,exist=lsimnl)

!--- read sim_nl
if (lsimnl) then
  open(nusimnl,file=sim_namelist,iostat=ios_simnl)
  read(nusimnl,sim_nl)
else
  return
endif


return
end subroutine sim_readnl


! ************************
! * SUBROUTINE SIM_CASES *
! ************************
subroutine sim_cases
use simmod
implicit none

integer :: jx
integer :: jy

real(8) :: gpvar(1:ngx,1:ngy)

select case(sim)
  !------------------------!
  ! initial value problems !
  !------------------------!

  !--- top-hat jet
  case("iv00")
     gpvar(:,:) = 0.0
     do jy = 1, ngy
        if ( jy .ge. ngy/2+1-iv00scl*(iv00w1+iv00w2) .and. & 
           jy .le. ngy/2-iv00scl*iv00w1 ) then
           gpvar(:,jy) = -iv00qmax
        endif
        if ( jy .ge. ngy/2+1+iv00scl*iv00w1 .and. & 
           jy .le. ngy/2+iv00scl*(iv00w1+iv00w2) ) then
           gpvar(:,jy) =  iv00qmax
        endif
     enddo
     call sim_wrtgp(gpvar,qcde,1)
  !-------------------!
  ! forced turbulence !
  !-------------------!

  !--- top-hat jet
  case("fd00")
     gpvar(:,:) = 0.0
     do jy = 1, ngy
        if ( jy .ge. ngy/2+1-fd00scl*(fd00w1+fd00w2) .and. & 
           jy .le. ngy/2-fd00scl*fd00w1 ) then
           gpvar(:,jy) = -fd00qmax
        endif
        if ( jy .ge. ngy/2+1+fd00scl*fd00w1 .and. & 
           jy .le. ngy/2+fd00scl*(fd00w1+fd00w2) ) then
           gpvar(:,jy) =  fd00qmax
        endif
     enddo
     call sim_wrtgp(gpvar,qfrccde,1)

  case default
end select

return
end subroutine sim_cases


! ************************
! * SUBROUTINE SIM_WRTGP *
! ************************
subroutine sim_wrtgp(gpvar,kcode,klev)
use simmod
implicit none

character(2), parameter :: gtp = "GP"
character(256)          :: fname
integer                 :: kcode,klev
integer                 :: yy,mo,dd,hh,mm,ii
real(8)                 :: gpvar(ngx,ngy)



!--- Build a header for service format
ihead(1) = kcode
ihead(2) = klev
ihead(3) = 10101
ihead(4) = 0 
ihead(5) = ngx
ihead(6) = ngy
ihead(7) = 0
ihead(8) = 0

!--- open file 
call sim_fname(gtp,kcode,fname)
open(nutmp,file=fname,form='unformatted')


!--- write ouput
write (nutmp) ihead
write (nutmp) gpvar(:,:)

close(nutmp)
return
end subroutine sim_wrtgp


! ************************
! * SUBROUTINE SIM_WRTSP *
! ************************
subroutine sim_wrtsp(spvar,kcode,klev)
use simmod
implicit none

character(2), parameter :: gtp = "SP"
character(256)          :: fname
integer :: kcode,klev
integer :: yy,mo,dd,hh,mm,ii
real(8) :: spvar(0:nfx,0:nfy)

ihead(1) = kcode
ihead(2) = klev
ihead(3) = 10101
ihead(4) = hh
ihead(5) = nfx+1
ihead(6) = nfy+1
ihead(7) = 0
ihead(8) = 0

!--- open file 
call sim_fname(gtp,kcode,fname)
open(nutmp,file=fname,form='unformatted')

!--- write data
write (nutmp) ihead
write (nutmp) spvar(:,:)

close(nutmp)
return
end subroutine sim_wrtsp


! **********************
! * SUBROUTINE SIMSTEP *
! **********************
subroutine simstep
use simmod
implicit none

return
end subroutine simstep


! **********************
! * SUBROUTINE SIMSTOP *
! **********************
subroutine simstop
use simmod
implicit none

return
end subroutine simstop


! ************************
! * SUBROUTINE SIM_FNAME *
! ************************
subroutine sim_fname(gtp,kcode,fname)
use simmod
implicit none

character(2)   :: gtp
character(256) :: fname
integer :: kcode

if (ngx < 100) then
  write(fname,'(a2,i2.2,"_var",i4.4,".srv")') gtp,ngx,kcode
elseif (ngx < 1000) then
  write(fname,'(a2,i3.3,"_var",i4.4,".srv")') gtp,ngx,kcode
else
  write(fname,'(a2,i4.4,"_var",i4.4,".srv")') gtp,ngx,kcode
endif


return
end subroutine sim_fname
