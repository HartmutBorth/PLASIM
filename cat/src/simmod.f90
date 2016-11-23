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
character (256) :: ysim = "djet01" ! type of predefined simulation
                                           ! options are:

                                 !----------------------------------!
                                 ! Initial Value Problems (decaying !
                                 ! flows)                           !
                                 !----------------------------------!
           
                                 ! "djet01"   decaying top hat jet

                                 !-----------------------! 
                                 ! Forced decaying Flows !
                                 !-----------------------! 
          
                                 ! "fjet01"   forced top hat jet 


!--- parameters of djet01 (initial top hat jet)
integer :: w1     = 8      ! half width of jet center (in grid points)
integer :: w2     = 4      ! width of vortex sheet (in grid points)
integer :: hscl   = 1      ! horizontal scale of jet
real(8) :: qmax   = 1.0    ! amplitude of vortex sheets


!--- sim_namelist parameters to overwrite cat_namelist parameters

!---------------------------------------------------------------!
! A sim_namelist parameter sXXXXX corresponds to a cat_namelist !
! parameter XXXXX. The sim_namelist parameters sXXXXX can only  !
! take the same values as the corresponding cat_namelist        !
! parameter.                                                    !
!---------------------------------------------------------------!

real(8) :: sdt     = -1.0 ! length of time step [s]

integer :: snsteps = -1   ! total number of time steps
integer :: sngui   = -1   ! time steps between GUI calls
integer :: snforc  = -1   ! type of forcing
integer :: snpert  = -1   ! type of perturbation
integer :: snpost  = -1   ! type of post processing

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
namelist /sim_nl/ ysim       ,                        &
                  qmax       ,w1      ,w2     ,hscl  ,& 
                  sdt        ,snsteps ,sngui  ,       &
                  snpert     ,snforc  ,snpost

!--- check if sim_namelist is present
inquire(file=sim_namelist,exist=lsimnl)

!--- read sim_nl
if (lsimnl) then
  open(nusimnl,file=sim_namelist,iostat=ios_simnl)
  read(nusimnl,sim_nl)
else
  return
endif

!--- overwrite cat_namelist parameters
if (sdt     .ge. 0)     dt = sdt
if (snsteps .ge. 0) nsteps = snsteps
if (sngui   .ge. 0)   ngui = sngui
if (snpert  .ge. 0)  npert = snpert
if (snforc  .ge. 0)  nforc = snforc
if (snpost  .ge. 0)  npost = snpost

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

select case(ysim)
  !------------------------!
  ! initial value problems !
  !------------------------!

  !--- top-hat jet
  case("djet01")
     gpvar(:,:) = 0.0
     do jy = 1, ngy
        if ( jy .ge. ngy/2+1-hscl*(w1+w2) .and. & 
           jy .le. ngy/2-hscl*w1 ) then
           gpvar(:,jy) = -qmax
        endif
        if ( jy .ge. ngy/2+1+hscl*w1 .and. & 
           jy .le. ngy/2+hscl*(w1+w2) ) then
           gpvar(:,jy) =  qmax
        endif
     enddo
     call sim_wrtgp(gpvar,qcde,1)

  !-------------------!
  ! forced turbulence !
  !-------------------!

  !--- top-hat jet
  case("fjet01")
     gpvar(:,:) = 0.0
     do jy = 1, ngy
        if ( jy .ge. ngy/2+1-hscl*(w1+w2) .and. & 
           jy .le. ngy/2-hscl*w1 ) then
           gpvar(:,jy) = -qmax
        endif
        if ( jy .ge. ngy/2+1+hscl*w1 .and. & 
           jy .le. ngy/2+hscl*(w1+w2) ) then
           gpvar(:,jy) =  qmax
        endif
     enddo
     call sim_wrtgp(gpvar,qfrccde,1)

  case("fjet02")
     gpvar(:,:) = 0.0
!    do jx = 1, ngx
!       if ( jx .ge. ngx/2+1-hscl*(w1+w2) .and. & 
!          jx .le. ngx/2-hscl*w1 ) then
!          gpvar(jx,:) = -qmax
!       endif
!       if ( jx .ge. ngx/2+1+hscl*w1 .and. & 
!          jx .le. ngx/2+hscl*(w1+w2) ) then
!          gpvar(jx,:) =  qmax
!       endif
!    enddo
     do jy = 1, ngy
        if ( jy .ge. ngy/2+1-hscl*(w1+w2) .and. & 
           jy .le. ngy/2-hscl*w1 ) then
           gpvar(:,jy) = -qmax
        endif
        if ( jy .ge. ngy/2+1+hscl*w1 .and. & 
           jy .le. ngy/2+hscl*(w1+w2) ) then
           gpvar(:,jy) =  qmax
        endif
     enddo
     gpvar = transpose(gpvar)
     call sim_wrtgp(gpvar,qfrccde,1)

  case("fjet03")
     gpvar(:,:) = 0.0
!     do jx = 1, ngx
!        do jy = 1,ngy
!        if ( jx .ge. ngx/2+1-hscl*(w1+w2) .and. & 
!           jx .le. ngx/2-hscl*w1 ) then
!           gpvar(jx,:) = -qmax
!        endif
!        if ( jx .ge. ngx/2+1+hscl*w1 .and. & 
!           jx .le. ngx/2+hscl*(w1+w2) ) then
!           gpvar(jx,:) =  qmax
!        endif
!        enddo
!     enddo
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
