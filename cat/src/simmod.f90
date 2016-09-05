module simmod
use catmod

!--- flags and switches
logical :: lsimnl = .false.    ! .true. if sim_namelist is present in run 
                               ! directory

integer :: ios_simnl = 1       ! 0 if <sim_namelist> is readable

!--- i/o units
integer, parameter :: nusimnl    = 110  ! sim namelist

!--- i/o file names
character (256) :: sim_namelist = "sim_namelist"

!--- temporary gridpoint fields (GP) 
real(8), allocatable :: sim_gtmp1(:,:)   ! temporary gridpoint field
real(8), allocatable :: sim_gtmp2(:,:)   ! temporary gridpoint field

!--- temporary fourier fields in real representation (F)
real(8), allocatable :: sim_ftmp1(:,:)   ! temporary gridpoint field
real(8), allocatable :: sim_ftmp2(:,:)   ! temporary gridpoint field

!--- predefined experimets
character (256) :: sim = "jet01" ! type of predefined simulation
                                 ! options are:
                                 ! "jet01"   top hat jets
                                 ! "jet02"   fourier jets
                                 ! "vor01"   gaussian vortices
                                 ! "dec01"   decaying turbulence
                                 ! "for01"   forced decaying turbulence


!--- parameters of jet01 (top hat jets)
real(8) :: jet01_amp    = 1.0    ! amplitude of vortex sheets
integer :: jet01_scl    = 1      ! scale
integer :: jet01_w1     = 8      ! half width of jet center (in grid points)
integer :: jet01_w2     = 4      ! width of vortex sheet (in grid points)

!--- parameters of jet02 (fourier jets)



!--- parameters of vor01 (gaussian vortices)



!--- parameters of dec01 (decaying turbulence)



!--- parameters of for01 (forced decaying turbulence)



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
namelist /sim_nl/ sim        ,                                      &
                  jet01_amp  ,jet01_scl  ,jet01_w1  ,jet01_w2


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


select case (sim)
  case ("jet01")
    allocate(sim_gtmp1(1:ngx,1:ngy)) ; sim_gtmp1(:,:) = 0.0
    do jy = 1, ngy
    enddo
    deallocate(sim_gtmp1)
  case ("vor01")
  case ("dec01")
  case ("for01")
  case default
end select

return
end subroutine sim_cases


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
subroutine sim_fname(kcode,gtp,fname)
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

fname = trim(fname)

return
end subroutine sim_fname


! ***********************
! * SUBROUTINE SIM_OPEN *
! ***********************
subroutine sim_open
use simmod
implicit none


return
end subroutine sim_open
