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
character (256) :: sim = "e00"   ! type of predefined simulation
                                 ! options are:
                                 ! "e00"   top hat jets
                                 ! "e01"   fourier jets


!--- parameters of e01 (top hat jets)
real(8) :: e00umax  = 1.0    ! amplitude of vortex sheets
integer :: e00w1    = 8      ! half width of jet center (in grid points)
integer :: e00w2    = 4      ! width of vortex sheet (in grid points)
integer :: e00scl   = 1      ! horizontal scale of jet

!--- parameters of e02 (fourier jets)



!--- parameters of e03 (gaussian vortices)


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
                  e00umax   ,e00w1     ,e00w2     ,e00scl    

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
  !--- top-hat jet
  case ("e00")
    allocate(sim_gtmp1(1:ngx,1:ngy))  
    sim_gtmp1(:,:) = 0.0
    do jy = 1, ngy
      if ( jy .ge. ngy/2+1-e00scl*(e00w1+e00w2) .and. & 
           jy .le. ngy/2-e00scl*e00w1 ) then
         sim_gtmp1(:,jy) =  e00umax
      endif
      if ( jy .ge. ngy/2+1+e00scl*e00w1 .and. & 
           jy .le. ngy/2+e00scl*(e00w1+e00w2) ) then
         sim_gtmp1(:,jy) = -e00umax  
      endif
    enddo
    deallocate(sim_gtmp1)
  !--- gaussian jet
  case ("e01")
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
