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

!--- parameters of e01 (top hat jets)
integer :: e00w1    = 8      ! half width of jet center (in grid points)
integer :: e00w2    = 4      ! width of vortex sheet (in grid points)
integer :: e00scl   = 1      ! horizontal scale of jet
real(8) :: e00umax  = 1.0    ! amplitude of vortex sheets

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
  case default
end select


deallocate(sim_gtmp1)
return
end subroutine sim_cases


! ************************
! * SUBROUTINE SIM_WRTGP *
! ************************
subroutine sim_wrtgp(ku,gpvar,kcode,klev)
use simmod
implicit none

integer :: ku,kcode,klev
integer :: yy,mo,dd,hh,mm,ii
real(8) :: gpvar(ngx,ngy)

! Build a header for service format

mm = mod(tstep,60)
ii = tstep / 60
hh = mod(ii,24)
ii = ii / 24
dd = mod(ii,30) + 1
ii = ii / 30
mo = mod(ii,12) + 1
yy = ii / 12

ihead(1) = kcode
ihead(2) = klev
ihead(3) = dd + 100 * mo + 10000 * yy
ihead(4) = mm + 100 * hh
ihead(5) = ngx
ihead(6) = ngy
ihead(7) = 0
ihead(8) = 0

write (ku) ihead
write (ku) gpvar(:,:)

return
end subroutine sim_wrtgp


! ************************
! * SUBROUTINE SIM_WRTSP *
! ************************
subroutine sim_wrtsp(ku,spvar,kcode,klev)
use simmod
implicit none

integer :: ku,kcode,klev
integer :: yy,mo,dd,hh,mm,ii
real(8) :: spvar(0:nfx,0:nfy)

mm = mod(tstep,60)
ii = tstep / 60
hh = mod(ii,24)
ii = ii / 24
dd = mod(ii,30) + 1
ii = ii / 30
mo = mod(ii,12) + 1
yy = ii / 12

ihead(1) = kcode
ihead(2) = klev
ihead(3) = dd + 100 * mo + 10000 * yy
ihead(4) = mm + 100 * hh
ihead(5) = nfx+1
ihead(6) = nfy+1
ihead(7) = 0
ihead(8) = 0

write (ku) ihead
write (ku) spvar(:,:)

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
