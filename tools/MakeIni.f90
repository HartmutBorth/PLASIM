program makeinit
implicit none

integer :: nout  = 20
integer :: kcode = 138

integer :: jx
integer :: jy
integer :: nx
integer :: ny

integer :: j1     = 20
integer :: j2     = 24
integer :: j3     = 39
integer :: j4     = 43
integer :: rr     = 1


real(8) :: pi
real(8) :: dx
real(8) :: dy 
real(8) :: amp   = 16.0d0

real(8), allocatable :: z(:,:)

character( 80) :: arg
character(256) :: fname


call get_command_argument(1,arg)
read(arg,*) nx
ny = nx

write (*,'("Allocate z(",i4,",",i4,")")') nx,ny

allocate(z(nx,ny))

pi = 4.0 * atan(1.0d0)
dx = 2.0d0 * pi / nx
dy = dx * 2.0d0

rr = nx/64

do jy = 1 , ny
   do jx = 1 , nx
!      z(jx,jy) = &
!               + 1.1 * cos((jx-1) * dx) &
!               + 0.1 * sin((jy-1) * dx)

!     z(jx,jy) = 0.0d0 
!
!     z(jx,jy) =  cos((jx-1) * dx +(jy-1) * dx)
!     z(jx,jy) =  1.0 * cos((jx-1)*dx)   +  3.0 * cos((jx-1)*dx*2.0) &
!               + 10.0 * sin((jy-1) * dy) + 20.0 * sin((jy-1) * dy * 2.0)



!     z(jx,jy) =  0.0d0
!     if (jx > rr*j1-(rr-1) .and. jx < rr*j2-(rr-1)) then
!        z(jx,jy) =  amp 
!     endif
!     if (jx > rr*j3-(rr-1) .and. jx < rr*j4-(rr-1)) then
!        z(jx,jy) = -amp
!     endif
!     z(jx,jy) = z(jx,jy)*(1 + 0.1*sin((jy-1)*dx))


      z(jx,jy) =  0.0d0
      if (jy > rr*j1-(rr-1) .and. jy < rr*j2-(rr-1)) then
         z(jx,jy) = -amp 
      endif
      if (jy > rr*j3-(rr-1) .and. jy < rr*j4-(rr-1)) then
         z(jx,jy) =  amp
      endif
      z(jx,jy) = z(jx,jy)*(1 + 0.1*sin((jx-1)*dx))


   enddo
enddo

!    z(nx/2,ny/2) = 2.0

call mk_fname(nx,kcode,fname,"GP")
open(nout,file=fname,form='unformatted')
write(nout) kcode,0,10101,0,nx,ny,0,0
write(nout) z
close(nout)
stop
end
      

! *###########################*
! * Functions and Subroutines *
! *###########################*


! ***********************
! * SUBROUTINE MK_FNAME *
! ***********************
subroutine mk_fname(kgx,kcode,fname,gtp)

character(2)   :: gtp
character(256) :: fname
integer :: kcode,kgx

if (kgx < 100) then
  write(fname,'(A2,I2.2,"_var",I4.4,".srv")') gtp,kgx,kcode
elseif (kgx < 1000) then
  write(fname,'(A2,I3.3,"_var",I4.4,".srv")') gtp,kgx,kcode
else
  write(fname,'(A2,I4.4,"_var",I4.4,".srv")') gtp,kgx,kcode
endif

fname = trim(fname)

return
end subroutine mk_fname
