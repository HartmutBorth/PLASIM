! *****************
! * MODULE FFTFOR *
! *****************

module fftfor
integer, parameter :: NRES = 10
integer :: nallowed(NRES) = [16,32,64,128,256,512,1024,2048,4096,8192]

! T3    - N16   : 8-2
! T10   - N32   : 8-2-2
! T21   - N64   : 8-4-2
! T42   - N128  : 8-4-4
! T85   - N256  : 8-4-4-2
! T170  - N512  : 8-4-4-4
! T341  - N1024 : 8-4-4-4-2
! T682  - N2048 : 8-4-4-4-4
! T1365 - N4096 : 8-4-4-4-4-2
! T2730 - N8192 : 8-4-4-4-4-4

integer :: lastn = 0 ! last used value for n
integer :: lsize = 0 ! last used size for w
integer :: nbit  = 0 ! bit # set for value of n
real(8),allocatable :: trigs(:)
real(8),allocatable :: w(:)
end module fftfor


! *******************
! * fourier_to_grid *
! *******************

subroutine fourier_to_grid(fc,gp,nfx,nfy,ngx,ngy)
implicit none
integer,intent(in ) :: nfx,nfy,ngx,ngy
real(8),intent(in ) :: fc(*)
real(8),intent(out) :: gp(*)

call fast_ftp(fc,gp,ngx)  ! matrix transposition & reordering
call fc2gp(gp,ngy,nfx+2)  ! 1-dim indirect fourier transformation
call fast_mtp(gp,ngx)     ! matrix transposition
call fc2gp(gp,ngx,ngy)    ! 1-dim indirect fourier transformation

return
end


! *******************
! * grid_to_fourier *
! *******************

subroutine grid_to_fourier(gp,fc,nfx,nfy,ngx,ngy)
implicit none
integer,intent(in ) :: nfx,nfy,ngx,ngy
real(8),intent(in ) :: gp(*)
real(8),intent(out) :: fc(*)

call gp2fc(gp,ngx,ngy)    ! 1-dim direct fourier transformation
call fast_mtp(gp,ngx)     ! matrix transposition
call gp2fc(gp,ngy,nfx+2)  ! 1-dim direct fourier transformation
call fast_gtp(gp,fc,ngx)  ! matrix transposition & reordering
 
return
end

! *********
! * fc2gp *
! *********

subroutine fc2gp(a,n,lot)
use fftfor
implicit none
real(8) :: a(*)
integer :: n,lot,j,la
logical :: l2,ld

if (n /= lastn .or. n * lot > lsize) call fftini(n,lot)
l2 = mod(nbit,2) == 0  ! start with ifft2s or ifft4s ?
ld = .true.            ! start ifft4m with w -> a

if (l2) then
   call ifft2s(a,w,trigs,n,lot)
   la = 2
else
   call ifft4s(a,w,trigs,n,lot)
   la = 4
endif

do j = 1 , nbit / 2
   if (ld) then
      call ifft4m(w,a,trigs,n,la,lot)
   else
      call ifft4m(a,w,trigs,n,la,lot)
   endif
   la = la * 4
   ld = .not. ld
enddo

if (ld) then
   call ifft8e(w,a,n,lot)
else
   call ifft8e(a,a,n,lot)
endif

return
end subroutine fc2gp


! *********
! * gp2fc *
! *********

subroutine gp2fc(a,n,lot)
use fftfor
implicit none
real(8) :: a(*)
integer :: n
integer :: lot
integer :: j
integer :: la
integer :: m4
logical :: ld 

if (n /= lastn .or. n * lot > lsize) call fftini(n,lot)

m4 = nbit / 2 ! # of calls to dfft4m

ld = mod(m4,2) == 0

if (ld) then
   call dfft8s(a,w,n,lot)
else
   call dfft8s(a,a,n,lot)
endif

la = n / 8
do j = 1 , m4
   la = la / 4
   if (ld) then
      call dfft4m(w,a,trigs,n,lot,la)
   else
      call dfft4m(a,w,trigs,n,lot,la)
   endif
   ld = .not. ld
enddo

if (la == 2) then
   call dfft2e(w,a,trigs,n,lot)
else
   call dfft4e(w,a,trigs,n,lot)
endif
return
end subroutine gp2fc


! **********
! * fftini *
! **********

subroutine fftini(n,lot)
use fftfor
implicit none
integer :: n,lot
logical labort

integer :: j,k,ibit
real(8) :: angle,del

! check for allowed values of n

labort = .true.
do j = 1 , NRES
   if (n == nallowed(j)) labort = .false.
enddo

if (labort) then
   write (*,*) '*** FFT does not support n = ',n,' ***'
   write (*,*) 'Following resolutions may be used:'
   write (*,*) '----------------------------------'
   do j = 1 , NRES
      write (*,1000) nallowed(j), nallowed(j)/2, nallowed(j)/3
   enddo
   stop
endif
 1000 format(' NLON=',I5,'  NLAT=',I5,'  NTRU=',I5)

! compute the bit # that is set in the integer n (power of 2)

ibit = n / 32 ! factors 8 and 4 set alreay
nbit = 0
do while (ibit > 0)
   ibit = ibit / 2
   nbit = nbit + 1
enddo  

! write (*,*) 'n = ',n,' nbit = ',nbit

! allocate trigs for new length n

if (allocated(trigs)) deallocate(trigs)
allocate(trigs(n))
lastn = n

! allocate w for new length n * lot

if (lot > 0) then
   lsize = n * lot
   if (allocated(w)) deallocate(w)
   allocate(w(lsize))
   ! write (*,*) 'allocate w',n,lot,lsize
endif

lastn = n

! compute values for trigs

del = 4.0 * asin(1.0) / n
do k=0,n/2-1
  angle = k * del
  trigs(2*k+1) = cos(angle)
  trigs(2*k+2) = sin(angle)
enddo
return
end subroutine fftini


