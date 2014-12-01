! ************************************************************************
! * module legsym - direct and indirect Legendre transformation routines *
! * using symmetric and antisymmetric fourier coefficients               *
! * E. Kirk 08-Sep-2010 tested for T21 - T682 resolutions 32 & 64 bit    *
! ************************************************************************

module legsym
! ************************
! * Legendre Polynomials *
! ************************

integer :: ntru
integer :: ntp1
integer :: ncsp
integer :: nesp
integer :: nlon
integer :: nlpp
integer :: nhpp
integer :: nlat
integer :: nlev
real    :: plavor

real   , allocatable :: qi(:,:) ! P(m,n) = Associated Legendre Polynomials
real   , allocatable :: qj(:,:) ! Q(m,n) = Used for d/d(mu)
real   , allocatable :: qc(:,:) ! P(m,n) * gwd              used in fc2sp
real   , allocatable :: qe(:,:) ! Q(mn,) * gwd / cos2       used in mktend
real   , allocatable :: qq(:,:) ! P(m,n) * gwd / cos2 * n * (n+1) / 2  "
real   , allocatable :: qu(:,:) ! P(m,n) / (n*(n+1)) * m    used in dv2uv
real   , allocatable :: qv(:,:) ! Q(m,n) / (n*(n+1))        used in dv2uv
real   , allocatable :: qm(:,:) ! P(m,n) * gwd / cos2 * m   used in mktend
end module legsym

! =================
! SUBROUTINE LEGINI
! =================

subroutine legini(klat,klpp,kesp,klev,vorpla,sid,gwd)
use legsym
implicit none

integer :: klat
integer :: klpp
integer :: kesp
integer :: klev

real          :: vorpla   ! planetary vorticity
real (kind=8) :: sid(*)   ! sin(phi)
real (kind=8) :: gwd(*)   ! Gaussian weight (phi)

integer :: jlat ! Latitude
integer :: lm
integer :: m
integer :: n

real (kind=8) :: amsq
real (kind=8) :: z1
real (kind=8) :: z2
real (kind=8) :: z3
real (kind=8) :: f1m
real (kind=8) :: f2m
real (kind=8) :: znn1
real (kind=8) :: zsin    ! sin
real (kind=8) :: zcsq    ! cos2
real (kind=8) :: zcos    ! cos
real (kind=8) :: zgwd    ! gw
real (kind=8) :: zgwdcsq ! gw / cos2
real (kind=8) :: zpli(kesp)
real (kind=8) :: zpld(kesp)

nlat = klat
nlpp = klpp
nhpp = klpp / 2
nlon = klat + klat
ntru = (nlon - 1) /3
ntp1 = ntru + 1
ncsp = ((ntru + 1) * (ntru + 2)) / 2
nesp = kesp
nlev = klev

plavor = vorpla

call legpar(ntp1,nlon,nhpp,nesp/2-1,plavor)

allocate(qi(ncsp,nhpp))
allocate(qj(ncsp,nhpp))
allocate(qc(ncsp,nhpp))
allocate(qe(ncsp,nhpp))
allocate(qm(ncsp,nhpp))
allocate(qq(ncsp,nhpp))
allocate(qu(ncsp,nhpp))
allocate(qv(ncsp,nhpp))

call legpola(qi,qj,qc,qe) ! copy address to assembler routines
call legpolb(qm,qq,qu,qv)

do jlat = 1 , nhpp

! set p(0,0) and p(0,1)

   zgwd    = gwd(jlat)            ! gaussian weight - from inigau
   zsin    = sid(jlat)            ! sin(phi) - from inigau
   zcsq    = 1.0_8 - zsin * zsin  ! cos(phi) squared
   zgwdcsq = zgwd / zcsq          ! weight / cos squared
   zcos    = sqrt(zcsq)           ! cos(phi)
   f1m     = sqrt(1.5_8)
   zpli(1) = sqrt(0.5_8)
   zpli(2) = f1m * zsin
   zpld(1) = 0.0
   lm      = 2

! loop over wavenumbers

   do m = 0 , ntru
      if (m > 0) then
         lm  = lm + 1
         f2m = -f1m * sqrt(zcsq / (m+m))
         f1m =  f2m * sqrt(m+m + 3.0_8)
         zpli(lm) = f2m
         if (lm < ncsp) then
            lm = lm + 1
            zpli(lm  ) =       f1m * zsin
            zpld(lm-1) =  -m * f2m * zsin
         endif ! (lm < ncsp)
      endif ! (m > 0)

      amsq = m * m

      do n = m+2 , ntru
         lm = lm + 1
         z1 = sqrt(((n-1)*(n-1) - amsq) / (4*(n-1)*(n-1)-1))
         z2 = zsin * zpli(lm-1) - z1 * zpli(lm-2)
         zpli(lm  ) = z2 * sqrt((4*n*n-1) / (n*n-amsq))
         zpld(lm-1) = (1-n) * z2 + n * z1 * zpli(lm-2)
      enddo ! n

      if (lm < ncsp) then ! mode (m,ntru)
         z3 = sqrt((ntru*ntru-amsq) / (4*ntru*ntru-1))
         zpld(lm)=-ntru*zsin*zpli(lm) + (ntru+ntru+1)*zpli(lm-1)*z3
      else                ! mode (ntru,ntru)
         zpld(lm)=-ntru*zsin*zpli(lm)
      endif
   enddo ! m

! Part 1: m=0 and n=0 (lm=1)

   qi(1,jlat) = zpli(1)
   qj(1,jlat) = 0.0
   qc(1,jlat) = zpli(1) * zgwd
   qu(1,jlat) = 0.0
   qv(1,jlat) = 0.0
   qe(1,jlat) = 0.0
   qq(1,jlat) = 0.0
   qm(1,jlat) = 0.0

! Part 2: m=0 and n>0

   do lm = 2 , ntp1 ! index for spectral mode lm = n + 1
      znn1 = 1.0_8 / ((lm-1) * lm)
      qi(lm,jlat) = zpli(lm)
      qj(lm,jlat) = zpld(lm)
      qc(lm,jlat) = zpli(lm) * zgwd
      qu(lm,jlat) = 0.0
      qv(lm,jlat) = zpld(lm) * znn1
      qe(lm,jlat) = zpld(lm) * zgwdcsq
      qq(lm,jlat) = zpli(lm) * zgwdcsq * (lm-1) * lm * 0.5_8
      qm(lm,jlat) = 0.0
   enddo ! lm

! Part 3: m>0 and n>0

   lm = ntp1
   do m = 1 , ntru
      do n = m , ntru
           lm = lm + 1
           znn1 = 1.0_8 / (n * (n+1))
           qi(lm,jlat) = zpli(lm)
           qj(lm,jlat) = zpld(lm)
           qc(lm,jlat) = zpli(lm) * zgwd
           qu(lm,jlat) = zpli(lm) * znn1 * m
           qv(lm,jlat) = zpld(lm) * znn1
           qe(lm,jlat) = zpld(lm) * zgwdcsq
           qq(lm,jlat) = zpli(lm) * zgwdcsq * n * (n+1) * 0.5_8
           qm(lm,jlat) = zpli(lm) * zgwdcsq * m
      enddo ! n
   enddo ! m
enddo ! jlat
return
end


! ==================
! SUBROUTINE REG2ALT
! ==================

subroutine reg2alt(pr,klev)
use legsym
implicit none

real :: pr(nlon,nlat,klev)
real :: pa(nlon,nlat,klev)

integer :: jlat
integer :: klev

do jlat = 1 , nlat / 2
  pa(:,2*jlat-1,:) = pr(:,jlat       ,:)
  pa(:,2*jlat  ,:) = pr(:,nlat-jlat+1,:)
enddo

pr = pa

return
end


! ==================
! SUBROUTINE ALT2REG
! ==================

subroutine alt2reg(pa,klev)
use legsym
implicit none

real :: pa(nlon,nlat,klev)
real :: pr(nlon,nlat,klev)

integer :: jlat
integer :: klev

do jlat = 1 , nlat / 2
  pr(:,jlat       ,:) = pa(:,2*jlat-1,:)
  pr(:,nlat-jlat+1,:) = pa(:,2*jlat  ,:)
enddo

pa = pr

return
end


! ================
! SUBROUTINE ALTCS
! ================

subroutine altcs(pcs)
use legsym
implicit none
real :: pcs(nlat,nlev)
real :: pal(nlat,nlev)
integer :: jlat

do jlat = 1 , nlat / 2
  pal(jlat       ,:) = pcs(2*jlat-1,:)
  pal(nlat-jlat+1,:) = pcs(2*jlat  ,:)
enddo

pcs = pal

return
end


! =================
! SUBROUTINE ALTLAT 
! =================

subroutine altlat(pr,klat)
implicit none
integer :: jlat
integer :: klat
real    :: pr(klat)  ! regular     grid
real    :: pa(klat)  ! alternating grid

do jlat = 1 , klat / 2
  pa(2*jlat-1) = pr(jlat       )
  pa(2*jlat  ) = pr(klat-jlat+1)
enddo

pr(:) = pa(:)

return
end

