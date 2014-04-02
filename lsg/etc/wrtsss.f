!     ==================================================================
!     ------------------------------------------------------------------
!
      subroutine wrtsss(pfield,mon)
      use lsgvar
      implicit none
!
!     ------------------------------------------------------------------
!
      integer :: mon,kpar,ktape,inum,kkk
      real (kind=8) :: anum
      real (kind=8) :: pfield(ienjen)
      real (kind=8) :: zhead(20),zprel(6)
!
      open (69,file="sssforo",form="formatted")
      rewind 69
!
      anum=10.
      kpar=-5
      ktape=71
      zhead(1)=real(kpar)     ! field code of the variable
      zhead(2)=-100.
      zhead(3)=real(ktape)
      zhead(4)=real(ien)
      zhead(5)=real(jen)
      zhead(6)=1.
      zhead(7)=0.3750000000e+01
      zhead(8)=0.9375000000e+02
      zhead(9)=0.2111030000e+06
      zhead(10)=0.8806140000e+06
!
      write (69,7000) anum
!
      inum=nint(anum)
      inum=min0(inum,20)
      inum=max0(inum,8)
!
!     Header.
!
      do kkk=1,inum
        write (69,7000) zhead(kkk)
      end do
!
!
!*    2.        Write data
!     ----------
!
      zprel(1)=real(kpar)     ! field code of the variable
      zprel(2)=real(mon)
      zprel(2)=1.
      zprel(3)=0.
      zprel(4)=0.
      zprel(5)=0.
      zprel(6)=0.
      write (69,7050) (zprel(kkk),kkk=1,6)
 7050 format (6e12.4)
      write (69,7100) (pfield(kkk),kkk=1,ienjen)
 7100 format (4e20.10)
      close (69)
 7000 format (e20.10)

      end subroutine wrtsss
