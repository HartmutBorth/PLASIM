!     ==================================================================
!     ------------------------------------------------------------------
!
      subroutine wrtsst(pfield,mon,nye)
      use lsgvar
      implicit none
!
!     ------------------------------------------------------------------


      integer :: mon,nye,kpar,ktape,inum,itim,kkk
      real (kind=8) :: anum
      real (kind=8) :: pfield(ienjen)
      real (kind=8) :: zhead(20),zprel(6)

      character(len=12) monthsst
!
!   ??? Attention - this is only valid for ntyear=12
      if (ntyear.ne.12) then
        write(*,*) ' This form of routine WRTSST runs is not valid!'
        write(*,*) '   - only suitable with ntyear=12'
        write(*,*) '     ntyear=',ntyear
        return
      endif
!
!
!     SST file for input to planet simulator, filename with year/month
      write (monthsst,"(a5,i5.5,i2.2)") "bsst_",nye,mon
!
      open (68,file=monthsst,form="formatted")
      rewind 68
!
      anum=10.
      kpar=-2
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

      anum=11.
      zhead(3)=12.
      zhead(9)=204065.
      call actdate(idat,itim)
      zhead(10)=real(idat)
      zhead(11)=real(itim)
      write(*,'(3a)') ' WRTSST: SST data written on file ',monthsst
!
      write (68,*) anum
!
      inum=nint(anum)
      inum=min0(inum,20)
      inum=max0(inum,8)
!
!     Header.
!
      write (68,'(6f6.0,2f6.3,3f8.0)') (zhead(kkk),kkk=1,inum)
!
!
!*    2.        Write data
!     ----------
!
      zprel(1)=real(kpar)     ! field code of the variable
      zprel(2)=real(mon)
      zprel(3)=0.
      zprel(4)=0.
      zprel(5)=0.
      zprel(6)=0.
      write (68,7050) (zprel(kkk),kkk=1,6)
 7050 format (6e12.4)
      write (68,7100) (pfield(kkk),kkk=1,ienjen)
 7100 format (4e20.10)
      close (68)
 7000 format (e20.10)
      return
      end subroutine wrtsst
