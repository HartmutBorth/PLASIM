      module mpimod
      use pumamod
      use mpi

      integer :: mpi_itype = MPI_INTEGER4
      integer :: mpi_rtype = MPI_REAL4
      integer :: mpi_ltype = MPI_LOGICAL

      end module mpimod
!
!     interface routines to MPI:
!


!     ================
!     SUBROUTINE MPBCI
!     ================

      subroutine mpbci(k) ! broadcast 1 integer
      use mpimod

      call mpi_bcast(k,1,mpi_itype,NROOT,myworld,mpinfo)

      return
      end subroutine mpbci

!     =================
!     SUBROUTINE MPBCIN
!     =================

      subroutine mpbcin(k,n) ! broadcast n integer
      use mpimod

      integer k(n)

      call mpi_bcast(k,n,mpi_itype,NROOT,myworld,mpinfo)

      return
      end subroutine mpbcin

!     ================
!     SUBROUTINE MPBCR
!     ================

      subroutine mpbcr(p) ! broadcast 1 real
      use mpimod

      call mpi_bcast(p,1,mpi_rtype,NROOT,myworld,mpinfo)

      return
      end subroutine mpbcr

!     =================
!     SUBROUTINE MPBCRN
!     =================

      subroutine mpbcrn(p,n) ! broadcast n real
      use mpimod

      real p(n)

      call mpi_bcast(p,n,mpi_rtype,NROOT,myworld,mpinfo)

      return
      end subroutine mpbcrn

!     ================
!     SUBROUTINE MPBCL
!     ================

      subroutine mpbcl(k) ! broadcast 1 logical
      use mpimod
      logical k

      call mpi_bcast(k,1,mpi_ltype,NROOT,myworld,mpinfo)

      return
      end subroutine mpbcl

!     =================
!     SUBROUTINE MPSCIN
!     =================

      subroutine mpscin(k,n) ! scatter n integer
      use mpimod

      integer k(*)

      call mpi_scatter(k,n,mpi_itype,k,n,mpi_itype                  &
     &                ,NROOT,myworld,mpinfo)

      return
      end subroutine mpscin

!     =================
!     SUBROUTINE MPSCRN
!     =================

      subroutine mpscrn(p,n) ! scatter n real
      use mpimod

      real p(*)

      call mpi_scatter(p,n,mpi_rtype,p,n,mpi_rtype,NROOT,myworld,mpinfo)

      return
      end subroutine mpscrn

!     =================
!     SUBROUTINE MPSCDN
!     =================

      subroutine mpscdn(p,n) ! scatter n double precision
      use mpimod

      real (kind=8) :: p(*)

      call mpi_scatter(p,n,MPI_REAL8,p,n,MPI_REAL8,NROOT,myworld,mpinfo)

      return
      end subroutine mpscdn

!     =================
!     SUBROUTINE MPSCGP
!     =================

      subroutine mpscgp(pf,pp,klev) ! scatter gridpoint fields
      use mpimod

      real pf(NUGP,klev)
      real pp(NHOR,klev)

      do jlev = 1 , klev
         call mpi_scatter(pf(1,jlev),NHOR,mpi_rtype,                     &
     &                    pp(1,jlev),NHOR,mpi_rtype,                     &
     &                    NROOT,myworld,mpinfo)
      enddo

      return
      end subroutine mpscgp

!     =================
!     SUBROUTINE MPGAGP
!     =================

      subroutine mpgagp(pf,pp,klev) ! gather gridpoint fields
      use mpimod

      real pf(NLON*NLAT,klev)
      real pp(NHOR,klev)

      do jlev = 1 , klev
         call mpi_gather(pp(1,jlev),NHOR,mpi_rtype,                      &
     &                   pf(1,jlev),NHOR,mpi_rtype,                      &
     &                   NROOT,myworld,mpinfo)
      enddo

      return
      end subroutine mpgagp

!     ===================
!     SUBROUTINE MPGALLGP
!     ===================

      subroutine mpgallgp(pf,pp,klev) ! gather gritpoint to all
      use mpimod

      real pf(NLON*NLAT,klev)
      real pp(NHOR,klev)

      do jlev = 1 , klev
         call mpi_allgather(pp(1,jlev),NHOR,mpi_rtype,                   &
     &                      pf(1,jlev),NHOR,mpi_rtype,                   &
     &                      myworld,mpinfo)
      enddo

      return
      end subroutine mpgallgp

!     =================
!     SUBROUTINE MPSCSP
!     =================

      subroutine mpscsp(pf,pp,klev) ! scatter spectral fields
      use mpimod

      real pf(NESP,klev)
      real pp(NSPP,klev)

      do jlev = 1 , klev
         call mpi_scatter(pf(1,jlev),NSPP,mpi_rtype                      &
     &                   ,pp(1,jlev),NSPP,mpi_rtype                      &
     &                   ,NROOT,myworld,mpinfo)
      enddo

      return
      end subroutine mpscsp

!     =================
!     SUBROUTINE MPGASP
!     =================

      subroutine mpgasp(pf,pp,klev) ! gather spectral fields
      use mpimod

      real pf(NESP,klev)
      real pp(NSPP,klev)

      do jlev = 1 , klev
         call mpi_gather(pp(1,jlev),NSPP,mpi_rtype                       &
     &                  ,pf(1,jlev),NSPP,mpi_rtype                       &
     &                  ,NROOT,myworld,mpinfo)
      enddo

      return
      end subroutine mpgasp

!     =================
!     SUBROUTINE MPGACS
!     =================

      subroutine mpgacs(pcs) ! gather cross sections
      use mpimod

      real pcs(NLAT,NLEV)

      do jlev = 1 , NLEV
         call mpi_gather(pcs(1,jlev),NLPP,mpi_rtype                      &
     &                  ,pcs(1,jlev),NLPP,mpi_rtype                      &
     &                  ,NROOT,myworld,mpinfo)

      enddo
      return
      end subroutine mpgacs

!     ===================
!     SUBROUTINE MPGALLSP
!     ===================

      subroutine mpgallsp(pf,pp,klev) ! gather spectral to all
      use mpimod

      real pf(NESP,klev)
      real pp(NSPP,klev)

      do jlev = 1 , klev
         call mpi_allgather(pp(1,jlev),NSPP,mpi_rtype                    &
     &                     ,pf(1,jlev),NSPP,mpi_rtype                    &
     &                     ,myworld,mpinfo)
      enddo

      return
      end subroutine mpgallsp

!     ================
!     SUBROUTINE MPSUM
!     ================

      subroutine mpsum(psp,klev) ! sum spectral fields
      use mpimod

      real psp(NESP,klev)
      real tmp(NESP,klev)

      call mpi_reduce(psp(1,1),tmp(1,1),NESP*klev,mpi_rtype,MPI_SUM     &
     &               ,NROOT,myworld,mpinfo)
      if (mypid == NROOT) psp = tmp

      return
      end subroutine mpsum

!     ==================
!     SUBROUTINE MPSUMSC
!     ==================

      subroutine mpsumsc(psf,psp,klev) ! sum & scatter spectral
      use mpimod

      real psf(NESP,klev)
      real psp(NSPP,klev)

      do jlev = 1 , klev
         call mpi_reduce_scatter(psf(1,jlev),psp(1,jlev),nscatsp        &
     &                          ,mpi_rtype,MPI_SUM,myworld,mpinfo)
      enddo

      return
      end subroutine mpsumsc

!     ====================
!     SUBROUTINE MPSUMR
!     ====================

      subroutine mpsumr(pr,kdim) ! sum kdim reals
      use mpimod

      real pr(kdim)
      real tmp(kdim)

      call mpi_reduce(pr,tmp,kdim,mpi_rtype,MPI_SUM,NROOT,myworld,mpinfo)
      if (mypid == NROOT) pr = tmp

      return
      end subroutine mpsumr

!     ====================
!     SUBROUTINE MPSUMBCR
!     ====================

      subroutine mpsumbcr(pr,kdim) ! sum & broadcast kdim reals
      use mpimod

      real :: pr(kdim)
      real :: tmp(kdim)

      call mpi_allreduce(pr,tmp,kdim,mpi_rtype,MPI_SUM,myworld,mpinfo)
      pr = tmp

      return
      end subroutine mpsumbcr

!     ==================
!     SUBROUTINE MPSTART
!     ==================

      subroutine mpstart ! initialization
      use mpimod
      integer :: itest = 0
      real    :: rtest = 0.0

      if (kind(itest) == 8) mpi_itype = MPI_INTEGER8
      if (kind(rtest) == 8) mpi_rtype = MPI_REAL8

      call mpi_init(mpinfo)
      myworld=MPI_COMM_WORLD

      call mpi_comm_size(myworld,npro ,mpinfo)
      call mpi_comm_rank(myworld,mypid,mpinfo)
      allocate(ympname(npro)) ; ympname(:) = ' '
      call mpi_get_processor_name(ympname(1),ilen,mpinfo)

      call mpi_gather(ympname,80,MPI_CHARACTER,   &
                      ympname,80,MPI_CHARACTER,   &
                      NROOT,myworld,mpinfo)

      return
      end subroutine mpstart

!     =================
!     SUBROUTINE MPSTOP
!     =================

      subroutine mpstop
      use mpimod

      call mpi_barrier(myworld,mpinfo)
      call mpi_finalize(mpinfo)

      return
      end subroutine mpstop

!     ===================
!     SUBROUTINE MPREADGP
!     ===================

      subroutine mpreadgp(ktape,p,kdim,klev)
      use mpimod

      real p(kdim,klev)
      real z(NLON*NLAT,klev)

      z = 0.0

      if (mypid == NROOT) read (ktape) z(1:NLON*NLAT,:)

      if (kdim == NHOR) then
       call mpscgp(z,p,klev)
      else
       if (mypid == NROOT) p = z
      endif

      return
      end subroutine mpreadgp

!     ====================
!     SUBROUTINE MPWRITEGP
!     ====================

      subroutine mpwritegp(ktape,p,kdim,klev)
      use mpimod

      real p(kdim,klev)
      real z(NLON*NLAT,klev)

      if (kdim == NHOR) then
       call mpgagp(z,p,klev)
       if (mypid == NROOT) write(ktape) z(1:NLON*NLAT,:)
      else
       if (mypid == NROOT) write(ktape) p(1:NLON*NLAT,:)
      endif

      return
      end subroutine mpwritegp


!     =====================
!     SUBROUTINE MPWRITEGPH
!     =====================

      subroutine mpwritegph(ktape,p,kdim,klev,ihead)
      use mpimod

      real p(kdim,klev)
      real z(NLON*NLAT,klev)
!
      real(kind=4) :: zp(kdim,klev)
      real(kind=4) :: zz(NLON*NLAT,klev)
!

      integer ihead(8)

      if (kdim == NHOR) then
       call mpgagp(z,p,klev)
       if (mypid == NROOT) then
          write(ktape) ihead
          zz(:,:)=z(:,:)
          write(ktape) zz(1:NLON*NLAT,:)
       endif
      else
       if (mypid == NROOT) then
          write(ktape) ihead
          zp(:,:)=p(:,:)
          write(ktape) zp(1:NLON*NLAT,:)
       endif
      endif

      return
      end subroutine mpwritegph

!     ===================
!     SUBROUTINE MPREADSP
!     ===================

      subroutine mpreadsp(ktape,p,kdim,klev)
      use mpimod

      real p(kdim,klev)
      real z(NESP,klev)

      z = 0.0
      if (mypid == NROOT) read(ktape) z(1:NRSP,:)
      if (kdim == NSPP) then
         call mpscsp(z,p,klev)
      else
         if (mypid == NROOT) p = z
      endif

      return
      end subroutine mpreadsp

!     ====================
!     SUBROUTINE MPWRITESP
!     ====================

      subroutine mpwritesp(ktape,p,kdim,klev)
      use mpimod

      real p(kdim,klev)
      real z(NESP,klev)

      if (kdim == NSPP) then
         call mpgasp(z,p,klev)
         if (mypid == NROOT) write(ktape) z(1:NRSP,:)
      else
         if (mypid == NROOT) write(ktape) p(1:NRSP,:)
      endif

      return
      end subroutine mpwritesp


!     ===================
!     SUBROUTINE MPI_INFO
!     ===================

      subroutine mpi_info(nprocess,npid)    ! get nproc and pid
      use mpimod

      myworld=MPI_COMM_WORLD

      call mpi_comm_size(myworld,nprocess,mpinfo)
      call mpi_comm_rank(myworld,npid,mpinfo)

      return
      end subroutine mpi_info


!     ==================
!     SUBROUTINE MPGETSP
!     ==================

      subroutine mpgetsp(yn,p,kdim,klev)
      use mpimod

      character (len=*) :: yn
      real :: p(kdim,klev)
      real :: z(NESP,klev)

      z(:,:) = 0.0
      if (mypid == NROOT) call get_restart_array(yn,z,NRSP,NESP,klev)
      call mpscsp(z,p,klev)

      return
      end subroutine mpgetsp


!     ==================
!     SUBROUTINE MPGETGP
!     ==================

      subroutine mpgetgp(yn,p,kdim,klev)
      use mpimod

      character (len=*) :: yn
      real :: p(kdim,klev)
      real :: z(NUGP,klev)

      if (mypid == NROOT) call get_restart_array(yn,z,NUGP,NUGP,klev)
      call mpscgp(z,p,klev)

      return
      end subroutine mpgetgp


!     ==================
!     SUBROUTINE MPPUTSP
!     ==================

      subroutine mpputsp(yn,p,kdim,klev)
      use mpimod

      character (len=*) :: yn
      real :: p(kdim,klev)
      real :: z(NESP,klev)

      call mpgasp(z,p,klev)
      if (mypid == NROOT) call put_restart_array(yn,z,NRSP,NESP,klev)

      return
      end subroutine mpputsp


!     ==================
!     SUBROUTINE MPPUTGP
!     ==================

      subroutine mpputgp(yn,p,kdim,klev)
      use mpimod

      character (len=*) :: yn
      real :: p(kdim,klev)
      real :: z(NUGP,klev)

      call mpgagp(z,p,klev)
      if (mypid == NROOT) call put_restart_array(yn,z,NUGP,NUGP,klev)

      return
      end subroutine mpputgp


!!     ===================
!!     SUBROUTINE MPSURFGP
!!     ===================
!
!      subroutine mpsurfgp(yn,p,kdim,klev)
!      use mpimod
!
!      character (len=*) :: yn
!      real :: p(kdim,klev)
!      real :: z(NUGP,klev)
!
!      if (mypid == NROOT) call get_surf_array(yn,z,NUGP,NUGP,klev,iread)
!      call mpbci(iread)
!      if (iread == 1) call mpscgp(z,p,klev)
!
!      return
!      end subroutine mpsurfgp
!
!
!!     =====================
!!     SUBROUTINE MPSURFYEAR
!!     =====================
!
!      subroutine mpsurfyear(yn,p,kdim,kmon)
!      use mpimod
!
!      character (len=*) :: yn
!      real :: p(kdim,kmon)
!      real :: z(NUGP,kmon)
!
!      if (mypid == NROOT) call get_surf_year(yn,z,NUGP,kmon,iread)
!      call mpbci(iread)
!      if (iread == 1) call mpscgp(z,p,kmon)
!
!      return
!      end subroutine mpsurfyear
!
!
!!     =====================
!!     SUBROUTINE MP3DYEAR
!!     =====================
!
!      subroutine mp3dyear(yn,p,kdim,klev,kmon)
!      use mpimod
!
!      character (len=*) :: yn
!      real :: p(kdim,klev,kmon)
!      real :: z(NUGP,klev,kmon)
!
!      if (mypid == NROOT) call get_3d_year(yn,z,NUGP,klev,kmon,iread)
!      call mpbci(iread)
!      if (iread == 1) call mpscgp(z,p,klev*kmon)
!
!      return
!      end subroutine mp3dyear


!     ===================
!     SUBROUTINE MPMAXVAL
!     ===================

      subroutine mpmaxval(p,kdim,klev,pmax)
      use mpimod

      real :: p(kdim,klev)

      zmax = maxval(p(:,:))
      call mpi_allreduce(zmax,pmax,1,mpi_rtype,MPI_MAX,myworld,mpinfo)

      return
      end subroutine mpmaxval


!     ===================
!     SUBROUTINE MPSUMVAL
!     ===================

      subroutine mpsumval(p,kdim,klev,psum)
      use mpimod

      real :: p(kdim,klev)

      zsum = sum(p(:,:))
      call mpi_allreduce(zsum,psum,1,mpi_rtype,MPI_SUM,myworld,mpinfo)

      return
      end subroutine mpsumval

!     Some dummy declarations that are use in multirun mode only

      subroutine mrdiff(p,d,n)
      real :: p(n)
      real :: d(n)
      return
      end

      subroutine mrsum(k) ! sum up 1 integer
      return
      end

      subroutine mrbci(k) ! broadcast 1 integer
      return
      end

      subroutine mrdimensions
      return
      end 


