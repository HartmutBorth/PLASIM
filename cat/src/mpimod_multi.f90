      module mpimod
      use catmod
      use mpi

      integer :: mpi_itype = MPI_INTEGER4
      integer :: mpi_rtype = MPI_REAL4
      integer :: mpi_ltype = MPI_LOGICAL

      end module mpimod


!     ==================
!     SUBROUTINE MPSTART
!     ==================

      subroutine mpstart ! initialization
      use mpimod
      integer :: itest = 0
      real    :: rtest = 0.0

      if (kind(itest) == 8) mpi_itype = MPI_INTEGER8
      if (kind(rtest) == 8) mpi_rtype = MPI_REAL8

      call mpi_init(mrinfo)
      mrworld=MPI_COMM_WORLD

      call mpi_comm_size(mrworld,mrnum,mrinfo)
      call mpi_comm_rank(mrworld,mrpid,mrinfo)
      allocate(ympname(mrnum)) ; ympname(:) = ' '  ! process names
      call mpi_get_processor_name(ympname(1),ilen,mrinfo)

      call mpi_gather(ympname,80,mpi_character,   &
                      ympname,80,mpi_character,   &
                      nroot,mrworld,mrinfo)
      return
      end subroutine mpstart


!     =======================
!     subroutine mrdimensions
!     =======================

      subroutine mrdimensions
      use mpimod

      allocate(mrtru(mrnum)) ! all truncations
      mrtru(1) = ntru
      call mpi_gather(mrtru,1,mpi_itype,   &
                      mrtru,1,mpi_itype,   &
                      nroot,mrworld,mrinfo)
      
      call mrbcin(mrtru,mrnum)  ! broadcast truncations
      mintru = minval(mrtru(1:mrnum))
      mrdim  = (mintru+1) * (mintru+2) 
      return
      end

!     =================
!     subroutine mpstop
!     =================

      subroutine mpstop
      use mpimod

      call mpi_barrier(mrworld,mrinfo)
      call mpi_finalize(mrinfo)

      return
      end subroutine mpstop


!     ================
!     subroutine mrsum
!     ================

      subroutine mrsum(k)
      use mpimod
      integer :: k,n
!     print *,'mrsum b[',mrpid,']',k
      call mpi_allreduce(k,n,1,mpi_itype,MPI_SUM,mrworld,mrinfo)
!     print *,'mrsum a[',mrpid,']',n
      k = n
      return
      end


!     =================
!     subroutine mrdiff
!     =================

      subroutine mrdiff(pa,pd,kesp,klev)
      use mpimod

      real :: pa(kesp,klev)
      real :: pd(kesp,klev)
      real :: px(mrdim)
      real :: py(mrdim)

      do jlev = 1 , klev
         if (kesp == mrdim) then
            px(:) = pa(:,jlev)
         else
            call mrtrunc(pa(:,jlev),ntru,px,mintru)
         endif
         call mpi_allreduce(px,py,mrdim,mpi_rtype,MPI_SUM,mrworld,mrinfo)
         py(:) = py(:) - 2.0 * px(:)
         if (kesp == mrdim) then
            pd(:,jlev) = py(:)
         else
            call mrexpand(py,mintru,pd(1,jlev),ntru)
         endif
      enddo
      return
      end


!     =================
!     subroutine mrbcin
!     =================

      subroutine mrbcin(k,n)
      use mpimod
      integer :: k(n)

      call mpi_bcast(k,n,mpi_itype,NROOT,mrworld,mrinfo)

      return
      end


!     ==================
!     subroutine mrtrunc
!     ==================

!     Reduce truncation from nthi to ntlo

      subroutine mrtrunc(sphi,nthi,splo,ntlo)
      real :: sphi(2,*)
      real :: splo(2,*)

      jl = 1
      jh = 1
      do jm = 0 , ntlo
         jn = ntlo - jm
         splo(:,jl:jl+jn) = sphi(:,jh:jh+jn)
         jl = jl + ntlo - jm + 1
         jh = jh + nthi - jm + 1
      enddo
      return
      end


!     ===================
!     subroutine mrexpand
!     ===================

!     Expand truncation from ntlo to nthi

      subroutine mrexpand(splo,ntlo,sphi,nthi)
      real :: sphi(2,*)
      real :: splo(2,*)

      sphi(:,1:nthi) = 0.0
      jl = 1
      jh = 1
      do jm = 0 , ntlo
         jn = ntlo - jm
         sphi(:,jh:jh+jn) = splo(:,jl:jl+jn)
         jl = jl + ntlo - jm + 1
         jh = jh + nthi - jm + 1
      enddo
      return
      end


      subroutine mpbci(k) ! broadcast 1 integer
      return
      end

      subroutine mpbcin(k,n) ! broadcast n integer
      integer :: k(n)
      return
      end

      subroutine mpbcr(p) ! broadcast 1 real
      return
      end

      subroutine mpbcrn(p,n) ! broadcast n real
      real :: p(n)
      return
      end

      subroutine mpbcl(k) ! broadcast 1 logical
      logical :: k
      return
      end

      subroutine mpscin(k,n) ! scatter n integer
      integer :: k(n)
      return
      end

      subroutine mpscrn(p,n) ! scatter n real
      real :: p(n)
      return
      end

      subroutine mpscdn(p,n) ! scatter n double precision
      real (kind=8) :: p(n)
      return
      end

      subroutine mpscsp(pf,pp,klev) ! scatter spectral fields
      use catmod
      real pf(nesp,klev)
      real pp(nspp,klev)
      pp(1:nspp,1:klev) = pf(1:nspp,1:klev)
      return
      end

      subroutine mpscgp(pf,pp,klev) ! scatter gridpoint fields
      use catmod
      real pf(nlon*nlat,klev)
      real pp(nhor,klev)
      pp(1:nhor,1:klev) = pf(1:nhor,1:klev)
      return
      end

      subroutine mpgasp(pf,pp,klev) ! gather spectral fields
      use catmod
      real pf(nesp,klev)
      real pp(nspp,klev)
      pf(1:nspp,1:klev) = pp(1:nspp,1:klev)
      return
      end

      subroutine mpgagp(pf,pp,klev) ! gather gridpoint fields
      use catmod
      real pf(nhor,klev)
      real pp(nhor,klev)
      pf = pp
      return
      end

      subroutine mpgacs(pcs) ! gather cross sections
      return
      end

      subroutine mpgallsp(pf,pp,klev) ! gather spectral to all
      use catmod
      real pf(nesp,klev)
      real pp(nspp,klev)
      pf(1:nspp,1:klev) = pp(1:nspp,1:klev)
      return
      end

      subroutine mpsum(psp,klev) ! sum spectral fields
      return
      end

      subroutine mpsumsc(psf,psp,klev) ! sum & scatter spectral
      use catmod
      real psf(nesp,klev)
      real psp(nspp,klev)
      psp(1:nspp,1:klev) = psf(1:nspp,1:klev)
      return
      end

      subroutine mpsumr(pr,kdim) ! sum kdim reals
      return
      end subroutine mpsumr

      subroutine mpsumbcr(pr,kdim) ! sum & broadcast kdim reals
      return
      end

      subroutine mpreadsp(ktape,p,kdim,klev)
      real p(kdim,klev)
      read (ktape) p
      return
      end

      subroutine mpreadgp(ktape,p,kdim,klev)
      real p(kdim,klev)
      read (ktape) p
      return
      end

      subroutine mpwritesp(ktape,p,kdim,klev)
      real p(kdim,klev)
      write (ktape) p
      return
      end

      subroutine mpwritegp(ktape,p,kdim,klev)
      real p(kdim,klev)
      write (ktape) p
      return
      end

      subroutine mpwritegph(ktape,p,kdim,klev,ihead)
      real :: p(kdim,klev)
      integer :: ihead(8)
      write (ktape) ihead
      write (ktape) p

      return
      end


      subroutine mpi_info(nprocess,pid)    ! get nproc and pid
      integer nprocess, pid
      nprocess = 1
      pid = 0
      return
      end subroutine mpi_info


      subroutine mpgetsp(yn,p,kdim,klev)
      character (len=*) :: yn
      real :: p(kdim,klev)
      call get_restart_array(yn,p,kdim,kdim,klev)
      return
      end subroutine mpgetsp


      subroutine mpgetgp(yn,p,kdim,klev)
      character (len=*) :: yn
      real :: p(kdim,klev)
      call get_restart_array(yn,p,kdim,kdim,klev)
      return
      end subroutine mpgetgp


      subroutine mpputsp(yn,p,kdim,klev)
      character (len=*) :: yn
      real :: p(kdim,klev)
      call put_restart_array(yn,p,kdim,kdim,klev)
      return
      end subroutine mpputsp


      subroutine mpputgp(yn,p,kdim,klev)
      character (len=*) :: yn
      real :: p(kdim,klev)
      call put_restart_array(yn,p,kdim,kdim,klev)
      return
      end subroutine mpputgp


!      subroutine mpsurfgp(yn,p,kdim,klev)
!      character (len=*) :: yn
!      real :: p(kdim,klev)
!      call get_surf_array(yn,p,kdim,kdim,klev,iread)
!      return
!      end subroutine mpsurfgp
!
!
!      subroutine mpsurfyear(yn,p,kdim,kmon)
!      character (len=*) :: yn
!      real :: p(kdim,kmon)
!      call get_surf_year(yn,p,kdim,kmon,iread)
!      return
!      end subroutine mpsurfyear
!
!
!      subroutine mp3dyear(yn,p,kdim,klev,kmon)
!      character (len=*) :: yn
!      real :: p(kdim,klev,kmon)
!      call get_3d_year(yn,p,kdim,klev,kmon,iread)
!      return
!      end subroutine mp3dyear


      subroutine mpmaxval(p,kdim,klev,pmax)
      real :: p(kdim,klev)
      pmax = maxval(p(:,:))
      return
      end subroutine mpmaxval


      subroutine mpsumval(p,kdim,klev,psum)
      real :: p(kdim,klev)
      psum = sum(p(:,:))
      return
      end subroutine mpsumval

