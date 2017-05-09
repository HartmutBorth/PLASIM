      module miscmod
      use pumamod
!
!     version identifier
!
      character(len=80) :: version = '03.03.2003 by Larry'
!
!     namelist parameters:
!
      integer :: nfixer = 1 ! switch for negative humidity fix (1/0 : on/off)
      integer :: nudge  = 0 ! switch for t-nudging (1/0 : on/off, 2: flxcorr)
      real    :: tnudget = 10. ! timescale for t-nudging (days)
!
!     global scalars
!
      real :: time4mis = 0.  ! CPU time needed for misc routines
      real :: time4fix = 0.  ! CPU time needed for fixer
!
!     global arrays
!
      real :: zgw(NHOR)          ! gaussian weights
      real :: dtnudge(NHOR,0:13) ! climatological temperature (annual cycle)
      real :: dfnudge(NHOR,0:13) ! climatological fluxes correction (annual cycle)
!
      end module miscmod

!     ==================
!     SUBROUTINE MISCINI
!     ==================

      subroutine miscini
      use miscmod
!
!     initialization
!
      namelist/miscmod_nl/nfixer,nudge,tnudget
!
      if (mypid == NROOT) then
         open(11,file=miscmod_namelist)
         read(11,miscmod_nl)
         close(11)
         write(nud,'(/," ***********************************************")')
         write(nud,'(" * MISCMOD ",a35," *")') trim(version)
         write(nud,'(" ***********************************************")')
         write(nud,'(" * Namelist MISCMOD_NL from <miscmod_namelist> *")')
         write(nud,'(" ***********************************************")')
         write(nud,miscmod_nl)
      endif

      call mpbci(nfixer)
      call mpbci(nudge)
      call mpbcr(tnudget)

      if (nudge == 1) then
         call mpsurfgp('dtnudge',dtnudge,NHOR,14)
      endif

      if (nudge == 2) then
         call mpsurfgp('dfnudge',dfnudge,NHOR,14)
      endif
!
!     gaussian weights
!
      do jlat=1,NLPP
       jhor1=(jlat-1)*NLON+1
       jhor2=jlat*NLON
       zgw(jhor1:jhor2)=gwd(jlat)
      enddo
!
!     nudging time in seconds
!
      tnudget=TWOPI*tnudget/ww
!
      return
      end subroutine miscini

!     ===================
!     SUBROUTINE MISCSTEP
!     ===================

      subroutine miscstep
      use miscmod

      if(ntime == 1) then
       call mksecond(zsec,0.)
       call mksecond(zsec1,0.)
      endif
      call fixer
      if(ntime == 1) then
       call mksecond(zsec1,zsec1)
       time4fix=time4fix+zsec1
      endif
!
      if(nudge == 1) then
       call mknudge
      elseif(nudge == 2) then
       call mkflcor
      endif
!
      if(ntime == 1) then
       call mksecond(zsec,zsec)
       time4mis=time4mis+zsec
      endif

      return
      end subroutine miscstep

!     ===================
!     SUBROUTINE MISCSTOP
!     ===================

      subroutine miscstop
      use miscmod

      if(mypid == NROOT .and. ntime == 1) then
       write(nud,*)'*********************************************'
       write(nud,*)' CPU usage in MISCSTEP (ROOT process only):  '
       write(nud,*)'    All routines : ',time4mis,' s'
       write(nud,*)'    Fixer        : ',time4fix,' s'
       write(nud,*)'*********************************************'
      endif

      return
      end subroutine miscstop


      subroutine fixer
      use miscmod
!
      real zqn(NHOR,NLEV)    ! humidity in and out
      real zneg(NHOR)        ! moisture needed
      real zpos(NHOR)        ! moisture available
!
      real zsum(2)           ! work array for mpi (needs to be an array)
!
!*    franks dbug
!
      if(ndiaggp > 0) then
       dgp3d(:,1:NLEV,12)=dq(:,1:NLEV)
      endif
!
      zqn(:,:)=dq(:,1:NLEV)
      zneg(:)=0.
      zpos(:)=0.
      do jhor=1,NHOR
       zneg(jhor)=-SUM(zqn(jhor,:)*dsigma(:),MASK=(zqn(jhor,:) < 0.))   &
     &           *dp(jhor)*zgw(jhor)
       zpos(jhor)=SUM(zqn(jhor,:)*dsigma(:),MASK=(zqn(jhor,:) > 0.))    &
     &           *dp(jhor)*zgw(jhor)
      enddo
!
!*    fix within grid point collumn
!
      if(SUM(zneg(:)) > 0.) call fix0d(zqn,zneg,zpos)
!
!*    fix within latitude
!
      if(SUM(zneg(:)) > 0.) call fix1d(zqn,zneg,zpos)
!
!*    fix globally
!
      zsum(1)=SUM(zneg(:))
      zsum(2)=SUM(zpos(:))
      call mpsumbcr(zsum,2)
      zzneg=zsum(1)
      zzpos=zsum(2)
      if(zzneg > 0.) then
       zfac=(zzpos-zzneg)/zzpos
       zqn(:,:)=AMAX1(zqn(:,:)*zfac,0.)
      endif
!
!*    tendencies
!
      do jlev=1,NLEV
       gqdt(:,jlev)=gqdt(:,jlev)                                        &
     &             +dp(:)*(zqn(:,jlev)-dq(:,jlev))/delt2/psurf
      enddo
!
!*    set new dq
!
      dq(:,1:NLEV)=zqn(:,:)
!
!*    franks dbug
!
      if(ndiaggp > 0) then
       do jlev=1,NLEV
        dgp3d(:,jlev,12)=(dq(:,jlev)-dgp3d(:,jlev,12))/deltsec2*dp(:)   &
     &                  *dsigma(jlev)/ga/1000.
       enddo
      endif
!
      return
      end subroutine fixer

!     ================
!     SUBROUTINE FIX0D
!     ================

      subroutine fix0d(pqn,pneg,ppos)
      use miscmod
!
      real pqn(NHOR,NLEV)
      real pneg(NHOR),ppos(NHOR)
!
      do jhor=1,NHOR
       if (pneg(jhor) > 0. .and. ppos(jhor) >= pneg(jhor)) then
        zfac=(ppos(jhor)-pneg(jhor))/ppos(jhor)
        pqn(jhor,:)=AMAX1(pqn(jhor,:)*zfac,0.)
        pneg(jhor)=0.
        ppos(jhor)=ppos(jhor)*zfac
       endif
      enddo
!
      return
      end subroutine fix0d

!     ================
!     SUBROUTINE FIX1D
!     ================

      subroutine fix1d(pqn,pneg,ppos)
      use miscmod
!
      real pqn(NHOR,NLEV)
      real pneg(NHOR),ppos(NHOR)
!
      do jlat=1,NLPP
       jhor1=(jlat-1)*NLON+1
       jhor2=jlat*NLON
       zneg=SUM(pneg(jhor1:jhor2))
       if(zneg > 0.) then
        zpos=SUM(ppos(jhor1:jhor2))
        if(zpos >= zneg) then
          zfac=(zpos-zneg)/zpos
          pqn(jhor1:jhor2,:)=AMAX1(pqn(jhor1:jhor2,:)*zfac,0.)
          pneg(jhor1:jhor2)=0.
          ppos(jhor1:jhor2)=ppos(jhor1:jhor2)*zfac
        endif
       endif
      enddo
!
      return
      end subroutine fix1d

!     ==================
!     SUBROUTINE MKNUDGE
!     ==================

      subroutine mknudge
      use miscmod
!
!     make nudging of upper level temperature
!
      real ztcl(NHOR)
      real zdtdt(NHOR)
!
!     get climatology
!
      call momint(nperpetual,nstep+1,imonth,jm,zw)
      ztcl(:) = (1.0 - zw) * dtnudge(:,imonth) + zw * dtnudge(:,jm)
!
!     tendencies due to nudging
!
      zdtdt(:)=(ztcl(:)-dt(:,1))/tnudget
!
!     add tendencies
!
      dtdt(:,1)=dtdt(:,1)+zdtdt(:)
!
!     franks diagnostic
!
      if(ndiaggp==1) then
       dgp3d(:,2:NLEV,2)=0.
       dgp3d(:,1,2)=zdtdt(:)
      endif
!
!     entropy/energy diagnostics
!
      if(nentropy > 0) then
       dentropy(:,6)=zdtdt(:)/dentrot(:,1)                              &
     &              *acpd*(1.+adv*dentroq(:,1))*dentrop(:)/ga*dsigma(1)
       if(nentro3d > 0) then
        dentro3d(:,2:NLEV,6)=0.
        dentro3d(:,1,6)=dentropy(:,6) 
       endif
      endif
      if(nenergy > 0) then
       denergy(:,6)=zdtdt(:)                                            &
     &              *acpd*(1.+adv*dq(:,1))*dp(:)/ga*dsigma(1)
       if(nener3d > 0) then
        dener3d(:,2:NLEV,6)=0.
        dener3d(:,1,6)=denergy(:,6)
       endif
      endif
!
      return
      end

!     ==================
!     SUBROUTINE MKFLCOR
!     ==================

      subroutine mkflcor
      use miscmod
      implicit none
      integer :: imonth   ! current month
      integer :: jm       ! next or previous month
      real    :: zw       ! interpolation weight
!
!     make flux correction for upper level temperature
!
      real :: zdtdt(NHOR)
!
!     get climatology
!
      call momint(nperpetual,nstep+1,imonth,jm,zw)
      zdtdt(:) = (1.0 - zw) * dfnudge(:,imonth) + zw * dfnudge(:,jm)
!
!     add tendencies due to flux correction
!
      dtdt(:,1)=dtdt(:,1)+zdtdt(:)
!
!     franks diagnostic
!
      if(ndiaggp==1) then
       dgp3d(:,2:NLEV,2)=0.
       dgp3d(:,1,2)=zdtdt(:)
      endif
!
      return
      end
