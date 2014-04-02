      module rainmod
      use pumamod
!
!     version identifier (date)
!
      character(len=80) :: version = '08.09.2006 by Larry'
!
!     namelist parameter:
!

      integer :: nprl   = 1  ! switch for large scale precip (1/0=yes/no)
      integer :: nprc   = 1  ! switch for large convective precip (1/0=yes/no)
      integer :: ndca   = 1  ! switch for dry convective adjustment (1/0=yes/no)
      integer :: nmoment= 0  ! switch for momentum mixing (1/0=yes/no)

      real :: rsoft    =  1.0 ! rel. humidity value for 'soft' convection
      real :: clwcrit1 = -0.1 ! 1st critical vertical velocity for clouds
      real :: clwcrit2 =  0.0 ! 2nd critical vertical velocity for clouds
                              ! not active if clwcrit1 < clwcrit2
      real :: tadjust  =  1.0 ! adustment time scale (in timesteps; >= 1.0!)
                              ! to enable slow adjustment like betts miller

      real :: rcrit(NLEV)    ! critical relative hum. for non conv. clouds

!
!     global scalras
!

      real :: clwfac    = 0. ! smothing for cloud suppression
      real :: time4rain = 0. ! CPU time needed for rainmod
      real :: time4prl  = 0. ! CPU time spent for large scale precip
      real :: time4prc  = 0. ! CPU time spent for convective precip
      real :: time4dca  = 0. ! CPU time spent dry convection
      real :: time4cl   = 0. ! CPU time spent for cloud formation

!
!     global arrays
!

      integer :: icclev(NHOR,NLEV)  ! flag for convective cloud layer
      integer :: icctot(NHOR)       ! total convective cloud layers

      end module rainmod

!     ==================
!     SUBROUTINE RAININI
!     ==================

      subroutine rainini
      use rainmod
!
      namelist/rainmod_nl/nprl,nprc,ndca,rcrit,nmoment,clwcrit1,clwcrit2   &
     &                ,rsoft,tadjust
!
!     set default rcrit
!
      rcrit(:)=MAX(0.85,MAX(sigma(:),1.-sigma(:)))
!
      if(mypid==NROOT) then 
         open(11,file=rainmod_namelist)
         read(11,rainmod_nl)
         close(11)
         write(nud,'(/," ***********************************************")')
         write(nud,'(" * RAINMOD ",a35," *")') trim(version)
         write(nud,'(" ***********************************************")')
         write(nud,'(" * Namelist RAINMOD_NL from <rainmod_namelist> *")')
         write(nud,'(" ***********************************************")')
         write(nud,rainmod_nl)
      endif
!
!     broadcast namelist parameter
!
      call mpbci(nprl)
      call mpbci(nprc)
      call mpbci(ndca)
      call mpbci(nmoment)
      call mpbcr(rsoft)
      call mpbcr(tadjust)
      call mpbcr(clwcrit1)
      call mpbcr(clwcrit2)
      call mpbcrn(rcrit,NLEV)
      call mpbcrn(rcrit,NLEV)
!
!     smoothing factor for cloud suppression
!
      if(clwcrit1 > clwcrit2) then
       clwfac=1./(clwcrit1-clwcrit2)
      elseif(clwcrit1 < clwcrit2) then
       clwfac=-1.
      else
       clwfac=1.
      endif
!
      return
      end subroutine rainini

!     ===================
!     SUBROUTINE RAINSTEP
!     ===================

      subroutine rainstep
      use rainmod
!
!**   performance estimate
!
      if(ntime==1) call mksecond(zsec,0.)
!
!
!**   0) preset precipitation
!

      dprl(:)=0.
      dprc(:)=0.
      dprs(:)=0.

!
!**   1) dry convective adjustment
!

      if(ntime ==1) call mksecond(zsec1,0.)
      if(ndca == 1) call mkdca
      if(ntime ==1) then
       call mksecond(zsec1,zsec1)
       time4dca=time4dca+zsec1
      endif

!
!**   2) moist convective adjustment
!

      if(ntime ==1) call mksecond(zsec1,0.)
      if(nprc == 1) call mca
      if(ntime ==1) then
       call mksecond(zsec1,zsec1)
       time4prc=time4prc+zsec1
      endif

!
!**   3) large scale precipitation
!

      if(ntime ==1) call mksecond(zsec1,0.)
      if(nprl == 1) call mklsp
      if(ntime ==1) then
       call mksecond(zsec1,zsec1)
       time4prl=time4prl+zsec1
      endif

!
!**   4) diagnose cloud cover
!

      if(ntime ==1) call mksecond(zsec1,0.)
      call mkclouds
      if(ntime ==1) then
       call mksecond(zsec1,zsec1)
       time4cl=time4cl+zsec1
       call mksecond(zsec,zsec)
       time4rain=time4rain+zsec
      endif
!
      return
      end subroutine rainstep

!     ===================
!     SUBROUTINE RAINSTOP
!     ===================

      subroutine rainstop
      use rainmod
!
      if(mypid == NROOT .and. ntime == 1) then
       write(nud,*)'*********************************************'
       write(nud,*)' CPU usage in RAINSTEP (ROOT process only):  '
       write(nud,*)'    All routines   : ',time4rain,' s'
       write(nud,*)'    Large scale p  : ',time4prl,' s'
       write(nud,*)'    Convective p   : ',time4prc,' s'
       write(nud,*)'    Cloud formation: ',time4cl,' s'
       write(nud,*)'    Dry convection : ',time4dca,' s'
       write(nud,*)'*********************************************'
      endif
!
      return
      end subroutine rainstop


!     ==================
!     SUBROUTINE MKDQTGP
!     ==================

      subroutine mkdqtgp
      use rainmod
!
      real zsqt(NESP,NLEV)
!
      call mpgallsp(zsqt,sqt,NLEV)
      call sp2fl(zsqt,dqt,NLEV)
      call fc2gp(dqt,NLON,NLPP*NLEV)
!
      do jlev=1,NLEV
       dqt(:,jlev)=dqt(:,jlev)*psurf*ww/dp(:)
      enddo
!
      return
      end subroutine mkdqtgp

!     ================
!     SUBROUTINE MKLSP
!     ================

      subroutine mklsp
      use rainmod

!
!     local arrays
!

      real zq(NHOR)          ! preliminary q at t+dt
      real zt(NHOR)          ! preliminary t at t+dt
      real zqn(NHOR)         ! q after condensation
      real ztn(NHOR)         ! t after condensation
      real zdqdt(NHOR)       ! local q-tendency
      real zdtdt(NHOR)       ! local t-tendency
      real zqsat(NHOR)       ! saturation humidity
      real zlcpe(NHOR)       ! L/CP
      real zcor(NHOR)        ! correction term for saturation humidity

!
!     dbuq arrays
!

      integer :: kmaxl(1),kminl(1)

!
!     allocatable arrays for dbug diagnostic
!

      real, allocatable :: zprf1(:)
      real, allocatable :: zprf2(:)
      real, allocatable :: zprf3(:)
      real, allocatable :: zprf4(:)
      real, allocatable :: zprf5(:)
      real, allocatable :: zprf6(:)

!
!*    make condensation
!

      do jlev=1,NLEV

!
!     preset tendencies
!

       zdqdt(:)=0.
       zdtdt(:)=0.

!
!     compute preliminary t and q at t+dt
!

       zq(:)=dq(:,jlev)+dqdt(:,jlev)*deltsec2
       zt(:)=dt(:,jlev)+dtdt(:,jlev)*deltsec2

!
!*    set l/cp depending on t
!

       where(zt(:)>TMELT)
        zlcpe(:)=ALV/(ACPD*(1.+ADV*zq(:)))
       elsewhere
        zlcpe(:)=ALS/(ACPD*(1.+ADV*zq(:)))
       endwhere

!
!     saturation humidity
!

       zqsat(:)=RDBRV*RA1*exp(RA2*(zt(:)-TMELT)/(zt(:)-RA4))            &
     &         /(dp(:)*sigma(jlev))
       zcor(:)=1./(1.-(1./RDBRV-1.)*zqsat(:))
       zqsat(:)=zqsat(:)*zcor(:)

!
!     look for supersaturation
!

       where(zq(:) > zqsat(:))

!
!     make condensation (2 iterations)
!

        zdqdt(:)=(zqsat(:)-zq(:))                                       &
     &          /(1.0+zlcpe(:)*RA2*(TMELT-RA4)                          &
     &           *zqsat(:)*zcor(:)/(zt(:)-RA4)**2)
        ztn(:)=zt(:)-zdqdt(:)*zlcpe(:)
        zqn(:)=zq(:)+zdqdt(:)
        zqsat(:)=RDBRV*RA1*exp(RA2*(ztn(:)-TMELT)/(ztn(:)-RA4))         &
     &          /(dp(:)*sigma(jlev))
        zcor(:)=1./(1.-(1./RDBRV-1.)*zqsat(:))
        zqsat(:)=zqsat(:)*zcor(:)
        zdqdt(:)=(zqsat(:)-zqn(:))                                      &
     &          /(1.0+zlcpe(:)*RA2*(TMELT-RA4)                          &
     &           *zqsat(:)*zcor(:)/(ztn(:)-RA4)**2)
        ztn(:)=ztn(:)-zdqdt(:)*zlcpe(:)
        zqn(:)=zqn(:)+zdqdt(:)

!
!     final tendencies
!

        zdqdt(:)=(zqn(:)-zq(:))/deltsec2
        zdtdt(:)=(ztn(:)-zt(:))/deltsec2

!
!     add tendencies
!

        dqdt(:,jlev)=dqdt(:,jlev)+zdqdt(:)
        dtdt(:,jlev)=dtdt(:,jlev)+zdtdt(:)

!
!     add precipitation
!

        dprl(:)=dprl(:)-AMIN1(zdqdt(:),0.0)*dsigma(jlev)

       endwhere

!
!     dbug print out
!

       if(nprint==2) then
        allocate(zprf1(NLON*NLAT))
        allocate(zprf2(NLON*NLAT))
        allocate(zprf3(NLON*NLAT))
        allocate(zprf4(NLON*NLAT))
        allocate(zprf5(NLON*NLAT))
        allocate(zprf6(NLON*NLAT))
        call mpgagp(zprf1,zq,1)
        call mpgagp(zprf2,zt,1)
        call mpgagp(zprf3,zqsat,1)
        call mpgagp(zprf4,zdqdt,1)
        call mpgagp(zprf5,zdtdt,1)
        call mpgagp(zprf6,dprl,1)
        if(mypid==NROOT) then
         if(jlev==1) write(nud,*)'In mklsp:'
         write(nud,*)'L= ',jlev,' zt= ',zprf2(nprhor)                        &
     &         ,' zq= ',zprf1(nprhor),' zqsat= ',zprf3(nprhor)
         write(nud,*)'L= ',jlev,' dtdt= ',zprf5(nprhor)                      &
     &         ,' dqdt= ',zprf4(nprhor)
         write(nud,*)'L= ',jlev,' prl = ',zprf6(nprhor)                      &
     &                                      *0.5*deltsec2*ntspd*1000.
        endif
        deallocate(zprf1)
        deallocate(zprf2)
        deallocate(zprf3)
        deallocate(zprf4)
        deallocate(zprf5)
        deallocate(zprf6)
       endif

!
!     franks dbug
!

       if(ndiaggp==1) then
        dgp3d(:,jlev,7)=zdtdt(:)
        dgp3d(:,jlev,15)=zdqdt(:)
       endif

      enddo

!
!     convert precipitation into m/s and set snow fall
!

      dprl(:)=dprl(:)*dp(:)/GA/1000.
      where(dt(:,NLEV) <= TMELT) dprs(:)=dprs(:)+dprl(:)

!
!     dbug print out
!

      if(nprint==1) then
       allocate(zprf1(NLON*NLAT))
       call mpgagp(zprf1,dprl,1)
       if(mypid==NROOT) then
        zzmax=MAXVAL(zprf1)*0.5*deltsec2*ntspd*1000.
        kmaxl=MAXLOC(zprf1)
        zzmin=MINVAL(zprf1)*0.5*deltsec2*ntspd*1000.
        kminl=MINLOC(zprf1)
        write(nud,*)'In mklsp:'
        write(nud,*)'MAX prl= ',zzmax,' NHOR= ',kmaxl
        write(nud,*)'MIN prl= ',zzmin,' NHOR= ',kminl
       endif
       deallocate(zprf1)
      endif
!
      return
      end subroutine mklsp

!     ==============
!     SUBROUTINE MCA
!     ==============

      subroutine mca
      use rainmod
!
      parameter(zeps=1.E-3)   ! max. moist static energy residuum
      parameter(niter=10)     ! max. no iterations
!
!     local arrays
!
      integer ilift(NHOR)    ! lifting level ( = cloud base)
      integer itop(NHOR)     ! cloud top level
      integer icon(NLEV)     ! flag
      integer icf(NHOR,NLEV) ! flag (1 = effected by convection)
      integer iiter(NHOR)    ! counter
      real zqe(NHOR,NLEV)    ! moisture of environmental air
      real zte(NHOR,NLEV)    ! temperatur of environmental air
      real ztnew(NHOR,NLEV)  ! new temperature
      real zqnew(NHOR,NLEV)  ! new moisture
      real zlcp(NHOR,NLEV)   ! help array for factors (L/cp)
      real zl(NHOR,NLEV)     ! help array for L
      real zdqdt(NHOR,NLEV)  ! moisture diff. (parcel-environment)
      real zdtdt(NHOR,NLEV)  ! temperature diff. (parcel-environment)
      real zqse(NHOR,NLEV)   ! sat. moisture of environmental air

!!
!     array for momentum mixing
!
      real zue(NHOR,NLEV)    ! u of environmental air
      real zve(NHOR,NLEV)    ! v of environmental air
      real zdudt(NHOR,NLEV)  ! u-momentum difference (cloud-environment)
      real zdvdt(NHOR,NLEV)  ! v-momentum difference (cloud-environment)
      real zup(NHOR)         ! reference (cloud) u
      real zvp(NHOR)         ! reference (cloud) v
      real zsums(NHOR)       ! sum(dsigma) for mixed levels
!!

!
!     dbuq arrays
!

      integer :: kmaxl(1),kminl(1)

!
!     allocatable arrays for dbug diagnostic
!

      real, allocatable :: zprf1(:,:)
      real, allocatable :: zprf2(:,:)
      real, allocatable :: zprf3(:)
      real, allocatable :: zprf4(:)
      real, allocatable :: zprf5(:)
      real, allocatable :: zprf6(:)
      real, allocatable :: zprf7(:)
      real, allocatable :: zprf8(:)
      real, allocatable :: zprf9(:)

!
!*    0) initialize some fields
!

      iiter(:)=0
      ilift(:)=-1
      itop(:)=-1
      icctot(:)=0
      icclev(:,:)=0
      icf(:,:)=0
      zdtdt(:,:)=0.
      zdqdt(:,:)=0.

!
!*    1) update q and t, compute qs
!

      zqe(:,:)=dq(:,1:NLEV)+dqdt(:,1:NLEV)*deltsec2
      zte(:,:)=dt(:,1:NLEV)+dtdt(:,1:NLEV)*deltsec2

      do jlev=1,NLEV
      do jhor=1,NHOR
       zes=RA1*exp(RA2*(zte(jhor,jlev)-TMELT)/(zte(jhor,jlev)-RA4))
       zqs=RDBRV*zes/(dp(jhor)*sigma(jlev))
       zcor=1./(1.-(1./RDBRV-1.)*zqs)
       zqse(jhor,jlev)=zqs*zcor
      enddo
      enddo

!
!*    2) set L and L/CP depending on T
!

      where(zte(:,:) < TMELT)
       zlcp(:,:)=ALS/(ACPD*(1.+ADV*zqe(:,:)))
       zl(:,:)=ALS
      elsewhere
       zlcp(:,:)=ALV/(ACPD*(1.+ADV*zqe(:,:)))
       zl(:,:)=ALV
      endwhere

!
!*    3) make convection
!

      do jlev1=NLEV,1,-1
       do jhor=1,NHOR
        icon(:)=0
        kiter=0
!
!     compute saturation humidity (or threshold humidity if rsoft < 1)
!
        zes=RA1*exp(RA2*(zte(jhor,jlev1)-TMELT)/(zte(jhor,jlev1)-RA4))
        zqs=RDBRV*zes/(dp(jhor)*sigma(jlev1))
        zcor=1./(1.-(1./RDBRV-1.)*zqs)
        zqs=zqs*zcor*rsoft
!
!     test if moisture exceeds threshold
!
        if(zqs <= zqe(jhor,jlev1)) then
         icon(jlev1)=1
!
!     compute moist adiabatic temperature gradient (exponent)
!
         zdefac=RA2*(TMELT-RA4)/(zte(jhor,jlev1)-RA4)**2
         zesdt=zdefac*zes
         zp=sigma(jlev1)*dp(jhor)
         zexp=GASCON/ACPD*(zp+RDBRV*zl(jhor,jlev1)                      &
     &                    *zes/GASCON/zte(jhor,jlev1))                  &
     &                   /(zp+RDBRV*zesdt*zl(jhor,jlev1)/ACPD)
!
!     moist static energy and its derivative
!
         ztn=zte(jhor,jlev1)
         ztnew(jhor,jlev1)=ztn
         zqnew(jhor,jlev1)=zqs
         zdqsdt=zdefac*zqs*zcor
         zenergy1=(ztn-zte(jhor,jlev1))*dsigma(jlev1)
         zenergy2=zlcp(jhor,jlev1)*(zqs-zqe(jhor,jlev1))*dsigma(jlev1)
         zdenergy=(1.+zlcp(jhor,jlev1)*zdqsdt)*dsigma(jlev1)
!
!     search depth of the convective column
!
         do jlev2=jlev1-1,1,-1
          zfac=(sigma(jlev2)/sigma(jlev1))**zexp
          ztn=zte(jhor,jlev1)*zfac
          if(ztn > zte(jhor,jlev2) .and. icon(jlev2+1) == 1             &
     &      .and. zqse(jhor,jlev2)*rsoft <= zqe(jhor,jlev2)) then
           icon(jlev2)=1
           zqs=RDBRV*RA1*exp(RA2*(ztn-TMELT)/(ztn-RA4))                 &
     &        /(dp(jhor)*sigma(jlev2))
           zcor=1./(1.-(1./RDBRV-1.)*zqs)
           zqs=zcor*zqs*rsoft
           ztnew(jhor,jlev2)=ztn
           zqnew(jhor,jlev2)=zqs
           zdqsdt=RA2*(TMELT-RA4)*zqs*zcor/((ztn-RA4)*(ztn-RA4))
           zenergy1=(ztn-zte(jhor,jlev2))*dsigma(jlev2)                 &
     &             +zenergy1
           zenergy2=zlcp(jhor,jlev2)*(zqs-zqe(jhor,jlev2))*dsigma(jlev2)&
     &             +zenergy2
           zdenergy=zdenergy                                            &
     &             +(1.+zlcp(jhor,jlev2)*zdqsdt)*dsigma(jlev2)*zfac
           itop(jhor)=jlev2
          endif
         enddo
!
!     if convection: compute new t(jlev1)
!
         if(SUM(icon(:)) > 1) then
          ilift(jhor)=jlev1
          zenergy=zenergy1+zenergy2
          do while(zenergy > zeps .and. kiter < niter)
           ztn=ztnew(jhor,jlev1)-zenergy/zdenergy
           zqs=RDBRV*RA1*exp(RA2*(ztn-TMELT)/(ztn-RA4))                 &
     &        /(dp(jhor)*sigma(jlev1))
           zcor=1./(1.-(1./RDBRV-1.)*zqs)
           zqs=zcor*zqs*rsoft
           ztnew(jhor,jlev1)=ztn
           zqnew(jhor,jlev1)=zqs
           zdqsdt=RA2*(TMELT-RA4)*zqs*zcor/((ztn-RA4)*(ztn-RA4))
           zenergy1=(ztn-zte(jhor,jlev1))*dsigma(jlev1)
           zenergy2=zlcp(jhor,jlev1)*(zqs-zqe(jhor,jlev1))*dsigma(jlev1)
           zdenergy=(1.+zlcp(jhor,jlev1)*zdqsdt)*dsigma(jlev1)
           do jlev2=jlev1-1,1,-1
            if(icon(jlev2) == 1) then
             zfac=(sigma(jlev2)/sigma(jlev1))**zexp
             ztn=ztnew(jhor,jlev1)*zfac
             zqs=RDBRV*RA1*exp(RA2*(ztn-TMELT)/(ztn-RA4))               &
     &          /(dp(jhor)*sigma(jlev2))
             zcor=1./(1.-(1./RDBRV-1.)*zqs)
             zqs=zcor*zqs*rsoft
             ztnew(jhor,jlev2)=ztn
             zqnew(jhor,jlev2)=zqs
             zdqsdt=RA2*(TMELT-RA4)*zqs*zcor/((ztn-RA4)*(ztn-RA4))
             zenergy1=(ztn-zte(jhor,jlev2))*dsigma(jlev2)               &
     &               +zenergy1
             zenergy2=zlcp(jhor,jlev2)                                  &
     &               *(zqs-zqe(jhor,jlev2))*dsigma(jlev2)               &
     &               +zenergy2
             zdenergy=zdenergy                                          &
     &       +(1.+zlcp(jhor,jlev2)*zdqsdt)*dsigma(jlev2)*zfac
            endif
           enddo
           kiter=kiter+1
           zenergy=zenergy1+zenergy2
          enddo
         else
          icon(jlev1)=0
         endif
        endif
        icf(jhor,:)=icf(jhor,:)+icon(:)
        iiter(jhor)=iiter(jhor)+kiter
        do jlev2=jlev1,1,-1
         if(icon(jlev2)==1 .and. zenergy2 <= 0.) then
          zdtdt(jhor,jlev2)=ztnew(jhor,jlev2)-zte(jhor,jlev2)           &
     &                     +zdtdt(jhor,jlev2)
          zdqdt(jhor,jlev2)=zqnew(jhor,jlev2)-zqe(jhor,jlev2)           &
     &                     +zdqdt(jhor,jlev2)
          icf(jhor,jlev2)=icf(jhor,jlev2)+icon(jlev2)
          zte(jhor,jlev2)=ztnew(jhor,jlev2)
          zqe(jhor,jlev2)=zqnew(jhor,jlev2)
         endif
        enddo
       enddo
      enddo

!
!*    4) update dqdt and dtdt, compute precip.
!

      do jlev=1,NLEV
!
!     update dqdt and dtdt
!
       where(icf(:,jlev) > 0)
        dqdt(:,jlev)=zdqdt(:,jlev)/(tadjust*deltsec2)                   &
     &              +dqdt(:,jlev)
        dtdt(:,jlev)=zdtdt(:,jlev)/(tadjust*deltsec2)                   &
     &              +dtdt(:,jlev)
!
!     calc. precipitation
!
        dprc(:)=dprc(:)-zdqdt(:,jlev)*dsigma(jlev)

       endwhere
!
!     flag and counter for clouds:
!

       where(icf(:,jlev) > 0)
        icclev(:,jlev)=1
        icctot(:)=icctot(:)+1
       endwhere
       if(sigma(jlev) < 0.4) then
        where(icclev(:,jlev) == 1 .and. itop(:) == jlev)
         icclev(:,jlev)=2
        endwhere
       endif
!
      enddo

!
!*    5) final precip. in m/s and snow fall
!

      dprc(:)=dprc(:)*dp(:)/GA/1000./(tadjust*deltsec2)
      where(dt(:,NLEV) <= TMELT) dprs(:)=dprs(:)+dprc(:)

!!
!
!*    6) momentum mixing
!

      if(nmoment == 1) then

!
!     6a) update u and v
!

       zue(:,:)=du(:,1:NLEV)+dudt(:,1:NLEV)*deltsec2
       zve(:,:)=dv(:,1:NLEV)+dvdt(:,1:NLEV)*deltsec2

!
!     6b) reference u and v
!

       zup(:)=0.
       zvp(:)=0.
       zsums(:)=0.
       do jlev=1,NLEV
        where(icf(:,jlev) > 0)
         zup(:)=zup(:)+zue(:,jlev)*dsigma(jlev)
         zvp(:)=zvp(:)+zve(:,jlev)*dsigma(jlev)
         zsums(:)=zsums(:)+dsigma(jlev)
        endwhere
       enddo
       where(zsums(:) > 0.)
        zup(:)=zup(:)/zsums(:)
        zvp(:)=zvp(:)/zsums(:)
       endwhere

!
!     6c) make and add tendencies
!

       do jlev=1,NLEV
        where(icf(:,jlev) > 0)
          zdudt(:,jlev)=zup(:)-zue(:,jlev)
          zdvdt(:,jlev)=zvp(:)-zve(:,jlev)
          dudt(:,jlev)=dudt(:,jlev)+zdudt(:,jlev)/(tadjust*deltsec2)
          dvdt(:,jlev)=dvdt(:,jlev)+zdvdt(:,jlev)/(tadjust*deltsec2)
        endwhere
       enddo
!!!
!     test printout
!
!       allocate(zprf1(NLON*NLAT,NLEV))
!       allocate(zprf2(NLON*NLAT,NLEV))
!       call mpgagp(zprf1,zdudt,NLEV)
!       call mpgagp(zprf2,zdvdt,NLEV)
!       if(mypid==NROOT) then
!        write(nud,*)'In Kuo momentum:'
!        write(nud,*)'MIN/MAX du= ',minval(zprf1),maxval(zprf1)
!        write(nud,*)'MIN/MAX dv= ',minval(zprf2),maxval(zprf2)
!       endif
!       deallocate(zprf1)
!       deallocate(zprf2)
!!
      endif
!!

!
!*    dbug diagnostic for nprint > 0 (see pumamod)
!

      if(nprint==1) then
       allocate(zprf6(NLON*NLAT))
       call mpgagp(zprf6,dprc,1)
       if(mypid==NROOT) then
        zzmax=MAXVAL(zprf6)*0.5*deltsec2*ntspd*1000.
        kmaxl=MAXLOC(zprf6)
        zzmin=MINVAL(zprf6)*0.5*deltsec2*ntspd*1000.
        kminl=MINLOC(zprf6)
        write(nud,*)'In MCA:'
        write(nud,*)'MAX prc= ',zzmax,' NHOR= ',kmaxl
        write(nud,*)'MIN prc= ',zzmin,' NHOR= ',kminl
       endif
       deallocate(zprf6)
      elseif(nprint==2) then
       allocate(zprf1(NLON*NLAT,NLEV))
       allocate(zprf2(NLON*NLAT,NLEV))
       allocate(zprf3(NLON*NLAT))
       allocate(zprf4(NLON*NLAT))
       allocate(zprf5(NLON*NLAT))
       allocate(zprf6(NLON*NLAT))
       allocate(zprf7(NHOR))
       allocate(zprf8(NHOR))
       allocate(zprf9(NHOR))
       zprf7(:)=real(ilift(:))
       zprf8(:)=real(itop(:))
       zprf9(:)=real(iiter(:))
       call mpgagp(zprf1,zdqdt,NLEV)
       call mpgagp(zprf2,zdtdt,NLEV)
       call mpgagp(zprf3,zprf7,1)
       call mpgagp(zprf4,zprf8,1)
       call mpgagp(zprf5,zprf9,1)
       call mpgagp(zprf6,dprc,1)
       if(mypid==NROOT) then
        write(nud,*)'In MCA:'
        write(nud,*)'ilift= ',zprf3(nprhor),' itop= ',zprf4(nprhor)
        write(nud,*)'No iterations= ',zprf5(nprhor)
        do jlev=1,NLEV
         if(zprf3(nprhor) >= jlev .and. zprf4(nprhor) <= jlev) then
          zzdt=zprf2(nprhor,jlev)/deltsec2
          zzdq=zprf1(nprhor,jlev)/deltsec2
          write(nud,*)'L= ',jlev,' dtdt= ',zzdt,' dqdt= ',zzdq
         endif
        enddo
        write(nud,*)'prc (mm/day)= ',zprf6(nprhor)*1000.*ntspd*deltsec2*0.5
       endif
       deallocate(zprf1)
       deallocate(zprf2)
       deallocate(zprf3)
       deallocate(zprf4)
       deallocate(zprf5)
       deallocate(zprf6)
       deallocate(zprf7)
       deallocate(zprf8)
       deallocate(zprf9)
      endif

!
!     franks dbug
!

      if(ndiaggp==1) then
       dgp3d(:,:,8)=zdtdt(:,:)/deltsec2
       dgp3d(:,:,16)=zdqdt(:,:)/deltsec2
      end if
!
      return
      end subroutine mca

!     ===================
!     SUBROUTINE MKCLOUDS
!     ===================

      subroutine mkclouds
      use rainmod
!
!     compute cloud properties (cover and liquid water content)
!

      parameter(zclmax=1.)
      parameter(zcca=0.245,zccb=0.125)
      parameter(zccmax=0.8,zccmin=0.05)

      real zq(NHOR),zt(NHOR),zqsat(NHOR),zrh(NHOR)
      real zrcrit(NHOR),zccconv(NHOR),zcctot(NHOR)
      real zwfac(NHOR),zcc(NHOR)
      real zzf(NHOR,NLEV),zzh(NHOR),zdh(NHOR)

!
!*    preset
!

      dql(:,:)=0.
      dcc(:,1:NLEV)=0.

!
!*    convective cloud cover for random overlab:
!

      zrfac=86400.*1000. ! convert m/s into mm/day
      where(dprc(:) > 0.)
       zcctot(:)=zcca+zccb*log(dprc(:)*zrfac)
       zcctot(:)=AMIN1(zccmax,AMAX1(zccmin,zcctot(:)))
       zccconv(:)=1.-(1.-zcctot(:))**(1./real(icctot(:)))
      end where

!
!*    cloud cover for the levels
!

      do jlev=1,NLEV
       zt(:)=dt(:,jlev)+dtdt(:,jlev)*deltsec2
       zq(:)=dq(:,jlev)+dqdt(:,jlev)*deltsec2
       zqsat(:)=RDBRV*RA1*exp(RA2*(zt(:)-TMELT)/(zt(:)-RA4))            &
     &         /(dp(:)*sigma(jlev))
       zqsat(:)=zqsat(:)/(1.-(1./RDBRV-1.)*zqsat(:))
       zrh(:)=zq(:)/zqsat(:)
!
!     convective clouds
!
!
!     a) cu/cb
!
       where(icclev(:,jlev) > 0)
        dcc(:,jlev)=zccconv(:)
       endwhere
!
!     b) ci
!
       where(icclev(:,jlev)==2 .and. zcctot(:) >= 0.3)
        dcc(:,jlev)=AMIN1(dcc(:,jlev)+2.*(zcctot(:)-0.3),zclmax)
       endwhere
!
!     reduce relative humidity
!
        zrh(:)=zrh(:)*(1.-dcc(:,jlev))
!
!     c) non convective clouds
!
       zrcrit(:)=rcrit(jlev)
       zwfac(:)=1.
       if(sigma(jlev) > 0.8 .and. clwfac > 0.)  then
        zwfac(:)=clwfac*AMAX1(0.,clwcrit1-dw(:,jlev))
        where(clwcrit2 > dw(:,jlev)) zwfac(:)=1.
       endif
       zcc(:)=zwfac(:)*AMAX1(0.,(zrh(:)-zrcrit(:))/(1.-zrcrit(:)))**2
       dcc(:,jlev)=AMIN1(dcc(:,jlev)+zcc(:),zclmax)

      enddo

!
!*    precipitable water (and z which is later used to compute dql)
!

      dqvi(:)=0.
      zzh(:)=0.
      do jlev=NLEV,2,-1
       dqvi(:)=dq(:,jlev)*dsigma(jlev)+dqvi(:)
       zdh(:)=-dt(:,jlev)*GASCON/GA*ALOG(sigmah(jlev-1)/sigmah(jlev))
       zzf(:,jlev)=zzh(:)+zdh(:)*0.5
       zzh(:)=zzh(:)+zdh(:)
      enddo
      dqvi(:)=dq(:,1)*dsigma(1)+dqvi(:)
      zdh(:)=-dt(:,1)*GASCON/GA*ALOG(sigma(1)/sigmah(1))*0.5
      zzf(:,1)=zzh(:)+zdh(:)

!
!     cloud liquid water (see CCM3 description: Kiehl et al. 1996 pp49+50)
!     and total cloud cover for random overlap
!

      dcc(:,NLEP)=1.
      dqvi(:)=dqvi(:)*dp(:)/GA
      zzh(:)=700.*ALOG(1.+dqvi(:))
      do jlev=1,NLEV
       where(zzh(:) > 0. .and. dcc(:,jlev) > 0.)
        dql(:,jlev)=0.00021*EXP(-zzf(:,jlev)/zzh(:))*GASCON*dt(:,jlev)  &
     &             /(sigma(jlev)*dp(:))
        dql(:,jlev)=MAX(dql(:,jlev),1.E-9)
       endwhere
       dcc(:,NLEP)=dcc(:,NLEP)*(1.-dcc(:,jlev))
      enddo
      dcc(:,NLEP)=1.-dcc(:,NLEP)
!
      return
      end subroutine mkclouds

!     ==============
!     SUBROUTINE DCA
!     ==============

      subroutine mkdca
      use rainmod
!
!     calculate t and q tendencies due to dry convection
!
      integer itop(NHOR),ibase(NHOR),icon(NHOR)
      real zdtdt(NHOR,NLEV)
      real zdqdt(NHOR,NLEV)
      real zt(NHOR,NLEV),zth(NHOR,NLEV)
      real zq(NHOR,NLEV),zqn(NHOR,NLEV)
      real ztht(NHOR)
      real zqt(NHOR)
      real zsum1(NHOR),zsum2(NHOR)
      real zsumq1(NHOR),zsumq2(NHOR)
      real zskap(NLEV)
!
      do jlev = 1,NLEV
       zskap(jlev)=sigma(jlev)**akap
       zt(:,jlev)=dt(:,jlev)+dtdt(:,jlev)*deltsec2
       zth(:,jlev)=zt(:,jlev)/zskap(jlev)
       zq(:,jlev)=dq(:,jlev)+dqdt(:,jlev)*deltsec2
       zqn(:,jlev)=zq(:,jlev)
      enddo
!
      jiter=0
 1000 continue
       icon(:)=0
       ibase(:)=0
       itop(:)=NLEV
       ztht(:)=0.
       zqt(:)=0.
       do jlev=1,NLEM
        jlep=jlev+1
        where(zth(:,jlev) < zth(:,jlep))
         ibase(:)=jlep
         itop(:)=jlev-1
         zsum1(:)=zskap(jlep)*dsigma(jlep)                              &
     &              +zskap(jlev)*dsigma(jlev)
         zsum2(:)=zskap(jlep)*dsigma(jlep)*zth(:,jlep)                  &
     &           +zskap(jlev)*dsigma(jlev)*zth(:,jlev)
         zsumq1(:)=dsigma(jlep)+dsigma(jlev)
         zsumq2(:)=dsigma(jlep)*zqn(:,jlep)+dsigma(jlev)*zqn(:,jlev)
         ztht(:)=zsum2(:)/zsum1(:)
         zqt(:)=zsumq2(:)/zsumq1(:)
         icon(:)=1
        endwhere
       enddo
       if(SUM(icon(:)) > 0) then
        do jlev=NLEM-1,1,-1
         where(itop(:) == jlev .and. zth(:,jlev) < ztht(:))
          itop(:)=jlev-1
          zsum1(:)=zsum1(:)+zskap(jlev)*dsigma(jlev)
          zsum2(:)=zsum2(:)+zskap(jlev)*dsigma(jlev)*zth(:,jlev)
          zsumq1(:)=zsumq1(:)+dsigma(jlev)
          zsumq2(:)=zsumq2(:)+dsigma(jlev)*zqn(:,jlev)
          ztht(:)=zsum2(:)/zsum1(:)
          zqt(:)=zsumq2(:)/zsumq1(:)
         endwhere
        enddo
        do jlev=1,NLEV
         where(jlev > itop(:) .and. jlev <= ibase(:))
          zth(:,jlev)=ztht(:)
          zqn(:,jlev)=zqt(:)
         endwhere
        enddo
        jiter=jiter+1
        if(jiter > 2*NLEV) goto 1001
        goto 1000
      endif
!
      do jlev=1,NLEV
       zdtdt(:,jlev)=(zth(:,jlev)*zskap(jlev)-zt(:,jlev))/deltsec2
       zdqdt(:,jlev)=(zqn(:,jlev)-zq(:,jlev))/deltsec2
       dtdt(:,jlev)=dtdt(:,jlev)+zdtdt(:,jlev)
       dqdt(:,jlev)=dqdt(:,jlev)+zdqdt(:,jlev)
      enddo
!
!     franks dbug
!
      if(ndiaggp==1) then
       dgp3d(:,1:NLEV,9)=zdtdt(:,1:NLEV)
       dgp3d(:,1:NLEV,17)=zdqdt(:,1:NLEV)
      endif
!
      return
 1001 continue
      end subroutine mkdca
