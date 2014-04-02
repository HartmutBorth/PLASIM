      module rainmod
      use pumamod
!
!     version identifier (date)
!
      character(len=80) :: version = '22.11.2002 by Larry'
!
!     namelist parameter:
!

      integer :: nprl   = 1  ! switch for large scale precip (1/0=yes/no)
      integer :: nprc   = 1  ! switch for large convective precip (1/0=yes/no)
      integer :: ndca   = 1  ! switch for dry convective adjustment (1/0=yes/no)
      integer :: nshallow = 0! switch for shallow convection (1/0=yes/no)
      integer :: nmoment  = 0! switch for momentum mixing (1/0=yes/no)

      real :: clwcrit1 = -0.1 ! 1st critical vertical velocity for clouds
      real :: clwcrit2 =  0.0 ! 2nd critical vertical velocity for clouds
                              ! not active if clwcrit1 < clwcrit2

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
      namelist/rainmod_nl/nprl,nprc,ndca,nshallow,nmoment,rcrit,clwcrit1   &
     &                ,clwcrit2
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
      call mpbci(nshallow)
      call mpbci(nmoment)
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
      if(nprc == 1) call mkbm
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

!     ===============
!     SUBROUTINE MKBM
!     ===============

      subroutine mkbm
      use rainmod
      parameter(zpdeep=70000.)   ! deep convection threshold
      parameter(zstab=0.85)      ! fraction of moist static stability
      parameter(zdpbase=-2500.)  ! saturation point difference at cloud base
      parameter(zdpfreez=-4000.) ! saturation point difference at freezing l.
      parameter(zdptop=-2000.)   ! saturation point difference at cloud top
      parameter(zdpsh=-5000.)    ! saturation point difference for shallow c.
      parameter(ztaud=7200.)     ! adjustment time scale for deep c.
      parameter(ztaus=14400.)    ! adjustment time scale for shallow c.
      parameter(zbeta=1.0)       ! beta

!
!     local arrays
!

      integer ilift(NHOR)    ! lifting level (= cloud base)
      integer ifreez(NHOR)   ! freezing level
      integer itop(NHOR)     ! cloud top level
      integer ideep(NHOR)    ! flag for deep convection
      integer ishallow(NHOR) ! flag for shallow convection
      integer icf(NHOR,NLEV) ! flag for convective layers
      real zqe(NHOR,NLEV)    ! moisture of environmental air
      real zte(NHOR,NLEV)    ! temperatur of environmental air
      real zdqrdt(NHOR,NLEV) ! dqref/dTref
      real zdqdt(NHOR,NLEV)  ! q-tendency
      real zdtdt(NHOR,NLEV)  ! t-tendency
      real zl(NHOR,NLEV)     ! latent heat of vapourisation/sublimation
      real zcp(NHOR,NLEV)    ! specific heat
      real zqp(NHOR,NLEV)    ! humidity of lifted parcel
      real ztp(NHOR,NLEV)    ! temperature of lifted parcel
      real ztref(NHOR,NLEV)  ! reference humidity profile
      real zqref(NHOR,NLEV)  ! reference temperature profile
      real zdpres(NHOR,NLEV) ! saturation pressure difference
      real zptop(NHOR)       ! cloud top pressure
      real zpfreez(NHOR)     ! freezing level pressure
      real zplift(NHOR)      ! lifting level pressure
      real zdtfreez(NHOR)    ! Tref-Tp at freezing level
      real zener(NHOR)       ! sum(cp (Tr-T) + L (qr-q))
      real zsumt(NHOR)       ! sum(tr-t)
      real zsumq(NHOR)       ! sum(qr-q)
      real zdsig(NHOR)       ! sum(dsigma)
      real zprc(NHOR)        ! precipitation
      real zspb(NHOR)        ! saturation pressure at cloud base
      real zslope(NHOR)      ! mixing line slope


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
!*    dbug arrays
!

      integer :: kmaxl(1), kminl(1)

!
!*    allocatable arrays for dbug diagnostic
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
      real, allocatable :: zprf10(:)
      real, allocatable :: zprf11(:)
      real, allocatable :: zprf12(:)

      zeps=0.01              ! a small number

!
!*    0) initialize some fields
!

      ilift(:)=-1
      itop(:)=-1
      ifreez(:)=-1
      ideep(:)=0
      ishallow(:)=0
      icctot(:)=0
!
      zptop(:)=1.E10
      zpfreez(:)=1.E10
      zplift(:)=1.E10
      zprc(:)=0.
!
      zdtdt(:,:)=0.
      zdqdt(:,:)=0.
      icclev(:,:)=0

!
!*    1) update q and t
!

      zqe(:,:)=dq(:,1:NLEV)+dqdt(:,1:NLEV)*deltsec2
      zte(:,:)=dt(:,1:NLEV)+dtdt(:,1:NLEV)*deltsec2

!
!     set L & CP depending on t
!

      zcp(:,:)=ACPD*(1.+ADV*zqe(:,:))
      where(zte(:,:) < TMELT)
       zl(:,:)=ALS
      elsewhere
       zl(:,:)=ALV
      endwhere

!
!*    2) search for lifting level
!

      do jlev=NLEM,2,-1
       jlep=jlev+1
       zfac=(sigma(jlev)/sigma(jlep))**akap
       do jhor=1,NHOR
        if(ilift(jhor) == -1) then

!
!     adiabatic uplift of parcel from below
!

         ztnew=zte(jhor,jlep)*zfac
         zqnew=zqe(jhor,jlep)

!
!     saturation humidity:
!

         zqsnl=RDBRV*RA1*exp(RA2*(ztnew-TMELT)/(ztnew-RA4))             &
     &        /(dp(jhor)*sigma(jlev))
         zcor=1./(1.-(1./RDBRV-1.)*zqsnl)
         zqsnl=zqsnl*zcor

!
!     condensation if supersaturated
!

         if(zqsnl < zqnew) then

!
!     update constants (ice/water phase)
!

          if(ztnew < TMELT) then
           zlcp=ALS/(ACPD*(1.+ADV*zqnew))
          else
           zlcp=ALV/(ACPD*(1.+ADV*zqnew))
          endif

!
!     temperature correction term
!

          zdqsdt=RA2*(TMELT-RA4)*zqsnl*zcor/((ztnew-RA4)*(ztnew-RA4))


!
!     update t and q
!

          zrain=(zqnew-zqsnl)/(1.+zlcp*zdqsdt)
          ztnew=ztnew+zlcp*zrain
          zqnew=zqnew-zrain

!
!     second iteration:
!

          zqsnl=RDBRV*RA1*exp(RA2*(ztnew-TMELT)/(ztnew-RA4))            &
     &         /(dp(jhor)*sigma(jlev))
          zcor=1./(1.-(1./RDBRV-1.)*zqsnl)
          zqsnl=zqsnl*zcor
          zdqsdt=RA2*(TMELT-RA4)*zqsnl*zcor/((ztnew-RA4)*(ztnew-RA4))
          zrain=(zqnew-zqsnl)/(1.+zlcp*zdqsdt)
          ztnew=ztnew+zlcp*zrain
          zqnew=zqnew-zrain

!
!     check for buyoncy
!

          if(ztnew > zte(jhor,jlev)) then

!
!     if convection: set lifting level and parcel q,t (zpq,ztp)
!

           ilift(jhor)=jlep
           itop(jhor)=jlev
           zplift(jhor)=dp(jhor)*sigma(jlep)
           zptop(jhor)=dp(jhor)*sigma(jlev)
           ztp(jhor,jlev)=ztnew
           zqp(jhor,jlev)=zqnew
          endif
         endif
        endif
       enddo
      enddo

!
!*    3) search for top level
!

      do jlev=NLEM,1,-1
       jlep=jlev+1
       zfac=(sigma(jlev)/sigma(jlep))**akap
       do jhor=1,NHOR
        if(jlep == itop(jhor)) then

!
!     3a) parcel uplift
!

         ztnew=ztp(jhor,jlep)*zfac
         zqnew=zqp(jhor,jlep)

!
!     saturation humidity
!

         zqsnl=RDBRV*RA1*exp(RA2*(ztnew-TMELT)/(ztnew-RA4))             &
     &        /(dp(jhor)*sigma(jlev))
         zcor=1./(1.-(1./RDBRV-1.)*zqsnl)
         zqsnl=zqsnl*zcor

!
!     3b) more condensation
!

         if(zqsnl < zqnew) then

!
!     update constants (ice/water phase)
!

          if(ztnew < TMELT) then
           zlcp=ALS/(ACPD*(1.+ADV*zqnew))
          else
           zlcp=ALV/(ACPD*(1.+ADV*zqnew))
          endif

!
!     temperature correction term
!

          zdqsdt=RA2*(TMELT-RA4)*zqsnl*zcor/((ztnew-RA4)*(ztnew-RA4))

!
!     update t,q
!

          zrain=(zqnew-zqsnl)/(1.+zlcp*zdqsdt)
          ztnew=ztnew+zlcp*zrain
          zqnew=zqnew-zrain

!
!     second iteration:
!

          zqsnl=RDBRV*RA1*exp(RA2*(ztnew-TMELT)/(ztnew-RA4))            &
     &         /(dp(jhor)*sigma(jlev))
          zcor=1./(1.-(1./RDBRV-1.)*zqsnl)
          zqsnl=zqsnl*zcor
          zdqsdt=RA2*(TMELT-RA4)*zqsnl*zcor/((ztnew-RA4)*(ztnew-RA4))
          zrain=(zqnew-zqsnl)/(1.+zlcp*zdqsdt)
          ztnew=ztnew+zlcp*zrain
          zqnew=zqnew-zrain

!
!     3c) check for buyoncy
!

          if(ztnew > zte(jhor,jlev)) then

!
!     3d) if further convection: update parcel q,t (zpq,ztp)
!

           itop(jhor)=jlev
           zptop(jhor)=dp(jhor)*sigma(jlev)
           ztp(jhor,jlev)=ztnew
           zqp(jhor,jlev)=zqnew

          endif
         endif
        endif
       enddo
      enddo

!
!     set flags for convetion types
!

      where(ilift(:) > 0 .and. zptop(:) <= zpdeep) ideep(:)=1
      where(ilift(:) > 0 .and. zptop(:) > zpdeep)  ishallow(:)=1

!
!*    4) deep convection
!
!
!     construct reference profiles
!
!     4a) temperature
!

      do jlev=NLEV,2,-1
       jlep=jlev+1
       jlem=jlev-1
       do jhor=1,NHOR
        if(ideep(jhor) == 1) then
         if(ilift(jhor) == jlev) then

!
!     set reference temperature at cloud base
!

          ztref(jhor,jlev)=zte(jhor,jlev)
          ztp(jhor,jlev)=zte(jhor,jlev)

         elseif(ilift(jhor) > jlev .and. jlev >= itop(jhor)) then

!
!     make 1st guess reference temperature profile
!

          if(jlev > ifreez(jhor)) then
           zfac=(sigma(jlev)/sigma(jlep))**akap
           ztref(jhor,jlev)=ztref(jhor,jlep)*zfac                       &
     &                     +zstab*(ztp(jhor,jlev)-ztp(jhor,jlep)*zfac)
           if(ztref(jhor,jlev) <= TMELT .or. jlev <= itop(jhor)) then
            ifreez(jhor)=jlev
            zpfreez(jhor)=dp(jhor)*sigma(jlev)+zeps
            zdtfreez(jhor)=ztref(jhor,jlev)-ztp(jhor,jlev)
           endif
          else
           zpres=dp(jhor)*sigma(jlev)
           zy=(zpfreez(jhor)-zpres)/(zpfreez(jhor)-zptop(jhor))
           ztref(jhor,jlev)=ztp(jhor,jlev)                              &
     &                     +zdtfreez(jhor)*(1.-zy*zy)
          endif
         endif
        endif
       enddo
      enddo

!
!     4b) moisture
!

      do jlev=NLEV,2,-1
       do jhor=1,NHOR
        if(ideep(jhor) == 1) then
         if(ilift(jhor) >= jlev .and. jlev >= itop(jhor)) then
          zpres=dp(jhor)*sigma(jlev)
          if(jlev > ifreez(jhor)) then
           zdpres(jhor,jlev)=((zplift(jhor)-zpres)*zdpfreez             &
     &                       +(zpres-zpfreez(jhor))*zdpbase)            &
     &                      /(zplift(jhor)-zpfreez(jhor))
          else
           zdpres(jhor,jlev)=((zpfreez(jhor)-zpres)*zdptop              &
     &                       +(zpres-zptop(jhor))*zdpfreez)             &
     &                      /(zpfreez(jhor)-zptop(jhor))
          endif
          zfac=(1.+zdpres(jhor,jlev)/zpres)**akap
          ztrp=ztref(jhor,jlev)*zfac
          zqref(jhor,jlev)=RDBRV*RA1*exp(RA2*(ztrp-TMELT)/(ztrp-RA4))   &
     &                    /(zpres+zdpres(jhor,jlev))
          zcor=1./(1.-(1./RDBRV-1.)*zqref(jhor,jlev))
          zqref(jhor,jlev)=zqref(jhor,jlev)*zcor
          zdqrdt(jhor,jlev)=zfac*RA2*(TMELT-RA4)*zqref(jhor,jlev)*zcor                       &
     &                     /((ztrp-RA4)*(ztrp-RA4))
         endif
        endif
       enddo
      enddo

!
!     4c) energy conservation
!

      do jiter=1,2
       zener(:)=0.
       zdsig(:)=0.
       do jlev=NLEV,2,-1
        do jhor=1,NHOR
         if(ideep(jhor) == 1) then
          if(ilift(jhor) >= jlev .and. jlev >= itop(jhor)) then
           zener(jhor)=zener(jhor)+dsigma(jlev)                         &
     &                *(zcp(jhor,jlev)*(ztref(jhor,jlev)-zte(jhor,jlev))&
     &                 +zl(jhor,jlev)*(zqref(jhor,jlev)-zqe(jhor,jlev)))
           zdsig(jhor)=zdsig(jhor)+dsigma(jlev)
          endif
         endif
        enddo
       enddo
       where(ideep(:)==1) zener(:)=zener(:)/zdsig(:)
       do jlev=NLEV,2,-1
        do jhor=1,NHOR
         if(ideep(jhor) == 1) then
          if(ilift(jhor) >= jlev .and. jlev >= itop(jhor)) then
           zpres=dp(jhor)*sigma(jlev)
           zfac=(1.+zdpres(jhor,jlev)/zpres)**akap
           zdener=zcp(jhor,jlev)+zl(jhor,jlev)*zdqrdt(jhor,jlev)
           ztref(jhor,jlev)=ztref(jhor,jlev)-zener(jhor)/zdener
           ztrp=ztref(jhor,jlev)*zfac
           zqref(jhor,jlev)=RDBRV*RA1*exp(RA2*(ztrp-TMELT)/(ztrp-RA4))  &
     &                     /(zpres+zdpres(jhor,jlev))
           zcor=1./(1.-(1./RDBRV-1.)*zqref(jhor,jlev))
           zqref(jhor,jlev)=zqref(jhor,jlev)*zcor
           zdqrdt(jhor,jlev)=zfac*RA2*(TMELT-RA4)*zqref(jhor,jlev)*zcor &
     &                      /((ztrp-RA4)*(ztrp-RA4))
          endif
         endif
        enddo
       enddo
      enddo

!
!     4d) make tendencies and precip.
!

      icf(:,:)=0.
      do jlev=NLEV,2,-1
       do jhor=1,NHOR
        if((ideep(jhor) == 1) .and.                                     &
     &  (ilift(jhor) >= jlev .and. jlev >= itop(jhor))) then
         zdtdt(jhor,jlev)=(ztref(jhor,jlev)-zte(jhor,jlev))/ztaud
         zdqdt(jhor,jlev)=(zqref(jhor,jlev)-zqe(jhor,jlev))/ztaud
         zprc(jhor)=zprc(jhor)-zdqdt(jhor,jlev)*dsigma(jlev)
         if(jlev == ilift(jhor)) then
          icf(jhor,jlev)=0
         else
          icf(jhor,jlev)=1
         endif
        endif
       enddo
      enddo

!
!     4e) switch points with negative precip to shallow convecton
!

      do jhor=1,NHOR
       if(ideep(jhor) > 0 .and. zprc(jhor) < 0.) then
        ishallow(jhor)=1
        ideep(jhor)=2
        zdtdt(jhor,:)=0.
        zdqdt(jhor,:)=0.
        ztref(jhor,:)=0.
        zqref(jhor,:)=0.
        icf(jhor,:)=0
       endif
      enddo

!
!     set new top level for switched points
!

      do jlev=NLEM,3,-1
       where(ideep(:) == 2 .and. sigma(jlev)*dp(:) > zpdeep)
        itop(:)=jlev
        zptop(:)=sigma(jlev)*dp(:)
       endwhere
      enddo
      where(ideep(:) == 2 .and. zptop(:) < zpdeep) ishallow(:)=2

!
!*    5) shallow convection
!

      if(nshallow == 1) then
       do jhor=1,NHOR
        if(ishallow(jhor) == 1) then
         jtop=itop(jhor)-2
         jbase=ilift(jhor)
         zzptop=dp(jhor)*sigma(jtop)
         zttop=zte(jhor,jtop)
         zqtop=zqe(jhor,jtop)
         zpbase=dp(jhor)*sigma(jbase)
         ztbase=zte(jhor,jbase)
         zqbase=zqe(jhor,jbase)

!
!        mean th, q and t for mixing line
!

         zqm=0.5*(zqtop+zqbase)
         zthm=0.5*(zttop*(100000./zzptop)**AKAP                         &
     &            +ztbase*(100000./zpbase)**AKAP)
         ztmix=zthm*(zzptop/100000.)**AKAP

!
!     saturation points
!
!     a) at profile base
!

         ze=zqbase*zpbase/(RDBRV+(1.-RDBRV)*zqbase)
         ztsp=55.+2840./(3.5*ALOG(ztbase)-ALOG(ze)-0.2)
         zpsbase=zpbase*(ztbase/ztsp)**AKAP
!
!     b) at profile top
!

         ze=zqtop*zzptop/(RDBRV+(1.-RDBRV)*zqtop)
         ztsp=55.+2840./(3.5*ALOG(ztmix)-ALOG(ze)-0.2)
         zpstop=zzptop*(ztmix/ztsp)**AKAP

!
!     slope of mixing line
!

         zdp=zpstop-zpsbase
         zslope(jhor)=(zthm-ztbase*(100000./zpbase)**AKAP)/zdp
        endif
       enddo
!
!     reference profile
!
       do jlev=NLEV,1,-1
        do jhor=1,NHOR
         if(ishallow(jhor) == 1) then
          if(jlev == ilift(jhor)) then
           zp=dp(jhor)*sigma(jlev)
           ztref(jhor,jlev)=zte(jhor,jlev)
           zspb(jhor)=zp+zdpsh
           zts=ztref(jhor,jlev)*(zspb(jhor)/zp)**AKAP
           zqref(jhor,jlev)=RDBRV*RA1*exp(RA2*(zts-TMELT)/(zts-RA4))    &
     &                     /zspb(jhor)
           zcor=1./(1.-(1./RDBRV-1.)*zqref(jhor,jlev))
           zqref(jhor,jlev)=zqref(jhor,jlev)*zcor
          endif
          if(jlev < ilift(jhor) .and. jlev >= itop(jhor)-1) then
           zp=dp(jhor)*sigma(jlev)
           zpl=dp(jhor)*sigma(jlev+1)
           ztref(jhor,jlev)=(ztref(jhor,jlev+1)*(100000./zpl)**akap     &
     &                      +zstab*zbeta*zslope(jhor)/(zp-zpl))         &
     &                     *(zp/100000.)**AKAP


           zsp=zspb(jhor)+zbeta*(zp-zspb(jhor)-zdpsh)
           zts=ztref(jhor,jlev)*(zsp/zp)**AKAP
           zqref(jhor,jlev)=RDBRV*RA1*exp(RA2*(zts-TMELT)/(zts-RA4))    &
     &                     /zsp
           zcor=1./(1.-(1./RDBRV-1.)*zqref(jhor,jlev))
           zqref(jhor,jlev)=zqref(jhor,jlev)*zcor
          endif
         endif
        enddo
       enddo

!
!     energy conservation
!

       zsumt(:)=0.
       zsumq(:)=0.
       zdsig(:)=0.
       do jlev=NLEV,1,-1
        do jhor=1,NHOR
         if(ishallow(jhor) == 1) then
          if(jlev <= ilift(jhor) .and. jlev >= itop(jhor)-1) then
           zsumt(jhor)=zsumt(jhor)                                      &
     &                +(ztref(jhor,jlev)-zte(jhor,jlev))*dsigma(jlev)
           zsumq(jhor)=zsumq(jhor)                                      &
     &                +(zqref(jhor,jlev)-zqe(jhor,jlev))*dsigma(jlev)
           zdsig(jhor)=zdsig(jhor)+dsigma(jlev)
          endif
         endif
        enddo
       enddo
       do jlev=NLEV,1,-1
        do jhor=1,NHOR
         if(ishallow(jhor) == 1) then
          if(jlev <= ilift(jhor) .and. jlev >= itop(jhor)-1) then
           ztref(jhor,jlev)=ztref(jhor,jlev)-zsumt(jhor)/zdsig(jhor)
           zqref(jhor,jlev)=zqref(jhor,jlev)-zsumq(jhor)/zdsig(jhor)
           if(zqref(jhor,jlev) < 0.) ishallow(jhor)=3
          endif
         endif
        enddo
       enddo

!
!      add tendencies
!

       do jlev=1,NLEV
        do jhor=1,NHOR
         if(ishallow(jhor) == 1) then
          if(jlev <= ilift(jhor) .and. jlev >= itop(jhor)-1) then
           zdtdt(jhor,jlev)=(ztref(jhor,jlev)-zte(jhor,jlev))/ztaus
           zdqdt(jhor,jlev)=(zqref(jhor,jlev)-zqe(jhor,jlev))/ztaus
           if(jlev == ilift(jhor)) then
            icf(jhor,jlev)=0
           else
            icf(jhor,jlev)=1
           endif
          endif
         endif
        enddo
       enddo

      endif

!
!*    6) make final dqdt, dtdt and precip.
!

!
!     dqdt and dtdt
!

      dqdt(:,1:NLEV)=zdqdt(:,:)+dqdt(:,1:NLEV)
      dtdt(:,1:NLEV)=zdtdt(:,:)+dtdt(:,1:NLEV)

!
!     precip. in m/s and snow fall
!

      dprc(:)=AMAX1(zprc(:),0.)*dp(:)/GA/1000.
      where(dt(:,NLEV) <= TMELT) dprs(:)=dprs(:)+dprc(:)

!
!     flag and counter for clouds:
!

      do jlev=1,NLEV
       where(icf(:,jlev) > 0)
        icclev(:,jlev)=1
        icctot(:)=icctot(:)+1
       endwhere
       if(sigma(jlev) < 0.4) then
        where(icclev(:,jlev) == 1 .and. itop(:) == jlev)
         icclev(:,jlev)=2
        endwhere
       endif
      enddo

!!
!
!*    7) momentum mixing
!

      if(nmoment == 1) then

!
!     7a) update u and v
!

       zue(:,:)=du(:,1:NLEV)+dudt(:,1:NLEV)*deltsec2
       zve(:,:)=dv(:,1:NLEV)+dvdt(:,1:NLEV)*deltsec2

!
!     7b) reference u and v
!

       zup(:)=0.
       zvp(:)=0.
       zsums(:)=0.
       do jlev=1,NLEV
        where(ilift(:) >= jlev .and. itop(:) <= jlev                    &
     &       .and. ideep(:) == 1)
         zup(:)=zup(:)+zue(:,jlev)*dsigma(jlev)
         zvp(:)=zvp(:)+zve(:,jlev)*dsigma(jlev)
         zsums(:)=zsums(:)+dsigma(jlev)
        endwhere
        if(nshallow==1) then
         where(ilift(:) >= jlev .and. jlev >= itop(:)-1                 &
     &        .and. ishallow(:) == 1)
          zup(:)=zup(:)+zue(:,jlev)*dsigma(jlev)
          zvp(:)=zvp(:)+zve(:,jlev)*dsigma(jlev)
          zsums(:)=zsums(:)+dsigma(jlev)
         endwhere
        endif
       enddo
       where(zsums(:) > 0.)
        zup(:)=zup(:)/zsums(:)
        zvp(:)=zvp(:)/zsums(:)
       endwhere

!
!     7c) make and add tendencies
!

       do jlev=1,NLEV
        where(ilift(:) >= jlev .and. itop(:) <= jlev                    &
     &       .and. ideep(:) == 1)
         zdudt(:,jlev)=(zup(:)-zue(:,jlev))/ztaud
         zdvdt(:,jlev)=(zvp(:)-zve(:,jlev))/ztaud
         dudt(:,jlev)=dudt(:,jlev)+zdudt(:,jlev)
         dvdt(:,jlev)=dvdt(:,jlev)+zdvdt(:,jlev)
        endwhere
        if(nshallow==1) then
         where(ilift(:) >= jlev .and. itop(:) <= jlev                   &
     &        .and. ishallow(:) == 1)
          zdudt(:,jlev)=(zup(:)-zue(:,jlev))/ztaus
          zdvdt(:,jlev)=(zvp(:)-zve(:,jlev))/ztaus
          dudt(:,jlev)=dudt(:,jlev)+zdudt(:,jlev)
          dvdt(:,jlev)=dvdt(:,jlev)+zdvdt(:,jlev)
         endwhere
        endif
       enddo
!!!
!     test printout
!
!       allocate(zprf1(NLON*NLAT,NLEV))
!!       allocate(zprf2(NLON*NLAT,NLEV))
!       call mpgagp(zprf1,zdudt,NLEV)
!       call mpgagp(zprf2,zdvdt,NLEV)
!       if(mypid==NROOT) then
!        write(nud,*)'In momentum:'
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
        write(nud,*)'In MKBM:'
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
       zprf7(:)=real(ideep(:))
       zprf8(:)=real(ishallow(:))
       call mpgagp(zprf3,zprf7,1)
       call mpgagp(zprf4,zprf8,1)
       if(mypid==NROOT) then
        write(nud,*)'In MKBM:'
        write(nud,*)'ideep= ',zprf3(nprhor),' ishallow= ',zprf4(nprhor)
       endif
       zprf7(:)=real(ilift(:))
       zprf8(:)=real(itop(:))
       call mpgagp(zprf3,zprf7,1)
       call mpgagp(zprf4,zprf8,1)
       call mpgagp(zprf1,zdqdt,NLEV)
       call mpgagp(zprf2,zdtdt,NLEV)
       call mpgagp(zprf6,dprc,1)
       if(mypid==NROOT) then
        write(nud,*)'ilift= ',zprf3(nprhor),' itop= ',zprf4(nprhor)
        do jlev=1,NLEV
          zzdt=zprf2(nprhor,jlev)/deltsec2
          zzdq=zprf1(nprhor,jlev)/deltsec2
          write(nud,*)'L= ',jlev,' dtdt= ',zzdt,' dqdt= ',zzdq
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

      return
      end subroutine mkbm

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
