      module rainmod
      use pumamod
!
!     version identifier (date)
!
      character(len=80) :: version = '22.11.2002 by Larry'
!
!     namelist parameter:
!

      integer :: kbetta = 1  ! switch for betta in kuo (1/0=yes/no)
      integer :: nprl   = 1  ! switch for large scale precip (1/0=yes/no)
      integer :: nprc   = 1  ! switch for large convective precip (1/0=yes/no)
      integer :: ndca   = 1  ! switch for dry convective adjustment (1/0=yes/no)
      integer :: ncsurf = 1  ! switch for conv. starts from surface (1/0=yes/no)
      integer :: nmoment= 0  ! switch for momentum mixing (1/0=yes/no)

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
      namelist/rainmod_nl/kbetta,nprl,nprc,ndca,ncsurf,nmoment             &
     &                ,rcrit,clwcrit1,clwcrit2
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
      call mpbci(kbetta)
      call mpbci(nprl)
      call mpbci(nprc)
      call mpbci(ndca)
      call mpbci(ncsurf)
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
!**   1) transfere adiabatic q-tendencies to grid point space
!        (needed for kuo scheme)
!

      if(nprc == 1) call mkdqtgp

!
!**   2) kuo convection
!

      if(ntime ==1) call mksecond(zsec1,0.)
      if(nprc == 1) call kuo
      if(ntime ==1) then
       call mksecond(zsec1,zsec1)
       time4prc=time4prc+zsec1
      endif

!
!**   3) dry convective adjustment
!

      if(ntime ==1) call mksecond(zsec1,0.)
      if(ndca == 1) call mkdca
      if(ntime ==1) then
       call mksecond(zsec1,zsec1)
       time4dca=time4dca+zsec1
      endif

!
!**   4) large scale precipitation
!

      if(ntime ==1) call mksecond(zsec1,0.)
      if(nprl == 1) call mklsp
      if(ntime ==1) then
       call mksecond(zsec1,zsec1)
       time4prl=time4prl+zsec1
      endif

!
!**   5) diagnose cloud cover
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
!     allocatable arrys for dbug diagnostic
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
         write(nud,*)'L= ',jlev,' prl (mm/day)= ',zprf6(nprhor)              &
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

      return
      end subroutine mklsp

!     ==============
!     SUBROUTINE KUO
!     ==============

      subroutine kuo
      use rainmod

!
!     local arrays
!

      integer ilift(NHOR)    ! lifting level (= cloud base + 1)
      integer itop(NHOR)     ! cloud top level
      real zqe(NHOR,NLEV)    ! moisture of environmental air
      real zte(NHOR,NLEV)    ! temperatur of environmental air
      real zdqtot(NHOR,NLEV) ! total moisture accession in resp. layer
      real zdqdt(NHOR,NLEV)  ! moisture diff. (parcel-environment)
      real zdtdt(NHOR,NLEV)  ! temperature diff. (parcel-environment)
      real zalcpe(NHOR,NLEV) ! help array for factors
      real zqp(NHOR)         ! humidity of lifted parcel
      real ztp(NHOR)         ! temperature of lifted parcel
      real zi(NHOR)          ! moisture supply in the collumn
      real zpt(NHOR)         ! energy surplus of the cloud (t-part)
      real zpq(NHOR)         ! energy surplus of the cloud (q-part)
      real zbetta(NHOR)      ! betta factor (new kuo scheme)
      real zpbetta(NHOR)     ! p-scaling for betta factor
      real zat(NHOR)         ! distribution factor (temperature)
      real zaq(NHOR)         ! distribution factor (moisture)
!
!     array for momentum mixing
!
      real zue(NHOR,NLEV)    ! u of environmental air
      real zve(NHOR,NLEV)    ! v of environmental air
      real zdudt(NHOR,NLEV)  ! u-momentum difference (cloud-environment)
      real zdvdt(NHOR,NLEV)  ! v-momentum difference (cloud-environment)
      real zau(NHOR)         ! distribution factor (momentum)
      real zup(NHOR)         ! reference (cloud) u
      real zvp(NHOR)         ! reference (cloud) v
      real zsums(NHOR)       ! sum(dsigma) for mixed levels

!
!     allocatable arrays for dbug diagnostic
!

      real, allocatable :: zprf1(:,:)
      real, allocatable :: zprf2(:,:)
      real, allocatable :: zprf3(:,:)
      real, allocatable :: zprf4(:,:)
      real, allocatable :: zprf5(:)
      real, allocatable :: zprf6(:)
      real, allocatable :: zprf7(:)
      real, allocatable :: zprf8(:)
      real, allocatable :: zprf9(:)
      real, allocatable :: zprf10(:)
      real, allocatable :: zprf11(:)
      real, allocatable :: zprf12(:)

!
!*    0) initialize some fields
!

      ilift(:)=-1
      itop(:)=-1
      icctot(:)=0
      zi(:)=0.
      zpt(:)=0.
      zpq(:)=0.
      zat(:)=0.
      zaq(:)=0.
      zbetta(:)=0.
      zpbetta(:)=0.
      zdtdt(:,:)=0.
      zdqdt(:,:)=0.
      icclev(:,:)=0

!
!     moisture accession (level and column)
!

      zdqtot(:,:)=(dqt(:,1:NLEV)+dqdt(:,1:NLEV))*deltsec2

!
!     for original kuo (whole column) use the next loop
!
!     do jlev=1,NLEV
!      zi(:)=zi(:)+zdqtot(:,jlev)*dsigma(jlev)
!     enddo

!
!     franks dbug
!

      if(ndiaggp==1) then
       dgp3d(:,1:NLEV,16)=-dqdt(:,1:NLEV)
       dgp3d(:,1:NLEV,18)=dqt(:,1:NLEV)
      end if

!
!*    1) update q and t
!

      zqe(:,:)=dq(:,1:NLEV)+dqdt(:,1:NLEV)*deltsec2
      zte(:,:)=dt(:,1:NLEV)+dtdt(:,1:NLEV)*deltsec2

!
!     set L/CP depending on t
!

      where(zte(:,:) < TMELT)
       zalcpe(:,:)=ALS/(ACPD*(1.+ADV*zqe(:,:)))
      elsewhere
       zalcpe(:,:)=ALV/(ACPD*(1.+ADV*zqe(:,:)))
      endwhere

!
!*    2) search for lifting level
!

      if(ncsurf == 1) then

!
!     a) start at surface (lifting level only if dry unstable)
!

      zfac=sigma(NLEV)**akap
      do jhor=1,NHOR
       ztnew=dt(jhor,NLEP)*zfac
       if(zdqtot(jhor,NLEV) > 0. .and. ztnew > zte(jhor,NLEV)) then

!
!     surface air specific humidity
!

        zqsnl=RDBRV*RA1*exp(RA2*(zte(jhor,NLEV)-TMELT)                  &
     &                          /(zte(jhor,NLEV)-RA4))                  &
     &       /(dp(jhor)*sigma(NLEV))
        zqsnl=zqsnl/(1.-(1./RDBRV-1.)*zqsnl)
        zqnew=dq(jhor,NLEP)*AMIN1(1.,zqe(jhor,NLEV)/zqsnl)

!
!     saturation humidity:
!

        zqsnl=RDBRV*RA1*exp(RA2*(ztnew-TMELT)/(ztnew-RA4))              &
     &       /(dp(jhor)*sigma(NLEV))
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

         zqsnl=RDBRV*RA1*exp(RA2*(ztnew-TMELT)/(ztnew-RA4))             &
     &        /(dp(jhor)*sigma(NLEV))
         zcor=1./(1.-(1./RDBRV-1.)*zqsnl)
         zqsnl=zqsnl*zcor
         zdqsdt=RA2*(TMELT-RA4)*zqsnl*zcor/((ztnew-RA4)*(ztnew-RA4))
         zrain=(zqnew-zqsnl)/(1.+zlcp*zdqsdt)
         ztnew=ztnew+zlcp*zrain
         zqnew=zqnew-zrain

!
!     if convection: set lifting level and parcel q,t (zpq,ztp)
!                    calculate prelim. dq,dt (zdqdt,zdtdt)
!                    update surplus of total energy (zp)
!

         ilift(jhor)=NLEP
         itop(jhor)=NLEV
         ztp(jhor)=ztnew
         zqp(jhor)=zqnew
         zqsa=RDBRV*RA1*exp(RA2*(zte(jhor,NLEV)-TMELT)                  &
     &                          /(zte(jhor,NLEV)-RA4))                  &
     &       /(dp(jhor)*sigma(NLEV))
         zqsa=zqsa/(1.-(1./RDBRV-1.)*zqsa)
         zdqdt(jhor,NLEV)=AMAX1(zqsa-zqe(jhor,NLEV),0.)
         zdtdt(jhor,NLEV)=ztp(jhor)-zte(jhor,NLEV)
         zpt(jhor)=zdtdt(jhor,NLEV)*dsigma(NLEV)/zalcpe(jhor,NLEV)
         zpq(jhor)=zdqdt(jhor,NLEV)*dsigma(NLEV)
         if(kbetta > 0 ) then
          zbetta(jhor)=AMIN1(1.,zqe(jhor,NLEV)/zqsa)*dsigma(NLEV)
          zpbetta(jhor)=dsigma(NLEV)
         endif

!
!     moisture accession (delete for original kuo)
!

         zi(jhor)=zdqtot(jhor,NLEV)*dsigma(NLEV)
        endif
       endif
      enddo

      endif

!
!     b) continue search for free atmosphere
!
      do jlev=NLEM,1,-1
       jlep=jlev+1
       zfac=(sigma(jlev)/sigma(jlep))**akap
       do jhor=1,NHOR
        if(ilift(jhor) == -1                                            &
!    &    .and. zi(jhor) > 0.                   ! use for original kuo
     &    .and. zdqtot(jhor,jlep) > 0.) then

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
!                    calculate prelim. dq,dt (zdqdt,zdtdt)
!                    update surplus of total energy (zp)
!

           ilift(jhor)=jlep
           itop(jhor)=jlev
           ztp(jhor)=ztnew
           zqp(jhor)=zqnew
           zqsa=RDBRV*RA1*exp(RA2*(zte(jhor,jlev)-TMELT)                &
     &                           /(zte(jhor,jlev)-RA4))                 &
     &         /(dp(jhor)*sigma(jlev))
           zqsa=zqsa/(1.-(1./RDBRV-1.)*zqsa)
           zdqdt(jhor,jlev)=AMAX1(zqsa-zqe(jhor,jlev),0.)
           zdtdt(jhor,jlev)=ztp(jhor)-zte(jhor,jlev)
           zpt(jhor)=zdtdt(jhor,jlev)*dsigma(jlev)/zalcpe(jhor,jlev)
           zpq(jhor)=zdqdt(jhor,jlev)*dsigma(jlev)
           if(kbetta > 0 ) then
            zbetta(jhor)=AMIN1(1.,zqe(jhor,jlev)/zqsa)*dsigma(jlev)
            zpbetta(jhor)=dsigma(jlev)
           endif

!
!          check if lifting level is supersaturated
!

           ztnew=zte(jhor,jlep)
           zqnew=zqe(jhor,jlep)
           zqsnl=RDBRV*RA1*exp(RA2*(zte(jhor,jlep)-TMELT)               &
     &                            /(zte(jhor,jlep)-RA4))                &
     &         /(dp(jhor)*sigma(jlep))
           zqsnl=zqsnl/(1.-(1./RDBRV-1.)*zqsnl)
           if(zqe(jhor,jlep) > zqsnl) then

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
     &          /(dp(jhor)*sigma(jlep))
           zcor=1./(1.-(1./RDBRV-1.)*zqsnl)
           zqsnl=zqsnl*zcor
           zdqsdt=RA2*(TMELT-RA4)*zqsnl*zcor/((ztnew-RA4)*(ztnew-RA4))
           zrain=(zqnew-zqsnl)/(1.+zlcp*zdqsdt)
           ztnew=ztnew+zlcp*zrain
           zqnew=zqnew-zrain
           zdqdt(jhor,jlep)=AMAX1(zqnew-zqe(jhor,jlep),0.)
           zdtdt(jhor,jlep)=ztnew-zte(jhor,jlep)
           zpt(jhor)=zdtdt(jhor,jlep)*dsigma(jlep)/zalcpe(jhor,jlep)    &
     &              +zpt(jhor)
           zpq(jhor)=zdqdt(jhor,jlep)*dsigma(jlep)                      &
     &              +zpq(jhor)
           if(kbetta > 0 ) then
            zbetta(jhor)=AMIN1(1.,zqe(jhor,jlep)/zqsnl)*dsigma(jlep)    &
     &                  +zbetta(jhor)
            zpbetta(jhor)=dsigma(jlep)+zpbetta(jhor)
           endif
           endif

!
!     moisture accession (delete for original kuo)
!

           zi(jhor)=zdqtot(jhor,jlev)*dsigma(jlev)                      &
     &             +zdqtot(jhor,jlep)*dsigma(jlep)
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

         ztnew=ztp(jhor)*zfac
         zqnew=zqp(jhor)

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
!                                calculate prelim. dq,dt (zdqdt,zdtdt)
!                                update surplus of total energy (zp)
!

           itop(jhor)=jlev
           ztp(jhor)=ztnew
           zqp(jhor)=zqnew
           zqsa=RDBRV*RA1*exp(RA2*(zte(jhor,jlev)-TMELT)                &
     &                            /(zte(jhor,jlev)-RA4))                &
     &         /(dp(jhor)*sigma(jlev))
           zqsa=zqsa/(1.-(1./RDBRV-1.)*zqsa)
           zdqdt(jhor,jlev)=AMAX1(zqsa-zqe(jhor,jlev),0.)
           zdtdt(jhor,jlev)=ztp(jhor)-zte(jhor,jlev)
           zpt(jhor)=zdtdt(jhor,jlev)*dsigma(jlev)/zalcpe(jhor,jlev)    &
     &              +zpt(jhor)
           zpq(jhor)=zdqdt(jhor,jlev)*dsigma(jlev)                      &
     &              +zpq(jhor)
           if(kbetta > 0 ) then
            zbetta(jhor)=AMIN1(1.,zqe(jhor,jlev)/zqsa)*dsigma(jlev)     &
     &                  +zbetta(jhor)
            zpbetta(jhor)=zpbetta(jhor)+dsigma(jlev)
           endif

!
!     moisture accession (delete for original kuo)
!

           zi(jhor)=zdqtot(jhor,jlev)*dsigma(jlev)                      &
     &             +zi(jhor)
          endif
         endif
        endif
       enddo
      enddo

!
!*    4) calc. moisture distribution factors
!

      zi(:)=AMAX1(0.,zi(:))

      if(kbetta==0) then

!
!     kuo without betta
!

       where(zpt(:)+zpq(:) > 0.)
        zat(:)=zi(:)/(zpt(:)+zpq(:))
        zaq(:)=zat(:)
       end where
      else

!
!     kuo with betta
!

       where(zpbetta(:) > 0.) zbetta(:)=(1.-zbetta(:)/zpbetta(:))**3
       where(zpt(:) > 0.) zat(:)=zi(:)*(1.-zbetta(:))/zpt(:)
       where(zpq(:) > 0.) zaq(:)=zi(:)*zbetta(:)/zpq(:)

      endif

!
!*    5) update dqdt and dtdt, compute precip.
!

      do jlev=1,NLEV

!
!     5a) in cloud: add tendencies due to convection and sub. advection
!         lifting level: sub. advection. (dtdt != 0 only if level supersat)
!

       where(ilift(:) >= jlev .and. itop(:) <= jlev .and. zi(:) > 0.)
         dqdt(:,jlev)=zaq(:)*zdqdt(:,jlev)/deltsec2                     &
     &               -dqt(:,jlev)
         dtdt(:,jlev)=zat(:)*zdtdt(:,jlev)/deltsec2                     &
     &               +dtdt(:,jlev)

!
!     5b) calc. precipitation
!

         dprc(:)=dprc(:)                                                &
     &          +zat(:)*zdtdt(:,jlev)*dsigma(jlev)/zalcpe(:,jlev)

       endwhere

!
!     5c) original kuo:
!         rest of column: subtrac all q-tendencies (i.e. reset dqdt)
!
!      where(zi(:) > 0 .and. (jlev >= ilift(:) .or. jlev < itop(:)))
!       dqdt(:,jlev)=-1.*dqt(:,jlev)
!      endwhere

!
!     5d) flag and counter for clouds:
!

       where(zat(:)*zdtdt(:,jlev) > 0.)
        icclev(:,jlev)=1
        icctot(:)=icctot(:)+1
       endwhere
       if(sigma(jlev) < 0.4) then
        where(icclev(:,jlev)==1 .and. itop(:)==jlev)
         icclev(:,jlev)=2
        endwhere
       endif

      enddo

!
!*    6) final precip. in m/s and snow fall
!

      zfac1=1./GA/deltsec2/1000.
      dprc(:)=dprc(:)*dp(:)*zfac1
      where(dt(:,NLEV) <= TMELT) dprs(:)=dprs(:)+dprc(:)

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
!     7a) make distribution factor
!

       zau(:)=0.
       where(zpt(:)+zpq(:) > 0.)
        zau(:)=AMIN1(1.,zi(:)/(zpt(:)+zpq(:)))
       end where

!
!     7b) reference u and v
!

       zup(:)=0.
       zvp(:)=0.
       zsums(:)=0.
       do jlev=1,NLEV
        where(ilift(:) >= jlev .and. itop(:) <= jlev .and. zi(:) > 0.)
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
!     7c) make and add tendencies
!

       do jlev=1,NLEV
        where(ilift(:) >= jlev .and. itop(:) <= jlev .and. zi(:) > 0.)
          zdudt(:,jlev)=zau(:)*(zup(:)-zue(:,jlev))
          zdvdt(:,jlev)=zau(:)*(zvp(:)-zve(:,jlev))
          dudt(:,jlev)=dudt(:,jlev)+zdudt(:,jlev)/deltsec2
          dvdt(:,jlev)=dvdt(:,jlev)+zdvdt(:,jlev)/deltsec2
        endwhere
       enddo

      endif

!
!*    dbug diagnostic for nprint==2 (see pumamod)
!

      if(nprint==2) then
       allocate(zprf1(NLON*NLAT,NLEV))
       allocate(zprf2(NLON*NLAT,NLEV))
       allocate(zprf3(NLON*NLAT,NLEV))
       allocate(zprf4(NLON*NLAT,NLEV))
       allocate(zprf5(NLON*NLAT))
       allocate(zprf6(NLON*NLAT))
       allocate(zprf7(NLON*NLAT))
       allocate(zprf8(NLON*NLAT))
       allocate(zprf9(NLON*NLAT))
       allocate(zprf10(NLON*NLAT))
       allocate(zprf11(NHOR))
       allocate(zprf12(NHOR))
       zprf11(:)=real(ilift(:))
       zprf12(:)=real(itop(:))
       call mpgagp(zprf1,zdqdt,NLEV)
       call mpgagp(zprf2,zdtdt,NLEV)
       call mpgagp(zprf3,zdqtot,NLEV)
       call mpgagp(zprf4,dqt,NLEV)
       call mpgagp(zprf5,zprf11,1)
       call mpgagp(zprf6,zprf12,1)
       call mpgagp(zprf7,dprc,1)
       call mpgagp(zprf8,zaq,1)
       call mpgagp(zprf9,zat,1)
       call mpgagp(zprf10,zi,1)
       deallocate(zprf11)
       deallocate(zprf12)
       if(mypid==NROOT) then
        write(nud,*)'In Kuo:'
        write(nud,*)'ilift= ',zprf5(nprhor),' itop= ',zprf6(nprhor)
        do jlev=1,NLEV
         write(nud,*)'L= ',jlev,' dqtot= ',zprf3(nprhor,jlev)                &
     &                    ,' dqt= ',zprf4(nprhor,jlev)
         if(zprf5(nprhor) >= jlev .and. zprf6(nprhor) <= jlev           &
     &     .and. zprf10(nprhor) > 0.) then
          zzdt=zprf9(nprhor)*zprf2(nprhor,jlev)/deltsec2
          zzdq=zprf8(nprhor)*zprf1(nprhor,jlev)/deltsec2
          write(nud,*)'L= ',jlev,' dtdt= ',zzdt,' dqdt= ',zzdq
         endif
        enddo
        write(nud,*)'prc (mm/day)= ',zprf7(nprhor)*1000.*ntspd*deltsec2*0.5
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
       deallocate(zprf10)
      endif

!
!     franks dbug
!

      if(ndiaggp==1) then
       do jlev=1,NLEV
        where(ilift(:) >= jlev .and. itop(:) <= jlev .and. zi(:) > 0.)
         dgp3d(:,jlev,8)=zat(:)*zdtdt(:,jlev)/deltsec2
         dgp3d(:,jlev,16)=zaq(:)*zdqdt(:,jlev)/deltsec2                 &
     &                   -dqt(:,jlev)
        elsewhere
         dgp3d(:,jlev,8)=0.
         dgp3d(:,jlev,16)=0.
        endwhere
       enddo
      end if

      return
      end subroutine kuo

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
