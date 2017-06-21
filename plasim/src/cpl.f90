      module cplmod
!
      parameter(nxa=64,nya=32)
      parameter(nxo=72,nyot=76,nyo=nyot-4)
      parameter(radea=6.371E6)
      parameter(nud=6)           ! output unit
!
      character(len=30) :: version='cplmod 20/06/08'
!
      real :: xa(nxa,nya)        ! atm. longitudes
      real :: xa1(nxa,nya)       ! atm. longitudes western boundary 
      real :: xa2(nxa,nya)       ! atm. longitudes eastern boundary
      real :: ya(nxa,nya)        ! atm. latitudes 
      real :: ya1(nxa,nya)       ! atm. latitudes boundary 1
      real :: ya2(nxa,nya)       ! atm. latitudes boundary 2
      real :: xos(nxo,nyo)       ! oce. longitudes (scalar)
      real :: xos1(nxo,nyo)      ! oce. longitudes western boundary (scalar)
      real :: xos2(nxo,nyo)      ! oce. longitudes eastern boundary (scalar)
      real :: xov(nxo,nyo)       ! oce. longitudes (vector)
      real :: xov1(nxo,nyo)      ! oce. longitudes western boundary (vector)
      real :: xov2(nxo,nyo)      ! oce. longitudes eastern boundary (vector)
      real :: yo(nxo,nyo)        ! oce. latitudes 
      real :: yo1(nxo,nyo)       ! oce. latitudes boundary 1
      real :: yo2(nxo,nyo)       ! oce. latitudes boundary 2
!
      real :: agw(nxa,nya)       ! atm. area weights
      real :: ogw(nxo,nyo)       ! oce. area weights
!
      real :: aslm(nxa,nya)   = 0. ! atm. land(0)-see(1)-mask
      real :: oslms(nxo,nyo)  = 0. ! oce. land(0)-see(1)-mask (scalar)
      real :: oslmv(nxo,nyo)  = 0. ! oce. land(0)-see(1)-mask (vector)

      real :: asst(nxa,nya)   = 0. ! atm. sst
      real :: aiced(nxa,nya)  = 0. ! atm. ice thickness
      real :: apme(nxa,nya)   = 0. ! atm. fresh water flux
      real :: ataux(nxa,nya)  = 0. ! atm. u-stress
      real :: atauy(nxa,nya)  = 0. ! atm. v-stress
      real :: ahfl(nxa,nya)   = 0. ! atm. heat flux
!
      real :: osst(nxo,nyo)   = 0. ! oce. sst
      real :: oiced(nxo,nyo)  = 0. ! oce. ice thickness
      real :: opme(nxo,nyo)   = 0. ! oce. fresh water flux
      real :: otaux(nxo,nyo)  = 0. ! oce. u-stress
      real :: otauy(nxo,nyo)  = 0. ! oce. v-stress
      real :: ohfl(nxo,nyo)   = 0. ! oce. heat flux
!
      real :: assta(nxa,nya)  = 0. ! atm. sst (accumulated)
      real :: aiceda(nxa,nya) = 0. ! atm. ice thickness (accumulated)
      real :: apmea(nxa,nya)  = 0. ! atm. frsh water flux (accumulated)
      real :: atauxa(nxa,nya) = 0. ! atm. u-stress (accumulated)
      real :: atauya(nxa,nya) = 0. ! atm. v-stress (accumulated)
      real :: ahfla(nxa,nya)  = 0. ! atm. heat flux (accumulated)
!
      real :: ossta(nxo,nyo)  = 0. ! oce. sst (accumulated)
      real :: oiceda(nxo,nyo) = 0. ! oce. ice thickness(accumulated)
      real :: opmea(nxo,nyo)  = 0. ! oce. fresh water flux (accumulated)
      real :: otauxa(nxo,nyo) = 0. ! oce. u-stress (accumulated)
      real :: otauya(nxo,nyo) = 0. ! oce. v-stress (accumulated)
      real :: ohfla(nxo,nyo)  = 0. ! oce. heat flux (accumulated)
!
      real :: flukosst(nxo,nyot,12)  = 0. ! sst flux correction
      real :: flukotaux(nxo,nyot,12) = 0. ! u-stress flux correction
      real :: flukotauy(nxo,nyot,12) = 0. ! v-stress flux correction
      real :: flukofresh(nxo,nyot,12)= 0. ! frsh water flux flux correction
      real :: flukoice(nxo,nyot,12)  = 0. ! ice flux correction
      real :: flukooheat(nxo,nyot,12)= 0. ! heat flux flux correction
!
      integer :: iroffi(nxa,nya) = 0  ! runoff mask input (nroffc=2,3)
      integer :: iroffo(nxa,nya) = 0  ! runoff mask output (nroffc=3)
!  
      integer :: ncoupling = 1 ! switch for coupling method 
      integer :: nfluko    = 0 ! switch for flux correction
      integer :: ngui      = 0 ! switch for gui
      integer :: nroffc    = 0 ! switch for runoff correction
!
      integer :: nssta = 0 ! counter for atm. sst accumulation
      integer :: nicea = 0 ! counter for atm. ice accumulation
      integer :: npmea = 0 ! counter for atm. fresh water flux accumulation
      integer :: ntaua = 0 ! counter for atm. stress accumulation
      integer :: nhfla = 0 ! counter for atm. heat flux accumulation
      integer :: nssto = 0 ! counter for oce. sst accumulation
      integer :: niceo = 0 ! counter for oce. ice accumulation
      integer :: npmeo = 0 ! counter for oce. fresh water flux accumulation
      integer :: ntauo = 0 ! counter for oce. stress accumulation 
      integer :: nhflo = 0 ! counter for oce. heat flux accumulation
!
      integer :: naomod = 0 ! atm-oce coupling interval
!
      integer :: ndatim(7) = 0 ! array containing time/date information 
      integer :: nyear     = 0 ! actual year
      integer :: nmonth    = 0 ! actual month
      integer :: nday      = 0 ! actual day
      integer :: nhour     = 0 ! actual hour
      integer :: nmin      = 0 ! actual minute
      integer :: nwday     = 0 ! actual day of the week
!
!     namelist parameter
!
      integer :: nprint    = 0 ! extended diagnostics
      integer :: ndbox     = 1 ! ocean diagnostic gp (i)
      integer :: ndboy     = 1 ! ocean diagnostic gp (j)
      integer :: nflsst    = 0 ! sst fluko
      integer :: nfltau    = 0 ! stress fluko
      integer :: nflice    = 0 ! ice fluko
      integer :: nflfresh  = 0 ! pme fluko
      integer :: nfloheat  = 0 ! heat flux (deep ocean) fluko
      integer :: nissta    = 0 ! choise for ssta interpolation
      integer :: niicea    = 0 ! choise for icea interpolation
      integer :: nipmea    = 0 ! choise for pmea interpolation
      integer :: nitaua    = 0 ! choise for taua interpolation
      integer :: nihfla    = 0 ! choise for hfla interpolation
      integer :: nissto    = 1 ! choise for ssto interpolation
      integer :: nihflo    = 1 ! choise for hflo interpolation
      integer :: ncssta    = 3 ! choise for ssta correction
      integer :: ncicea    = 1 ! choise for icea correction
      integer :: ncpmea    = 3 ! choise for pmea correction
      integer :: nctaua    = 3 ! choise for taua correction
      integer :: nchfla    = 3 ! choise for hfla correction
      integer :: ncssto    = 3 ! choise for ssto correction
      integer :: nchflo    = 3 ! choise for hflo correction
!
      real :: cfacsst   = 1. ! coupling factor sst
      real :: cfactau   = 1. ! coupling factor stress
      real :: cfacice   = 1. ! coupling factor ice
      real :: cfacfresh = 1. ! coupling factor fresh water
!
      character(len=80) :: flukofile='fluko_data.txt'! file for flux corrction
      character(len=80) :: runoffmap='runoffmap.txt' ! file for runoff map
!
      end module cplmod
!
!     ==================
!     subroutine CLSGINI
!     ==================
!
      subroutine clsgini(kdatim,ktspd,kaomod,kcpl,kgui,pslm,kxa,kya)
      use cplmod
!
      real :: pslm(nxa,nya)
      integer :: kdatim(7)
      integer :: ktspd
      integer :: kaomod
      integer :: kcpl
      integer :: kxa,kya
      logical :: lopen
      character(len=10) :: chform
!
!     dbug
!
      real :: zduma(nxa,nya) = 1
      real :: zdumo(nxo,nyo) = 1
      real :: zslmaos(nxo,nyo) 
      real :: zslmaov(nxo,nyo)
      real :: zslmoa(nxa,nya)
!
      character(len=8) :: cform
!
      namelist/cplpar/ndbox,ndboy,nprint,nissta,niicea,nipmea,nitaua    &
     &               ,nihfla,nissto,nihflo,nflsst,nfltau,nflice,nflfresh&
     &               ,nfloheat,ncssta,ncicea,ncpmea,nctaua,nchflo,ncssto&
     &               ,cfacsst,cfacice,cfactau,cfacfresh,nroffc          &
     &               ,flukofile,runoffmap
!
!     check atm dimensions
!
      if(kxa /= nxa .or. kya /= nya) then
       write(nud,*)'!ERROR! wrong Plasim dimensions in cplmod'
       stop
      endif
!
!     read namelist
!
      do jt=30,99
       ntape=jt
       inquire(unit=ntape,opened=lopen)
       if(.not. lopen) exit
      enddo
      open(ntape,file='cpl_namelist',form='formatted')
       read(ntape,cplpar)
       write(nud,'(/," ****************************************")')
       write(nud,'(" * CPLMOD ",a30," *")') trim(version)
       write(nud,'(" ******************************************")')
       write(nud,'(" * Namelist CPLPAR from <cpl_namelist>    *")')
       write(nud,'(" ******************************************")')
       write(nud,cplpar)
      close(ntape)
!
!     initialize coupling method from kcpl (= nlsg in oceanmod)
!     ncoupling=1 (atmospheric sst and oceanic hfl)
!     ncoupling=2 (atmospheric hfl and oceanic sst)
!     
      ncoupling=kcpl
      ngui=kgui
!
!     initialize flux correction
!
      nfluko=nflsst+nfltau+nflice+nflfresh+nfloheat
      if(nfluko > 0) then
       call cflukoini
      else
       write(nud,*) ' '
       write(nud,*) 'No flux correction chosen in CLSGINI'
       write(nud,*) ' '
      endif
!
!     initialize runoff correction
!
      if(nroffc > 1) then
       write(chform,'(A1,I6.6,A3)') '(',nxa,'I1)'
       open(ntape,file=trim(runoffmap),form='formatted')
       read(ntape,chform,end=9000) iroffi(:,:)
       if(nroffc == 3) then
        read(ntape,*,end=9000)
        read(ntape,chform,end=9000) iroffo(:,:)
       endif
       close(ntape)
      endif
!
!     make plasim and lsg grids
!
      call cinigrids
!
!     put atm land sea mask to common
!
      call cputlsma(pslm,nxa,nya)
!
!     put atm time to common
!
      call cputtime(kdatim)
!
!     initialize LSG
!
      naomod=kaomod
      ndcoup=kaomod/ktspd
      call lsgini(ndcoup,ncoupling,ndbox,ndboy)
      if(ngui > 0) call do_lsg_visual
!
!     dbug print out
!
      if(nprint > 0) then
       zgr=180./(2.*asin(1.))
       write(nud,*) ' '
       write(nud,*) 'End of CLSGINI: '
       write(cform,'(A1,I2.2,A3)') '(',nxa,'I1)'
       write(nud,*) 'Atmosphere Mask'
       write(nud,cform) nint(aslm)
       write(cform,'(A1,I2.2,A3)') '(',nxo,'I1)'
       write(nud,*) 'Ocean Scalar Mask'
       write(nud,cform) nint(oslms)
       write(nud,*) 'Ocean Vector MasK'
       write(nud,cform) nint(oslmv)
       write(nud,*) ' '
       call cglm(aslm,zduma,zgs,zgm,zgw,agw,nxa,nya)
       write(nud,*) 'Atm. ocean area (sum/mean) + global area : '         &
     &            ,zgs*radea,zgm,zgw*radea
       call cglm(oslms,zdumo,zgs,zgm,zgw,ogw,nxo,nyo)
       write(nud,*) 'Oce. ocean area scalar (sum/mean) + global area : '  &
     &            ,zgs*radea,zgm,zgw*radea
       call cglm(oslmv,zdumo,zgs,zgm,zgw,ogw,nxo,nyo)
       write(nud,*) 'Oce. ocean area vector (sum/mean) + global area : '  &
     &            ,zgs*radea,zgm,zgw*radea
!
!      check mask overlaps:
!
       write(nud,*) 'Check mask overlaps for area interpolation'
       write(nud,*) ' '
       call cinter(aslm,zduma,zslmaos,zdumo,agw,ogw,xa1,xa2,xos1,xos2   &
     &            ,ya1,ya2,yo1,yo2,nxa,nxo,nya,nyo,0,nprint)
       kmiss=0
       do j=1,nyo
       do i=1,nxo
        zdumo(i,j)=0.
        if(zslmaos(i,j) < 1.E-9 .and. oslms(i,j) > 0.5) then
         kmiss=kmiss+1
         zdumo(i,j)=1.
        endif
       enddo
       enddo
       write(nud,*) kmiss                                                 &
     &           ,' Oce Scalar points must be extrapolated'  
       if(kmiss > 0) then
        write(nud,*) 'These points are:'
        do j=1,nyo
        do i=1,nxo
         if(zdumo(i,j) > 0.5) then
          write(nud,*) 'i,j= ',i,j,' Lon,Lat= ',xos(i,j)*zgr,yo(i,j)*zgr
         endif
        enddo
        enddo
        write(nud,*) 'Map of missed points'
        write(cform,'(A1,I2.2,A3)') '(',nxo,'I1)'
        write(nud,cform) nint(zdumo)
       endif
       zdumo=1.
       call cinter(aslm,zduma,zslmaov,zdumo,agw,ogw,xa1,xa2,xov1,xov2   &
     &            ,ya1,ya2,yo1,yo2,nxa,nxo,nya,nyo,0,nprint)
       kmiss=0
       do j=1,nyo
       do i=1,nxo
        zdumo(i,j)=0.
        if(zslmaov(i,j) < 1.E-9 .and. oslmv(i,j) > 0.5) then
         kmiss=kmiss+1
         zdumo(i,j)=1.
        endif
       enddo
       enddo
       write(nud,*) kmiss                                                 &
     &           ,' Oce Vector points must be extrapolated:'
       if(kmiss > 0) then
        write(nud,*) 'These points are:'
        do j=1,nyo
        do i=1,nxo
         if(zdumo(i,j) > 0.5) then
          write(nud,*) 'i,j= ',i,j,' Lon,Lat= ',xov(i,j)*zgr,yo(i,j)*zgr
         endif
        enddo
        enddo
        write(nud,*) 'Map of missed points'
        write(cform,'(A1,I2.2,A3)') '(',nxo,'I1)'
        write(nud,cform) nint(zdumo)
       endif
       zdumo=1.
       call cinter(oslms,zdumo,zslmoa,zduma,ogw,agw,xos1,xos2,xa1,xa2   &
     &            ,yo1,yo2,ya1,ya2,nxo,nxa,nyo,nya,0,nprint)
       kmiss=0
       do j=1,nya
       do i=1,nxa
        zduma(i,j)=0.
        if(zslmoa(i,j) < 1.E-9 .and. aslm(i,j) > 0.5) then
         kmiss=kmiss+1
         zduma(i,j)=1.
        endif
       enddo
       enddo
       write(nud,*) kmiss                                                 &
     &           ,' Atm points must be extrapolated:'
       if(kmiss > 0) then
        write(nud,*) 'These points are:'
        do j=1,nya
        do i=1,nxa
         if(zduma(i,j) > 0.5) then
          write(nud,*) 'i,j= ',i,j,' Lon,Lat= ',xa(i,j)*zgr,ya(i,j)*zgr
         endif
        enddo
        enddo
        write(nud,*) 'Map of missed points'
        write(cform,'(A1,I2.2,A3)') '(',nxa,'I1)'
        write(nud,cform) nint(zduma)
       endif
       write(nud,*)
!
!     check mask overlaps for bi-linear interpolation:
!
       write(nud,*) 'Check for bi-linear interpolation'
       write(nud,*)
       zduma=1.
       zdumo=1.
       call cinter2(aslm,zduma,zslmaos,zdumo,agw,ogw,xa,xos             &
     &             ,ya,yo,nxa,nxo,nya,nyo,0,nprint)
       kmiss=0
       do j=1,nyo
       do i=1,nxo
        zdumo(i,j)=0.
        if(zslmaos(i,j) < 1.E-9 .and. oslms(i,j) > 0.5) then
         kmiss=kmiss+1
         zdumo(i,j)=1.
        endif
       enddo
       enddo
       write(nud,*) kmiss                                                 &
     &           ,' Oce Scalar points must be extrapolated'  
       if(kmiss > 0) then
        write(nud,*) 'These points are:'
        do j=1,nyo
        do i=1,nxo
         if(zdumo(i,j) > 0.5) then
          write(nud,*) 'i,j= ',i,j,' Lon,Lat= ',xos(i,j)*zgr,yo(i,j)*zgr
         endif
        enddo
        enddo
        write(nud,*) 'Map of missed points'
        write(cform,'(A1,I2.2,A3)') '(',nxo,'I1)'
        write(nud,cform) nint(zdumo)
       endif
       zdumo=1.
       call cinter2(aslm,zduma,zslmaov,zdumo,agw,ogw,xa,xov             &
     &             ,ya,yo,nxa,nxo,nya,nyo,0,nprint)
       kmiss=0
       do j=1,nyo
       do i=1,nxo
        zdumo(i,j)=0.
        if(zslmaov(i,j) < 1.E-9 .and. oslmv(i,j) > 0.5) then
         kmiss=kmiss+1
         zdumo(i,j)=1.
        endif
       enddo
       enddo
       write(nud,*) kmiss                                                 &
     &           ,' Oce Vector points must be extrapolated:'
       if(kmiss > 0) then
        write(nud,*) 'These points are:'
        do j=1,nyo
        do i=1,nxo
         if(zdumo(i,j) > 0.5) then
          write(nud,*) 'i,j= ',i,j,' Lon,Lat= ',xov(i,j)*zgr,yo(i,j)*zgr
         endif
        enddo
        enddo
        write(nud,*) 'Map of missed points'
        write(cform,'(A1,I2.2,A3)') '(',nxo,'I1)'
        write(nud,cform) nint(zdumo)
       endif
       zdumo=1.
       call cinter2(oslms,zdumo,zslmoa,zduma,ogw,agw,xos,xa             &
     &             ,yo,ya,nxo,nxa,nyo,nya,0,nprint)
       kmiss=0
       do j=1,nya
       do i=1,nxa
        zduma(i,j)=0.
        if(zslmoa(i,j) < 1.E-9 .and. aslm(i,j) > 0.5) then
         kmiss=kmiss+1
         zduma(i,j)=1.
        endif
       enddo
       enddo
       write(nud,*) kmiss                                                 &
     &           ,' Atm points must be extrapolated:'
       if(kmiss > 0) then
        write(nud,*) 'These points are:'
        do j=1,nya
        do i=1,nxa
         if(zduma(i,j) > 0.5) then
          write(nud,*) 'i,j= ',i,j,' Lon,Lat= ',xa(i,j)*zgr,ya(i,j)*zgr
         endif
        enddo
        enddo
        write(nud,*) 'Map of missed points'
        write(cform,'(A1,I2.2,A3)') '(',nxa,'I1)'
        write(nud,cform) nint(zduma)
       endif
       write(nud,*)
      endif
!
      return
!
 9000 continue
      write(nud,*) 'ERROR: nroffc set to ',nroffc,' but file '            &
     &           ,trim(runoffmap)                                       &
     &           ,'not correctly given!'
      stop
!
      end subroutine clsgini
!
!     ===================
!     subroutine CLSGSTEP
!     ===================
!
      subroutine clsgstep(kdatim,kstep,psst,ptaux,ptauy,ppme,proff,pice &
     &                   ,pheat,pfldo)
      use cplmod
!
      real :: psst(nxa,nya)      ! atm sst (input)
      real :: ptaux(nxa,nya)     ! atm u-stress (input)
      real :: ptauy(nxa,nya)     ! atm v-stress (input)
      real :: ppme(nxa,nya)      ! atm p-e (input)
      real :: proff(nxa,nya)     ! atm runoff (input)
      real :: pice(nxa,nya)      ! atm ice thickness (incl. snow) (input)
      real :: pheat(nxa,nya)     ! atm heat flux (input; not used yet)
      real :: pfldo(nxa,nya)     ! atm deep ocan heat flux (output)
!
      real :: zfresh(nxa,nya)    ! atm fresh water flux (p-e+runoff)
!
      integer :: kdatim(7)       ! date and time
      integer :: kstep           ! current atm time step
!
!     put all atm input to common (+ accumulate)
!
      if(nroffc > 0) call croffc(proff)
!
      zfresh(:,:)=ppme(:,:)+proff(:,:)
      call cputtime(kdatim)
      call cputtaua(ptaux,ptauy)
      call cputpmea(zfresh)
      if(ncoupling == 2) call cputhfla(pheat)
!
!     run lsg if its time
!
      if (mod(kstep,naomod) == naomod-1) then
!
       if(nprint > 0) then
        write(nud,*) 'in cpl: make lsg step'
       endif
!
       if(ncoupling == 1) then
        call cputssta(psst)
       endif
       call cputicea(pice)
       call lsgstep
       if(ngui > 0) call do_lsg_visual
!
!      get the new deep ocean heat flux
!
       if(ncoupling == 1) then
        call cgethfla(pfldo)
       elseif(ncoupling == 2) then
        call cgetssta(psst)
       endif
      endif
!
      return
      end subroutine clsgstep
!
!     ===================
!     subroutine CLSGSTOP
!     ===================
!
      subroutine clsgstop
      use cplmod
!
!     finalize lsg
!
      call lsgstop
!
      return
      end subroutine clsgstop
!
!     ====================
!     subroutine CINIGRIDS
!     ====================
!
      subroutine cinigrids
      use cplmod
!
      real (kind=8) :: zsi(nya),zgw(nya)
!
      zpi=2.*asin(1.)
!
!     atmospheric grid
!
      call inigau(nya,zsi,zgw)
      zdxa=2.*zpi/real(nxa)
      do i=1,nxa
      do j=1,nya
       xa(i,j)=zdxa*(i-1)
       ya(i,j)=asin(zsi(j))
      enddo
      enddo
      do i=2,nxa
       xa1(i,:)=0.5*(xa(i,:)+xa(i-1,:))
      enddo
      xa1(1,:)=xa1(2,:)-zdxa
      do i=1,nxa-1
       xa2(i,:)=0.5*(xa(i,:)+xa(i+1,:))
      enddo
      xa2(nxa,:)=xa2(nxa-1,:)+zdxa
      zgw(:)=2.*zgw(:)/sum(zgw(:))
      zsinm=sin(zpi*0.5)
      ya1(:,1)=asin(zsinm)
      do j=2,nya
       zsin=zsinm-zgw(j-1)
       ya1(:,j)=asin(zsin)
       ya2(:,j-1)=asin(zsin)
       zsinm=zsin
      enddo
      ya2(:,nya)=-zpi*0.5
!
!     oceanic grid
!
      call cmklonlate(xos,yo,nxo,nyo)
      xov(:,:)=xos(:,:)+zpi/real(nxo)
      where(xos(:,:) < 0.) xos(:,:)=xos(:,:)+2.*zpi
      where(xov(:,:) < 0.) xov(:,:)=xov(:,:)+2.*zpi
      zdxo=2.*zpi/real(nxo)
      xos1(:,:)=xos(:,:)-0.5*zdxo
      xos2(:,:)=xos(:,:)+0.5*zdxo
      xov1(:,:)=xov(:,:)-0.5*zdxo
      xov2(:,:)=xov(:,:)+0.5*zdxo
      zdyo=zpi/real(nyo)
      yo1(:,:)=yo(:,:)-zdyo*0.5
      yo2(:,:)=yo(:,:)+zdyo*0.5
!
!     weights
!
      do j=1,nya
      do i=1,nxa
       agw(i,j)=zdxa*abs(sin(ya1(i,j))-sin(ya2(i,j)))
      enddo
      enddo
      do j=1,nxo
      do i=1,nyo
       ogw(i,j)=zdxo*abs(sin(yo1(i,j))-sin(yo2(i,j)))
      enddo
      enddo
!
      return
      end
!
!     =====================
!     subroutine CMKLONLATE
!     =====================
!
      subroutine cmklonlate(x,y,nlon,nlat)
      parameter(radea=6.371E6)
      real :: x(nlon,nlat)
      real :: y(nlon,nlat)
      real :: rlat(nlat)
!
      delta=360./real(nlon)
!
      pi=2.*asin(1.)
!
!     lon und lat berechnen(umrechnung vom E Gitter in geog. koord.) 
!
      do j=1,nlat
       y(:,j)=90./real(nlat)-180.*real(j-nlat/2)/real(nlat)
      enddo	 
!
      rlonlimit=180.-delta/100. ! lower the limit for lon with respect to delta
      do j=1,nlat
       do i=1,nlon
        x(i,j)=(((i-(j/2.))/(nlon/2.))*180.)-delta 
!
!     lonlimit zur Vermeidung von Fehlern durch Ungenauigkeiten
!
        if(x(i,j).gt.rlonlimit) x(i,j)=-(360.-x(i,j))
        if(x(i,j).lt.-rlonlimit) x(i,j)=-x(i,j)
       enddo
      enddo	 

      x(:,:)=x(:,:)*pi/180.
      y(:,:)=y(:,:)*pi/180.
!
      return
      end
!
!     ===================
!     subroutine CPUTTIME
!     ===================
!
      subroutine cputtime(kdatim)
      use cplmod
!
!     put time to common
!
      integer :: kdatim(7)
!
      ndatim(:)=kdatim(:)
      nyear=ndatim(1)
      nmonth=ndatim(2)
      nday=ndatim(3)
      nhour=ndatim(4)
      nmin=ndatim(5)
      nwday=ndatim(6)
!      
      return
      end
!
!     ===================
!     subroutine CGETTIME
!     ===================
!
      subroutine cgettime(kdatim)
      use cplmod
!
!     get time from common
!
      integer :: kdatim(7)
!
      kdatim(:)=ndatim(:)
!
      if(nprint > 0) then
       write(nud,*) ' '
       write(nud,*) 'In CGETTIME, datim= ',ndatim(:)
       write(nud,*) ' '
      endif
!       
      return
      end
!
!     ===================
!     subroutine CPUTLSMA
!     ===================
!
      subroutine cputlsma(plsm,kx,ky)
      use cplmod
!
!     put atm land see mask to common
!
      real :: plsm(nxa,nya)      ! land see mask (input)
!
!     check dimensions
!
      if(nxa /= kx .or. nya /= ky) then
       write(nud,*) '!ERROR! inconsistent atm. dimensions in cpl'
       stop
      endif
!
!     copy land see mask and convert 0/1 to 1/0
!
      aslm(:,:)=real(nint(1.-plsm(:,:)))
!
      return
      end
!
!     ===================
!     subroutine CPUTLSMO
!     ===================
!
      subroutine cputlsmo(plsms,plsmv,kx,ky)
      use cplmod
!
!     put oce land see mask to common
!
      real :: plsms(nxo,nyot)   ! land see mask scalar (input)
      real :: plsmv(nxo,nyot)   ! land see mask vector (input)
!
!     check dimensions
!
      if(nxo /= kx .or. nyot /= ky) then
       write(nud,*) '!ERROR! inconsistent oce. dimensions in cpl'
       stop
      endif
!
!     copy land see mask and convert 0/1 to 1/0
!
      oslms(1:nxo,1:nyo)=real(nint(plsms(1:nxo,3:nyot-2)))
      oslmv(1:nxo,1:nyo)=real(nint(plsmv(1:nxo,3:nyot-2)))
!
      return
      end
!
!     ===================
!     subroutine CPUTSSTA
!     ===================
!
      subroutine cputssta(psst)
      use cplmod
!
!     put atm sst to common (accumulated)
!
      real :: psst(nxa,nya) ! sst (input)
!
!     accumulate
!
      assta(:,:)=assta(:,:)+psst(:,:)
      nssta=nssta+1
!
      return
      end
!
!     ===================
!     subroutine CPUTICEA
!     ===================
!
      subroutine cputicea(piced)
      use cplmod
!
!     put atm ice to common (accumulated)
!
      real :: piced(nxa,nya) ! ice thickness (input)
!
!     accumulate
!
      aiceda(:,:)=aiceda(:,:)+piced(:,:)
      nicea=nicea+1
!
      return
      end
!
!     ===================
!     subroutine CPUTPMEA
!     ===================
!
      subroutine cputpmea(ppme)
      use cplmod
!
!     put atm fresh water flux to common (accumulated)
!
      real :: ppme(nxa,nya) ! fresh water flux (input)
!
!     accumulate
!
      apmea(:,:)=apmea(:,:)+ppme(:,:)
      npmea=npmea+1
!
      return
      end
!
!     ===================
!     subroutine CPUTTAUA
!     ===================
!
      subroutine cputtaua(ptaux,ptauy)
      use cplmod
!
!     put atm wind stress to common (accumulated)
!
      real :: ptaux(nxa,nya) ! u-stress (input)
      real :: ptauy(nxa,nya) ! v-stress (input)
!
!     accumulate
!
      atauxa(:,:)=atauxa(:,:)+ptaux(:,:)
      atauya(:,:)=atauya(:,:)+ptauy(:,:)
      ntaua=ntaua+1
!
      return
      end
!
!     ===================
!     subroutine CPUTHFLA
!     ===================
!
      subroutine cputhfla(phfl)
      use cplmod
!
!     put atm heat flux to common (accumulated)
!
      real :: phfl(nxa,nya) ! heat flux (input)
!
!     accumulate
!
      ahfla(:,:)=ahfla(:,:)+phfl(:,:)
      nhfla=nhfla+1
!
      return
      end
!
!     ===================
!     subroutine CPUTHFLO
!     ===================
!
      subroutine cputhflo(phfl)
      use cplmod
!
!     put oce heat flux to common (accumulated)
!
      real :: phfl(nxo,nyot) ! heat flux (input)
!
!     accumulate
!
      ohfla(:,:)=ohfla(:,:)+phfl(1:nxo,3:nyot-2)
      nhflo=nhflo+1
!
      return
      end
!
!     ===================
!     subroutine CPUTSSTO
!     ===================
!
      subroutine cputssto(psst)
      use cplmod
!
!     put oce sst to common (accumulated)
!
      real :: psst(nxo,nyot) ! heat flux (input)
!
!     accumulate
!
      ossta(:,:)=ossta(:,:)+psst(1:nxo,3:nyot-2)
      nssto=nssto+1
!
      return
      end
!
!     ===================
!     subroutine CGETSSTO
!     ===================
!
      subroutine cgetssto(psst)
      use cplmod
      parameter(rkelvin=273.16)
!
!     get interpolated sst from common 
!
      real :: psst(nxo,nyot) ! sst (output)
!
!     get accumulated atm sst
!
      if(nssta == 0) then
       write(nud,*) '!ERROR! no accumulated atm sst found'
       stop
      endif
      asst(:,:)=assta(:,:)/real(nssta)
!
!     interpolate
!
      if(nissta== 1) then
       call cinter(asst,aslm,osst,oslms,agw,ogw,xa1,xa2,xos1,xos2       &
     &            ,ya1,ya2,yo1,yo2,nxa,nxo,nya,nyo,ncssta,nprint)
      else
       call cinter2(asst,aslm,osst,oslms,agw,ogw,xa,xos                 &
     &             ,ya,yo,nxa,nxo,nya,nyo,ncssta,nprint)
      endif
!
!     copy to output (convert into C)
!
      psst(1:nxo,3:nyot-2)=osst(:,:)-rkelvin
!
!     dbug printout
!
      if(nprint > 0) then
       write(nud,*) ' '
       write(nud,*) 'In CGETSSTO: '
       write(nud,*) 'sst from ',nssta,' accumulated values'
       call cglm(asst,aslm,zgs,zgm,zgw,agw,nxa,nya)
       write(nud,*) 'Atm (sum,mean,max,min): '                            &
     &            ,zgs,zgm,MAXVAL(asst,MASK=(aslm > 0.5))               &
     &            ,MINVAL(asst,MASK=(aslm > 0.5))
       call cglm(osst,oslms,zgs,zgm,zgw,ogw,nxo,nyo)
       write(nud,*) 'Oce (sum,mean,max,min): '                            &
     &            ,zgs,zgm,MAXVAL(osst,MASK=(oslms > 0.5))              &
     &            ,MINVAL(osst,MASK=(oslms > 0.5))
       write(nud,*) ' '
      endif
!
!     reset accumulated field
!
      assta(:,:)=0.
      nssta=0
!
      return
      end
!
!     ===================
!     subroutine CGETICEO
!     ===================
!
      subroutine cgeticeo(piced)
      use cplmod
      parameter(rhoi=1000.)
      parameter(rho0=1030.)
      real :: zpiced(nxo,nyot-4)
!
!     get interpolated ice from common 
!
      real :: piced(nxo,nyot) ! ice thickness (output)
!
!     get accumulated atm ice
!
      if(nicea == 0) then
       write(nud,*) '!ERROR! no accumulated atm ice found'
       stop
      endif
      aiced(:,:)=aiceda(:,:)/real(nicea)
!
!     interpolate
!
      if(niicea == 1) then
       call cinter(aiced,aslm,oiced,oslms,agw,ogw,xa1,xa2,xos1,xos2     &
     &            ,ya1,ya2,yo1,yo2,nxa,nxo,nya,nyo,ncicea,nprint)
      else
       call cinter2(aiced,aslm,oiced,oslms,agw,ogw,xa,xos               &
     &             ,ya,yo,nxa,nxo,nya,nyo,ncicea,nprint)
      endif
!
!     copy to output
!
      piced(1:nxo,3:nyot-2)=oiced(:,:)*rhoi/rho0
!
!     dbug printout
!
      if(nprint > 0) then
       write(nud,*) ' '
       write(nud,*) 'In CGETICEO: '
       write(nud,*) 'ice from ',nicea,' accumulated values'
       call cglm(aiced,aslm,zgs,zgm,zgw,agw,nxa,nya)
       write(nud,*) 'Atm (sum,mean,max,min): '                            &
     &            ,zgs,zgm,MAXVAL(aiced,MASK=(aslm > 0.5))              &
     &            ,MINVAL(aiced,MASK=(aslm > 0.5))
       call cglm(oiced,oslms,zgs,zgm,zgw,ogw,nxo,nyo)
       write(nud,*) 'Oce uncorrected (sum,mean,max,min): '                &
     &            ,zgs,zgm,MAXVAL(oiced,MASK=(oslms > 0.5))             &
     &            ,MINVAL(oiced,MASK=(oslms > 0.5))
       zpiced(:,:) = piced(:,3:nyot-2)
       call cglm(zpiced,oslms,zgs,zgm,zgw,ogw,nxo,nyo)  
       write(nud,*) 'Oce corrected (sum,mean,max,min): '                  &
     &        ,zgs,zgm,MAXVAL(zpiced(:,:),MASK=(oslms > 0.5))           &
     &        ,MINVAL(zpiced(:,:),MASK=(oslms > 0.5))
       write(nud,*) ' '
      endif
!
!     reset accumulated field
!
      aiceda(:,:)=0.
      nicea=0
!
      return
      end
!
!     ===================
!     subroutine CGETPMEO
!     ===================
!
      subroutine cgetpmeo(ppme)
      use cplmod
      parameter(rhop=1000.)
      parameter(rho0=1030.)
!
!     get interpolated pme from common 
!
      real :: ppme(nxo,nyot) ! pme (output)
!
!     get accumulated atm pme
!
      if(npmea == 0) then
       write(nud,*) '!ERROR! no accumulated atm pme found'
       stop
      endif
      apme(:,:)=apmea(:,:)/real(npmea)
!
!     interpolate
!
      if(nipmea == 1) then
       call cinter(apme,aslm,opme,oslms,agw,ogw,xa1,xa2,xos1,xos2       &
     &            ,ya1,ya2,yo1,yo2,nxa,nxo,nya,nyo,ncpmea,nprint)
      else
       call cinter2(apme,aslm,opme,oslms,agw,ogw,xa,xos                 &
     &             ,ya,yo,nxa,nxo,nya,nyo,ncpmea,nprint)
      endif
!
!     copy to output and correct for density differences
!
      ppme(1:nxo,3:nyot-2)=opme(:,:)*rhop/rho0
!
!     dbug printout
!
      if(nprint > 0) then
       write(nud,*) ' '
       write(nud,*) 'In CGETPMEO: '
       write(nud,*) 'pme from ',npmea,' accumulated values'
       call cglm(apme,aslm,zgs,zgm,zgw,agw,nxa,nya)
       write(nud,*) 'Atm (sum,mean,max,min): '                            &
     &            ,zgs,zgm,MAXVAL(apme,MASK=(aslm > 0.5))               &
     &            ,MINVAL(apme,MASK=(aslm > 0.5))
       call cglm(opme,oslms,zgs,zgm,zgw,ogw,nxo,nyo)
       write(nud,*) 'Oce (sum,mean,max,min): '                            &
     &            ,zgs,zgm,MAXVAL(opme,MASK=(oslms > 0.5))              &
     &            ,MINVAL(opme,MASK=(oslms > 0.5))
       write(nud,*) ' '
      endif
!
!     reset accumulated field
!
      apmea(:,:)=0.
      npmea=0
!
      return
      end
!
!     ===================
!     subroutine CGETTAUO
!     ===================
!
      subroutine cgettauo(ptaux,ptauy)
      use cplmod
!
!     get interpolated stress from common 
!
      real :: ptaux(nxo,nyot) ! u-stress (output)
      real :: ptauy(nxo,nyot) ! v-stress (output)
!
!     get accumulated atm stress
!
      if(ntaua == 0) then
       write(nud,*) '!ERROR! no accumulated atm stress found'
       stop
      endif
      ataux(:,:)=atauxa(:,:)/real(ntaua)
      atauy(:,:)=atauya(:,:)/real(ntaua)
!
!     interpolate
!
      if(nitaua == 1) then
       call cinter(ataux,aslm,otaux,oslmv,agw,ogw,xa1,xa2,xov1,xov2     &
     &            ,ya1,ya2,yo1,yo2,nxa,nxo,nya,nyo,nctaua,nprint)
       call cinter(atauy,aslm,otauy,oslmv,agw,ogw,xa1,xa2,xov1,xov2     &
     &            ,ya1,ya2,yo1,yo2,nxa,nxo,nya,nyo,nctaua,nprint)
      else
       call cinter2(ataux,aslm,otaux,oslmv,agw,ogw,xa,xov               &
     &             ,ya,yo,nxa,nxo,nya,nyo,nctaua,nprint)
       call cinter2(atauy,aslm,otauy,oslmv,agw,ogw,xa,xov               &
     &             ,ya,yo,nxa,nxo,nya,nyo,nctaua,nprint)
      endif
!
!     copy to output
!
      ptaux(1:nxo,3:nyot-2)=otaux(:,:)
      ptauy(1:nxo,3:nyot-2)=otauy(:,:)
!
!     dbug printout
!
      if(nprint > 0) then
       write(nud,*) ' '
       write(nud,*) 'In CGETTAUO: '
       write(nud,*) 'stress from ',ntaua,' accumulated values'
       call cglm(ataux,aslm,zgs,zgm,zgw,agw,nxa,nya)
       write(nud,*) 'Atm u (sum,mean,max,min): '                          &
     &            ,zgs,zgm,MAXVAL(ataux,MASK=(aslm > 0.5))              &
     &            ,MINVAL(ataux,MASK=(aslm > 0.5))
       call cglm(atauy,aslm,zgs,zgm,zgw,agw,nxa,nya)
       write(nud,*) 'Atm v (sum,mean,max,min): '                          &
     &            ,zgs,zgm,MAXVAL(atauy,MASK=(aslm > 0.5))              &
     &            ,MINVAL(atauy,MASK=(aslm > 0.5))
       call cglm(otaux,oslmv,zgs,zgm,zgw,ogw,nxo,nyo)
       write(nud,*) 'Oce u (sum,mean,max,min): '                          &
     &            ,zgs,zgm,MAXVAL(otaux,MASK=(oslmv > 0.5))             &
     &            ,MINVAL(otaux,MASK=(oslmv > 0.5))
       call cglm(otauy,oslmv,zgs,zgm,zgw,ogw,nxo,nyo)
       write(nud,*) 'Oce v (sum,mean,max,min): '                          &
     &            ,zgs,zgm,MAXVAL(otauy,MASK=(oslmv > 0.5))             &
     &            ,MINVAL(otauy,MASK=(oslmv > 0.5))
       write(nud,*) ' '
      endif
!
!     reset accumulated field
!
      atauxa(:,:)=0.
      atauya(:,:)=0.
      ntaua=0
!
      return
      end
!
!     ===================
!     subroutine CGETHFLO
!     ===================
!
      subroutine cgethflo(phfl)
      use cplmod
!
!     get interpolated atm hfl from common 
!
      real :: phfl(nxo,nyot) ! hfl (output)
!
!     get accumulated atm hfl
!
      if(nhfla == 0) then
       write(nud,*) '!ERROR! no accumulated atm hfl found'
       stop
      endif
      ahfl(:,:)=ahfla(:,:)/real(nhfla)
!
!     interpolate
!
      if(nihfla == 1) then
       call cinter(ahfl,aslm,ohfl,oslms,agw,ogw,xa1,xa2,xos1,xos2       &
     &            ,ya1,ya2,yo1,yo2,nxa,nxo,nya,nyo,nchfla,nprint)
      else
       call cinter2(ahfl,aslm,ohfl,oslms,agw,ogw,xa,xos                 &
     &             ,ya,yo,nxa,nxo,nya,nyo,nchfla,nprint)
      endif
!
!     copy to output
!
      phfl(1:nxo,3:nyot-2)=ohfl(:,:)
!
!     dbug printout
!
      if(nprint > 0) then
       write(nud,*) ' '
       write(nud,*) 'In CGETHFLO: '
       write(nud,*) 'hfl from ',nhfla,' accumulated values'
       call cglm(ahfl,aslm,zgs,zgm,zgw,agw,nxa,nya)
       write(nud,*) 'Atm (sum,mean,max,min): '                            &
     &            ,zgs,zgm,MAXVAL(ahfl,MASK=(aslm > 0.5))               &
     &            ,MINVAL(ahfl,MASK=(aslm > 0.5))
       call cglm(ohfl,oslms,zgs,zgm,zgw,ogw,nxo,nyo)
       write(nud,*) 'Oce (sum,mean,max,min): '                            &
     &            ,zgs,zgm,MAXVAL(ohfl,MASK=(oslms > 0.5))              &
     &            ,MINVAL(ohfl,MASK=(oslms > 0.5))
       write(nud,*) ' '
      endif
!
!     reset accumulated field
!
      ahfla(:,:)=0.
      nhfla=0
!
      return
      end
!
!     ===================
!     subroutine CGETHFLA
!     ===================
!
      subroutine cgethfla(phfl)
      use cplmod
!
!     get interpolated oce hfl from common 
!
      real :: phfl(nxa,nya) ! hfl (output)
!
!     get accumulated oce hfl
!
      if(nhflo == 0) then
       write(nud,*) '!ERROR! no accumulated oce hfl found'
       stop
      endif
      ohfl(:,:)=ohfla(:,:)/real(nhflo)
!
!     interpolate
!
      if(nihflo == 1) then
       call cinter(ohfl,oslms,ahfl,aslm,ogw,agw,xos1,xos2,xa1,xa2       &
     &            ,yo1,yo2,ya1,ya2,nxo,nxa,nyo,nya,nchflo,nprint)
      else
       call cinter2(ohfl,oslms,ahfl,aslm,ogw,agw,xos,xa                 &
     &             ,yo,ya,nxo,nxa,nyo,nya,nchflo,nprint)
      endif
!
!     copy to output
!
      phfl(:,:)=ahfl(:,:)
!
!     dbug printout
!
      if(nprint > 0) then
       write(nud,*) ' '
       write(nud,*) 'In CGETHFLA: '
       write(nud,*) 'hfl from ',nhflo,' accumulated values'
       call cglm(ohfl,oslms,zgs,zgm,zgw,ogw,nxo,nyo)
       write(nud,*) 'Oce (sum,mean,max,min): '                            &
     &            ,zgs,zgm,MAXVAL(ohfl,MASK=(oslms > 0.5))              &
     &            ,MINVAL(ohfl,MASK=(oslms > 0.5))
       call cglm(ahfl,aslm,zgs,zgm,zgw,agw,nxa,nya)
       write(nud,*) 'Atm (sum,mean,max,min): '                            &
     &            ,zgs,zgm,MAXVAL(ahfl,MASK=(aslm > 0.5))               &
     &            ,MINVAL(ahfl,MASK=(aslm > 0.5))
       write(nud,*) ' '
      endif
!
!     reset accumulated field
!
      ohfla(:,:)=0.
      nhflo=0
!
      return
      end
!
!     ===================
!     subroutine CGETSSTA
!     ===================
!
      subroutine cgetssta(psst)
      use cplmod
      parameter(tfreeze=271.25)
      parameter(rkelvin=273.16)
!
!     get interpolated oce sst from common 
!
      real :: psst(nxa,nya) ! sst (output)
!
!     get accumulated oce sst
!
      if(nssto == 0) then
       write(nud,*) '!ERROR! no accumulated oce sst found'
       stop
      endif
      osst(:,:)=ossta(:,:)/real(nssto)+rkelvin
!
!     interpolate
!
      if (nissto == 1) then
       call cinter(osst,oslms,asst,aslm,ogw,agw,xos1,xos2,xa1,xa2       &
     &            ,yo1,yo2,ya1,ya2,nxo,nxa,nyo,nya,ncssto,nprint)
      else
       call cinter2(osst,oslms,asst,aslm,ogw,agw,xos,xa                 &
     &             ,yo,ya,nxo,nxa,nyo,nya,ncssto,nprint)
      endif
!
!     copy to output
!
      psst(:,:)=asst(:,:)
!
!     dbug printout
!
      if(nprint > 0) then
       write(nud,*) ' '
       write(nud,*) 'In CGETSSTA: '
       write(nud,*) 'sst from ',nssto,' accumulated values'
       call cglm(osst,oslms,zgs,zgm,zgw,ogw,nxo,nyo)
       write(nud,*) 'Oce (sum,mean,max,min): '                            &
     &            ,zgs,zgm,MAXVAL(osst,MASK=(oslms > 0.5))              &
     &            ,MINVAL(osst,MASK=(oslms > 0.5))
       call cglm(asst,aslm,zgs,zgm,zgw,agw,nxa,nya)
       write(nud,*) 'Atm (sum,mean,max,min): '                            &
     &            ,zgs,zgm,MAXVAL(asst,MASK=(aslm > 0.5))               &
     &            ,MINVAL(asst,MASK=(aslm > 0.5))
       write(nud,*) ' '
      endif
!
!     reset accumulated field
!
      ossta(:,:)=0.
      nssto=0
!
      return
      end
!
!     =================
!     subroutine CINTER
!     =================
!
      subroutine cinter(pfi,pslmi,pfo,pslmo,pgwi,pgwo,pxi1,pxi2,pxo1    &
     &                 ,pxo2,pyi1,pyi2,pyo1,pyo2,nxi,nxo,nyi,nyo,nc,npr)
!
      parameter(zerr=-9.E09)
!
      real :: fi(nxi,nyi)      ! input field
      real :: slmi(nxi,nyi)    ! input mask (0/1)
      real :: fo(nxo,nyo)      ! output field
      real :: slmo(nxo,nyo)    ! output mask (0/1)
      real :: xi1(nxi,nyi)     ! left x's for input
      real :: xi2(nxi,nyi)     ! right x's for input
      real :: xo1(nxo,nyo)     ! left x's for output
      real :: xo2(nxo,nyo)     ! right x's for output
      real :: yi1(nxi,nyi)     ! southern y's for input
      real :: yi2(nxi,nyi)     ! northern y's for input
      real :: yo1(nxo,nyo)     ! southern y's for output
      real :: yo2(nxo,nyo)     ! northern y's for output
!
      real :: pfi(nxi,nyi)
      real :: pslmi(nxi,nyi)
      real :: pfo(nxo,nyo)
      real :: pslmo(nxo,nyo)
      real :: pxi1(nxi,nyi)
      real :: pxi2(nxi,nyi)
      real :: pxo1(nxo,nyo)
      real :: pxo2(nxo,nyo)
      real :: pyi1(nxi,nyi)
      real :: pyi2(nxi,nyi)
      real :: pyo1(nxo,nyo)
      real :: pyo2(nxo,nyo)
      real :: pgwi(nxi,nyi)
      real :: pgwo(nxo,nyo)
!
      pi=2.*asin(1.)
!
      fi(:,:)=pfi(:,:)
      slmi(:,:)=pslmi(:,:)
      if(pyi1(1,1) < pyi2(1,1)) then
       yi1(:,:)=pyi1(:,:)
       yi2(:,:)=pyi2(:,:)
      else
       yi1(:,:)=pyi2(:,:)
       yi2(:,:)=pyi1(:,:)
      endif
      xi1(:,:)=pxi1(:,:)
      xi2(:,:)=pxi2(:,:)
      if (pyo1(1,1) < pyo2(1,1)) then
       yo1(:,:)=pyo1(:,:)
       yo2(:,:)=pyo2(:,:)
      else
       yo1(:,:)=pyo2(:,:)
       yo2(:,:)=pyo1(:,:)
      endif
      xo1(:,:)=pxo1(:,:)
      xo2(:,:)=pxo2(:,:)
      slmo(:,:)=pslmo(:,:)

      kmiss=0
      do jo=1,nyo
      do io=1,nxo
       zgw=0.
       zf=0.
       fo(io,jo)=zerr
       do ji=1,nyi
       do ii=1,nxi
        xol=xo1(io,jo)
        xor=xo2(io,jo)
        xil=xi1(ii,ji)
        xir=xi2(ii,ji)
        if(ABS(xol-xil) > pi) then
         if(xol > xil) then
          xor=xor-2.*pi
          xol=xol-2.*pi
         else
          xir=xir-2.*pi
          xil=xil-2.*pi
         endif
        endif
!
        if((xol <= xil .and. xor > xil)                   &
     &  .and. (yo1(io,jo) <= yi1(ii,ji) .and. yo2(io,jo) > yi1(ii,ji))) &
     &  then
         zweight=(MIN(xor,xir)-xil)                &
     &          *ABS(sin(MIN(yo2(io,jo),yi2(ii,ji)))-sin(yi1(ii,ji)))
         zf=zf+zweight*fi(ii,ji)*slmi(ii,ji) 
         zgw=zgw+zweight*slmi(ii,ji)
        endif
!
        if((xol > xil .and. xol < xir)                    &
     &  .and. (yo1(io,jo) <= yi1(ii,ji) .and. yo2(io,jo) > yi1(ii,ji))) &
     &  then
         zweight=(MIN(xor,xir)-xol)                &
     &          *ABS(sin(MIN(yo2(io,jo),yi2(ii,ji)))-sin(yi1(ii,ji)))
         zf=zf+zweight*fi(ii,ji)*slmi(ii,ji) 
         zgw=zgw+zweight*slmi(ii,ji)
        endif
!
        if((xol <= xil .and. xor > xil)                   &
     &  .and. (yo1(io,jo) > yi1(ii,ji) .and. yo1(io,jo) < yi2(ii,ji)))  &
     &  then
         zweight=(MIN(xor,xir)-xil)                &
     &          *ABS(sin(MIN(yo2(io,jo),yi2(ii,ji)))-sin(yo1(io,jo)))
         zf=zf+zweight*fi(ii,ji)*slmi(ii,ji) 
         zgw=zgw+zweight*slmi(ii,ji)
        endif
!
        if((xol > xil .and. xol < xir)                    &
     &  .and. (yo1(io,jo) > yi1(ii,ji) .and. yo1(io,jo) < yi2(ii,ji)))  &
     &  then
         zweight=(MIN(xor,xir)-xol)                &
     &          *ABS(sin(MIN(yo2(io,jo),yi2(ii,ji)))-sin(yo1(io,jo)))
         zf=zf+zweight*fi(ii,ji)*slmi(ii,ji) 
         zgw=zgw+zweight*slmi(ii,ji)
        endif
       enddo
       enddo
       if(zgw > 0.) fo(io,jo)=zf/zgw
       if(fo(io,jo) == zerr .and. slmo(io,jo) > 0.) then 
        kmiss=kmiss+1
       else
        fo(io,jo)=fo(io,jo)*slmo(io,jo)
       endif
      enddo
      enddo
!
      if(npr > 0) then
       write(nud,*) ' '
       write(nud,*) 'In CINTER: ',kmiss,' points to be extrapolated'
      endif
!
      if(kmiss > 0) call cfill(fo,slmo,nxo,nyo,zerr,npr)
!
      if(nc > 0) then
       call conserve(fi,fo,slmi,slmo,pgwi,pgwo,nxi,nxo,nyi,nyo,nc,npr)
      endif
!
      pfo(:,:)=fo(:,:)
!
      return
      end subroutine cinter
!
!     ==================
!     subroutine CINTER2
!     ==================
!
      subroutine cinter2(pfi,pslmi,pfo,pslmo,pgwi,pgwo,pxi,pxo          &
     &                  ,pyi,pyo,nxi,nxo,nyi,nyo,nc,npr)
!
      parameter(zerr=-9.E09)
!
      real :: fi(0:nxi,0:nyi+1)    ! input field
      real :: zzfi(nxi,nyi)
      real :: slmi(0:nxi,0:nyi+1)  ! input mask (0/1)
      real :: zzslmi(nxi,nyi)
      real :: fo(nxo,nyo)          ! output field
      real :: slmo(nxo,nyo)        ! output mask (0/1)
      real :: xi(0:nxi,0:nyi+1)    ! x's for input
      real :: xo(nxo,nyo)          ! x's for output
      real :: yi(0:nxi,0:nyi+1)    ! y's for input
      real :: yo(nxo,nyo)          ! y's for output
      real :: gwi(0:nxi,0:nyi+1)   ! input weights
!
      real :: pfi(nxi,nyi)
      real :: pslmi(nxi,nyi)
      real :: pfo(nxo,nyo)
      real :: pslmo(nxo,nyo)
      real :: pxi(nxi,nyi)
      real :: pxo(nxo,nyo)
      real :: pyi(nxi,nyi)
      real :: pyo(nxo,nyo)
      real :: pgwi(nxi,nyi)
      real :: pgwo(nxo,nyo)
!
      pi=2.*asin(1.)
!
      fi(1:nxi,1:nyi)=pfi(:,:)
      fi(0,1:nyi)=pfi(nxi,:)
      fi(:,0)=fi(:,1)
      fi(:,nyi+1)=fi(:,nyi)
      slmi(1:nxi,1:nyi)=pslmi(:,:)
      slmi(0,1:nyi)=pslmi(nxi,:)
      slmi(:,0)=slmi(:,1)
      slmi(:,nyi+1)=slmi(:,nyi)
      gwi(1:nxi,1:nyi)=pgwi(:,:)
      gwi(0,1:nyi)=pgwi(nxi,:)
      gwi(:,0)=gwi(:,1)
      gwi(:,nyi+1)=gwi(:,nyi)
      yi(1:nxi,1:nyi)=pyi(:,:)
      yi(0,1:nyi)=pyi(nxi,:)
      if(pyi(1,1) < pyi(1,nyi)) then
       yi(:,0)=-pi*0.5
       yi(:,nyi+1)=pi*0.5
      else
       yi(:,0)=pi*0.5
       yi(:,nyi+1)=-pi*0.5
      endif
      xi(1:nxi,1:nyi)=pxi(:,:)
      xi(0,1:nyi)=pxi(nxi,:)-2.*pi
      xi(:,0)=xi(:,1)
      xi(:,nyi+1)=xi(:,nyi)
      yo(:,:)=pyo(:,:)
      xo(:,:)=pxo(:,:)
      xo(:,:)=pxo(:,:)
      slmo(:,:)=pslmo(:,:)
!
      kmiss=0
      do jo=1,nyo
      do io=1,nxo
       zgw=0.
       zf=0.
       zxo=xo(io,jo)
       zyo=yo(io,jo)
       fo(io,jo)=zerr
       do ji=1,nyi+1
       do ii=1,nxi
        zxi=xi(ii,ji)
        zxil=xi(ii-1,ji)
        if(zxil > zxi) zxil=zxil-2.*pi
        if(ABS(zxil-zxi) > pi) then
         if(zxil > zxi) then
          zxil=zxil-2.*pi
         else
          zxi=zxi-2.*pi
         endif
        endif
        if(ABS(zxo-zxi) > pi) then
         if(zxo > zxi) then
          zxo=zxo-2.*pi
         else
          zxi=zxi-2.*pi
          zxil=zxil-2.*pi
         endif
        endif
        if(((zyo <= yi(ii,ji) .and. zyo > yi(ii,ji-1))                  &
     &      .or.(zyo > yi(ii,ji) .and. zyo <= yi(ii,ji-1)))             &
     &     .and.(zxo <= zxi .and. zxo > zxil)) then
         zgwx1=abs(zxo-zxi)*slmi(ii-1,ji)
         zgwx2=abs(zxo-zxil)*slmi(ii,ji)
         zgwx=zgwx1+zgwx2
         if(zgwx > 0.) then
          zf1=(zgwx1*fi(ii-1,ji)+zgwx2*fi(ii,ji))/zgwx
          zgwy1=abs(zyo-yi(ii,ji-1))                                    &
     &         *(gwi(ii,ji)*slmi(ii,ji)+gwi(ii-1,ji)*slmi(ii-1,ji))
         else
          zf1=0.
          zgwy1=0.
         endif
         zxi=xi(ii,ji-1)
         zxil=xi(ii-1,ji-1)
         if(zxil > zxi) zxil=zxil-2.*pi
         if(ABS(zxil-zxi) > pi) then
          if(zxil > zxi) then
           zxil=zxil-2.*pi
          else
           zxi=zxi-2.*pi
          endif
         endif
         if(ABS(zxo-zxi) > pi) then
          if(zxo > zxi) then
           zxo=zxo-2.*pi
          else
           zxi=zxi-2.*pi
           zxil=zxil-2.*pi
          endif
         endif
         zgwx1=abs(zxo-zxi)*slmi(ii-1,ji-1)
         zgwx2=abs(zxo-zxil)*slmi(ii,ji-1)
         zgwx=zgwx1+zgwx2
         if(zgwx > 0.) then
          zf2=(zgwx1*fi(ii-1,ji-1)+zgwx2*fi(ii,ji-1))/zgwx
          zgwy2=abs(zyo-yi(ii,ji))                                      &
     &         *(gwi(ii,ji-1)*slmi(ii,ji-1)                             &
     &          +gwi(ii-1,ji-1)*slmi(ii-1,ji-1))
         else
          zf2=0.
          zgwy2=0.
         endif
         zgwy=zgwy1+zgwy2
         if(zgwy > 0.) fo(io,jo)=(zf1*zgwy1+zf2*zgwy2)/zgwy
        endif
       enddo
       enddo
       if(fo(io,jo) == zerr .and. slmo(io,jo) > 0.) then 
        kmiss=kmiss+1
       else
        fo(io,jo)=fo(io,jo)*slmo(io,jo)
       endif
      enddo
      enddo
!
      if(npr > 0) then
       write(nud,*) ' '
       write(nud,*) 'In CINTER2: ',kmiss,' points to be extrapolated'
      endif
!
      if(kmiss > 0) call cfill(fo,slmo,nxo,nyo,zerr,npr)
!
      if(nc > 0) then
       zzfi(:,:) = fi(1:nxi,1:nyi)
       zzslmi(:,:) = slmi(1:nxi,1:nyi)
       call conserve(zzfi,fo,zzslmi,slmo,pgwi,pgwo,nxi,nxo,nyi,nyo,nc   &
     &               ,npr)
      endif
!
      pfo(:,:)=fo(:,:)
!
      return
      end subroutine cinter2
!
!     ================
!     subroutine CFILL
!     ================
!
      subroutine cfill(pf,pm,nx,ny,perr,npr)
!
!     complete the field by a simple extrapolation
!
      real :: pf(nx,ny),pm(nx,ny)
      real :: ppf(0:nx+1,0:ny+1),ppm(0:nx+1,0:ny+1)
!
      ppf(1:nx,1:ny)=pf(:,:)
      ppf(0,1:ny)=pf(nx,:)
      ppf(nx+1,1:ny)=pf(1,:)
      ppf(:,0)=ppf(:,1)
      ppf(:,ny+1)=ppf(:,ny)
      ppm(1:nx,1:ny)=pm(:,:)
      ppm(0,1:ny)=pm(nx,:)
      ppm(nx+1,1:ny)=pm(1,:)
      ppm(:,0)=ppm(:,1)
      ppm(:,ny+1)=ppm(:,ny)
!
      kmiss=1
      do while(kmiss > 0)
       kmiss=0
       do j=1,ny
       do i=1,nx
        if(ppf(i,j) == perr .or. ppm(i,j) < 0.5) then
         zf=0.
         zgw=0.
         if(ppf(i-1,j) /= perr) then
          zf=zf+ppf(i-1,j)*ppm(i-1,j)
          zgw=zgw+ppm(i-1,j)
         endif
         if(ppf(i+1,j) /= perr) then
          zf=zf+ppf(i+1,j)*ppm(i+1,j)
          zgw=zgw+ppm(i+1,j)
         endif
         if(ppf(i,j-1) /= perr) then
          zf=zf+ppf(i,j-1)*ppm(i,j-1)
          zgw=zgw+ppm(i,j-1)
         endif
         if(ppf(i,j+1) /= perr) then
          zf=zf+ppf(i,j+1)*ppm(i,j+1)
          zgw=zgw+ppm(i,j+1)
         endif
         if(zgw > 0.) then
          ppf(i,j)=zf/zgw
          ppm(i,j)=1.
         endif
         if(i==1) ppf(nx+1,j)=ppf(i,j)
         if(i==nx) ppf(0,j)=ppf(i,j)
         if(j==1) ppf(i,0)=ppf(i,j)
         if(j==ny) ppf(i,ny+1)=ppf(i,j)
        endif
        if(ppf(i,j) == perr) kmiss=kmiss+1
       enddo
       enddo
       if(npr > 0) then 
        write(nud,*) 'In CFILL: ',kmiss,' points still missed'
       endif
      enddo
!
      pf(:,:)=ppf(1:nx,1:ny)
!
      return
      end
!
!     ===================
!     subroutine CONSERVE
!     ===================
!
      subroutine conserve(pfi,pfo,pmi,pmo,pgwi,pgwo                     &
     &                   ,nxi,nxo,nyi,nyo,nc,npr)
!
!     correct global mean (multiplicative or additive depending on nc)
!
!     nc=1 : multiply f to give sum(fi*gwi)=sum(fon*gwo) with fon=f*fo
!     nc=2 : multiply f to give sum(fi*gwi)/sum(gwi)=sum(fon*gwo)/sum(gwo)
!     nc=3 : add a to give sum(fi*gwi)=sum(fon*gwo) with fon=fo+a
!     nc=4 : add a to give sum(fi*gwi)/sum(gwi)=sum(fon*gwo)/sum(gwo)
!
      real :: pfi(nxi,nyi)
      real :: pfo(nxo,nyo)
      real :: pmi(nxi,nyi)
      real :: pmo(nxo,nyo)
      real :: pgwi(nxi,nyi)
      real :: pgwo(nxo,nyo)
!
      zgwi=0.
      zfi=0.
      do j=1,nyi
      do i=1,nxi
       zgw=pgwi(i,j)*pmi(i,j)
       zfi=zfi+pfi(i,j)*zgw
       zgwi=zgwi+zgw
      enddo
      enddo
      zgwo=0.
      zfo=0.
      do j=1,nyo
      do i=1,nxo
       zgw=pgwo(i,j)*pmo(i,j)
       zfo=zfo+pfo(i,j)*zgw
       zgwo=zgwo+zgw
      enddo
      enddo
!
      if(nc < 3) then
       if(zfo == 0.) then
        if(zfi /= 0.) write(nud,*) 'WARNING (CPL): No Conservation (zfo=0.)'
        zfac1=1.
        zfac2=1.
       else
        zfac1=zfi/zfo
        zfac2=zfi/zfo*zgwo/zgwi
       endif
      endif
      zadd1=(zfi-zfo)/zgwo
      zadd2=zfi/zgwi-zfo/zgwo
!
      if(nc==1) pfo(:,:)=pfo(:,:)*zfac1
      if(nc==2) pfo(:,:)=pfo(:,:)*zfac2
      if(nc==3) pfo(:,:)=pfo(:,:)+zadd1*pmo(:,:)
      if(nc==4) pfo(:,:)=pfo(:,:)+zadd2*pmo(:,:)
!
      if(npr > 0) then
       if(nc==1) write(nud,*) 'In CONSERVE: nc= ',nc,' Fac= ',zfac1
       if(nc==2) write(nud,*) 'In CONSERVE: nc= ',nc,' Fac= ',zfac2
       if(nc==3) write(nud,*) 'In CONSERVE: nc= ',nc,' Add= ',zadd1
       if(nc==4) write(nud,*) 'In CONSERVE: nc= ',nc,' Add= ',zadd2
      endif
!
      return
      end
!
!     ===============
!     subroutine CGLM
!     ===============
!
      subroutine cglm(pf,pm,pgs,pgm,pgw,pgwi,nx,ny)
!
!     compute global means
!
      real :: pf(nx,ny)
      real :: pm(nx,ny)
      real :: pgwi(nx,ny)
!
      real :: zgwd,zgsd,zgw
!      
      zgwd=0.
      zgsd=0.
      do j=1,ny
      do i=1,nx
       zgw=pgwi(i,j)*pm(i,j)
       zgsd=zgsd+pf(i,j)*zgw
       zgwd=zgwd+zgw
      enddo
      enddo
      if(zgwd > 0.) then
       pgm=zgsd/zgwd
      else
       write(nud,*)'!WARNING! no sea points found in cglm2'
      endif
      pgs=zgsd
      pgw=zgwd
!
      return
      end subroutine cglm
!
!     ====================
!     subroutine CFLUKOINI
!     ====================
!
      subroutine cflukoini
      use cplmod
!
      integer ih(8)
      real :: zin(nxo,nyot) = 0.
      logical lopen
      character (len=20) :: yformat = "(8E12.6)"
!
      do jt=30,99
       ntape=jt
       inquire(unit=ntape,opened=lopen)
       if(.not. lopen) exit
      enddo
!
      open(ntape,file=trim(flukofile),form='formatted')
      jflfresh=0
      jflsst=0
      jfltaux=0
      jfltauy=0
      jflice=0
      jfloheat=0
      do
       read(ntape,'(8I10)',end=1001) ih
       read(ntape,yformat) zin(:,3:74)
       jmon=ih(3)
       if(jmon >= 100) jmon=(jmon-(jmon/10000)*10000)/100  ! oroginal service
       if(ih(1)==75 .or. ih(1)== -82) then
        flukofresh(:,:,jmon)=zin(:,:)
        jflfresh=jflfresh+1
       elseif(ih(1)==73 .or. ih(1)== -83) then
        flukotaux(:,:,jmon)=zin(:,:)
        jfltaux=jfltaux+1
       elseif(ih(1)==74 .or. ih(1)== -84) then
        flukotauy(:,:,jmon)=zin(:,:)
        jfltauy=jfltauy+1
       elseif(ih(1)==72 .or. ih(1)== -80) then
        flukosst(:,:,jmon)=zin(:,:)
        jflsst=jflsst+1
       elseif(ih(1)==78) then
        flukooheat(:,:,jmon)=zin(:,:)
        jfloheat=jfloheat+1
       elseif(ih(1)==79 .or. ih(1)== -81) then
        flukoice(:,:,jmon)=zin(:,:)
        jflice=jflice+1
       endif
      enddo
 1001 continue
      close(ntape)
!
      write(nud,*)' '
      write(nud,*)'Reading Flux corrections'
      write(nud,*)jflsst,' months sst correction found in ',trim(flukofile)
      write(nud,*)jfltaux,' months taux correction found in ',trim(flukofile)
      write(nud,*)jfltauy,' months tauy correction found in ',trim(flukofile)
      write(nud,*)jflice,' months ice correction found in ',trim(flukofile)
      write(nud,*)jflfresh,' months pme correction found in ',trim(flukofile)
      write(nud,*)jfloheat,' months heat correction found in ',trim(flukofile)
      write(nud,*)' '
!
      return
!
      end subroutine cflukoini
!
!     ==================
!     subroutine CFLUKOA
!     ==================
!
      subroutine cflukoa(psst,ptaux,ptauy,pice,pfresh)
      use cplmod
!
      real :: psst(nxo,nyot),ptaux(nxo,nyot),ptauy(nxo,nyot)
      real :: pice(nxo,nyot),pfresh(nxo,nyot)
!
      jmon=ndatim(2)
!
      if(nflsst==1) then
       psst(:,:)=psst(:,:)+flukosst(:,:,jmon)
      endif       
!
      if(nflice==1) then
       pice(:,:)=pice(:,:)+flukoice(:,:,jmon)
      endif       
!
      if(nfltau==1) then
       ptaux(:,:)=ptaux(:,:)+flukotaux(:,:,jmon)
       ptauy(:,:)=ptauy(:,:)+flukotauy(:,:,jmon)
      endif       
!
      if(nflfresh==1) then
       pfresh(:,:)=pfresh(:,:)+flukofresh(:,:,jmon)
      endif       
!
      return
      end subroutine cflukoa
!
!     ==================
!     subroutine CFLUKOB
!     ==================
!
      subroutine cflukob(poheat)
      use cplmod
!
      real :: poheat(nxo,nyot)
!
      jmon=ndatim(2)
!
      if(nfloheat==1) then
       poheat(:,:)=poheat(:,:)+flukooheat(:,:,jmon)
      endif       
!
      return
      end subroutine cflukob
!
!     ================
!     subroutine CPFAC
!     ================
!
      subroutine cpfac(pfsst,pftau,pffresh,pfice)
      use cplmod
!
      real :: pfsst,pftau,pffresh,pfice
!
      pfsst=cfacsst 
      pftau=cfactau
      pffresh=cfacfresh
      pfice=cfacice
!
      return
      end subroutine cpfac
!
!     =================
!     subroutine CROFFC
!     =================
!
      subroutine croffc(proff)
      use cplmod
!
      real :: proff(nxa,nya)
      real :: zroff(nxa,nya)
!     
      if(nroffc == 1) then
!
!     substitute lokal runoff by global mean
!
       call cglm(proff,aslm,zgs,zgm,zgw,agw,nxa,nya)
       proff(:,:)=zgm*aslm(:,:)
      endif
!
      if(nroffc == 2) then
!
!     substitute lokal runoff by area averages
!     where iroffi==0 local runoff is kept
!
!!     kroff=0
       nroffa=maxval(iroffi(:,:))
       do ja=1,nroffa
        zroffa=SUM(proff(:,:)*agw(:,:)*aslm(:,:)                        &
     &            ,mask=(ja == iroffi(:,:)))
        zgwa=SUM(agw(:,:)*aslm(:,:),mask=(ja == iroffi(:,:))) 
!!      kroff=kroff+SUM(nint(aslm(:,:)),mask=(ja == iroffi(:,:)))
        if(zgwa <= 0.) then
         zgwa=-1.
         zroffa=0.
        endif
        where(ja == iroffi(:,:))
         proff(:,:)=zroffa/zgwa*aslm(:,:)
        endwhere      
       enddo
!!     klsm=SUM(nint(aslm(:,:)))
!!     if(kroff /= klsm) then
!!      write(nud,*)'WARNING: runoff corrected only for ',kroff,' from '     &
!!   &        ,klsm,' gridpoints!'
!!     endif
      elseif(nroffc == 3) then
       zroff(:,:)=0.
!
!     redistribute lokal runoff using area averages
!     where iroffi==0 local runoff + area average is used
!
       where(iroffi(:,:) == 0) 
        zroff(:,:)=proff(:,:)
       endwhere
       nroffa=maxval(iroffi(:,:))
!!     kroffi=0
       do ja=1,nroffa
        zroffa=SUM(proff(:,:)*agw(:,:)*aslm(:,:)                        &
     &            ,mask=(ja == iroffi(:,:)))
        zgwa=SUM(agw(:,:)*aslm(:,:),mask=(ja == iroffo(:,:))) 
!!      kroffi=kroffi+SUM(nint(aslm(:,:)),mask=(ja == iroffi(:,:)))
        if(zgwa <= 0. .and. zroffa > 0.) then
         write(nud,*)'ERROR: zero area for runoff output no ',ja
         stop
        endif
        if(zgwa > 0.) then
         where(ja == iroffo(:,:))
          zroff(:,:)=zroffa/zgwa*aslm(:,:)+zroff(:,:)
         endwhere    
        endif  
       enddo
!!     klsm=SUM(nint(aslm(:,:)))
!!     if(kroffi /= klsm) then
!!      write(nud,*)'ERROR: not all runoff redistributed!',klsm,kroffi,nroffa 
!!      stop    
!!     endif
       proff(:,:)=zroff(:,:)
      endif
!
      return
      end subroutine croffc

