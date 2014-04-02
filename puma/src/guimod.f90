!     ***************************************************
!     * GUIMOD - Graphical User Interface Routines      *
!     * 21-Sep-2006 - Edilbert Kirk                     *
!     ***************************************************
!     * This file contains all interface routines for   *
!     * communication between model (PUMA or PLASIM)    *
!     * and the GUI routines in file "pumax.c".         *
!     * This file is identical for both models, if you  *
!     * make changes, make sure, that these either      *
!     * affect both models in a proper way or use the   *
!     * if (model == PUMA)    ...  endif                *
!     * if (model == PLASIM)  ...  endif                *
!     * statements. After changes copy the new version, *
!     * so that ../puma/src/guimod.f90  and             *
!     *       ../plasim/src/guimod.f90 are identical.   *
!     ***************************************************

      subroutine guistart
      use pumamod

      if (ngui == 0) return

      call initgui(model,nguidbg,NLAT,mrpid,mrnum,trim(yplanet)//char(0))

      return
      end subroutine guistart

!     ==================
!     SUBROUTINE GUISTOP
!     ==================

      subroutine guistop
      use pumamod

      if (mypid == NROOT .and. ngui > 0) call guiclose

      return
      end subroutine guistop

!     =======================
!     SUBROUTINE GUISTEP_PUMA
!     =======================

      subroutine guistep_puma
      use pumamod

      interface 
         integer(kind=4) function iguistep(parc,idatim)
            real   (kind=4), intent(inout) :: parc(*)
            integer(kind=4), intent(in)    :: idatim(6)
         end function iguistep
      end interface

      integer (kind=4) idatim(6)

      nsyncold = nsync
      call mrsum(nsyncold)

      if (mypid == NROOT) then
         parc(1)   = disp
         parc(2)   = dtep * CT;
         parc(3)   = dtns * CT;
         parc(4)   = nsync
         parc(5)   = syncstr
         crap(:)   = parc(:)
         idatim(:) = ndatim(:)
         nshutdown = iguistep(parc,idatim)    ! GUI event handler
         if (parc(1) /= crap(1)) call change_disp(1)
         if (parc(2) /= crap(2)) call change_dtep(2)
         if (parc(3) /= crap(3)) call change_dtns(3)
         if (parc(4) /= crap(4)) call change_nsync(4)
         if (parc(5) /= crap(5)) call change_syncstr(5)
      endif
   
      call mrsum(nshutdown)   ! Any instance can signal shutdown
      call mpbci(nshutdown)
      call mpbcr(disp)

      nsyncnew = nsync
      call mrsum(nsyncnew)    ! Any instance can switch NSYNC
      if (nsyncnew > nsyncold) nsync = 1 
      if (nsyncnew < nsyncold) nsync = 0 
   
      call mpbcl(ldtep)
      if (ldtep) then
         call mpbcr(dtep)            ! broadcast changed value
         call mpscsp(sr1,srp1,NLEV)  ! scatter   changed array
         ldtep = .false.
      endif 
   
      call mpbcl(ldtns)
      if (ldtns) then
         call mpbcr(dtns)            ! broadcast changed value
         call mpscsp(sr2,srp2,NLEV)  ! scatter   changed array
         ldtns = .false.
      endif 
   
      return
      end subroutine guistep_puma

!     ================
!     SUBROUTINE GUIPS
!     ================

      subroutine guips(f,pmean)
      use pumamod
      real :: f(NLON,NLAT)
      real :: z(NLON,NLAT)
      real (kind=4) :: x(NLON+1,NLAT)

      if (ngui == 0) return
      if (mypid == NROOT) then
         z(:,:) = f(:,:)
         call alt2reg(z,1)
         mlon = NLON/2
         pm   = pmean * 0.01 ! [hPa]
         x(     1:mlon  ,:) = z(mlon+1:NLON  ,:) * pm
         x(mlon+1:NLON+1,:) = z(     1:mlon+1,:) * pm
         call guiput("GP" // char(0) ,x,NLON+1,NLAT,1)
      endif
      return
      end

!     =================
!     SUBROUTINE GUIHOR
!     =================

      subroutine guihor(yname,f,klev,pm,pa)
      use pumamod
      character (len=*) :: yname
      real :: f(NLON,NLPP,klev)
      real :: z(NLON,NLAT,klev)
      real (kind=4) :: x(NLON+1,NLAT,klev)

      if (ngui == 0) return

      ! Incoming array f stores longitudes from 0 deg to (360 - delta lambda)
      ! GUI gets rotated array x stored from -180 deg to +180 deg

      call mpgagp(z,f,klev)
      if (mypid == NROOT) then
         mlon = NLON/2
         x(     1:mlon  ,:,:) = z(mlon+1:NLON  ,:,:) * pm + pa
         x(mlon+1:NLON+1,:,:) = z(     1:mlon+1,:,:) * pm + pa
         call guiput(yname,x,NLON+1,NLAT,klev)
      endif
      return
      end

!     ================
!     SUBROUTINE GUIGV
!     ================

      subroutine guigv(yname,f)
      use pumamod
      character (len=*) :: yname
      real :: f(NLON,NLPP,NLEV)
      real :: z(NLON,NLAT,NLEV)
      real (kind=4) :: x(NLON+1,NLAT,NLEV)

      if (ngui == 0) return
      
      call mpgagp(z,f,NLEV)
      if (model == PUMA) call alt2reg(z,NLEV)
      if (mypid == NROOT) then
         mlon = NLON/2
         do jlat = 1 , NLAT
            x(     1:mlon  ,jlat,:) = z(mlon+1:NLON  ,jlat,:) * CV * rcs(jlat)
            x(mlon+1:NLON+1,jlat,:) = z(     1:mlon+1,jlat,:) * CV * rcs(jlat)
         enddo
         call guiput(yname,x,NLON+1,NLAT,NLEV)
      endif
      return
      end

!     ================
!     SUBROUTINE GUIGT
!     ================

      subroutine guigt(f)
      use pumamod
      real :: f(NLON,NLPP,NLEV)
      real :: z(NLON,NLAT,NLEV)
      real (kind=4) :: x(NLON+1,NLAT,NLEV)

      if (ngui == 0) return
      
      call mpgagp(z,f,NLEV)
      if (model == PUMA) call alt2reg(z,NLEV)
      if (mypid == NROOT) then
         mlon = NLON/2
         do jlon = 1 , mlon
            do jlat = 1 , NLAT
               x(jlon,jlat,:) =  (z(jlon+mlon,jlat,:) + t0(:))*CT - 273.16
            enddo
         enddo
         do jlon = mlon+1,NLON+1
            do jlat = 1 , NLAT
               x(jlon,jlat,:) =  (z(jlon-mlon,jlat,:) + t0(:))*CT - 273.16
            enddo
         enddo
         call guiput("GT" // char(0),x,NLON+1,NLAT,NLEV)
      endif
      return
      end

!     ===================
!     SUBROUTINE GUIGVCOL
!     ===================

      subroutine guigvcol(yname,f,klon)
      use pumamod
      character (len=*) :: yname
      real :: f(NLON,NLPP,NLEV)
      real :: z(NLON,NLAT,NLEV)
      real (kind=4) :: x(NLEV,NLAT)

      if (ngui == 0) return

      call mpgagp(z,f,NLEV)
      if (mypid == NROOT) then
         do jlat = 1 , NLAT
          do jlev = 1 , NLEV
            x(jlev,jlat) = z(klon,jlat,jlev) * CV * rcs(jlat)
          enddo
         enddo
         call guiput(yname,x,NLEV,NLAT,1)
      endif
      return
      end

!     ===================
!     SUBROUTINE GUIGTCOL
!     ===================

      subroutine guigtcol(f,klon)
      use pumamod
      real :: f(NLON,NLPP,NLEV)
      real :: z(NLON,NLAT,NLEV)
      real (kind=4) :: x(NLEV,NLAT)

      if (ngui == 0) return

      call mpgagp(z,f,NLEV)
      if (mypid == NROOT) then
         do jlat = 1 , NLAT
          do jlev = 1 , NLEV
            x(jlev,jlat) = z(klon,jlat,jlev) - TMELT
          enddo
         enddo
         call guiput("GTCOL" // char(0),x,NLEV,NLAT,1)
      endif
      return
      end

!     ====================
!     SUBROUTINE GUID3DCOL
!     ====================

      subroutine guid3dcol(yname,f,klon,klev,pm,pa)
      use pumamod
      character (len=*) :: yname
      real :: f(NLON,NLPP,klev)
      real :: z(NLON,NLAT,klev)
      real (kind=4) :: x(NLEV,NLAT)

      if (ngui == 0) return

      call mpgagp(z,f,klev)
      if (mypid == NROOT) then
         do jlat = 1 , NLAT
          do jlev = 1 , NLEV
            x(jlev,jlat) = z(klon,jlat,jlev)*pm + pa
          enddo
         enddo
         call guiput(yname,x,NLEV,NLAT,1)
      endif
      return
      end

!     =======================
!     SUBROUTINE CHANGE_NSYNC
!     =======================

      subroutine change_nsync(k)
      use pumamod
      write (*,7000) nstep,'NSYNC',crap(k),parc(k)
      nsync = parc(k) + 0.001
      return
 7000 format('Step',i8,': User changed ',a,' from ',f6.2,' to ',f6.2)
      end

!     =========================
!     SUBROUTINE CHANGE_SYNCSTR
!     =========================

      subroutine change_syncstr(k)
      use pumamod
      write (*,7000) nstep,'SYNCSTR',crap(k),parc(k)
      syncstr = parc(k)
      return
 7000 format('Step',i8,': User changed ',a,' from ',f6.2,' to ',f6.2)
      end

!     ======================
!     SUBROUTINE CHANGE_DISP
!     ======================

      subroutine change_disp(k)
      use pumamod
      write (*,7000) nstep,'DISP  ',crap(k),parc(k)
      disp = parc(k)
      return
 7000 format('Step',i8,': User changed ',a,' from ',f6.2,' to ',f6.2)
      end

!     ======================
!     SUBROUTINE CHANGE_DTEP
!     ======================

      subroutine change_dtep(k)
      use pumamod
      write (*,7000) nstep,'DTEP  ',crap(k),parc(k)
      dtep = parc(k) / CT;
      zttrop = tgr-dtrop*ALR
      ztps   = (zttrop/tgr)**(GA/(ALR*GASCON))
      do jlev = 1 , NLEV
         zfac = sin(0.5*PI*(sigma(jlev)-ztps)/(1.-ztps))
         if (zfac < 0.0) zfac = 0.0
         sr1(5,jlev) = -2.0/3.0 * sqrt(0.4) * dtep * zfac
      enddo
      ldtep = .true.
      return
 7000 format('Step',i8,': User changed ',a,' from ',f7.2,' to ',f7.2)
      end

!     ======================
!     SUBROUTINE CHANGE_DTNS
!     ======================

      subroutine change_dtns(k)
      use pumamod
      write (*,7000) nstep,'DTNS  ',crap(k),parc(k)
      dtns = parc(k) / CT;
      zttrop = tgr-dtrop*ALR
      ztps   = (zttrop/tgr)**(GA/(ALR*GASCON))
      do jlev = 1 , NLEV
         zfac = sin(0.5*PI*(sigma(jlev)-ztps)/(1.-ztps))
         if (zfac < 0.0) zfac = 0.0
         sr2(3,jlev) = (1.0 / sqrt(6.0)) * dtns * zfac
      enddo
      ldtns = .true.
      return
 7000 format('Step',i8,': User changed ',a,' from ',f7.2,' to ',f7.2)
      end


!     ========================
!     SUBROUTINE CHANGE_SELLON
!     ========================

      subroutine change_sellon(k)
      use pumamod
      sellon  = nint(parc(k) * NLON / 360.0 + 1.0)
      if (sellon <    1) sellon = 1
      if (sellon > NLON) sellon = NLON
      return
      end
