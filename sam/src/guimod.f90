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

      call initgui(model,nguidbg,NLAT,-1,1)

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

      if (mypid == NROOT) then
         parc(1)   = disp
         parc(2)   = 0
         parc(3)   = 0
         parc(4)   = 0
         parc(5)   = 0
         crap(:)   = parc(:)
         idatim(:) = ndatim(:)
         nshutdown = iguistep(parc,idatim)    ! GUI event handler
         if (parc(1) /= crap(1)) call change_disp(1)
      endif
   
      call mpbci(nshutdown)
      call mpbcr(disp)
   
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
      call alt2reg(z,klev)
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
      real :: f(NLON,NLPP)
      real :: z(NLON,NLAT)
      real (kind=4) :: x(NLON+1,NLAT)

      if (ngui == 0) return
      
      call mpgagp(z,f,1)
      if (mypid == NROOT) then
         call alt2reg(z,1)
         mlon = NLON/2
         do jlat = 1 , NLAT
            x(     1:mlon  ,jlat) = z(mlon+1:NLON  ,jlat) * CV * rcs(jlat)
            x(mlon+1:NLON+1,jlat) = z(     1:mlon+1,jlat) * CV * rcs(jlat)
         enddo
         call guiput(yname,x,NLON+1,NLAT,NLEV)
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
