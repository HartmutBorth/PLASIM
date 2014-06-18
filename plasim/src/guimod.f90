! *********************************************
! * GUIMOD - Graphical User Interface Routines*
! * Planet Simulater Version                  *
! * 21-Oct-2010 - Edilbert Kirk               *
! *********************************************
! * This file contains interface routines for *
! * communication between the model PLASIM    *
! * and the GUI routines in file "pumax.c".   *
! * <pumax.c> is the same for all models, but *
! * each model has its own version of         *
! * of  <guimod.f90>                          *
! *********************************************


! ===================
! SUBROUTINE GUISTART
! ===================

subroutine guistart
use pumamod
character (80) :: ypl

if (naqua /= 0) then
   ypl = "Aqua " // trim(yplanet) // char(0)
else
   ypl = trim(yplanet) // char(0)
endif

if (ngui == 0) return

call initgui(model,nguidbg,NLAT,mrpid,mrnum,ypl)

return
end subroutine guistart


! ==================
! SUBROUTINE GUISTOP
! ==================

subroutine guistop
use pumamod

if (mypid == NROOT .and. ngui > 0) call guiclose

return
end subroutine guistop


! =========================
! SUBROUTINE GUISTEP_PLASIM
! =========================

subroutine guistep_plasim
use radmod

interface 
   integer(kind=4) function iguistep(p2,k3)
      real   (kind=4), intent(inout) :: p2(*)
      integer(kind=4), intent(in)    :: k3(6)
   end function iguistep
end interface

integer(kind=4) :: idatim(7)

nsyncold = nsync
call mrsum(nsyncold)

if (mypid == NROOT) then
   parc(1)   = co2 
   parc(2)   = gsol0
   parc(3)   = (sellon-1) * 360.0 / NLON
   parc(4)   = nsync
   parc(5)   = syncstr

   crap(:)   = parc(:)
   idatim(:) = ndatim(:)
   nshutdown = iguistep(parc,idatim)    ! GUI event handler
   if (parc(1) /= crap(1)) call change_co2(1)
   if (parc(2) /= crap(2)) call change_gsol0(2)
   if (parc(3) /= crap(3)) call change_sellon(3)
   if (parc(4) /= crap(4)) call change_nsync(4)
   if (parc(5) /= crap(5)) call change_syncstr(5)
endif
   
call mrsum(nshutdown)   ! Any instance can signal shutdown
call mpbci(nshutdown)
call mpbcr(co2)
call mpbcr(gsol0)
call mpbci(sellon)
nsyncnew = nsync
call mrsum(nsyncnew)    ! Any instance can switch NSYNC
if (nsyncnew > nsyncold) nsync = 1
if (nsyncnew < nsyncold) nsync = 0
   
return
end subroutine guistep_plasim


! ================
! SUBROUTINE GUIPS
! ================

subroutine guips(f)
use pumamod
real :: f(NLON,NLAT)
real (kind=4) :: x(NLON+1,NLAT)

if (ngui == 0) return
if (mypid == NROOT) then
   mlon = NLON / 2
   x(     1:mlon  ,:) = f(mlon+1:NLON  ,:) * PSURF * 0.01 ! [hPa]
   x(mlon+1:NLON+1,:) = f(     1:mlon+1,:) * PSURF * 0.01 ! [hPa]
   call guiput("GP" // char(0) ,x,NLON+1,NLAT,1)
endif
return
end


! =================
! SUBROUTINE GUIHOR
! =================

subroutine guihor(yname,f,klev,pm,pa)
use pumamod
character (len=*) :: yname
real :: f(NLON,NLPP,klev)
real :: z(NLON,NLAT,klev)
real (kind=4) :: x(NLON+1,NLAT,klev)

if (ngui == 0) return

call mpgagp(z,f,klev)
if (mypid == NROOT) then
   mlon = NLON / 2
   x(     1:mlon  ,:,:) = z(mlon+1:NLON  ,:,:) * pm + pa
   x(mlon+1:NLON+1,:,:) = z(     1:mlon+1,:,:) * pm + pa
   call guiput(yname,x,NLON+1,NLAT,klev)
endif
return
end


! ================
! SUBROUTINE GUIGV
! ================

subroutine guigv(yname,f)
use pumamod
character (len=*) :: yname
real :: f(NLON,NLPP,NLEV)
real :: z(NLON,NLAT,NLEV)
real (kind=4) :: x(NLON+1,NLAT,NLEV)

if (ngui == 0) return

call mpgagp(z,f,NLEV)
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


! ===================
! SUBROUTINE GUIGVCOL
! ===================

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


! ===================
! SUBROUTINE GUIGTCOL
! ===================

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


! ====================
! SUBROUTINE GUID3DCOL
! ====================

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

! =======================
! SUBROUTINE CHANGE_NSYNC
! =======================

subroutine change_nsync(k)
use pumamod
write(nud,7000) nstep,'NSYNC',crap(k),parc(k)
if (parc(k) > crap(k)) then
   nsync = 1
else
   nsync = 0
endif
return
 7000 format('Step',i8,': User changed ',a,' from ',f6.2,' to ',f6.2)
end

! =========================
! SUBROUTINE CHANGE_SYNCSTR
! =========================

subroutine change_syncstr(k)
use pumamod
write (*,7000) nstep,'SYNCSTR',crap(k),parc(k)
syncstr = parc(k)
return
 7000 format('Step',i8,': User changed ',a,' from ',f6.2,' to ',f6.2)
end 

! ======================
! SUBROUTINE CHANGE_DISP
! ======================

subroutine change_disp(k)
use pumamod
write(nud,7000) nstep,'DISP  ',crap(k),parc(k)
disp = parc(k)
return
 7000 format('Step',i8,': User changed ',a,' from ',f6.2,' to ',f6.2)
end


! ======================
! SUBROUTINE CHANGE_DTEP
! ======================

subroutine change_dtep(k)
use pumamod
write(nud,7000) nstep,'DTEP  ',crap(k),parc(k)
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


! ======================
! SUBROUTINE CHANGE_DTNS
! ======================

subroutine change_dtns(k)
use pumamod
write(nud,7000) nstep,'DTNS  ',crap(k),parc(k)
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


! =====================
! SUBROUTINE CHANGE_CO2
! =====================

subroutine change_co2(k)
use pumamod
write(nud,7000) nstep,'CO2   ',crap(k),parc(k)
co2  = parc(k)
return
 7000 format('Step',i8,': User changed ',a,' from ',f7.2,' to ',f7.2)
end


! =======================
! SUBROUTINE CHANGE_GSOL0
! =======================

subroutine change_gsol0(k)
use pumamod
use radmod
write(nud,7000) nstep,'GSOL0 ',crap(k),parc(k)
gsol0  = parc(k)
return
 7000 format('Step',i8,': User changed ',a,' from ',f7.2,' to ',f7.2)
end


! ======================
! SUBROUTINE CHANGE_DAWN
! ======================

subroutine change_dawn(k)
use radmod
write(nud,7000) nstep,'DAWN ',crap(k),parc(k)
dawn  = parc(k)
return
 7000 format('Step',i8,': User changed ',a,' from ',f7.2,' to ',f7.2)
end



! ========================
! SUBROUTINE CHANGE_SELLON
! ========================

subroutine change_sellon(k)
use pumamod
sellon  = nint(parc(k) * NLON / 360.0 + 1.0)
if (sellon <    1) sellon = 1
if (sellon > NLON) sellon = NLON
return
end


! ====================
! SUBROUTINE GUIHORLSG
! ====================

subroutine guihorlsg(yname,f,mask,nxoce,nyoce,klev,pm,pa)
use pumamod
real (kind=4) :: pm,pa
real (kind=4) :: f(nxoce,nyoce,klev)
real (kind=4) :: x(nxoce+1,nyoce,klev)
real (kind=4) :: mask(nxoce,nyoce,klev)
integer :: klev
character (len=*) :: yname
logical :: even

if (ngui == 0) return
if (mypid == NROOT) then
 x(:,:,:) = 0.0
 do jlev = 1 , klev
  even = .false.
  do jlat = 1 , nyoce
   do jlon = 1 , nxoce + 1
    ilon = jlon + nxoce / 2 ! rotate from (0:360) to (-180:180)
    if (ilon > nxoce) ilon = ilon - nxoce
    ilon = ilon + (jlat-1) / 2 ! correct arakawa grid shift
    if (ilon > nxoce) ilon = ilon - nxoce
    ilp1 = ilon + 1
    if (ilp1 > nxoce) ilp1 = 1
    if (mask(ilon,jlat,jlev) >= 1.0) &
     x(jlon,jlat,jlev) = f(ilon,jlat,jlev) * pm + pa
    if (even.and.mask(ilon,jlat,jlev)>=1.0.and.mask(ilp1,jlat,jlev)>=1.0) &
     x(jlon,jlat,jlev) = 0.5*(f(ilon,jlat,jlev)+f(ilp1,jlat,jlev))*pm+pa
   enddo ! jlon
   even = .not. even
  enddo ! jlat
 enddo ! jlev
 call guiput(yname,x,nxoce+1,nyoce,klev)
endif
return
end
          
