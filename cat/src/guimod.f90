! ========
! GUISTART
! ========

subroutine guistart
use catmod
implicit none

integer :: model   =  4 ! CAT 
integer :: nlat    = 32 ! dummy value
integer :: mrpid   = -1 ! # process ID
integer :: mrnum   =  1 ! # of instances

character (80) :: yplanet = "CAT" ! dummy

if (ngui == 0) return

call initgui(model,nguidbg,nlat,mrpid,mrnum,trim(yplanet)//char(0))

return
end subroutine guistart

! =======
! GUISTOP
! =======

subroutine guistop
implicit none

integer :: ngui = 1 ! set by caller

if (ngui > 0) call guiclose

return
end subroutine guistop

! ============
! GUISTEP_CAT
! ============

subroutine guistep_cat
use catmod
implicit none

interface 
   integer(kind=4) function iguistep(parc,idatim)
      real   (kind=4), intent(inout) :: parc(*)
      integer(kind=4), intent(in)    :: idatim(6)
   end function iguistep
end interface

integer (kind=4) idatim(6)

parc(1)   = ngui
! parc(2)   = qmax * 1.0e7
parc(3)   = 3.0
parc(4)   = 4.0
parc(5)   = 5.0

idatim(:) = 0
idatim(1) = tstep
nshutdown = iguistep(parc,idatim)    ! GUI event handler

if (parc(1) >= 1.0) ngui = parc(1)
! qmax = parc(2) * 1.0e-7
   
return
end subroutine guistep_cat
