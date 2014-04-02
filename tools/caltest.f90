!     =============
!     MODULE CALMOD
!     =============

      module calmod
      implicit none
      integer :: mondays(0:12) = (/0,31,28,31,30,31,30,31,31,30,31,30,31/)
      integer :: monaccu(0:12)
      integer :: ny400d = 400 * 365 + 97
      integer :: ny100d = 100 * 365 + 24
      integer :: ny004d =   4 * 365 +  1
      integer :: ny001d =       365
      
      end module calmod


!     =================
!     FUNCTION NWEEKDAY
!     =================

      integer function nweekday(kday)
      nweekday = mod(kday+6,7)
      return
      end

!     ===================
!     SUBROUTINE STEP2CAL
!     ===================

      subroutine step2cal(kstep,ktspd,kyea,kmon,kday,khou,kmin)
      use calmod
      implicit none
      integer, intent(IN ) :: kstep ! time step since simulation start
      integer, intent(IN ) :: ktspd ! time steps per day
      integer, intent(OUT) :: kyea  ! current year   of simulation
      integer, intent(OUT) :: kmon  ! current month  of simulation
      integer, intent(OUT) :: kday  ! current day    of simulation
      integer, intent(OUT) :: khou  ! current hour   of simulation
      integer, intent(OUT) :: kmin  ! current minute of simulation

      integer :: idall
      integer :: istp
      integer :: iy400,id400
      integer :: iy100,id100
      integer :: iy004,id004
      integer :: iy001,id001
      integer :: jmon

      logical :: leap

      idall = kstep   / ktspd

      iy400 = idall   / ny400d          ! segment [of 400 years]
      id400 = mod(idall,ny400d)

      if (id400 <= ny100d) then         ! century year is leap year
         iy100 =         0              ! century in segment [0]
         id100 =     id400
         iy004 =     id100 / ny004d     ! tetrade in century [0..24]
         id004 = mod(id100 , ny004d)
         leap  = (id004 <= ny001d)
         if (leap) then
            iy001 =     0               ! year in tetrade [0]
            id001 = id004
         else
            iy001 =    (id004-1)/ny001d ! year in tetrade [1,2,3]
            id001 = mod(id004-1, ny001d)
         endif
      else                              ! century year is not leap year
         iy100 =    (id400-1)/ny100d    ! century in segment [1,2,3]
         id100 = mod(id400-1, ny100d)
         if (id100 < ny004d-1) then
            iy004 = 0                   ! tetrade in century [0]
            id004 = id100
            leap  = .false.
            iy001 =     id004/ny001d    ! year in tetrade [1,2,3]
            id001 = mod(id004,ny001d)
         else
            iy004 = (id100+1)/ny004d    ! tetrade in century [1..24]
            id004 = mod(id100+1,ny004d)
            leap  = (id004 <= ny001d)
            if (leap) then
               iy001 =     0            ! year in tetrade [0]
               id001 = id004
            else
               iy001 =    (id004-1)/ny001d
               id001 = mod(id004-1, ny001d)
            endif
         endif
      endif

      kyea  = iy400 * 400 + iy100 * 100 + iy004 * 4 + iy001
!     print *,"yea ",kyea,iy400,iy100,iy004,iy001,id004,ny001d

      monaccu(0) = mondays(0)
      monaccu(1) = mondays(1)
      monaccu(2) = mondays(1) + mondays(2)
      if (leap) monaccu(2) = monaccu(2) + 1
      do jmon = 3 , 12
         monaccu(jmon) = monaccu(jmon-1) + mondays(jmon)
      enddo
      kmon = 1
      id001 = id001 + 1          
      do while (id001 > monaccu(kmon))
         kmon = kmon + 1
      enddo
      kday = id001 - monaccu(kmon-1)

      istp = mod(kstep,ktspd)
      kmin = (istp * 1440) / ktspd
      khou = kmin / 60
      kmin = mod(kmin,60)

      return
      end subroutine step2cal

!     ===================
!     SUBROUTINE CAL2STEP
!     ===================

      subroutine cal2step(kstep,ktspd,kyea,kmon,kday,khou,kmin)
      use calmod
      implicit none
      integer, intent(OUT) :: kstep ! time step since simulation start
      integer, intent(IN ) :: ktspd ! time steps per day
      integer, intent(IN ) :: kyea  ! current year   of simulation
      integer, intent(IN ) :: kmon  ! current month  of simulation
      integer, intent(IN ) :: kday  ! current day    of simulation
      integer, intent(IN ) :: khou  ! current hour   of simulation
      integer, intent(IN ) :: kmin  ! current minute of simulation

      integer :: idall
      integer :: ilp
      integer :: iy400,id400
      integer :: iy100,id100
      integer :: iy004,id004
      integer :: jmon

      logical :: leap

      iy400 = kyea   /  400    ! segment [400]
      id400 = mod(kyea ,400)   ! year in segment [0..399]
      iy100 = id400   / 100    ! century [0,1,2,3]
      id100 = mod(id400,100)   ! year in century [0..99]
      iy004 = id100   /   4    ! tetrade [0..24]
      id004 = mod(id100,  4)   ! year in tetrade [0,1,2,3]

      leap  = (id004 == 0 .and. (id100 /= 0 .or. id400 == 0))

      ilp = -1
      if (id004 > 0) ilp = ilp + 1
      if (iy100 > 0 .and. id100 == 0 ) ilp = ilp + 1

      monaccu(0) = mondays(0)
      monaccu(1) = mondays(1)
      monaccu(2) = mondays(1) + mondays(2)
      if (leap) monaccu(2) = monaccu(2) + 1
      do jmon = 3 , 12
         monaccu(jmon) = monaccu(jmon-1) + mondays(jmon)
      enddo

      idall = iy400 * ny400d + iy100 * ny100d + iy004 * ny004d &
            + id004 * ny001d + monaccu(kmon-1)+ kday + ilp
      kstep = ktspd * idall + (ktspd * (khou * 60 + kmin)) / 1440

      return
      end subroutine cal2step

!     =====================
!     SUBROUTINE STEP2CAL30
!     =====================

      subroutine step2cal30(kstep,ktspd,kyea,kmon,kday,khou,kmin)
      use calmod
      implicit none
      integer, intent(IN ) :: kstep ! time step since simulation start
      integer, intent(IN ) :: ktspd ! time steps per day
      integer, intent(OUT) :: kyea  ! current year   of simulation
      integer, intent(OUT) :: kmon  ! current month  of simulation
      integer, intent(OUT) :: kday  ! current day    of simulation
      integer, intent(OUT) :: khou  ! current hour   of simulation
      integer, intent(OUT) :: kmin  ! current minute of simulation

      integer :: idall
      integer :: istp

      idall = kstep / ktspd
      kyea  = idall / 360
      idall = mod(idall,360)
      kmon  = idall / 30 + 1
      kday  = mod(idall, 30) + 1
      istp  = mod(kstep,ktspd)
      kmin  = (istp * 1440) / ktspd
      khou  = kmin / 60
      kmin  = mod(kmin,60)

      return
      end subroutine step2cal30

!     =====================
!     SUBROUTINE CAL2STEP30
!     =====================

      subroutine cal2step30(kstep,ktspd,kyea,kmon,kday,khou,kmin)
      use calmod
      implicit none
      integer, intent(OUT) :: kstep ! time step since simulation start
      integer, intent(IN ) :: ktspd ! time steps per day
      integer, intent(IN ) :: kyea  ! current year   of simulation
      integer, intent(IN ) :: kmon  ! current month  of simulation
      integer, intent(IN ) :: kday  ! current day    of simulation
      integer, intent(IN ) :: khou  ! current hour   of simulation
      integer, intent(IN ) :: kmin  ! current minute of simulation

      kstep = ktspd * (kyea * 360 + (kmon-1) * 30 + kday - 1) &
            +(ktspd * (khou *  60 + kmin)) / 1440

      return
      end subroutine cal2step30

      program caltest
      use calmod
      parameter (nstep=100)
      do istep = 1 , nstep
         call step2cal(istep,32,iyy,imm,idd,ihh,imi)
         call cal2step(jstep,32,iyy,imm,idd,ihh,imi)
         if (imm==2 .and. idd > 28 .and. ihh == 0 .and. imi == 0) &
         write (*,1000) istep,jstep,istep/32+1,iyy,imm,idd,ihh,imi
!        if (iyy == 3 .and. imm == 12) &
!        write (*,1000) istep,jstep,istep/32+1,iyy,imm,idd,ihh,imi
!        if (istep > 1215470) &
!        write (*,1000) istep,jstep,istep/32+1,iyy,imm,idd,ihh,imi
         if (istep /= jstep) then
         write (*,1000) istep,jstep,istep/32+1,iyy,imm,idd,ihh,imi
         stop
         endif
      enddo
 1000 format(2i10,i8,i5,i3,i3,i3,i3)
      print *,'monaccu'
      do j = 0,12
         print '(i2,i4)',j,monaccu(j)
      enddo
      print *,'monaccu'
      print *,'enter date1 (YYYY MM DD)'
      read (*,'(i4,i3,i3)') iy1,im1,id1
      print *,'enter date2 (YYYY MM DD)'
      read (*,'(i4,i3,i3)') iy2,im2,id2

      call cal2step(istep1,32,iy1,im1,id1,0,0)
      call cal2step(istep2,32,iy2,im2,id2,0,0)
      idd = (istep2-istep1)/32
      idy = idd / 365
      print *,idd,' Days  or ',idy,' Years'
      print *,'Weekday is ',nweekday(istep1/32)
      print *,'Weekday is ',nweekday(istep2/32)
      stop
      end
