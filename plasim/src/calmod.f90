!     =============
!     MODULE CALMOD
!     =============

      module calmod
      implicit none
      integer :: mondays(0:12) = (/0,31,28,31,30,31,30,31,31,30,31,30,31/)
      integer :: mona365(0:12) = (/0,31,59,90,120,151,181,212,243,273,304,334,365/)
      integer :: monaccu(0:12)
      integer :: ny400d = 400 * 365 + 97
      integer :: ny100d = 100 * 365 + 24
      integer :: ny004d =   4 * 365 +  1
      integer :: ny001d =       365
      integer :: nud        =  6

!     These values are copied from pumamod in subroutine calini

      integer :: n_days_per_month =  30
      integer :: n_days_per_year  = 360
      integer :: n_start_step     =   0
      integer :: ntspd            =   0
      real    :: solar_day        = 86400.0 ! [sec]

      end module calmod


!     =================
!     SUBROUTINE CALINI
!     =================

      subroutine calini(k_days_per_month,k_days_per_year,k_start_step &
                       ,ktspd,psolday,kpid)
      use calmod

      n_days_per_month = k_days_per_month
      n_days_per_year  = k_days_per_year
      n_start_step     = k_start_step
      ntspd            = ktspd
      solar_day        = psolday

      if (kpid == 0) then
         write(nud,1050)
         write(nud,1060)
         write(nud,1050)
         write(nud,1000) n_days_per_year
         write(nud,1010) n_days_per_month
         write(nud,1030) n_start_step
         write(nud,1040) ntspd
         write(nud,1050)
         endif
      return
 1060 format(" * Calmod initialization   *")
 1000 format(" * Days per year:    ",i5," *")
 1010 format(" * Days per month:   ",i5," *")
 1030 format(" * Start step:",i12," *")
 1040 format(" * Timesteps per day:",i5," *")
 1050 format(" ***************************")
      end subroutine calini

!     ====================
!     SUBROUTINE YDAY2MMDD
!     ====================

      subroutine yday2mmdd(kyday,kmon,kday)
      use calmod
      if (n_days_per_year == 365) then
         kmon = 1
         do while (kyday > mona365(kmon))
            kmon = kmon + 1
         enddo
         kday = kyday - monaccu(kmon-1)
      else
         kmon = (kyday-1) / n_days_per_month
         kday = kyday - n_days_per_month * kmon
      endif
      return
      end

!     =================
!     FUNCTION NWEEKDAY
!     =================

      integer function nweekday(kday)
      nweekday = mod(kday+5,7)
      return
      end


!     ===================
!     FUNCTION NDAYOFYEAR
!     ===================

      integer function ndayofyear(kstep)
      use calmod
      integer :: idatim(7)

      if (n_days_per_year == 365) then
         call step2cal(kstep,ntspd,idatim)
         ndayofyear = idatim(3) + monaccu(idatim(2)-1)
      else
         call step2cal30(kstep,ntspd,idatim)
         ndayofyear = idatim(3) + n_days_per_month * (idatim(2)-1)
      endif

      return
      end


!     ===================
!     SUBROUTINE STEP2CAL
!     ===================

      subroutine step2cal(kstep,ktspd,kdatim)
      use calmod
      implicit none
      integer, intent(IN ) :: kstep     ! time step since simulation start
      integer, intent(IN ) :: ktspd     ! time steps per day
      integer, intent(OUT) :: kdatim(7) ! year,month,day,hour,min,weekday,leapyear

      integer :: iyea  ! current year   of simulation
      integer :: imon  ! current month  of simulation
      integer :: iday  ! current day    of simulation
      integer :: ihou  ! current hour   of simulation
      integer :: imin  ! current minute of simulation
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

      iyea  = iy400 * 400 + iy100 * 100 + iy004 * 4 + iy001

      monaccu(0) = mondays(0)
      monaccu(1) = mondays(1)
      monaccu(2) = mondays(1) + mondays(2)
      if (leap) monaccu(2) = monaccu(2) + 1
      do jmon = 3 , 12
         monaccu(jmon) = monaccu(jmon-1) + mondays(jmon)
      enddo
      imon = 1
      id001 = id001 + 1
      do while (id001 > monaccu(imon))
         imon = imon + 1
      enddo
      iday = id001 - monaccu(imon-1)

      istp = mod(kstep,ktspd)
      imin = (istp * 1440) / ktspd
      ihou = imin / 60
      imin = mod(imin,60)

      kdatim(1) = iyea
      kdatim(2) = imon
      kdatim(3) = iday
      kdatim(4) = ihou
      kdatim(5) = imin
      kdatim(6) = mod(kstep/ktspd+5,7) ! day of week
      if (leap) then
         kdatim(7) = 1
      else
         kdatim(7) = 0
      endif

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

      if (n_days_per_year /= 365) then ! simplified calendar
         kstep = ktspd * (kyea    * n_days_per_year  &
                       + (kmon-1) * n_days_per_month &
                       +  kday-1)
      return
      endif
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

      subroutine step2cal30(kstep,ktspd,kdatim)
      use calmod
      implicit none
      integer, intent(IN ) :: kstep     ! time step since simulation start
      integer, intent(IN ) :: ktspd     ! time steps per day
      integer, intent(OUT) :: kdatim(7) ! year,month,day,hour,min,weekday,leapyear

      integer :: iyea  ! current year   of simulation
      integer :: imon  ! current month  of simulation
      integer :: iday  ! current day    of simulation
      integer :: ihou  ! current hour   of simulation
      integer :: imin  ! current minute of simulation
      integer :: idall
      integer :: istp

      idall = kstep / ktspd
      iyea  = idall / n_days_per_year
      idall = mod(idall,n_days_per_year)
      imon  = idall / n_days_per_month + 1
      iday  = mod(idall,n_days_per_month) + 1
      istp  = mod(kstep,ktspd)
      imin  = (istp * solar_day) / (ktspd * 60)
      ihou  = imin / 60
      imin  = mod(imin,60)

      kdatim(1) = iyea
      kdatim(2) = imon
      kdatim(3) = iday
      kdatim(4) = ihou
      kdatim(5) = imin
      kdatim(6) = 0 ! day of week
      kdatim(7) = 0 ! leap year

      return
      end subroutine step2cal30

!     =================
!     SUBROUTINE NTOMIN
!     =================

       subroutine ntomin(kstep,kmin,khou,kday,kmon,kyea)
       use calmod
       implicit none

       integer, intent(in ) :: kstep
       integer, intent(out) :: kmin,khou,kday,kmon,kyea
       integer :: idatim(7)

       if (n_days_per_year == 365) then
          call step2cal(kstep,ntspd,idatim)
       else
          call step2cal30(kstep,ntspd,idatim)
       endif
       kyea = idatim(1)
       kmon = idatim(2)
       kday = idatim(3)
       khou = idatim(4)
       kmin = idatim(5)
       return
       end

!     =================
!     SUBROUTINE NTODAT
!     =================

      subroutine ntodat(istep,datch)
      character(len=18) datch
      character(len=3) mona(12)
      data mona /'Jan','Feb','Mar','Apr','May','Jun',                   &
     &           'Jul','Aug','Sep','Oct','Nov','Dec'/
      call ntomin(istep,imin,ihou,iday,imon,iyea)
      write (datch,20030) iday,mona(imon),iyea,ihou,imin
20030 format(i2,'-',a3,'-',i4.4,2x,i2,':',i2.2)
      end


!     =================
!     SUBROUTINE DTODAT
!     =================

      subroutine dtodat(id,datch)
      integer :: id(6)
      character(len=18) datch
      character(len=3) mona(12)
      data mona /'Jan','Feb','Mar','Apr','May','Jun',                   &
     &           'Jul','Aug','Sep','Oct','Nov','Dec'/
      write (datch,20030) id(3),mona(id(2)),id(1),id(4),id(5)
20030 format(i2,'-',a3,'-',i4.4,2x,i2,':',i2.2)
      end


!     =================
!     SUBROUTINE MOMINT
!     =================

!     Compute month indices and weights for time interpolation from
!     monthly mean data to current timestep

      subroutine momint(kperp,kstep,kmona,kmonb,pweight)
      use calmod
      implicit none
      integer, intent(in ) :: kperp   ! perpetual mode ?
      integer, intent(in ) :: kstep   ! target step
      integer, intent(out) :: kmona   ! current month (1-12)
      integer, intent(out) :: kmonb   ! next or previous month (0-13)
      real   , intent(out) :: pweight ! interpolation weight

      integer :: idatim(7) ! date time array
      integer :: iday
      integer :: ihour
      integer :: imin
      integer :: jmonb  ! next or previous month (1-12)

      real    :: zday   ! fractional day (including hour & minute)
      real    :: zdpma  ! days per month a
      real    :: zdpmb  ! days per month b
      real    :: zmeda  ! median day of the month a
      real    :: zmedb  ! median day of the month b

!     convert time step to date / time

      idatim(:) = 0
      if (kperp > 0) then                     ! perpetual date
         call yday2mmdd(kperp,idatim(2),idatim(3))
      else if (n_days_per_year == 365) then   ! real calendar
         call step2cal(kstep,ntspd,idatim)
      else                                    ! simple calendar
         call step2cal30(kstep,ntspd,idatim)
      endif

      kmona = idatim(2)
      iday  = idatim(3)
      ihour = idatim(4)
      imin  = idatim(5)

!     set fractional day

      zday = iday + ((ihour * 60.0 + imin) * 60.0) / solar_day

!     compute median of month a

      zdpma = n_days_per_month
      if (n_days_per_year == 365) then
         zdpma = mondays(kmona)
         if (kmona == 2) zdpma = zdpma + idatim(7) ! leap year
      endif
      zmeda = 0.5 * (zdpma + 1.0) ! median day a

!     define neighbour month

      if (zday > zmeda) then
         kmonb = kmona + 1 !     next month (maybe 13)
      else
         kmonb = kmona - 1 ! previous month (maybe  0)
      endif

!     compute median of month b

      zdpmb = n_days_per_month
      if (n_days_per_year == 365) then
         jmonb = mod(kmonb+11,12) + 1 ! convert month (0-13) -> (1-12)
         zdpmb = mondays(jmonb)
         if (jmonb == 2) zdpmb = zdpmb + idatim(7) ! leap year
      endif
      zmedb = 0.5 * (zdpmb + 1.0) ! median day b

!     compute weight

      pweight = abs(zday - zmeda) / (zmeda + zmedb - 1.0)

      return
      end subroutine momint

