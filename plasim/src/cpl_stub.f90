!     ==================
!     subroutine CLSGINI
!     ==================
!
      subroutine clsgini(kdatim,ktspd,kaomod,pslm,kxa,kya)
!
      real :: pslm(kxa,kya)
      integer :: kdatim(7)
      integer :: ktspd
      integer :: kaomod
      integer :: kxa,kya
!
      return
      end subroutine clsgini
!
!     ===================
!     subroutine CLSGSTEP
!     ===================
!
      subroutine clsgstep(kdatim,kstep,psst,ptaux,ptauy,pfresh,pice     &
     &                   ,pheat,pfldo)
!
      real :: psst(*)      ! atm sst (input)
      real :: ptaux(*)     ! atm u-stress (input)
      real :: ptauy(*)     ! atm v-stress (input)
      real :: pfresh(*)    ! atm fresh water flux (input)
      real :: pice(*)      ! atm ice thickness (incl. snow) (input)
      real :: pheat(*)     ! atm heat flux (input; not used yet)
      real :: pfldo(*)     ! atm deep ocan heat flux (output)
      integer :: kdatim(7) ! date and time
      integer :: kstep     ! current atm time step
!
      return
      end subroutine clsgstep
!
!     ===================
!     subroutine CLSGSTOP
!     ===================
!
      subroutine clsgstop
!
      return
      end subroutine clsgstop
