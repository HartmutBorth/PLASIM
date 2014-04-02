!     ==================================================================
!     ------------------------------------------------------------------
!
      subroutine datfrnt(ntt,kyear,kdate)
      use lsgvar
      implicit none
!
!     ------------------------------------------------------------------
!
!**** *datfrnt*
!
!     by E. Kirk 22-Feb-2008
!
!     Purpose.
!     --------
!     *datfrnt* transfers date and time information copied from
!     Planet Simulator to LSG
!
!     Input
!     -----
      integer, intent(in) :: ntt ! LSG timestep (ignored)
!
!     mdatim(6) = Year,Month,Day,Hour,Minute,Weekday
!
!     Output.
!     -------
      integer, intent(out) :: kyear ! Simulation year
      integer, intent(out) :: kdate ! MMDD
!
!*    1.        Define output values
!     ------------------------------
      kyear = mdatim(1)
      kdate = -(100 * mdatim(2) + mdatim(3))
      return
      end subroutine datfrnt
