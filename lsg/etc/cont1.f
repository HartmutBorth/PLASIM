      subroutine cont1
      use lsgvar
      implicit none
!
!     ------------------------------------------------------------------
!
!**** *cont1*.
!
!     by E. Maier-Reimer.
!     Last modified by U. Mikolajewicz 9/87.
!
!     Purpose.
!     --------
!     *cont1* computes the vertical velocities *w* from the total
!     velocities *utot* and *vtot* using the equation of continuity.
!
!**   Input.
!     ------
!     utot,vtot   horizontal velocities ( common/lsgfie/).
!
!     Output.
!     -------
!     w           vertical velocities   (common/lsgfie/).
!     trup        sum of vertical transports up (m**3/s).
!     trdo         "  "     "          "     down   "   .
!                      both in common /lsgdia/.
!
!     Interface.
!     ----------
!     *call* *cont1*.
!
!     ------------------------------------------------------------
!
!     Declaration of local variables.
!     -------------------------------
      integer :: i,j,k,il
      real (kind=8) :: hordiv,wups
!
!*    1.        Initialisation.
!     ---------------
      w(:,:,:) = 0.0
!
!*    2.        Computation of the vertical velocity from bottom upward.
!     --------------------------------------------------------
!
!     Determination of the values of the neighbour cell.
!
      do k=ken-1,1,-1
        do j=3,jen-2
          do i=1,ien
            il=i-1
            if (il<1) il=ien
            hordiv=utot(il,j,k+1)*delta(il,j,k+1)*dphi-utot(i,j,k+1)    &
     &            *delta(i,j,k+1)*dphi+vtot(i,j+1,k+1)*delta(i,j+1,k+1) &
     &            *dl(j+1)-vtot(il,j-1,k+1)*delta(il,j-1,k+1)*dl(j-1)
            w(i,j,k)=w(i,j,k+1)+hordiv/(dl(j)*dphi)*wet(i,j,k+1)
          end do
        end do
      end do
!
!*    3.        Computation of control parameters.
!     ----------------------------------
!     The parameters *trup* and *trdo* describe the sum of
!     the transports up and down.
!
      do k=1,ken
        trup(k)=0.
        trdo(k)=0.
        do j=1,jen
          do i=1,ien
            wups=sign(0.5_8,w(i,j,k))
            trup(k)=trup(k)+w(i,j,k)*(wups+abs(wups))*dlh(i,j)
            trdo(k)=trdo(k)+w(i,j,k)*(abs(wups)-wups)*dlh(i,j)
          end do
        end do
        trup(k)=trup(k)*dphi
        trdo(k)=trdo(k)*dphi
      end do
      return
      end subroutine cont1
