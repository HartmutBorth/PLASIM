!     ==================================================================
!     ------------------------------------------------------------------
!
      subroutine diva
      use lsgvar
      implicit none
!
!     ------------------------------------------------------------------
!
!**** *diva*.
!
!     by E. Maier-Reimer.
!     Last modified by U. Mikolajewicz 9/87.
!
!     Purpose.
!     --------
!     *diva* controls the computation of diagnostic variables
!     and computes total velocities (*utot* and *vtot*).
!     Dependent on *nsve* different subroutines are called.
!
!
!     *nsve*=2  (exact solution for barotropic velocities).
!     --------
!      *diva* calls  *press*     computes normalized pressure *p*.
!                    *uvtrop2*   computes barotropic velocities *ub*,
!                                *vb* and surface elevation *zeta*.
!                    *geost2*    computes baroclinic velocities.
!                    *cont1*     computes vertical velocities *w*.
!
!**   Input.
!     ------
!     rhdif     pot density-difference in the vertical direction.
!     taux      ) windstress divided by density.
!     tauy      )
!     old velocities and surface elevation.
!     all via common /lsgfie/.
!
!     Output.
!     -------
!     utot      ) total velocities.
!     vtot      )
!     w         vertical velocities.
!     ub        ) barotropic velocities.
!     vb        )
!     zeta      surface elevation.
!     zetado    time derivative of *zeta*.
!     psi       barotropic stream function (pure diagnostic).
!
!    -------------------------------------------------------------------
!
!     Declaration of local constants/variables.
!     -----------------------------------------
!
      integer :: i,j,k
!
!     fields for predictor-corrector scheme
!
      real (kind=8) :: ubold(ien,jen)
      real (kind=8) :: vbold(ien,jen)
      real (kind=8) :: zold(ien,jen)
!
!
!*    1.        Set initial constants and variables.
!     ------------------------------------
!
      zold(:,:) = zeta(:,:)
!
!*    2.        Computation of normalized pressure.
!     -----------------------------------
!
      call press
!
!*    3.        Computation of velocities.
!     --------------------------
!
!*    3.2.      *nsve*=2
!     --------
!
!     Compute barotropic velocities.
!
      call uvtrop2
!
!     Compute baroclinic velocities.
!
      call geost2
!
!     predictor-corrector for sea ice/sea level
!
      zeta(:,:) = zold(:,:)
!
      call uvtrop2
!
      ubold(:,:) = 0.0
      vbold(:,:) = 0.0

      do k=1,ken
        do j=3,jen-2
          do i=1,ien
            ubold(i,j)=ubold(i,j)+utot(i,j,k)*delta(i,j,k)
            vbold(i,j)=vbold(i,j)+vtot(i,j,k)*delta(i,j,k)
          end do
        end do
      end do
      do j=3,jen-2
        do i=1,ien
          if (depth(i,j)<0.5) cycle
          ubold(i,j)=ubold(i,j)/depth(i,j)
          vbold(i,j)=vbold(i,j)/depth(i,j)
        end do
      end do
      do k=1,ken
        do j=3,jen-2
          do i=1,ien
            utot(i,j,k)=(utot(i,j,k)+ub(i,j)-ubold(i,j))*wetvec(i,j,k)
            vtot(i,j,k)=(vtot(i,j,k)+vb(i,j)-vbold(i,j))*wetvec(i,j,k)
          end do
        end do
      end do
!
!
!*    4.        Compute vertical velocities.
!     ----------------------------
!
      call cont1
      return
      end subroutine diva
