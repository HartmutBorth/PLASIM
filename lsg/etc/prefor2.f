      subroutine prefor2(b)
      use lsgvar
      implicit none
!
!     ------------------------------------------------------------------
!
!**** *prefor2*.
!
!     by E. Maier-Reimer.
!     optimsed by E. Kirk (19-Jan-2009)
!
!**   Purpose.
!     --------
!     *prefor2* computes the right hand side of the equation for the
!     barotropic velocities, when *nsve*=2.
!
!     Input.
!     ------
!     zeta       surface elevation
!     p          norm pressure
!     taux       windstress/density
!     tauy       windstress/density
!
!     Output.
!     -------
!     b          array containing the value of the right hand side of
!                the equation for the barotropic velocities.
!
!     Interface.
!     ----------
!     *call* *prefor2*(*b*)
!                 *b* an array with dimension b(matrx)
!
!     ------------------------------------------------------------------
!
!     Parameter declaration.
!     ----------------------
!
      real (kind=8) :: b(2,matrx/2)
!
!     ------------------------------------------------------------------
!
!
!     Declaration of local variables.
!     -------------------------------
!
      integer :: i,j,k,l
      real (kind=8) :: zx,zy,zdli
!
!*    1.        Computation of the terms.
!     -------------------------
!
!     Computation of the terms containing windstress and surfaceslope
!     and of the horizontal pressure differences in the entire water mass.
!
      do j = 3 , jen-2
        zdli = dli(j)
        do i = 1 , ien
          if (numx(i,j) > 0) then
            l = mod(i,ien) + 1
            zx=taux(i,j)+g*zdli*depth(i,j)*(zeta(i,j  )-zeta(l,j  ))
            zy=tauy(i,j)+g*dpin*depth(i,j)*(zeta(l,j+1)-zeta(i,j-1))
            do k = 1 , ken
              zx = zx - zdli * (p(l,j  ,k)-p(i,j  ,k)) * delta(i,j,k)
              zy = zy - dpin * (p(i,j-1,k)-p(l,j+1,k)) * delta(i,j,k)
            enddo
            b(1,numx(i,j)) = zx / depth(i,j)
            b(2,numx(i,j)) = zy / depth(i,j)
          endif
        enddo
      enddo

      return
      end subroutine prefor2
