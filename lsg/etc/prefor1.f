!     ==================================================================
!     ------------------------------------------------------------------
!
      subroutine prefor1(b)
      use lsgvar
      implicit none
!
!     ------------------------------------------------------------------
!
!**** *prefor1*.
!
!     by E. Maier-Reimer.
!     Last modified by U. Mikolajewicz 9/87.
!
!**   Purpose.
!     --------
!     *prefor1* computes the right hand side of the equation for the
!     barotropic velocities, when *nsve*=2.
!
!     Input.
!     ------
!     zeta       surface elevation         (common /lsgsur/ ).
!     p          norm pressure.            (common /lsgpre/ ).
!     taux       windstress/density.       (common /lsgfie/ ).
!     tauy                                 (common /lsgfie/ ).
!
!     Output.
!     -------
!     b          array containing the value of the right hand side of
!                the equation for the barotropic velocities.
!                parameter.
!
!     Interface.
!     ----------
!     *call* *prefor1*(*b*)
!                 *b* an array with dimension b(matrx)
!
!     ------------------------------------------------------------------
!
!     Parameter declaration.
!     ----------------------
!
      real (kind=8) :: b(matrx)
!
!     ------------------------------------------------------------------
!
!
!     Declaration of local variables.
!     -------------------------------
!
      integer :: i,j,k,l,m1

      real (kind=8) :: b1(ien,jen,2),b2(matrx),zetar(ien,jen)
      real (kind=8) :: pr(ien,jen,ken)
      real (kind=8) :: pra(ien*jen*ken),zetara(ien*jen)
      real (kind=8) :: b1a(ien*jen*2)
      real (kind=8) :: b1b(ien*jen,2)
      real (kind=8) :: prb(ien*jen,ken)

      equivalence (pra,pr)
      equivalence (zetar,zetara)
      equivalence (b1,b1a)
      equivalence (b1,b1b)
      equivalence (prb,pr)
!
!*    1.         Set initial values.
!     -------------------
!
!     *zetar*  next surface elevation to the right.
!
      l=ien*jen-1
      do i=1,l
        zetara(i)=zetaa(i+1)
      end do
      do j=1,jen
        zetar(ien,j)=zeta(1,j)
      end do
!
!     *pr* norm pressure the right.
!
      l=ien*jen*ken-1
      do i=2,l+1
        pra(i-1)=pb(i)
      end do
      do k=1,ken
        do j=1,jen
          pr(ien,j,k)=p(1,j,k)
        end do
      end do
!
!*    2.        Computation of the terms.
!     -------------------------
      l=ien*jen*2
      do i=1,l
        b1a(i)=0.
      end do
!
!     Computation of the terms containing windstress and surfaceslope.
!
      l=ien*(jen-2)
      do i=ien+1,ien+l
        b1b(i,1)=tauxa(i)+g*dliha(i)*(zetaa(i)-zetara(i))*deptha(i)
        b1b(i,2)=tauya(i)+g*dpin*(zetara(i+ien)-zetaa(i-ien))*deptha(i)
      end do
!
!     Computation of the horizontal pressure differences in the entire
!     water mass.
!
      do k=1,ken
        do i=ien+1,ien+l
          b1b(i,1)=b1b(i,1)-dliha(i)*(prb(i,k)-pa(i,k))*deltaa(i,k)
          b1b(i,2)=b1b(i,2)-dpin*(pa(i-ien,k)-prb(i+ien,k))*deltaa(i,k)
        end do
      end do
!
!     Scaling with the depth.
!
      do i=ien+1,ien+l
        if (deptha(i)>0) then
          b1b(i,1)=b1b(i,1)/deptha(i)
          b1b(i,2)=b1b(i,2)/deptha(i)
        end if
      end do
!
!     3.         Construction of the output array *b*.
!     -------------------------------------
!
      l=ien*jen
      m1=matrx/2
      j=-1
      do i=1,l
        if (numxa(i)>0.) then
          b(numxa(i)*2-1)=b1b(i,1)
          b(numxa(i)*2)=0.
        end if
      end do
      j=-1
      do i=1,l
        if (numxa(i)>0.) then
          b2(2*numxa(i)-1)=0.
          b2(2*numxa(i))=b1b(i,2)
        end if
      end do
      do i=1,matrx
        b(i)=b(i)+b2(i)
      end do
!      do k=1,jen
!         write (80,'("numx(:,",i4,")=",i5)') k
!         write (80,'(16I5)') numx(:,k)
!      enddo
!      write (81,'(8e10.3)') b

!
      end subroutine prefor1
