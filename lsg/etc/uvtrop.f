!     ==================================================================
!     ------------------------------------------------------------------
!
      subroutine uvtrop(b1)
      use lsgvar
      implicit none
!
!     ------------------------------------------------------------------
!
!**** *uvtrop*.
!
!     by E. Maier-Reimer.
!     last modified by U. Mikolajewicz 9/87.
!     optimized by E. Kirk 2/2008
!
!     Purpose.
!     --------
!     *uvtrop* computes barotropic velocities and stream function
!     according to the equation for the barotropic velocities.
!     *uvtrop* computes the surface elevation *zeta* and its
!     time derivative *zetado* from the depth-integrated equation
!     of continuity.
!
!**   Input.
!     ------
!     b1       right hand side of the equation (from *prefor1*).
!     skal     scaling factors.           )
!     elim     elimination factors.       )   via common /lsgmat/.
!     trisys   tridiagonalsystem.         )
!     zeta     old surface elevation.       (common /lsgsur/).
!
!     Output.
!     -------
!     ub       ) barotropic velocities  (common /lsgfie/).
!     vb       )
!     psi      barotropic stream function.
!     zeta     surface elevation.       (common /lsgsur/).
!     zetado   time derivative of *zeta*.(common /lsgsur/).
!
!     Interface.
!     ----------
!     *call* *uvtrop(b1)*.
!             parameter b1, dimension b1(matrx).
!
!
!     Parameter declaration.
!     ----------------------
      real (kind=8) :: b1(matrx)
!
!
!     Declaration of local variables.
!     -------------------------------
      integer,save :: ncall = 0
      integer :: i,j,k,l,jem2,itgar,itgrr,matrx1

      real (kind=8) :: xr,psi1,dti,fact,aree,areo,zetame,zetamo
      real (kind=8) :: apafac
      real (kind=8) :: b(matot),x(matot)
      real (kind=8) :: ubla(ienjen),vbla(ienjen)
      real (kind=8) :: zetas(ienjen)
!
!
!
!*    1.        Set initial values and constants.
!     ---------------------------------
!
      jem2=jen-2
      itgar=33
      matrx1=matrx-1
      itgrr=31
      mindi=0
      mindj=0
!
      x(:) = 0.0
!
!
!     *   2.        Comput. of the right side of the triangular matrix.
!     ---------------------------------------------------
!
!     Scaling.
!
      b(1:matrx)  = b1(:) * skal(:) * dt
      b(matrx+1:) = 0.0
!
!     Elimination.
!
      do j=1,matrx1
        b(j+1:j+kb) = b(j+1:j+kb) - elim(:,j) * b(j)
      end do
!
!
!*    3.        Computation of the barotropic velocities.
!     -----------------------------------------
!
!     Solving the matrix equation.
!
      x(matrx)=b(matrx)/trisys(1,matrx)
      do l=matrx1,1,-1
        xr = dot_product(x(l+1:l+kb),trisys(2:kb+1,l))
        x(l)=(b(l)-xr)/trisys(1,l)
      end do
!
!
!     Computation of *ub* and *vb*.
!
      where (numxa(:) > 0)
         uba(:) = x(2*numxa(:)-1)
         vba(:) = x(2*numxa(:)  )
      else where
         uba(:) = 0.0
         vba(:) = 0.0
      end where
!
!
!*    4.        Computation of the barotropic stream function.
!     ----------------------------------------------
!     with zero at Antarctica.
      psi(:,jem2:jen) = 0.0

      do j=jem2,1,-1
        psi(ien,j)=psi(1,j+2)+depth(ien,j+1)*dphi*ub(ien,j+1)*1.e-6
        do i=1,ien-1
          psi(i,j)=psi(i+1,j+2)+depth(i,j+1)*dphi*ub(i,j+1)*1.e-6
        end do
      end do
!
      psimax=0.0
      do j=3,jen-2
        do i=1,ien
          psi1=abs(psi(i,j))
          if (psi1<psimax) cycle
          psimax=psi1
          mindi=i
          mindj=j
        end do
      end do
!
!
!*    5.        New surface elevation.
!     ----------------------
!
!     Defines arrays containing velocities of next grid points.
!
      do i=1,ienjen,ien
         ubla(i) = uba(i+ien-1)
         vbla(i) = vba(i+ien-1)
         ubla(i+1:i+ien-1) = uba(i:i+ien-2)
         vbla(i+1:i+ien-1) = vba(i:i+ien-2)
      end do
!
!     Computes surface elevation.
!
      zetas(:) = zetaa(:) ! Save zetaa for later use

      l=ien*(jen-2)
      do i=ien+1,ien+l
        zetaa(i)=zetaa(i)+dt*dliha(i)                                   &
     &           *(ubla(i)*depthla(i)-uba(i)*deptha(i)                  &
     &           +dpin*(vba(i+ien)*deptha(i+ien)*dlha(i+ien)-vbla(i-ien)&
     &           *depthla(i-ien)*dlha(i-ien)))
      end do
!
!     Surface elevation at the poles.
!
      do i=1,ien
        zeta(i,1)=0.
        zeta(i,2)=0.
        zeta(i,jen)=0.
        zeta(i,jen-1)=0.
      end do
!
      dti=1./dt
!
!     Computes time derivative of surface elevation.
!
      do i=ien+1,ien+l
        fact=1.
        zetadoa(i)=dti*(zetaa(i)-zetas(i))*fact
        zetadoa(i)=dmax1(zetadoa(i),-dw(1)/4.)
        zetadoa(i)=dmin1(zetadoa(i),dw(1)/4.)
      end do
!
!     Compensation for cutoff errors of zeta.
!
      aree=0.
      areo=0.
      zetame=0.
      zetamo=0.
      do j=2,jen-1,1
        do i=1,ien
          zetamo=zetamo+wet(i,j,1)*dlh(i,j)*dphi*zeta(i,j)
          areo=areo+wet(i,j,1)*dlh(i,j)*dphi
          zetame=zetame+wet(i,j+1,1)*dlh(i,j+1)*dphi*zeta(i,j+1)
          aree=aree+wet(i,j+1,1)*dlh(i,j+1)*dphi
        end do
      end do
      apafac=1./120.
      zetame=apafac*zetame/aree
      zetamo=apafac*zetamo/areo
      print *,"ZETAM0",zetamo,zetame
      do j=3,jen-2,2
        do i=1,ien
          zeta(i,j)=wet(i,j,1)*(zeta(i,j)-zetamo)
          zeta(i,j+1)=wet(i,j+1,1)*(zeta(i,j+1)-zetame)
        end do
      end do
      !if (ncall == 0) then
      write (87,'(4e20.9)') ub
      write (87,'(4e20.9)') vb
      write (87,'(4e20.9)') psi
      write (87,'(4e20.9)') zeta
      write (87,'(4e20.9)') zetado
      !ncall = 1
      !endif
  
      return
      end subroutine uvtrop
