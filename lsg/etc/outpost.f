!     ==================================================================
!     ------------------------------------------------------------------
!
      subroutine outpost
      use lsgvar
      implicit none
!
!     ------------------------------------------------------------------
!
!**** *outpost*
!
!     by U. Mikolajewicz 12/87.
!
!     Purpose.
!     --------
!     *outpost writes output data on tapeunit *nopost*.
!     This file is further processed in the postprocessor.
!     The file is defined during the run.
!     The time for the next output is read from file and the
!     name of the next output file is generated.
!
!     Input.
!     ------
!     common blocks /lsgfie/ and /lsgsur/
!     filpnew   file for output to be written on.
!
!     Output.
!     -------
!     File with output data.
!     filpnew   new name for the next file with output data.
!     ntnout    number of time step to write next output file.
!
!     Interface.
!     ----------
!     *call* *outpost*
!
!     Calls subroutines *datfrnt*.
!     ------------------------------------------------------------------
!
!
!     Dimension of local variables.
!     -----------------------------
      integer :: ji,jr,k,iyy,idd,jlev,jl,i,j
      real (kind=8) :: zero,zdum,zlen,zprel,zcode,zlev
      integer :: iddr(512)
      real (kind=8) :: zddr(512)
!
!*    1.        Actualize *iddr*.
!     -----------------
      zero=0.
      zdum=zero
!
!     *iddr* and *zddr* are header fields for output file.
!
      do ji=1,512
        iddr(ji)=nddr(ji)
        zddr(ji)=oddr(ji)
      end do
      iddr(404)=0
      iddr(22)=6*ken+8
      iddr(16)=6*ken+8
      call datfrnt(nt,iddr(12),iddr(11))
!
!     Pot. temperature in Kelvin.
!
      do jr=1,ken
        iddr(22+2*jr-1)=-2
      end do
      do jr=1,ken
!
!       Salinity.
!
        iddr(22+2*ken+2*jr-1)=-5
        iddr(22+2*ken+2*jr)=nint(du(jr))
!
!       Velocities.
!
        iddr(22+4*ken+2*jr-1)=-3
        iddr(22+4*ken+2*jr)=nint(du(jr))
        iddr(22+6*ken+2*jr-1)=-4
        iddr(22+6*ken+2*jr)=nint(du(jr))
        iddr(22+8*ken+2*jr-1)=-7
        iddr(22+8*ken+2*jr)=nint(dw(jr))
      end do
!
!     Barotropic velocities.
!
      iddr(22+10*ken+1)=-37
      iddr(22+10*ken+2)=-100
      iddr(22+10*ken+3)=-38
      iddr(22+10*ken+4)=-100
!
!     Surface elevation.
!
      iddr(22+10*ken+5)=-1
      iddr(22+10*ken+6)=-100
!
!     Ice thickness.
!
      iddr(22+10*ken+7)=-13
      iddr(22+10*ken+8)=-100
!
!     Topography in vector-points.
!
      iddr(22+10*ken+9)=-99
      iddr(22+10*ken+10)=-100
!
!     Topography in scalar-points.
!
      iddr(22+10*ken+11)=-98
      iddr(22+10*ken+12)=-100
!
!     Fresh water fluxes due to Newtonian coupling.
!
      iddr(22+10*ken+13)=-67
      iddr(22+10*ken+14)=-100
!
!     Heat fluxes due to Newtonian coupling.
!
      iddr(22+10*ken+15)=-68
      iddr(22+10*ken+16)=-100
!
!     Convective adjustment.
!
      iddr(22+10*ken+17)=-69
      iddr(22+10*ken+18)=-100
!
!     Horizontal stream function
!
      iddr(22+10*ken+19)=-27
      iddr(22+10*ken+20)=-100
!
!     Freshwater flux
!
      iddr(22+10*ken+21)=-65
      iddr(22+10*ken+22)=-100
!
      do k=2,ken
        if (iddr(22+10*ken+22+(k-1)*2)>=300) cycle
        iddr(22+10*ken+21+(k-1)*2)=-69
        iddr(22+10*ken+22+(k-1)*2)=nint(dw(k-1))
      end do
!
!
!*    2.        Write output on file *filpnew*.
!     -------------------------------
!
!     Write *iddr* and *zddr* on file.
!
      iyy=iddr(12)
      idd=iddr(11)

!
!  Coupled with PlaSim: name without year to be edited by coupled script
      write (filpnew,'(a)') 'LSG_outpost'
 9885 format (' OUTPOST: standard output ',a,' on file ',a,             &
     &        ' YYYYY -MMDD :',i5.2,i6.4,'  +')
!#ifdef 1
!      write (no6,9885) '(real*4)',filpnew,iyy,idd,ksum
!#else
      write (no6,9885) ' ',filpnew,iyy,idd
!#endif

      open (nopost,file=filpnew,form="unformatted",position="append")




      write (nopost) iddr
      write (nopost) zddr
      zlen=ien*jen
      zprel=6.
!
!     Write temperature.
!
      zcode=-2
      do jlev=1,ken
        write (nopost) zprel,zlen,zcode,zddr(10+jlev-1),zdum,zdum
        write (nopost) (((t(jl,jr,jlev)+tkelvin),jl=1,ien),jr=1,jen)
      end do
!
!     Salinity.
!
      zcode=-5
      do jlev=1,ken
        write (nopost) zprel,zlen,zcode,zddr(10+jlev-1),zdum,zdum
        write (nopost) ((s(jl,jr,jlev),jl=1,ien),jr=1,jen)
      end do
!
!     Velocities.
!
      zcode=-3
      do jlev=1,ken
        write (nopost) zprel,zlen,zcode,zddr(10+jlev-1),zdum,zdum
        write (nopost) ((utot(jl,jr,jlev),jl=1,ien),jr=1,jen)
      end do
!
      zcode=-4
      do jlev=1,ken
        write (nopost) zprel,zlen,zcode,zddr(10+jlev-1),zdum,zdum
        write (nopost) ((vtot(jl,jr,jlev),jl=1,ien),jr=1,jen)
      end do
!
      zcode=-7
      do jlev=1,ken
        write (nopost) zprel,zlen,zcode,dw(jlev),zdum,zdum
        write (nopost) ((w(jl,jr,jlev),jl=1,ien),jr=1,jen)
      end do
!
!     Barotropic velocities.
!
      zlev=-100.
      zcode=-37.
      write (nopost) zprel,zlen,zcode,zlev,zdum,zdum
      write (nopost) ((ub(jl,jr),jl=1,ien),jr=1,jen)
!
      zcode=-38.
      write (nopost) zprel,zlen,zcode,zlev,zdum,zdum
      write (nopost) ((vb(jl,jr),jl=1,ien),jr=1,jen)
!
!     Surface elevation and ice thickness.
!
      zcode=-1.
      write (nopost) zprel,zlen,zcode,zlev,zdum,zdum
      write (nopost) ((zeta(jl,jr),jl=1,ien),jr=1,jen)
!
      zcode=-13.
      write (nopost) zprel,zlen,zcode,zlev,zdum,zdum
      write (nopost) ((sice(jl,jr),jl=1,ien),jr=1,jen)
!
!     Horizontal barotropic stream function
      zcode=-27.
      write (nopost) zprel,zlen,zcode,zlev,zdum,zdum
      write (nopost) ((psi(jl,jr),jl=1,ien),jr=1,jen)
!
!     Store topography in vector-points.
!
      zcode=-99.
      write (nopost) zprel,zlen,zcode,zlev,zdum,zdum
      write (nopost) ((depth(jl,jr),jl=1,ien),jr=1,jen)
!
!     Store topography in scalar-points.
!
      zcode=-98.
      write (nopost) zprel,zlen,zcode,zlev,zdum,zdum
      write (nopost) ((depp(jl,jr),jl=1,ien),jr=1,jen)
!
!     Store *fluwat*
!
      zcode=-65.
      write (nopost) zprel,zlen,zcode,zlev,zdum,zdum
      write (nopost) fluwat
!
!     Store *flukwat*.
!
      zcode=-67.
      write (nopost) zprel,zlen,zcode,zlev,zdum,zdum
      write (nopost) flukwat
!
!     Store *flukhea*.
!
      zcode=-68.
      write (nopost) zprel,zlen,zcode,zlev,zdum,zdum
      write (nopost) flukhea
!
!     Store convect ice adjustments.
!
      zcode=-66.
      zlev=-100.
      write (nopost) zprel,zlen,zcode,zlev,zdum,zdum
      write (nopost) ((convad(jl,jr,1),jl=1,ien),jr=1,jen)
      zcode=-69.
      do k=2,ken
        zlev=dw(k-1)
        write (nopost) zprel,zlen,zcode,zlev,zdum,zdum
        write (nopost) ((convad(jl,jr,k),jl=1,ien),jr=1,jen)
      end do
!
!     Store *taux*.
!
      zcode=52.
      zlev=-100.
      write (nopost) zprel,zlen,zcode,zlev,zdum,zdum
      write (nopost) ((taux(i,j)*rhonul,i=1,ien),j=1,jen)
!
!     Store *tauy*.
!
      zcode=53.
      zlev=-100.
      write (nopost) zprel,zlen,zcode,zlev,zdum,zdum
      write (nopost) ((tauy(i,j)*rhonul,i=1,ien),j=1,jen)
!
!     Store *tbound*.
!
      zcode=-92.
      zlev=-100.
      write (nopost) zprel,zlen,zcode,zlev,zdum,zdum
      write (nopost) ((tbound(i,j)+tkelvin,i=1,ien),j=1,jen)
!
!     Store *fluhea*.
!
      zcode=-18.
      write (nopost) zprel,zlen,zcode,zlev,zdum,zdum
      write (nopost) fluhea
!
!*    3.        Store output on file *filpnew*.
!     -------------------------------
      close (nopost)
      return
      end subroutine outpost
