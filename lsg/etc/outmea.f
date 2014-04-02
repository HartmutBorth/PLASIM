!     ------------------------------------------------------------------
!
      subroutine outmea(pt,ps,putot,pvtot,pzeta,psice,pflukhea,pflukwat,&
     &                  pub,pvb,pw,pconvad,ptaux,ptauy,ptbound,pfluhea, &
     &                  ppsi,pfluwat,ksum,knt)
      use lsgvar
      implicit none
      integer :: ksum,knt
!     ------------------------------------------------------------------
!
!**** *outmea*
!
!     by U. Mikolajewicz 12/87.
!
!     Purpose.
!     --------
!     *outmea writes output data on tapeunit *nopost*.
!     This file is further processed in the postprocessor.
!     The file is defined during the run.
!     The time for the next output is read from file and the
!     name of the next output file is generated.
!
!     Input.
!     ------
!     common blocks /lsgfie/ and /lsgsur/
!
!     Output.
!     -------
!     file with output data.
!
!     Interface.
!     ----------
!     *call* *outmea*
!
!     Calls subroutines *datfrnt*.
!     ------------------------------------------------------------------
!
!
!     Dimension of parameters.
!     -----------------------
      real (kind=8) :: pt(ien,jen,ken)
      real (kind=8) :: ps(ien,jen,ken)
      real (kind=8) :: putot(ien,jen,ken)
      real (kind=8) :: pvtot(ien,jen,ken)
      real (kind=8) :: pw(ien,jen,ken)
      real (kind=8) :: pub(ien,jen)
      real (kind=8) :: pvb(ien,jen)
      real (kind=8) :: pzeta(ien,jen)
      real (kind=8) :: psice(ien,jen)
      real (kind=8) :: pflukhea(ien,jen)
      real (kind=8) :: pflukwat(ien,jen)
      real (kind=8) :: pconvad(ien,jen,ken)
      real (kind=8) :: ptaux(ien,jen)
      real (kind=8) :: ptauy(ien,jen)
      real (kind=8) :: ptbound(ien,jen)
      real (kind=8) :: pfluhea(ien,jen)
      real (kind=8) :: ppsi(ien,jen)
      real (kind=8) :: pfluwat(ien,jen)


!
!     Declaration of local variables.
!     -----------------------------
      character(len=8) filemea
      integer :: i,j,k,ji,iyy,idd,iyear,khelp,i3dfields,i2dfields
      integer :: jr,jlev
      real (kind=8) :: zero,zdum,zcode
!
!     32bit variables for file <filemea>
!
      integer (kind=4)  :: jddr(512)
      real (kind=4)     :: yddr(512)
      real (kind=4)     :: ylen,yprel,ycode,ylev,ydum
      real (kind=4)     :: yhelp(ien,jen)
!
!     arrays of default kind
!
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
      iddr(11)=-9999
      iddr(12)=0
      iddr(508)=ksum
      call datfrnt(nt,iyy,idd)
      call datfrnt(nt-ksum+1,iddr(510),iddr(509))
      iddr(512)=iyy
      iddr(511)=idd

!     iyear=iyy+1  !  wrong year
      iyear=iyy
      khelp=mod(iyear,100000)

!
!  Coupled with PlaSim: name without year to be edited by coupled script
      write (filemea,'(a)') 'LSG_outm'
!
!  V2.2 2005/09/14
!  number of fields stored in iddr(16):
!  6 fields more than by Uwe M.,
!   where 8 horzontal fields were stored (iddr(16)=6*ken+8):
!   (2 more + 4 forcing fields):
!   code -27: ppsi
!   code -65: pfluwat
!   (evtl. tracers here)
!   code  52: ptaux
!   code  53: ptauy
!   code  92: ptbound (ts or ta, depending on nsmix)
!   code  18: pfluhea (=pflukheat, depending on mix)
!
      i3dfields=6
      i2dfields=8+6
      iddr(16)=i3dfields*ken+i2dfields
      iddr(22)=i3dfields*ken+i2dfields
!
 9880 format (  ' OUTMEA: output average ',a,' on file ',a,             &
     &        ' YYYYY -MMDD :',i5.2,i6.4,'  +'/                         &
     &          ' OUTMEA: output average of',i6,' timesteps')

      write (no6,9880) '(real*4)',filemea,iyy,idd,ksum
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
!     Fresh water fluxes due to newtonian coupling.
!
      iddr(22+10*ken+13)=-67
      iddr(22+10*ken+14)=-100
!
!     Heat fluxes due to newtonian coupling.
!
      iddr(22+10*ken+15)=-68
      iddr(22+10*ken+16)=-100
!
!     Convective adjustment.
!
      iddr(22+10*ken+17)=-69
      iddr(22+10*ken+18)=-100
!
!     Horizontal stream function.
!
      iddr(22+10*ken+19)=-27
      iddr(22+10*ken+20)=-100
!
!     Freshwater flux.
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
!*    2.        Write output on file *filemea*.
!     -------------------------------
!
!     Write *iddr* and *zddr* on file.
!
!     write (no6,*) " open mean output file ",filemea
      open (nopost,file=filemea,access="sequential",form="unformatted")
      rewind nopost

      jddr(:)=iddr(:)
      yddr(:)=zddr(:)
      write (nopost) jddr
      write (nopost) yddr
      ylen=ien*jen
      yprel=6.
!
!     Write temperature.
!
      zcode=-2
      do jlev=1,ken
        ycode=-2.
        write (nopost) yprel,ylen,ycode,yddr(10+jlev-1),ydum,ydum
        yhelp(:,:)=pt(:,:,jlev)+tkelvin
        write (nopost) yhelp(:,:)
      end do
!
!     Salinity.
!
      zcode=-5
      do jlev=1,ken
        ycode=-5.
        write (nopost) yprel,ylen,ycode,yddr(10+jlev-1),ydum,ydum
        yhelp(:,:)=ps(:,:,jlev)
        write (nopost) yhelp(:,:)
      end do
!
!     Velocities.
!
      zcode=-3
      do jlev=1,ken
        ycode=-3.
        write (nopost) yprel,ylen,ycode,yddr(10+jlev-1),ydum,ydum
        yhelp(:,:)=putot(:,:,jlev)
        write (nopost) yhelp(:,:)
      end do
!
      zcode=-4
      do jlev=1,ken
        ycode=-4.
        write (nopost) yprel,ylen,ycode,yddr(10+jlev-1),ydum,ydum
        yhelp(:,:)=pvtot(:,:,jlev)
        write (nopost) yhelp(:,:)
      end do
!
      zcode=-7
      do jlev=1,ken
        ycode=-7.
        write (nopost) yprel,ylen,ycode,yddr(10+jlev-1),ydum,ydum
        yhelp(:,:)=pw(:,:,jlev)
        write (nopost) yhelp
      end do
!
!     Now 2-d fields:
!
!  2-d fields output in real*4
      ylev=-100.
      ycode=-37.
      write (nopost) yprel,ylen,ycode,ylev,ydum,ydum
      yhelp(:,:)=pub(:,:)
      write (nopost) yhelp
      ycode=-38.
      write (nopost) yprel,ylen,ycode,ylev,ydum,ydum
      yhelp(:,:)=pvb(:,:)
      write (nopost) yhelp
!
!     Surface elevation and ice thickness.
!
      ycode=-1.
      write (nopost) yprel,ylen,ycode,ylev,ydum,ydum
      yhelp(:,:)=pzeta(:,:)
      write (nopost) yhelp
!
      ycode=-13.
      write (nopost) yprel,ylen,ycode,ylev,ydum,ydum
      yhelp(:,:)=psice(:,:)
      write (nopost) yhelp
!
!     Horizontal stream function.
!
      ycode=-27.
      write (nopost) yprel,ylen,ycode,ylev,ydum,ydum
      yhelp(:,:)=ppsi(:,:)
      write (nopost) yhelp
!
!     Store topography in vector-points.
!
      ycode=-99.
      write (nopost) yprel,ylen,ycode,ylev,ydum,ydum
      yhelp(:,:)=depth(:,:)
      write (nopost) yhelp
!
!     Store topography in scalar-points.
!
      ycode=-98.
      write (nopost) yprel,ylen,ycode,ylev,ydum,ydum
      yhelp(:,:)=depp(:,:)
      write (nopost) yhelp
!
!     Store *pfluwat*.
!
      ycode=-65.
      write (nopost) yprel,ylen,ycode,ylev,ydum,ydum
      yhelp(:,:)=pfluwat(:,:)
      write (nopost) yhelp
!
!     Store *pflukwat*.
!
      ycode=-67.
      write (nopost) yprel,ylen,ycode,ylev,ydum,ydum
      yhelp(:,:)=pflukwat(:,:)
      write (nopost) yhelp
!
!     Store *pflukhea*.
!
      ycode=-68.
      write (nopost) yprel,ylen,ycode,ylev,ydum,ydum
      yhelp(:,:)=pflukhea(:,:)
      write (nopost) yhelp
!
!     Store convective adjustments.
!
      ycode=-66.
      ylev=-100.
      write (nopost) yprel,ylen,ycode,ylev,ydum,ydum
      yhelp(:,:)=pconvad(:,:,1)
      write (nopost) yhelp
      ycode=-69.
      do k=2,ken
        ylev=dw(k-1)
        write (nopost) yprel,ylen,ycode,ylev,ydum,ydum
        yhelp(:,:)=pconvad(:,:,k)
        write (nopost) yhelp
      end do
!
!     Store *ptaux*
!
      ycode=52.
      ylev=-100.
      write (nopost) yprel,ylen,ycode,ylev,ydum,ydum
      yhelp(:,:)=ptaux(:,:)*rhonul
      write (nopost) yhelp
!
!     Store *ptauy*
!
      ycode=53.
      ylev=-100.
      write (nopost) yprel,ylen,ycode,ylev,ydum,ydum
      yhelp(:,:)=ptauy(:,:)*rhonul
      write (nopost) yhelp
!
!     Store *ptbound*
!
      ycode=-92.
      ylev=-100.
      write (nopost) yprel,ylen,ycode,ylev,ydum,ydum
      yhelp(:,:)=ptbound(:,:)+tkelvin
      write (nopost) yhelp
!
!     Store *pfluhea*.
!
      ycode=-18.
      write (nopost) yprel,ylen,ycode,ylev,ydum,ydum
      yhelp(:,:)=pfluhea(:,:)
      write (nopost) yhelp

!
!*    3.        Store output on file *filemea*.
!     -------------------------------
      close (nopost)
      return
      end subroutine outmea
