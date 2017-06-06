!     =================
!     SUBROUTINE OUTINI
!     =================

      subroutine outini
      use pumamod

!     Output has always 32 bit precision

      integer (kind=4) :: ihead(8)   ! header of first data set
      real    (kind=4) :: zsig(NUGP) ! first block contains settings

      call ntomin(nstep,nmin,nhour,nday,nmonth,nyear)

      ihead(1) = 333  ! ID for PUMA/PLASIM parameter block
      ihead(2) = 0
      ihead(3) = nday + 100 * nmonth + 10000 * nyear
      ihead(4) = 0
      ihead(5) = NLON
      ihead(6) = NLAT
      ihead(7) = NLEV
      ihead(8) = NTRU

!     The first data block with Code = 333 is used for transferring
!     planet properties to the postprocessor "burn"

      zsig(:)      = 0.0       ! initialize
      zsig(1:NLEV) = sigmah(:) ! vertical coordinate table
      zsig(NLEV+1) = n_days_per_year

      open  (40,file=plasim_output,form='unformatted')
      write (40) ihead(:)
      write (40) zsig(:)

      return
      end

!     ==================
!     SUBROUTINE WRITEGP
!     ==================

      subroutine writegp(kunit,pf,kcode,klev)
      use pumamod
      real :: pf(NHOR)
      real :: zf(NUGP)
      integer(kind=4) :: ihead(8)
      integer(kind=4) :: la(NPGP)
      real(kind=4) :: zmin
      real(kind=4) :: zsca
      real(kind=4) :: zzf(NUGP)

      istep = nstep
      call ntomin(istep,nmin,nhour,nday,nmonth,nyear)

      ihead(1) = kcode
      ihead(2) = klev
      ihead(3) = nday + 100 * nmonth + 10000 * nyear
      ihead(4) = nmin + 100 * nhour
      ihead(5) = NLON
      ihead(6) = NLAT
      ihead(7) = nstep - nstep1
      ihead(8) = n_days_per_year

      call mpgagp(zf,pf,1)

      if (mypid == NROOT) then
         write (kunit) ihead
         zzf(:) = zf(:)
         write (kunit) zzf
      endif

      return
      end

!     ==================
!     SUBROUTINE WRITESP
!     ==================

      subroutine writesp(kunit,pf,kcode,klev,pscale,poff)
      use pumamod
      real :: pf(NRSP)
      integer(kind=4) :: ihead(8)
      integer(kind=4) :: la(NTP1+1:NCSP)
      real(kind=4) :: zf(NRSP)
      real(kind=4) :: za(NTP1+1)

      istep = nstep
      call ntomin(istep,nmin,nhour,nday,nmonth,nyear)

      ihead(1) = kcode
      ihead(2) = klev
      ihead(3) = nday + 100 * nmonth + 10000 * nyear
      ihead(4) = nmin + 100 * nhour
      ihead(5) = NRSP
      ihead(6) = 1
      ihead(7) = nstep - nstep1
      ihead(8) = n_days_per_year

!     normalize ECHAM compatible and scale to physical dimensions

      zf(:) = pf(:) * spnorm(1:NRSP) * pscale
      zf(1) = zf(1) + poff ! Add offset if necessary
      write (kunit) ihead
      write (kunit) zf

      return
      end

!     ================
!     SUBROUTINE OUTSP
!     ================

      subroutine outsp
      use pumamod

!     ************
!     * orograpy *
!     ************

      call writesp(40,so,129,0,CV*CV,0.)

!     ************
!     * pressure *
!     ************

      call writesp(40,sp,152,0,1.0,log(psurf))

!     ***************
!     * temperature *
!     ***************

      do jlev = 1 , NLEV
         call writesp(40,st(1,jlev),130,jlev,ct,t0(jlev) * ct)
      enddo

!     *********************
!     * specific humidity *
!     *********************

      if (nqspec == 1) then
         do jlev = 1 , NLEV
            call writesp(40,sqout(1,jlev),133,jlev,1.0,0.0)
         enddo
      endif

!     **************
!     * divergence *
!     **************

      do jlev = 1 , NLEV
         call writesp(40,sd(1,jlev),155,jlev,ww,0.0)
      enddo

!     *************
!     * vorticity *
!     *************

      do jlev = 1 , NLEV
         zsave = sz(3,jlev)
         sz(3,jlev) = sz(3,jlev) - plavor
         call writesp(40,sz(1,jlev),138,jlev,ww,0.0)
         sz(3,jlev) = zsave
      enddo

      return
      end


!     ================
!     SUBROUTINE OUTGP
!     ================

      subroutine outgp
      use pumamod

!     *********************
!     * specific humidity *
!     *********************

      if (nqspec == 0) then ! Semi Langrangian advection active
         do jlev = 1 , NLEV
            call writegp(40,dq(1,jlev),133,jlev)
         enddo
      endif

!     **********************************
!     * mixed-layer depth (from ocean) *
!     **********************************

      call writegp(40,dmld,110,0)

!     ***********************
!     * surface temperature *
!     ***********************

      call writegp(40,dt(1,NLEP),139,0)

!     ****************
!     * soil wetness *
!     ****************

      call writegp(40,dwatc,140,0)

!     **************
!     * snow depth *
!     **************

      call writegp(40,dsnow,141,0)

!     **********************
!     * large scale precip *
!     **********************

      aprl(:)=aprl(:)/real(naccuout)
      call writegp(40,aprl,142,0)

!     *********************
!     * convective precip *
!     *********************

      aprc(:)=aprc(:)/real(naccuout)
      call writegp(40,aprc,143,0)

!     *************
!     * snow fall *
!     *************

      aprs(:)=aprs(:)/real(naccuout)
      call writegp(40,aprs,144,0)

!     **********************
!     * sensible heat flux *
!     **********************

      ashfl(:)=ashfl(:)/real(naccuout)
      call writegp(40,ashfl,146,0)

!     ********************
!     * latent heat flux *
!     ********************

      alhfl(:)=alhfl(:)/real(naccuout)
      call writegp(40,alhfl,147,0)

!     ************************
!     * liquid water content *
!     ************************

      do jlev = 1 , NLEV
         call writegp(40,dql(1,jlev),161,jlev)
      enddo

!     *************
!     * u-star**3 *
!     *************

      call writegp(40,dust3,159,0)

!     **********
!     * runoff *
!     **********

      aroff(:)=aroff(:)/real(naccuout)
      call writegp(40,aroff,160,0)

!     ***************
!     * cloud cover *
!     ***************

      do jlev = 1 , NLEV
        call writegp(40,dcc(1,jlev),162,jlev)
      enddo
      acc(:)=acc(:)/real(naccuout)
      call writegp(40,acc,164,0)

!     ***************************
!     * surface air temperature *
!     ***************************

      atsa(:)=atsa(:)/real(naccuout)
      call writegp(40,atsa,167,0)

!     ******************************
!     * surface temperature (accu) *
!     ******************************

      ats0(:)=ats0(:)/real(naccuout)
      call writegp(40,ats0,169,0)

!     *************************
!     * deep soil temperature *
!     *************************

      call writegp(40,dtd5,170,0)

!     *****************
!     * land sea mask *
!     *****************

      call writegp(40,dls,172,0)

!     *********************
!     * surface roughness *
!     *********************

      call writegp(40,dz0,173,0)

!     **********
!     * albedo *
!     **********

      call writegp(40,dalb,175,0)

!     ***************************
!     * surface solar radiation *
!     ***************************

      assol(:)=assol(:)/real(naccuout)
      call writegp(40,assol,176,0)

!     *****************************
!     * surface thermal radiation *
!     *****************************

      asthr(:)=asthr(:)/real(naccuout)
      call writegp(40,asthr,177,0)

!     ***********************
!     * top solar radiation *
!     ***********************

      atsol(:)=atsol(:)/real(naccuout)
      call writegp(40,atsol,178,0)

!     *************************
!     * top thermal radiation *
!     *************************

      atthr(:)=atthr(:)/real(naccuout)
      call writegp(40,atthr,179,0)

!     ************
!     * u-stress *
!     ************

      ataux(:)=ataux(:)/real(naccuout)
      call writegp(40,ataux,180,0)

!     *************
!     * v- stress *
!     *************

      atauy(:)=atauy(:)/real(naccuout)
      call writegp(40,atauy,181,0)

!     ***************
!     * evaporation *
!     ***************

      aevap(:)=aevap(:)/real(naccuout)
      call writegp(40,aevap,182,0)

!     *********************
!     * soil temperature *
!     *********************

      call writegp(40,dtsoil,183,0)

!     ***********************************
!     * maximum surface air temperature *
!     ***********************************

      call writegp(40,atsama,201,0)

!     ***********************************
!     * minimum surface air temperature *
!     ***********************************

      call writegp(40,atsami,202,0)

!     ********************
!     * top solar upward *
!     ********************

      atsolu(:)=atsolu(:)/real(naccuout)
      call writegp(40,atsolu,203,0)

!     ************************
!     * surface solar upward *
!     ************************

      assolu(:)=assolu(:)/real(naccuout)
      call writegp(40,assolu,204,0)

!     **************************
!     * surface thermal upward *
!     **************************

      asthru(:)=asthru(:)/real(naccuout)
      call writegp(40,asthru,205,0)

!     *******************************
!     * soil temperatures level 2-4 *
!     *******************************

      call writegp(40,dtd2,207,0)
      call writegp(40,dtd3,208,0)
      call writegp(40,dtd4,209,0)

!     *****************
!     * sea ice cover *
!     *****************

      call writegp(40,dicec,210,0)

!     *********************
!     * sea ice thickness *
!     *********************

      call writegp(40,diced,211,0)

!     ****************
!     * forest cover *
!     ****************

      call writegp(40,dforest,212,0)

!     *************
!     * snow melt *
!     *************

      asmelt(:)=asmelt(:)/real(naccuout)
      call writegp(40,asmelt,218,0)

!     *********************
!     * snow depth change *
!     *********************

      asndch(:)=asndch(:)/real(naccuout)
      call writegp(40,asndch,221,0)

!     ******************
!     * field capacity *
!     ******************

      call writegp(40,dwmax,229,0)

!     *****************************************
!     * vertical integrated specific humidity *
!     *****************************************

      aqvi(:)=aqvi(:)/real(naccuout)
      call writegp(40,aqvi,230,0)

!     ****************
!     * glacier mask *
!     ****************

      call writegp(40,dglac,232,0)

!     *********************
!     ***   S I M B A   ***
!     *********************

      if (nveg > 0) call vegout

      return
      end

!     ==================
!     SUBROUTINE OUTDIAG
!     ==================

      subroutine outdiag
      use pumamod

!     *****************************************
!     * 2-D diagnostic arrays, if switched on *
!     *****************************************

      if(ndiagsp2d > 0 .and. mypid == NROOT) then
       do jdiag=1,ndiagsp2d
        jcode=50+jdiag
        call writesp(40,dsp2d(1,jdiag),jcode,0,1.,0.0)
       enddo
      end if

!     *****************************************
!     * 3-D diagnostic arrays, if switched on *
!     *****************************************

      if(ndiagsp3d > 0 .and. mypid == NROOT) then
       do jdiag=1,ndiagsp3d
        jcode=60+jdiag
        do jlev=1,NLEV
         call writesp(40,dsp3d(1,jlev,jdiag),jcode,jlev,1.,0.0)
        enddo
       enddo
      end if

!     *****************************************
!     * 2-D diagnostic arrays, if switched on *
!     *****************************************

      if(ndiaggp2d > 0) then
       do jdiag=1,ndiaggp2d
        jcode=jdiag
        call writegp(40,dgp2d(1,jdiag),jcode,0)
       enddo
      end if

!     *****************************************
!     * 3-D diagnostic arrays, if switched on *
!     *****************************************

      if(ndiaggp3d > 0) then
       do jdiag=1,ndiaggp3d
        jcode=20+jdiag
        do jlev=1,NLEV
         call writegp(40,dgp3d(1,jlev,jdiag),jcode,jlev)
        enddo
       enddo
      end if

!     ************************************************
!     * cloud forcing (clear sky fluxes) diagnostics *
!     ************************************************

      if(ndiagcf > 0) then
       call writegp(40,dclforc(1,1),101,0)
       call writegp(40,dclforc(1,2),102,0)
       call writegp(40,dclforc(1,3),103,0)
       call writegp(40,dclforc(1,4),104,0)
       call writegp(40,dclforc(1,5),105,0)
       call writegp(40,dclforc(1,6),106,0)
       call writegp(40,dclforc(1,7),107,0)
      end if

!     **************************************
!     * entropy diagnostics if switched on *
!     **************************************

      if(nentropy > 0) then
       do jdiag=1,36
        jcode=319+jdiag
        if(jcode == 333) cycle                      !333 is reserved
        call writegp(40,dentropy(1,jdiag),jcode,0)
       enddo
      end if
      if(nentro3d > 0) then
       do jdiag=1,23
        jcode=419+jdiag
        do jlev=1,NLEV
         call writegp(40,dentro3d(1,jlev,jdiag),jcode,jlev)
        enddo
       enddo
      end if

!     *************************************
!     * energy diagnostics if switched on *
!     *************************************

      if(nenergy > 0) then
       do jdiag=1,28
        jcode=359+jdiag
        call writegp(40,denergy(1,jdiag),jcode,0)
       enddo
      end if
      if(nener3d > 0) then
       do jdiag=1,28
        jcode=459+jdiag
        do jlev=1,NLEV
         call writegp(40,dener3d(1,jlev,jdiag),jcode,jlev)
        enddo
       enddo
      end if
!
      return
      end

!     ===================
!     SUBROUTINE OUTRESET
!     ===================

      subroutine outreset
      use pumamod
!
!     reset accumulated arrays and counter
!
      aprl(:)=0.
      aprc(:)=0.
      aprs(:)=0.
      aevap(:)=0.
      ashfl(:)=0.
      alhfl(:)=0.
      acc(:)=0.
      assol(:)=0.
      asthr(:)=0.
      atsol(:)=0.
      atthr(:)=0.
      assolu(:)=0.
      asthru(:)=0.
      atsolu(:)=0.
      ataux(:)=0.
      atauy(:)=0.
      aroff(:)=0.
      asmelt(:)=0.
      asndch(:)=0.
      aqvi(:)=0.
      atsa(:)=0.
      ats0(:)=0.
      atsami(:)=1.E10
      atsama(:)=0.

      naccuout=0

!     ************************************************
!     * cloud forcing (clear sky fluxes) diagnostics *
!     ************************************************

      if(ndiagcf > 0) then
       dclforc(:,:)=0.
      end if
!
      return
      end

!     ==================
!     SUBROUTINE OUTACCU
!     ==================

      subroutine outaccu
      use pumamod
!
!     accumulate diagnostic arrays
!
      aprl(:)=aprl(:)+dprl(:)
      aprc(:)=aprc(:)+dprc(:)
      aprs(:)=aprs(:)+dprs(:)
      aevap(:)=aevap(:)+devap(:)
      ashfl(:)=ashfl(:)+dshfl(:)
      alhfl(:)=alhfl(:)+dlhfl(:)
      acc(:)=acc(:)+dcc(:,NLEP)
      assol(:)=assol(:)+dswfl(:,NLEP)
      asthr(:)=asthr(:)+dlwfl(:,NLEP)
      atsol(:)=atsol(:)+dswfl(:,1)
      atthr(:)=atthr(:)+dlwfl(:,1)
      assolu(:)=assolu(:)+dfu(:,NLEP)
      asthru(:)=asthru(:)+dftu(:,NLEP)
      atsolu(:)=atsolu(:)+dfu(:,1)
      ataux(:)=ataux(:)+dtaux(:)
      atauy(:)=atauy(:)+dtauy(:)
      aroff(:)=aroff(:)+drunoff(:)
      asmelt(:)=asmelt(:)+dsmelt(:)
      asndch(:)=asndch(:)+dsndch(:)
      aqvi(:)=aqvi(:)+dqvi(:)
      atsa(:)=atsa(:)+dtsa(:)
      ats0(:)=ats0(:)+dt(:,NLEP)
      atsami(:)=AMIN1(atsami(:),dtsa(:))
      atsama(:)=AMAX1(atsama(:),dtsa(:))

      naccuout=naccuout+1
!
      return
      end
