      module pumamod

!     ****************************************************************
!     *                      Puma Pre Processor                      *
!     ****************************************************************
!     *                 E. Kirk & T. Kunz & F. Lunkeit               *
!     *                   Meteorologisches Institut                  *
!     *                      Universitaet Hamburg                    *
!     *                        Bundesstrasse 55                      *
!     *                          20146 HAMBURG                       *
!     * 20-Apr-2009                 GERMANY                          *
!     ****************************************************************

!     ****************************************************************
!     * Insert your own code for modification of initial data to:    *
!     * subroutine modify_orography                                  *
!     * subroutine modify_ground_temperature                         *
!     ****************************************************************

!     ****************************************************************
!     * PUMA in its default setup can run without initial data       *
!     * The default setup is an aqua planet with no orography and    *
!     * zonally symmetric forcing (newtonian cooling with Tr)        *
!     * The atmosphere starts at rest with no horizontal gradients   *
!     *                                                              *
!     * This preprocessor program performs following tasks:          *
!     * 1) Prepare a realistic orography (T21 or T42)                *
!     * 2) Enable user modification of this orography                *
!     * 3) Enable user modification of the ground temperature field  *
!     * 4) Adjust vertical profiles of Restoration Temperature       *
!     * 5) Adjust the mean value of surface pressure                 *
!     * 6) Build an initial Ps field adjusted to orography           *
!     * 7) Setup Yoden-profiles                                      *
!     *                                                              *
!     *   Inputfile: <Naaa_surf_0129.sra> with aaa=032,048,064, ...  *
!     *              <Naaa_surf_0139.sra> : Ts anomalies             *
!     * Outputfiles: <Naaa_surf_0129.sra> : topography [m2/s2]       *
!     *              <Naaa_surf_0134.sra> : Surface pressure [hPa]   *
!     *              <Naaa_surf_0121.sra> : Constant part of Tr      *
!     *              <Naaa_surf_0122.sra> : Variable part of Tr      *
!     *              <Naaa_surf_0123.sra> : Damping time scales      *
!     *                                                              *
!     * The outputfiles contain topography, surface pressure,        *
!     * the part of Tr, that is constant in time and the part pf Tr  *
!     * that can be modulated by an annual cycle.                    *
!     *                                                              *
!     * All files are written formatted, such avoiding the problems  *
!     * assigned to big endian and little endian machines            *
!     ****************************************************************


! ****************************************************************
! * The horizontal resolution of PUMA by the number of latitudes *
! * nlat is read from file "resolution_namelist"                 *
! ****************************************************************
integer :: nlat = 32

! example values:  32,  48,  64, 128,  192,  256,  512,  1024
! truncation:     T21, T31, T42, T85, T127, T170, T341,  T682


! *****************************************************************
! * The number of sigma levels of PUMA are modified after reading *  
! * file <resolution_namelist>.                                   *     
! *****************************************************************
integer :: nlev =   10 ! Levels


! *****************************************************************!
! * Grid related paramters, which are reset after reading the file !
! * <resolution_namelist>. All parameters are initialized for the  !
! * T21 truncation                                                 !
! *****************************************************************!
integer :: nlon =   64 ! Longitudes = 2 * latitudes
integer :: ntru =   21 ! (nlon-1) / 3
integer :: nlpp =   32 ! Latitudes per process
integer :: nhor = 2048 ! Horizontal part
integer :: nlem =    9 ! Levels - 1
integer :: nlep =   11 ! Levels + 1
integer :: nlsq =  100 ! Levels squared
integer :: ntp1 =   22 ! ntru + 1
integer :: nrsp =  506 ! (ntru+1) * (ntru+2)
integer :: ncsp =  253 ! nrsp / 2
integer :: nspp =  506 ! nodes per process
integer :: nesp =  506 ! number of extended modes
integer :: nzom =   44 ! Number of zonal modes

integer :: nlah =   16 ! Half of latitudes
integer :: nhpp =   16 ! Half latitudes per process

!     ****************************************************************
!     * Don't touch the following parameter definitions !            *
!     ****************************************************************

      parameter(AKAP   = 0.286)            ! Kappa
      parameter(ALR    = 0.0065)           ! Lapse rate
      parameter(EZ     = 1.632993161855452D0) ! ez = 1 / sqrt(3/8)
      parameter(GA     = 9.81)             ! Gravity
      parameter(GASCON = 287.0)            ! Gas constant
      parameter(PI     = 3.141592653589793D0) ! Pi
      parameter(TWOPI  = PI + PI)             ! 2 Pi
      parameter(PLARAD_EARTH = 6371000.0)        ! Planet radius
      parameter(SID_DAY_EARTH= 86164.)      ! Siderial day Earth 23h 56m 04s
      parameter(PNU    = 0.02)             ! Time filter
      parameter(PNU21  = 1.0 - 2.0*PNU)    ! Time filter 2
      parameter(PSURF  = 101100.0)         ! Surface pressure [Pa]
      parameter(WW     = 0.00007292)       ! Rotation speed [1/sec]
      parameter(CV     = PLARAD_EARTH * WW)      ! cv
      parameter(CT     = CV*CV/GASCON)     ! ct
      parameter(OSCAR  = CV*CV/GA)         ! Scale Orography


! *************
! * filenames *
! *************
character (256) :: resolution_namelist = "resolution_namelist"
character (256) :: puma_namelist       = "puma_namelist"
character (256) :: ppp_puma_txt        = "ppp-puma.txt"

! *****************************************************************
! * For multiruns the instance number is appended to the filename *
! * e.g.: puma_namelist_1 puma_diag_1 etc. for instance # 1       *
! *****************************************************************










!     **************************
!     * Global Integer Scalars *
!     **************************

      integer :: kick     =  1 ! kick > 1 initializes eddy generation
      integer :: mpstep   =  0 ! PUMA
      integer :: nafter   =  0 ! PUMA
      integer :: nwpd     =  1 ! PUMA
      integer :: ncoeff   =  0 ! number of modes to print
      integer :: ndel     =  6 ! ndel
      integer :: ndiag    = 12 ! write diagnostics interval
      integer :: nkits    =  3 ! number of initial timesteps
      integer :: newsr    =  0 ! force sr recalculatio
      integer :: nlevt    =  9 ! tropospheric levels (set_vertical_grid)

      integer :: ngui     =  0 ! PUMA variable
      integer :: nguidbg  =  0 ! PUMA variable
      integer :: nqspec   =  1 ! PLASIM variable
      integer :: nrun     =  0 ! PUMA variable
      integer :: nstep    =  0 ! current timestep
      integer :: nstop    =  0 ! finishing timestep
      integer :: nsync    =  0 ! PUMA variable
      integer :: ntspd    = 24 ! number of timesteps per day
      integer :: ncu      =  0 ! check unit (debug output)
      integer :: noro     =  0 ! 1: read orography
                               ! 2: orography is computed (sine wave)
                               ! 3: orography is computed (gauss)

      integer :: norox    =  0 ! x-wavenumber of idealized orography
      integer :: nreverse =  0 ! t-gradient reversal in stratosphere; 0:no 1:yes
      integer :: ncorrect =  0 ! correct tr due to orography (0:no, 1:yes)
      integer :: nsym     =  0 ! produces total symmetric initial conditions 

      integer :: nsrv     =  1 ! 1: write gridpoint fields for diagnostics
      integer :: lon1oro  =  0 ! Define rectangle for anomaly
      integer :: lon2oro  =  0 ! 
      integer :: lat1oro  =  0 ! 
      integer :: lat2oro  =  0 ! 
      integer :: ntgr     =  0 ! 1: read ground temperature
      integer :: lon1tgr  =  0 ! Define rectangle for anomaly
      integer :: lon2tgr  =  0 ! 
      integer :: lat1tgr  =  0 ! 
      integer :: lat2tgr  =  0 ! 
      integer :: nstrato  =  0 ! 1: Torben's stratosphere forcing
      integer :: nsponge  =  0 ! Switch for sponge layer
      integer :: nhelsua  =  0 ! 1: Set up Held & Suarez T_R field
!                                   instead of original PUMA T_R field
!                                2: Set up Held & Suarez T_R field
!                                   instead of original PUMA T_R field
!                                   AND use latitudinally varying
!                                   heating timescale in PUMA (H&Z(94)),
!                                   irrelevant for PumaPreProcessor (ppp)
!                                3: Use latitudinally varying
!                                   heating timescale in PUMA (H&Z(94)),
!                                   irrelevant for PumaPreProcessor (ppp)
      integer  :: ntestgp = 0 ! switch 0/1 produce a second set of restoration 
                              ! restoration temperatures and damping time scales
      integer  :: nvg     = 0 ! type of vertical grid
                              ! 0 = linear
                              ! 1 = Scinocca & Haynes
                              ! 2 = Polvani & Kushner

      integer :: nyoden   = 0 ! > 0 Read yoden profile t0(:) and dt(:)
      integer :: npackgp  = 1 ! used in PUMA
      integer :: npacksp  = 1 ! used in PUMA
      integer :: noutput  = 0 ! used in PUMA
      integer :: noutsrv  = 0 ! used in PUMA


!     These three predifined Yoden profiles may be selected by setting
!     NYODEN= 1, 3 or 5 and nlev=20

      real :: t0yod1(20)  = (/224.14,213.91,211.36,212.95,217.56 &
                             ,224.07,231.24,237.73,243.54,248.87 &
                             ,253.41,257.83,261.72,265.60,268.96 &
                             ,272.33,275.36,278.34,281.06,283.74/)
      real :: dtyod1(20)  = (/  0.00,  0.00,  2.00,  9.19, 19.13 &
                             , 28.26, 34.91, 40.37, 45.10, 49.19 &
                             , 52.05, 54.69, 56.20, 57.70, 58.31 &
                             , 58.91, 59.26, 59.56, 59.76, 59.94/)
      real :: t0yod3(20)  = (/265.14,254.91,246.36,240.95,237.56 &
                             ,234.65,235.24,237.73,243.54,248.87 &
                             ,253.41,257.83,261.72,265.60,268.96 &
                             ,272.33,275.36,278.34,281.06,283.74/)
      real :: dtyod3(20)  = (/  0.00, -3.05,-20.32,-26.19,-28.13 &
                             , 28.26, 34.91, 40.37, 45.10, 49.19 &
                             , 52.05, 54.69, 56.20, 57.70, 58.31 &
                             , 58.91, 59.26, 59.56, 59.76, 59.94/)
      real :: t0yod5(20)  = (/305.00,298.00,292.00,286.00,280.00 &
                             ,274.00,268.00,261.00,256.54,254.87 &
                             ,253.41,257.83,261.72,265.60,268.96 &
                             ,272.33,275.36,278.34,281.06,283.74/)
      real :: dtyod5(20)  = (/  0.00, -3.05,-21.32,-27.19,-29.13 &
                             ,-32.26,-34.91,-40.37,-45.10,-49.19 &
                             , 52.05, 54.69, 56.20, 57.70, 58.31 &
                             , 58.91, 59.26, 59.56, 59.76, 59.94/)
      real :: t0yod7(20)  = (/265.14,254.91,246.36,240.95,237.56 &
                             ,234.65,235.24,237.73,243.54,248.87 &
                             ,253.41,257.83,261.72,265.60,268.96 &
                             ,272.33,275.36,278.34,281.06,283.74/)
      real :: dtyod7(20)  = (/  0.00,-25.05,-64.32,-58.19,-22.13 &
                             , 28.26, 34.91, 40.37, 45.10, 49.19 &
                             , 52.05, 54.69, 56.20, 57.70, 58.31 &
                             , 58.91, 59.26, 59.56, 59.76, 60.00/)
      real :: t0yod8(20)  = (/265.14,254.91,246.36,240.95,237.56 &
                             ,234.65,235.24,237.73,243.54,248.87 &
                             ,253.41,257.83,261.72,265.60,268.96 &
                             ,272.33,275.36,278.34,281.06,283.74/)
      real :: dtyod8(20)  = (/  0.00,-25.05,-64.32,-58.19,-22.13 &
                             , 18.26, 24.91, 30.37, 35.10, 39.19 &
                             , 42.05, 44.69, 46.20, 47.70, 48.31 &
                             , 48.91, 49.26, 49.56, 49.76, 50.00/)
      real :: t0yod9(20)  = (/265.14,254.91,246.36,240.95,237.56 &
                             ,234.65,235.24,237.73,243.54,248.87 &
                             ,253.41,257.83,261.72,265.60,268.96 &
                             ,272.33,275.36,278.34,281.06,283.74/)
      real :: dtyod9(20)  = (/  0.00,  0.00,  0.00,  0.00,  0.00 &
                             , 28.26, 34.91, 40.37, 45.10, 49.19 &
                             , 52.05, 54.69, 56.20, 57.70, 58.31 &
                             , 58.91, 59.26, 59.56, 59.76, 60.00/)

      integer :: nfls     = 1  ! =1: lower stratospheric forcing is applied,
!                                see parameters flsp0, flsdp, flsamp and
!                                flsoff below

!     ***********************
!     * Global Real Scalars *
!     ***********************

      real :: delt             ! 2 pi / ntspd timestep interval
      real :: delt2            ! 2 * delt
!
!     namelist parameter for yoden setup
!     (for arrays to and dt see namelist arrays below)

      real :: dtep  =    60.0  ! delta T equator <-> pole  [K]
      real :: dtns  =     0.0  ! delta T   north <-> south [K]
      real :: dtrop = 12000.0  ! Tropopause height [m]
      real :: dttrp =     2.0  ! Tropopause smoothing [K]
      real :: dtzz  =    10.0  ! delta(Theta)/H additional lapserate in
      real :: syncstr=    0.0  ! PUMA variable
!                                Held & Suarez T_R field
      real :: ttp   =   200.0  ! Tropopause temperature in
!                                Held & Suarez T_R field
      real :: plavor=    EZ    ! planetary vorticity
      real :: rotspd =     1.0  ! rotation speed 1.0 = normal Earth rotation
      real :: sigmax=  6.0e-7  ! sigma for top half level
      real :: tdiss =     0.25 ! diffusion time scale [days]
      real :: tac   =       0. ! length of annual cycle [days] (0 = no cycle)
      real :: pac   =       0. ! phase of the annual cycle [days]
      real :: tgr   =   288.0  ! Ground Temperature in mean profile [K]
      real :: oroano=     0.0  ! Orography anomaly in [gpm]
      real :: tgrano=     0.0  ! Ground temperature anomaly in [K]
!!
      real :: alrs  =   -0.0000! stratospheric lapse rate [K/m]
      real :: horo  =    0.0   ! height of the idealized orography [m]
      real :: orofac=    1.0   ! factor to scale the orograpy
      real :: ttropo=   250.0  ! temp. at tropopause
!!
      real :: dorox  =    0.0   ! lon. base point for gauss mountain [dec]
      real :: doroy  =    0.0   ! lat. base point for gauss mountain [dec]
      real :: doroxs =    0.0   ! gauss mountain scale in lon.-dir. [dec]
      real :: doroys =    0.0   ! gauss mountain scale in lat.-dir. [dec]

      real :: tauta  =    40.0  ! heating timescale far from surface [days]
      real :: tauts  =     4.0  ! heating timescale close to surface [days] 
      real :: sid_day=  SID_DAY_EARTH ! siderial day [sec]



      real :: zsigb  =     0.7  ! sigma_b for Held&Suarez frict. and heat. 
                                ! time scale 



!     * parameter for stratospheric forcing *

      real :: alrpv =   0.002  ! Vertical lapse rate of stratospheric
!                                polar vortex forcing restoration
!                                temperature field. It corresponds to
!                                'gamma' in Polvani & Kushner (2002).
!                                But alrpv is in [K/m], thus,
!                                alrpv=0.002 corresponds to 'gamma=2'
!                                in P&K (2002).
      real :: radpv = 50.      ! Radius of stratospheric polar vortex
!                                forcing [deg latitude]
      real :: edgepv = 10.     ! Width of edge of stratospheric polar
!                                vortex forcing [deg latitude].
!                                If edgepv=0., then no polar vortex is
!                                set up.
      real :: pmaxpv = 10000.  ! Lower boundary of stratospheric polar
!                                vortex forcing, max. pressure [Pa].
!                                If pmaxpv=0. then pmaxpv is set to the
!                                pressure at tropopause height, specified
!                                by dtrop, according to the US standard
!                                atmosphere (USSA) vertical profile, used
!                                for contruction of the restoration
!                                temperature field. With the standard
!                                setting (dtrop=11000m, ALR=0.0065K/m,
!                                tgr=288.15K) this gives pmaxpv=22632Pa.
!                                ++++NOTE: pmaxpv should be within the
!                                second USSA layer, that is the interval
!                                from 54749Pa to 22632Pa (in case of
!                                standard parameter setting).++++

      real :: flsp0 = 10000.   ! pressure of max. lower stratospheric forcing
!                                field [Pa]
      real :: flsdp =  9000.   ! half of vertical extension of lower
!                                stratospheric forcing field [Pa]
      real :: flsamp =  15.    ! amplitude of lower stratospheric forcing
!                                field [K]
      real :: flsoff =  -5.    ! offset of lower stratospheric forcing
!                                field [K]


!     *******************
!     * Namelist Arrays *
!     *******************
!     workaround for older fortran versions, where allocatable arrays
!     are not allowed in namelists


      integer, parameter :: MAXLEV   = 100
      integer, parameter :: MAXSELZW =  42
      integer, parameter :: MAXSELSP = ((MAXSELZW+1) * (MAXSELZW+2)) / 2
      integer :: nselect(0:MAXSELZW) = 1 ! NSELECT can be used up tp T42
      integer :: nspecsel(MAXSELSP)  = 1
      integer :: ndl(MAXLEV)         = 0
      real    :: restim(MAXLEV)      = 15.0
      real    :: sigmah(MAXLEV)      = 0.0
      real    :: t0k(MAXLEV)         = 250.0
      real    :: t0(MAXLEV)          = 250.0
      real    :: tfrc(MAXLEV)        = 0.0
      real    :: dt(MAXLEV)          = 0.0

!     **************************
!     * Global Spectral Arrays *
!     **************************

      real, allocatable :: sor(:)    ! Spectral Orography
      real, allocatable :: ssp(:)    ! Spectral surface pressure
      real, allocatable :: stg(:)    ! Spectral ground temperature
      real, allocatable :: sep(:)    ! Spectral equator-pole gradient
      real, allocatable :: sns(:)    ! Spectral north-south  gradient
      real, allocatable :: spnorm(:) ! ECHAM -> PUMA normalization factors
      real, allocatable :: sr1(:,:)  ! Constant part of Tr
      real, allocatable :: sr2(:,:)  ! Variable part of Tr

!     ***************************
!     * Global Gridpoint Arrays *
!     ***************************

      real, allocatable :: gor(:,:)     ! Orography
      real, allocatable :: gsp(:,:)     ! Surface pressure
      real, allocatable :: gtg(:,:)     ! Ground temperature
      real, allocatable :: gan(:,:)     ! Ground temperature anomaly
      real, allocatable :: gep(:,:)     ! Equator-pole gradient
      real, allocatable :: gns(:,:)     ! North-south  gradient
      real, allocatable :: gtc(:,:,:)   ! Restoration Temperature
      real, allocatable :: gtv(:,:,:)   ! Restoration Temperature NS mode
      real, allocatable :: gtdamp(:,:,:)! Reciprocal of damping time scale for T 
      real, allocatable :: gaf(:,:,:)   ! Anomaly factors
      real, allocatable :: gra(:)       ! Gradient factors

!     *******************
!     * Latitude Arrays *
!     *******************

      character (3), allocatable :: chlat(:) ! label for latitudes

      real (kind=8), allocatable :: sid(:)   ! sin(phi)
      real (kind=8), allocatable :: gwd(:)   ! Gaussian weight
      real (kind=8), allocatable :: csq(:)   ! cos(phi)^2
      real (kind=8), allocatable :: dla(:)   ! phi
 

!     ****************
!     * Level Arrays *
!     ****************

      real, allocatable :: t0d(:)   ! vertical t0k gradient
      real, allocatable :: dsigma(:)
      real, allocatable :: rdsig(:)
      real, allocatable :: sigma(:)
      real, allocatable :: tkp(:)

!     **********************
!     * Dummy declarations *
!     **********************

      real :: gp(2)
      real :: sp(2)
      real :: gpj(2)
      real :: gu(2,2)
      real :: gv(2,2)
      real :: sd(2,2)
      real :: st(2,2)
      real :: sq(2,2)
      real :: sz(2,2)
      real :: gd(2,2)
      real :: gq(2,2)
      real :: gt(2,2)
      real :: gz(2,2)

      end module pumamod

      program ppp

!     open file ppp-puma interface information
      if (mypid == NROOT) then 
         open(95,file="ppp-puma.txt",form='formatted')  
      endif

      call resolution
      call prolog
      call gridpoint

!     close file ppp-puma interface information
      if (mypid == NROOT) then 
         close(95)             
      endif

      stop
      end

 ! ******************************
 ! * SUBROUTINE ALLOCATE_ARRAYS *
 ! ******************************

 subroutine allocate_arrays
 use pumamod

 !---  Global Spectral Arrays 
 allocate(sor(nesp))       ; sor(:)      = 0.0 ! Spectral Orography
 allocate(ssp(nesp))       ; ssp(:)      = 0.0 ! Spectral surface pressure
 allocate(stg(nesp))       ; stg(:)      = 0.0 ! Spectral ground temperature
 allocate(sep(nesp))       ; sep(:)      = 0.0 ! Spectral equator-pole gradient
 allocate(sns(nesp))       ; sns(:)      = 0.0 ! Spectral north-south  gradient
 allocate(spnorm(nesp))    ; spnorm(:)   = 0.0 ! ECHAM -> PUMA normalization
 allocate(sr1(nesp,nlev))  ; sr1(:,:)    = 0.0 ! Constant part of Tr
 allocate(sr2(nesp,nlev))  ; sr2(:,:)    = 0.0 ! Variable part of Tr

 !---  Global Gridpoint Arrays
 allocate(gor(nlon,nlat))        ; gor(:,:)      =  0.0 ! Orography
 allocate(gsp(nlon,nlat))        ; gsp(:,:)      =  0.0 ! Surface pressure
 allocate(gtg(nlon,nlat))        ; gtg(:,:)      =  0.0 ! Ground temperature
 allocate(gan(nlon,nlat))        ; gan(:,:)      =  0.0 ! Ground temperature anomaly
 allocate(gep(nlon,nlat))        ; gep(:,:)      =  0.0 ! Equator-pole gradient
 allocate(gns(nlon,nlat))        ; gns(:,:)      =  0.0 ! North-south  gradient
 allocate(gtc(nlon,nlat,nlev))   ; gtc(:,:,:)    =  0.0 ! Restoration Temperature
 allocate(gtv(nlon,nlat,nlev))   ; gtv(:,:,:)    =  0.0 ! Rest. Temperature NS mode
 allocate(gtdamp(nlon,nlat,nlev)); gtdamp(:,:,:) =  0.0 ! Reciprocal of damping time scale for T
 allocate(gaf(nlon,nlat,nlev))   ; gaf(:,:,:)    =  0.0 ! Anomaly factors

 !---  Latitude Arrays
 allocate(chlat(nlat))  ; chlat(:) = "   "  ! label for latitudes
 allocate(sid(nlat))    ; sid(:)   = 0.0    ! sin(phi)
 allocate(gwd(nlat))    ; gwd(:)   = 0.0    ! Gaussian weight
 allocate(csq(nlat))    ; csq(:)   = 0.0    ! cos(phi)^2
 allocate(dla(nlat))    ; dla(:)   = 0.0    ! phi

 !---  Level Arrays
 allocate(gra(nlev))    ; gra(:)    = 0.0    ! Gradient factors
 allocate(t0d(nlev))    ; t0d(:)    = 0.0    ! vertical t0k gradient
 allocate(dsigma(nlev)) ; dsigma(:) = 0.0
 allocate(rdsig(nlev))  ; rdsig(:)  = 0.0
 allocate(sigma(nlev))  ; sigma(:)  = 0.0
 allocate(tkp(nlev))    ; tkp(:)    = 0.0

 return
 end subroutine allocate_arrays

!     ===========================
!     SUBROUTINE MODIFY_OROGRAPHY
!     ===========================

      subroutine modify_orography(por)
      use pumamod
      real :: por(nlon,nlat)

!     Array <por> contains the orography on a Gaussian grid
!     the units are that of a geopotential [m2/s2] (gpm * g)
!     You may modify the orography here with your own code
!     The new orography is spectrally fitted after this routine
!     A gridpoint representation is written in service format
!     to file <puma_oro_tnn.srv> with nn = spectral truncation.

      if (noro == 2) call mkoro(por)
      if (noro == 3) call mkorog(por)

!     Rectangular anomaly

!     if (oroano /= 0.0 .and. lat1oro > 0 .and. lat2oro <= nlat &
!                       .and. lon1oro > 0 .and. lon2oro <= nlon &
!                       .and. lat1oro < lat2oro .and. lon1oro < lon2oro) then
!        do jlat = lat1oro , lat2oro
!           do jlon = lon1oro , lon2oro
!              por(jlon,jlat) = por(jlon,jlat) + oroano * GA
!              if (por(jlon,jlat) < 0.0) por(jlon,jlat) = 0.0
!           enddo
!        enddo
!     endif

!     Elliptic anomaly

      if (oroano /= 0.0 .and. lat1oro > 0 .and. lat2oro <= nlat &
                        .and. lon1oro > 0 .and. lon2oro <= nlon &
                        .and. lat1oro < lat2oro .and. lon1oro < lon2oro) then
         x0 = (lon1oro + lon2oro) * 0.5
         y0 = (lat1oro + lat2oro) * 0.5
         xf = PI / (lon2oro - lon1oro)
         yf = PI / (lat2oro - lat1oro)
         do jlat = lat1oro , lat2oro
            yb = (jlat - y0) * yf
            do jlon = lon1oro , lon2oro
               xa = (jlon - x0) * xf
               cx = cos(xa)
               cy = cos(yb)
               if (cx > 0.0 .and. cy > 0.0) then
                  por(jlon,jlat) = por(jlon,jlat) + oroano * GA * cx * cy
               endif
               if (por(jlon,jlat) < 0.0) por(jlon,jlat) = 0.0
            enddo
         enddo
      endif

      return
      end

!     ====================================
!     SUBROUTINE MODIFY_GROUND_TEMPERATURE
!     ====================================

      subroutine modify_ground_temperature(ptgr)
      use pumamod
      real :: ptgr(nlon,nlat)

!     Array <ptgr> contains the ground temperature on a Gaussian grid.
!     The units are [K]
!     You may modify the ground temperature here with your own code.
!     The new ground temperature is used to construct the temperature
!     profile of the restoration temperature on each grid point.
!     A gridpoint representation is written in service format
!     to file <puma_gtgr_tnn.srv> with nn = spectral truncation.

      if (tgrano /= 0.0 .and. lat1tgr > 0 .and. lat2tgr <= nlat &
                        .and. lon1tgr > 0 .and. lon2tgr <= nlon &
                        .and. lat1tgr < lat2tgr .and. lon1tgr < lon2tgr) then
         do jlat = lat1tgr , lat2tgr
            do jlon = lon1tgr , lon2tgr
               gan(jlon,jlat) = gan(jlon,jlat) + tgrano
               ptgr(jlon,jlat) = ptgr(jlon,jlat) + tgrano
               if (ptgr(jlon,jlat) < 0.0) ptgr(jlon,jlat) = 0.0
            enddo
         enddo
      endif

      return
      end


!     ================
!     SUBROUTINE MKORO
!     ================

      subroutine mkoro(por)
      use pumamod
!
      real por(nlon,nlat)
!
      zscale=horo*GA*0.5
!
      por(:,:)=0.
!
      do jlat=1,nlat/2
       zlat=dla(jlat)
       zfacy=(sin(2.*zlat))**2
       do jlon=1,nlon
        if(norox > 0) then
         zfacx=(1.+cos(norox*real(jlon)*TWOPI/real(nlon)))
        else
         zfacx=1.
        endif
        por(jlon,jlat)=zscale*zfacx*zfacy
       enddo
      enddo
!
      if(nsym == 1) then
       do jlat=1,nlat/2
        j2=nlat+1-jlat
        por(:,j2)=por(:,jlat)
       enddo
      endif
!
      return
      end

!     =================
!     SUBROUTINE MKOROG
!     =================

      subroutine mkorog(por)
      use pumamod
!
      real por(nlon,nlat)
!
      zscale=horo*GA
!
      zlon0=dorox*PI/180.
      zlat0=doroy*PI/180.
      zlons=(180./(doroxs*PI))**2
      zlats=(180./(doroys*PI))**2
!
      do jlat=1,nlat
       zlat=dla(jlat)
!!     zcos2=cos(zlat)**2
       zcos2=1.
       zdlat2=(zlat-zlat0)**2
       do jlon=1,nlon
        zlon=TWOPI*real(jlon-1)/real(nlon)
        zdlon=abs(zlon-zlon0)
        if(zdlon > PI) zdlon=TWOPI-zdlon
        zdlon2=zdlon**2
        por(jlon,jlat)=zscale*EXP(-zlons*zcos2*zdlon2-zlats*zdlat2)
       enddo
      enddo
!
      if(nsym == 1 .and. zlat0 > 0.) then
       do jlat=1,nlat/2
        j2=nlat+1-jlat
        por(:,j2)=por(:,jlat)
       enddo
      elseif(nsym == 1 .and. zlat0 < 0.) then
       do jlat=1,nlat/2
        j2=nlat+1-jlat
        por(:,jlat)=por(:,j2)
       enddo
      endif
!
      return
      end

!     =================
!     SUBROUTINE PROLOG
!     =================

      subroutine prolog
      use pumamod
      integer :: ival(1)

      call allocate_arrays

      call printparameter
      call inigau(nlat,sid,gwd)
      call inilat
      call legpri
      call readnl
      call initpm
      call legini(nlat,nlpp,nesp,nlev,plavor,sid,gwd)
      call initfd


      if (mypid == NROOT) then
         ival(1) = nlat
         call ppp_write_i('NLAT',1,ival)
         ival(1) = nlev
         call ppp_write_i('NLEV',1,ival)
         ival(1) = nhelsua
         call ppp_write_i('NHELSUA',1,ival)
         call ppp_write_r('SIGMH',NLEV,sigmah)
      endif

      return
      end

!     =================
!     SUBROUTINE INITFD
!     =================

      subroutine initfd
      use pumamod

      dimension zrmean(nlev)

      zfmode0  = sqrt(2.0)
      zfmode1  = 1.0 / sqrt(6.0)
      zfmode2  = -2.0 / 3.0 * sqrt(0.4)

      stg(1)   = zfmode0 * tgr  ! Ground temperature  [K]
      sns(3)   = zfmode1 * dtns ! North-South gradient [K]
      sep(5)   = zfmode2 * dtep ! Equator-Pole gradient [K]

!     Find sigma at dtrop

      zttrop = tgr - dtrop * ALR
      ztps   = (zttrop/tgr)**(GA/(ALR*GASCON))

!     The North-South and Equator-Pole gradients are defined on z=0.
!     gra() modifies the gradient from full at z=0 to zero at tropopause
!     PUMA aquaplanet compatibility mode (sine function used)

      do jlev = 1 , nlev
         gra(jlev) = sin(0.5 * PI * (sigma(jlev) - ztps) / (1.0-ztps))
         if (gra(jlev) < 0.0) gra(jlev) = 0.0
      enddo

      return
      end

!     ==========================
!     SUBROUTINE ANOMALY_FACTORS
!     ==========================

      subroutine anomaly_factors
      use pumamod

      do jlat = 1 , nlat
      do jlon = 1 , nlon

!     Find sigma at dtrop

      zttrop = tgr - dtrop * ALR
      ztps   = (zttrop/gtg(jlon,jlat))**(GA/(ALR*GASCON))

!     The North-South and Equator-Pole gradients are defined on z=0.
!     gaf() modifies the gradient from full at z=0 to zero at tropopause
!     PUMA aquaplanet compatibility mode (sine function used)

      if(nreverse == 0) then
       do jlev = 1 , nlev
        gaf(jlon,jlat,jlev) = sin(0.5 * PI * (sigma(jlev) - ztps) / (1.0-ztps))
        if (gaf(jlon,jlat,jlev) < 0.0) gaf(jlon,jlat,jlev) = 0.0
       enddo
      else
       do jlev = 1 , nlev
         gaf(jlon,jlat,jlev) = sin(0.5 * PI * (sigma(jlev) - ztps) / (1.0-ztps))
         if (sigma(jlev) < ztps)                                        &
     &   gaf(jlon,jlat,jlev) = sin(PI*(sigma(jlev)-ztps)/(ztps-sigma(1)))
       enddo
      endif

      enddo
      enddo

      return
      end


!     ========================
!     SUBROUTINE READ_GAN_GRID
!     ========================

      subroutine read_gan_grid(kread)
      use pumamod

      logical :: lexist
      integer :: ihead(8)
      character(20)  :: ynumber
      character(256) :: yfilename

!     Read temperature anomalies for PUMA
!     Use formatted Service style

      kread = 0

      if (nlat < 1000) then
         write(ynumber,'(I3.3)') nlat
      else
         write(ynumber,'(I4.4)') nlat
      endif

      yfilename = "N" // trim(adjustl(ynumber)) // "_surf_0139.sra"
      inquire(file=yfilename,exist=lexist)
      if (lexist) then
         write(*,*) ' Reading anomaly temperature from file <',trim(yfilename),'>'
         open (65,file=yfilename,form='formatted')
         read (65,*) ihead(:)
         read (65,*) gan(:,:)
         close(65)
         kread = 1
      else
         write(*,*) ' No anomaly temperature file'
      endif
      return
      end subroutine read_gan_grid


!     ========================
!     SUBROUTINE READ_ORO_GRID
!     ========================

      subroutine read_oro_grid(kread)
      use pumamod

      logical :: lexist
      integer :: ihead(8)
      character(20)  :: ynumber
      character(256) :: yfilename

!     Read orography for PUMA
!     Use formatted Service style

      kread = 0

      if (nlat < 1000) then
         write(ynumber,'(I3.3)') nlat
      else
         write(ynumber,'(I4.4)') nlat
      endif

      yfilename = "N" // trim(adjustl(ynumber)) // "_surf_0129.sra"
      inquire(file=yfilename,exist=lexist)
      if (lexist) then
         write(*,*) ' Reading orography from file <',trim(yfilename),'>'
         open (65,file=yfilename,form='formatted')
         read (65,*) ihead(:)
         read (65,*) gor(:,:)
         close(65)
         kread = 1
      else
         gor(:,:) = 0.0
         write(*,*) ' No orography file - starting with zero orography'
      endif
      return
      end subroutine read_oro_grid


!     ====================
!     SUBROUTINE WRITE_ORO
!     ====================

      subroutine write_oro
      use pumamod

      dimension itime(8)
      dimension ihead(8)

      character(20)  :: ynumber
      character(256) :: yfilename

      call date_and_time(values=itime)

!     Write orography for PUMA
!     Use formatted Service style

      if (nlat < 1000) then
         write(ynumber,'(I3.3)') nlat
      else
         write(ynumber,'(I4.4)') nlat
      endif

      yfilename = "N" // trim(adjustl(ynumber)) // "_surf_0129.sra"
      open(60,file=yfilename,form='formatted')
      ihead(1) = 129   ! code for orography
      ihead(2) =   0   ! level
      ihead(3) = itime(1) * 10000 + itime(2) * 100 + itime(3) ! YYYYMMDD
      ihead(4) = itime(5) * 100   + itime(6)                  ! HHMM
      ihead(5) = nlon  ! 1. dimension
      ihead(6) = nlat  ! 2. dimension
      ihead(7) =    1
      ihead(8) =    0

      write (60,'(8I10)') ihead(:)
      write (60,'(8F10.3)') gor(:,:)

      close(60)

      return
      end subroutine write_oro


!     ===================
!     SUBROUTINE WRITE_PS
!     ===================

      subroutine write_ps
      use pumamod

      dimension itime(8)
      dimension ihead(8)

      character(20)  :: ynumber
      character(256) :: yfilename

      real :: zpres(nlon,nlat)

      call date_and_time(values=itime)

!     Write surface pressure for PUMA
!     Use formatted Service style

      if (nlat < 1000) then
         write(ynumber,'(I3.3)') nlat
      else
         write(ynumber,'(I4.4)') nlat
      endif
      yfilename = "N" // trim(adjustl(ynumber)) // "_surf_0134.sra"
      open(60,file=yfilename,form='formatted')
      ihead(1) = 134   ! code for surface pressure [hPa]
      ihead(2) =   0   ! level
      ihead(3) = itime(1) * 10000 + itime(2) * 100 + itime(3) ! YYYYMMDD
      ihead(4) = itime(5) * 100   + itime(6)                  ! HHMM
      ihead(5) = nlon  ! 1. dimension
      ihead(6) = nlat  ! 2. dimension
      ihead(7) =    1
      ihead(8) =    0

      zpres(:,:) = 0.01 * PSURF * exp(gsp(:,:)) ! Store as [hPa]

      write (60,'(8I10)') ihead(:)
      write (60,'(8F10.4)') zpres(:,:)

      close(60)

      return
      end subroutine write_ps


!     ====================
!     SUBROUTINE WRITE_GTC
!     ====================

      subroutine write_gtc
      use pumamod

      dimension itime(8)
      dimension ihead(8)

      character(20)  :: ynumber
      character(256) :: yfilename

      call date_and_time(values=itime)

!     Write constant part of Tr
!     Use formatted Service style

      if (nlat < 1000) then
         write(ynumber,'(I3.3)') nlat
      else
         write(ynumber,'(I4.4)') nlat
      endif
      yfilename = "N" // trim(adjustl(ynumber)) // "_surf_0121.sra"

      open(60,file=yfilename,form='formatted')
      ihead(1) = 121   ! code for Tr const
      ihead(2) =   0   ! level
      ihead(3) = itime(1) * 10000 + itime(2) * 100 + itime(3) ! YYYYMMDD
      ihead(4) = itime(5) * 100   + itime(6)                  ! HHMM
      ihead(5) = nlon  ! 1. dimension
      ihead(6) = nlat  ! 2. dimension
      ihead(7) =    1
      ihead(8) =    0

      do jlev = 1 , nlev
         ihead(2) = jlev
         write (60,'(8I10)') ihead(:)
         write (60,'(8F10.4)') gtc(:,:,jlev)
      enddo

      close(60)

      return
      end subroutine write_gtc



!     ====================
!     SUBROUTINE WRITE_GTV
!     ====================

      subroutine write_gtv
      use pumamod

      dimension itime(8)
      dimension ihead(8)

      character(20)  :: ynumber
      character(256) :: yfilename

      call date_and_time(values=itime)

!     Write variable part of Tr
!     Use formatted Service style

      if (nlat < 1000) then
         write(ynumber,'(I3.3)') nlat
      else
         write(ynumber,'(I4.4)') nlat
      endif
      yfilename = "N" // trim(adjustl(ynumber)) // "_surf_0122.sra"
      open(60,file=yfilename,form='formatted')
      ihead(1) = 122   ! code for Tr variable
      ihead(2) =   0   ! level
      ihead(3) = itime(1) * 10000 + itime(2) * 100 + itime(3) ! YYYYMMDD
      ihead(4) = itime(5) * 100   + itime(6)                  ! HHMM
      ihead(5) = nlon  ! 1. dimension
      ihead(6) = nlat  ! 2. dimension
      ihead(7) =    1
      ihead(8) =    0

      do jlev = 1 , nlev
         ihead(2) = jlev
         write (60,'(8I10)') ihead(:)
         write (60,'(8F10.4)') gtv(:,:,jlev)
      enddo

      close(60)

      return
      end subroutine write_gtv


!     ========================
!     SUBROUTINE WRITE_VARGP2D
!     ========================

      subroutine write_vargp2D(zgp,kcode)
      use pumamod

      dimension itime(8)
      dimension ihead(8)

      character(20)  :: ynumber
      character(256) :: yfilename
      real :: zgp(nlon,nlat)


      call date_and_time(values=itime)

!     produce file name to be written
      if (NLAT < 1000) then
         write(yfilename,'("N",I3.3,"_surf_",I4.4,".sra")') NLAT,kcode
      else
         write(yfilename,'("N",I4.4,"_surf_",I4.4,".sra")') NLAT,kcode
      endif



      open(60,file=yfilename,form='formatted')
      ihead(1) = kcode  ! code for reciprocal of damping time scale 
      ihead(2) =   0    ! level
      ihead(3) = itime(1) * 10000 + itime(2) * 100 + itime(3) ! YYYYMMDD
      ihead(4) = itime(5) * 100   + itime(6)                  ! HHMM
      ihead(5) = nlon  ! 1. dimension
      ihead(6) = nlat  ! 2. dimension
      ihead(7) =    1
      ihead(8) =    0

      select case(kcode)
      case(129,134)
         ihead(2) = 0 
         write (60,'(8I10)') ihead(:)
         write (60,'(8(X,E16.9))') zgp(:,:)
      end select

      close(60)

      return
      end subroutine write_vargp2D

!     ========================
!     SUBROUTINE WRITE_VARGP3D
!     ========================

      subroutine write_vargp3D(zgp,kcode,klev)
      use pumamod

      dimension itime(8)
      dimension ihead(8)

      character(20)  :: ynumber
      character(256) :: yfilename
      real :: zgp(nlon,nlat,klev)


      call date_and_time(values=itime)

!     produce file name to be written
      if (NLAT < 1000) then
         write(yfilename,'("N",I3.3,"_surf_",I4.4,".sra")') NLAT,kcode
      else
         write(yfilename,'("N",I4.4,"_surf_",I4.4,".sra")') NLAT,kcode
      endif



      open(60,file=yfilename,form='formatted')
      ihead(1) = kcode  ! code for reciprocal of damping time scale 
      ihead(2) =   0    ! level
      ihead(3) = itime(1) * 10000 + itime(2) * 100 + itime(3) ! YYYYMMDD
      ihead(4) = itime(5) * 100   + itime(6)                  ! HHMM
      ihead(5) = nlon  ! 1. dimension
      ihead(6) = nlat  ! 2. dimension
      ihead(7) =    1
      ihead(8) =    0

      select case(kcode)
      case(121,122,123,124,125,126)
         do jlev = 1 , klev
            ihead(2) = jlev
            write (60,'(8I10)') ihead(:)
            write (60,'(8(X,E16.9))') zgp(:,:,jlev)
         enddo
      end select

      close(60)

      return
      end subroutine write_vargp3D


!     ====================
!     SUBROUTINE GRIDPOINT
!     ====================

      subroutine gridpoint
      use pumamod
!!
      dimension ihead(8)
      dimension zprof(nlev)
      dimension zzprof(nlon,nlat,nlev)
      character(256) :: yfilename,yfohead,yfodata,ymessage
      logical :: exist

!     Read orography

      gor(:,:) = 0.0
      if (noro > 0) then
         call read_oro_grid(iread)
      endif

!     Read ground anomaly temperature

      call read_gan_grid(iread)

      call sp2fc(sns,gns)
      call sp2fc(sep,gep)
      call sp2fc(stg,gtg)

      call alt2reg(gns,1)
      call alt2reg(gep,1)
      call alt2reg(gtg,1)


      call fc2gp(gns,nlon,nlat)
      call fc2gp(gep,nlon,nlat)
      call fc2gp(gtg,nlon,nlat)

      gtg(:,:) = gtg(:,:) + gan(:,:)

      call modify_orography(gor) ! User interface

      if (nyoden /= 0) then
         call yoden
      else
!        compute ground temperature on orography surface

         gtg(:,:) = gtg(:,:) - ALR * gor(:,:) / GA ! [gpm]

         call modify_ground_temperature(gtg) ! User interface

         call anomaly_factors       ! Compute factors for NS & EP

!        Compute vertical profile for each column

         do jlat = 1 , nlat
         do jlon = 1 , nlon
            call tprofile(gtg(jlon,jlat),zprof,gor(jlon,jlat)/GA)
            gtc(jlon,jlat,:) = zprof(:)
         enddo
         enddo

!        Modify Restoration Temperature with EP mode

         do jlev = 1 , nlev
            gtc(:,:,jlev) = gtc(:,:,jlev) + gaf(:,:,jlev) * gep(:,:)
         enddo

!        Compute vertical profile for each column including variable NS mode

         do jlev = 1 , nlev
            gtv(:,:,jlev) = gaf(:,:,jlev) * gns(:,:)
         enddo

      endif ! (nyoden == 1)

!     Initialize surface pressure (LnPs)

      do jlat = 1 , nlat
      do jlon = 1 , nlon
         gsp(jlon,jlat) = -gor(jlon,jlat) / (GASCON*tgr)
      enddo
      enddo

!     if (nhelsua == 1 .or. nhelsua == 2) then
      if (nhelsua > 0) then
         call heldsuarez
         gor(:,:)   = 0.0
         gsp(:,:)   = 0.0
      endif


      if (nstrato == 1) then
         call setzt2 ! Torben's forcing initialisation
      endif

      call write_oro
      call write_ps
      call write_gtc
      call write_gtv
      call write_vargp3D(gtdamp,123,nlev)
      if (ntestgp == 1) then
         call write_vargp3D(gtc,124,nlev)
         call write_vargp3D(gtv,125,nlev)
         call write_vargp3D(gtdamp,126,nlev)
      endif

      call printprofile

      return
      end

!     =====================
!     SUBROUTINE HELDSUAREZ
!     =====================

      subroutine heldsuarez
      use pumamod

!     Set up the restoration temperature field according to that given
!     in Held & Suarez (1994, Bul. Amer. Meteor. Soc.).
!     Only difference: There is an offset of 1./3. added to the sine of
!     latitude. The reason is that the namelist parameter TGR still
!     represents the global mean of surface restoration temperature.
!
!     ==> Set TGR=295. in order to get exactly the same restoration
!         temperature field as in Held & Suarez (1994).
!
!     ==> DTNS in H&S (1994) is the equator-pole difference while
!         PUMA uses DTNS for pole-pole difference. Therefore we use
!         0.5 * dtns in this subroutine.


!     Produce restoration temperature

      real :: zp0,z3
      real :: zdtc,zdtv,zsip

      zp0 = 100000.
      z3  = 1.0 / 3.0
      if (nhelsua == 1 .or. nhelsua == 2) then
         do jlev=1,nlev
            do jlat=1,nlat
               zsip = sigma(jlev)*PSURF/zp0
               zdtc = dtep * (sid(jlat)**2-z3) + dtzz * log(zsip) * csq(jlat)
               zdtv = zdtc - 0.5 * dtns * sid(jlat)
               gtc(:,jlat,jlev) = max(ttp, (tgr - zdtc) * zsip**AKAP)
               gtv(:,jlat,jlev) = max(ttp, (tgr - zdtv) * zsip**AKAP)
            enddo
         enddo
         gtv(:,:,:) = gtv(:,:,:) - gtc(:,:,:)
      endif

!     Produce reciprocal of damping time scale for T [1/sec]

!     tauta = 1.0 / (TWOPI * tauta)
!     tauts = 1.0 / (TWOPI * tauts)

      rtauta_dim = 1.0 /(tauta*sid_day)
      rtauts_dim = 1.0 /(tauts*sid_day)
      if (nhelsua == 2 .or. nhelsua == 3) then
         do jlev=1,nlev
            if (sigma(jlev) > zsigb) then
               do jlat = 1,nlat
                  gtdamp(1:nlon,jlat,jlev) = rtauta_dim                 &
     &       + (rtauts_dim - rtauta_dim) * ((sigma(jlev) - zsigb) / (1. - zsigb)) &
     &                                          * (1. - sid(jlat)**2)**2
               enddo
            else
               gtdamp(:,:,jlev) = rtauta_dim
            endif
         enddo
      endif

      return
      end

!     =====================
!     SUBROUTINE RESOLUTION
!     =====================

      subroutine resolution
      use pumamod
      logical :: lex
      namelist /res/ nlat, nlev

      nlat = 32
      nlev = 10

      inquire(file=resolution_namelist,exist=lex)
      if (.not. lex) then
         resolution_namelist = trim(resolution_namelist) // "_00"
         inquire(file=resolution_namelist,exist=lex)
      endif

      if (lex) then
         open(14,file=resolution_namelist)
         read(14,res)
         close(14)
      endif

      nlem = nlev - 1
      nlep = nlev + 1
      nlsq = nlev * nlev

      nlon = nlat + nlat ! Longitudes
      nlah = nlat / 2
      nlpp = nlat
      nhpp = nlah
      nhor = nlon * nlpp

      ntru = (nlon - 1) / 3
      ntp1 = ntru + 1
      nzom = ntp1 + ntp1
      nrsp = (ntru + 1) * (ntru + 2)
      ncsp = nrsp / 2
      nspp = nrsp
      nesp = nspp

      return
      end

!     =================
!     SUBROUTINE READNL
!     =================

      subroutine readnl
      use pumamod

!     This namelist must be identical to namelist ppp_nl in puma.f90

       namelist /ppp_nl/ &
         alpha   , alrpv   , alrs    , disp    &
       , dorox   , doroxs  , doroy   , doroys  , dt      , dtep    &
       , dtns    , dtrop   , dttrp   , dtzz    , dvdiff  , edgepv  &
       , epsync  &
       , flsp0   , flsdp   , flsamp  , flsoff  , horo    , kick    &
       , lat1oro , lat1tgr , lat2oro , lat2tgr &
       , lon1oro , lon1tgr , lon2oro , lon2tgr &
       , mpstep  &
       , nafter  , ncoeff  , ncorrect, ncu     , ndel    , ndiag   &
       , ndl     , nextout , newsr   , nfls    , ngui    , nguidbg &
       , nhelsua , nsync    &
       , ntestgp , nkits   &
       , nlevt   , nmonths , noro    , norox   , noutput , noutsrv &
       , npackgp , npacksp , nreverse&
       , nruido  , nrun    , nselect , nspecsel, nsponge , nsrv    &
       , nstep   , nstop   , nstrato , ntgr    , nsym   , ntspd   &
       , nvg     , nwpd    , nwspini , nyears  , nyoden  &
       , oroano  , orofac  &
       , pac     , pmaxpv  , pspon   &
       , radpv   , restim  , rotspd  &
       , sigmah  , sigmax  , sponk   , syncstr &
       , t0      , t0k     , tac     , tauta   , tauts   &
       , tdiss   , tfrc    , tgr     , tgrano  , ttp     &
       , nenergy , nentropy, ndheat

      nselect(:)     = 1
      nspecsel(:)    = 1
      ndl(:)         = 0
      restim(:)      = 15.0
      sigmah(:)      = 0.0
      tfrc(1:nlev)   = (/ (0.0,i=1,nlem), 1.0 /)
      t0k(:)         = 250.0
      t0(:)          = 250.0
      dt(:)          = 0.0

      open(13,file='ppp_namelist')
      read (13,ppp_nl)

!     Use predefined Yoden profile ?

      if (nlev == 20 .and. nyoden > 0) then
!        noro  =   0    ! Don't read orography
!        norox =   2    ! Make idealized orography
!        horo  = 500.0  ! Height of idealized orography
         alrs  = -0.001 ! Stratospheric lapse rate
         nreverse = 0   ! No T-gradient reversal at tropopause
         ncorrect = 0   ! No T-correction due to orography
         idim = min(20,nlev)
         if (nyoden == 1) then
            t0(1:idim) = t0yod1(1:idim)
            dt(1:idim) = dtyod1(1:idim)
         else if (nyoden == 3) then
            t0(1:idim) = t0yod3(1:idim)
            dt(1:idim) = dtyod3(1:idim)
         else if (nyoden == 5) then
            t0(1:idim) = t0yod5(1:idim)
            dt(1:idim) = dtyod5(1:idim)
         else if (nyoden == 7) then
            t0(1:idim) = t0yod7(1:idim)
            dt(1:idim) = dtyod7(1:idim)
         else if (nyoden == 8) then
            t0(1:idim) = t0yod8(1:idim)
            dt(1:idim) = dtyod8(1:idim)
         else if (nyoden == 9) then
            t0(1:idim) = t0yod9(1:idim)
            dt(1:idim) = dtyod9(1:idim)
         endif
      endif
     
      close(13)
      write(*,ppp_nl)

      return
      end

!     =====================
!     * SET VERTICAL GRID *
!     =====================

      subroutine set_vertical_grid

      use pumamod

      if (sigmah(nlev) /= 0.0) return ! Already read in from namelist INP

      if (nvg == 1) then              ! Scinocca & Haynes sigma levels

         if (nlevt >= nlev) then      ! Security check for 'nlevt'
            write (*,*) '*** ERROR *** nlevt >= nlev'
            write (*,*) 'Number of levels (nlev): ',nlev
            write (*,*) 'Number of tropospheric levels (nlevt): ',nlevt
         endif   

!     troposphere: linear spacing in sigma
!     stratosphere: linear spacing in log(sigma)
!     after (see their Appendix):
!     Scinocca, J. F. and P. H. Haynes (1998): Dynamical forcing of
!        stratospheric planetary waves by tropospheric baroclinic eddies.
!        J. Atmos. Sci., 55 (14), 2361-2392

!     Here, zsigtran is set to sigma at dtrop (tropopause height for
!     construction of restoration temperature field). If tgr=288.15K,
!     ALR=0.0065K/km and dtrop=11.km, then zsigtran=0.223 (=0.1 in
!     Scinocca and Haynes (1998)).
!     A smoothing of the transition between linear and logarithmic
!     spacing, as noted in Scinocca and Haynes (1998), is not yet
!     implemented.

         zsigtran = (1. - ALR * dtrop / tgr)**(GA/(GASCON*ALR))
         zsigmin = 1. - (1. - zsigtran) / real(nlevt)

         do jlev=1,nlev
            if (jlev == 1) then
               sigmah(jlev) = 0.000001
               sigmah(jlev) = sigmax
            elseif (jlev > 1 .and. jlev < nlev - nlevt) then
               sigmah(jlev) = exp((log(sigmax) - log(zsigtran))         &
     &             / real(nlev - nlevt - 1) * real(nlev - nlevt - jlev) &
     &             + log(zsigtran))
            elseif (jlev >= nlev - nlevt .and. jlev < nlev - 1) then
               sigmah(jlev) = (zsigtran - zsigmin) / real(nlevt - 1)    &
     &                        * real(nlev - 1 - jlev) + zsigmin
            elseif (jlev == nlev - 1) then
               sigmah(jlev) = zsigmin
            elseif (jlev == nlev) then
               sigmah(jlev) = 1.
            endif
         enddo
         return 
      endif ! (nvg == 1)

      if (nvg == 2) then   ! Polvani & Kushner sigma levels
         inl = int(real(nlev)/(1.0 - sigmax**(1.0/5.0)))
         do jlev=1,nlev
            sigmah(jlev) = (real(jlev + inl - nlev) / real(inl))**5
         enddo
         return
      endif ! (nvg == 2)

!     Default: equidistant sigma levels

      if (nvg == 0) then
         do jlev = 1 , nlev
            sigmah(jlev) = real(jlev) / real(nlev)
         enddo
      endif ! (nvg == 0)

      return
      end

!     =================
!     SUBROUTINE INITPM
!     =================

      subroutine initpm
      use pumamod

      real (kind=8) radea,zakk,zzakk
      radea  = PLARAD_EARTH        ! Planet radius in high precision
      plavor = EZ * rotspd   ! Planetary vorticity

!     *************************************************************
!     * carries out all initialisation of model prior to running. *
!     * major sections identified with comments.                  *
!     * this s/r sets the model parameters and all resolution     *
!     * dependent quantities.                                     *
!     *************************************************************

!     *********************
!     * set vertical grid *
!     *********************

      call set_vertical_grid

      dsigma(1     ) = sigmah(1)
      dsigma(2:nlev) = sigmah(2:nlev) - sigmah(1:NLEM)

      rdsig(:) = 0.5 / dsigma(:)

      sigma(1     ) = 0.5 * sigmah(1)
      sigma(2:nlev) = 0.5 * (sigmah(1:NLEM) + sigmah(2:nlev))

!     annual cycle period and phase in timesteps

      if (tac > 0.0) tac = TWOPI / (ntspd * tac)
      pac = pac * ntspd

!     compute internal diffusion parameter

      jdelh = ndel/2
      if (tdiss > 0.0) then
         zakk = WW*(radea**ndel)/(TWOPI*tdiss*((ntru*(ntru+1.))**jdelh))
      else
         zakk = 0.0
      endif
      zzakk = zakk / (WW*(radea**ndel))

!     set coefficients which depend on wavenumber

      zrsq2 = 1.0 / sqrt(2.0)

      jr=-1
      do jm=0,ntru
         do jn=jm,ntru
            jr=jr+2
            ji=jr+1
            spnorm(jr)=zrsq2
            spnorm(ji)=zrsq2
         enddo
         zrsq2=-zrsq2
      enddo

      return
      end

      subroutine printparameter
      use pumamod
      
      print 8000
      print 8050
      print 8000
      print 8010,nlev
      print 8020,ntru
      print 8030,nlat
      print 8040,nlon
      print 8000
      print 8120
      return
 8000 format(' *****************************************')
 8010 format(' * nlev = ',i6,'   Number of levels      *')
 8020 format(' * ntru = ',i6,'   Triangular truncation *')
 8030 format(' * nlat = ',i6,'   Number of latitudes   *')
 8040 format(' * nlon = ',i6,'   Number of longitues   *')
 8050 format(' *       PPP - Puma Pre Processor        *')
 8120 format(/)
      end


!     ===================
!     SUBROUTINE GPROFILE
!     ===================
      subroutine gprofile(ptgr,prgrad,pgpm)
      use pumamod

!     *************************************************************
!     * Set up the restoration temperature profiles for gradient  *
!     * modes DTNS - mode[0,1] and DTEP - mode[0,2]               *
!     * The lapse rate of ALR K/m is assumed under the tropopause *
!     * and zero above. The tropopause is defined by <dtrop>.     *
!     * The profile is a sine wave with 0 at tropopause sigma and *
!     * 1 at sigma = 1.                                           *
!     ************************************************************* 

      dimension prgrad(nlev)

      ztpheight = dtrop - pgpm      ! Tropopause height over ground
      ztptemp   = tgr - ALR * dtrop ! Tropopause temperature
      ztps = (ztptemp/ptgr)**(GA/(ALR*GASCON)) ! Tropoause sigma

      do jlev = 1 , nlev
         prgrad(jlev) = sin(0.5*PI*(sigma(jlev)-ztps)/(1.0-ztps))
         if (sigma(jlev) < ztps) prgrad(jlev) = 0.0
      enddo

      return
      end

!     ===================
!     SUBROUTINE TPROFILE
!     ===================
      subroutine tprofile(ptgr,prof,pgpm)
      use pumamod

!     *************************************************************
!     * Set up the restoration temperature profile for one column *
!     * The temperature at sigma = 1 is <ptgr>, entered in kelvin *
!     * The lapse rate of ALR K/m is assumed under the tropopause *
!     * and zero above. The tropopause is defined by <ztpheight>. *
!     * The smoothing ot the tropopause depends on <dttrp>.       *
!     ************************************************************* 

      dimension prof(nlev)  ! Resulting temperature profile [K]

      zsigprev  = 1.0               ! sigma value
      ztprev    = ptgr              ! Temperature [K]
      zzprev    = 0.0               ! Height      [m]
      ztpheight = dtrop - pgpm      ! Tropopause height over ground
      ztptemp   = tgr - ALR * dtrop ! Tropopause temperature
      zalr      = (ptgr - ztptemp) / ztpheight

      do jlev = nlev , 1 , -1   ! from bottom to top of atmosphere
         zlogsig = GASCON / GA * log(zsigprev / sigma(jlev))
         zzp     = zzprev + ztprev * zlogsig
         ztp=ztptemp+sqrt((.5*zalr*(zzp-ztpheight))**2+dttrp**2)
         ztp=ztp-.5*zalr*(zzp-ztpheight)
         ztpm=.5*(ztprev+ztp)

         zzpp = zzprev + ztpm * zlogsig
         ztpp=ztptemp+sqrt((.5*zalr*(zzpp-ztpheight))**2+dttrp**2)
         ztpp=ztpp-.5*zalr*(zzpp-ztpheight)

         prof(jlev)=ztpp
         zzprev=zzprev + 0.5 * (ztpp+ztprev) * zlogsig
         ztprev=ztpp
         zsigprev=sigma(jlev)
      enddo
      return
      end

!     ======================
!     SUBROUTINE ppp_write_i
!     ======================

      subroutine ppp_write_i(yvarname,nvals,ivals)
      use pumamod
  
      character(*)   :: yvarname
      integer        :: nvals
      integer        :: ivals(nvals)

      if (mypid == NROOT) then 
         write(95,'("[",A,"]")') trim(yvarname)
         write(95,'(I4)') nvals
         write(95,'(I6)') ivals(:)
      endif
      return
      end

!     ======================
!     SUBROUTINE ppp_write_r 
!     ======================

      subroutine ppp_write_r(yvarname,nvals,pvals)
      use pumamod

      character(*)   :: yvarname
      integer        :: nvals
      real           :: pvals(nvals)

      if (mypid == NROOT) then
         write(95,'("[",A,"]")') trim(yvarname)
         write(95,'(I4)') nvals
         write(95,'(E14.7)') pvals(:)
      endif
      return
      end

!     ================
!     SUBROUTINE yoden
!     ================

      subroutine yoden
      use pumamod
!
      do jlev=1,nlev
       do jlat=1,nlat/2
        zlat=dla(jlat)
        ztres=t0(jlev)+dt(jlev)/2.*(cos(2.*zlat)-1./3.)
        do jlon=1,nlon
         gtc(jlon,jlat,jlev)=ztres
        enddo
       enddo
      enddo
!
      do jlat=1,nlat/2
       j2=nlat+1-jlat
       gtc(:,j2,:)=gtc(:,jlat,:)
      enddo

      write (*,*) ' Computed Yoden profile',nyoden
!
      return
      end
!
!     =======================
!     SUBROUTINE PRINTPROFILE
!     =======================

      subroutine printprofile
      use pumamod

!     **********************************
!     * write out vertical information *
!     **********************************

      dimension ztr(nlev+1)

      ztr(nlev+1) = tgr

      write(*,9001)
      write(*,9002)
      write(*,9003)
      write(*,9002)

      do jlev=1,nlev
         ztr(jlev) = sum(gtc(:,:,jlev)) / (nlon * nlat)
      enddo

      do jlev=1,nlev
      write(*,9004) jlev,sigma(jlev),ztr(jlev),ztr(jlev+1)-ztr(jlev),gra(jlev)
      enddo

      write(*,9002)
      write(*,9001)
      return
 9001 format(/)
 9002 format(1x,45('*'))
 9003 format(' * Lv *     Sigma Restor-T  Delta-T    Vfact *')
 9004 format(' *',i3,' * ',4f9.3,' *')
      end

!     =====================
!     * SUBROUTINE LEGPRI *
!     =====================

      subroutine legpri
      use pumamod

      write (*,231)
      write (*,232)
      write (*,233)
      write (*,232)
      do 14 jlat = 1 , nlat
         zalat = dla(jlat)*180.0/PI
         write (*,234) jlat,zalat,csq(jlat),gwd(jlat)
   14 continue
      write (*,232)
      write (*,231)
      return
  231 format(/)
  232 format(1x,36('*'))
  233 format(' * No *   Lat *       csq    weight *')
  234 format(' *',i3,' *',f6.1,' *',2f10.4,' *')
      end


!     =================
!     SUBROUTINE INILAT
!     =================

      subroutine inilat
      use pumamod
      character(1) :: ch

      ch = 'N'
      do jlat = 1 , nlat
         csq(jlat) = 1.0 - sid(jlat) * sid(jlat)
         dla(jlat) = asin(sid(jlat))
      enddo
      do jlat = 1 , nlat/2
         ideg = nint(180.0/PI * asin(sid(jlat)))
         write(chlat(jlat),'(i2,a1)') ideg,'N'
         write(chlat(nlat+1-jlat),'(i2,a1)') ideg,'S'
      enddo
      return
      end


!     =================
!     SUBROUTINE SETZT2
!     =================
      subroutine setzt2
      use pumamod

!                             US standard atmosphere (1976):
      parameter(INL = 7)      ! number of defined layers
      dimension zzus(0:INL)   ! height of interfaces between layers
      dimension zlus(INL)     ! temperature lapse rates of layers
      dimension zpus(0:INL)   ! pressure at interfaces between layers
      dimension ztus(0:INL)   ! temperature at interfaces between layers

      dimension ztrs(nlev)  ! Mean profile
      dimension ztpv(nlev)  ! Vertical profile of stratospheric polar
!                             vortex forcing
      dimension zdtep(nlat)
      dimension zdtns(nlat)
      dimension zff(nlev)
      dimension zfw(nlat,nlev)
      dimension zfsph(nlat,nlev)
      dimension zqc(nlat,nlev)
      dimension zphi(nlat)
      dimension zgr1(nlon,nlat,nlev)
      dimension zgr2(nlon,nlat,nlev)
      dimension zfls(nlat,nlev)

      real :: zp
      real :: zsigtp
      real :: pref
      real :: zpmaxsph

      sr1(:,:) = 0.0 ! NESP,nlev
      sr2(:,:) = 0.0 ! NESP,nlev

!     1. Mean vertical profile (MVP), approx. US standard atmosphere

      zzus(0) = 0.
      zzus(1) = dtrop     ! US standard atmosphere: zzus(1) = 11000.
      zzus(2) = 20000.
      zzus(3) = 32000.
      zzus(4) = 47000.
      zzus(5) = 51000.
      zzus(6) = 71000.
      zzus(7) = 84852.
      zlus(1) = ALR      ! US standard atmosphere: zlus(1) = 0.0065
      zlus(2) = 0.0
      zlus(3) = -0.001
      zlus(4) = -0.0028
      zlus(5) = 0.0
      zlus(6) = 0.0028
      zlus(7) = 0.002

! calculation of pressure and temperature at layer interfaces

      zpus(0) = PSURF ! US standard atmosphere: zpus(0) = 1013.25 hPa
      ztus(0) = tgr   ! US standard atmosphere: ztus(0) = 288.15 K

      do ji=1,INL
         ztus(ji) = ztus(ji-1) - zlus(ji) * (zzus(ji) - zzus(ji-1))
         if (zlus(ji) == 0.) then
            zpus(ji) = zpus(ji-1) * exp(-GA * (zzus(ji) - zzus(ji-1))   &
     &                  / (GASCON * ztus(ji-1)))
         else
            zpus(ji) = zpus(ji-1)                                       &
     &                 * (ztus(ji) / ztus(ji-1))**(GA/(GASCON*zlus(ji)))
         endif
      enddo

! calculation of temperature on given sigma full levels, sigma(1:nlev)
      do jlev=nlev,1,-1
         zp = sigma(jlev)*PSURF
         if (zp <= zpus(0) .and. zp > zpus(1)) then
            ztrs(jlev) = ztus(0) * (zp / zpus(0))**(GASCON*zlus(1)/GA)
         elseif (zp <= zpus(1) .and. zp > zpus(2)) then
            ztrs(jlev) = ztus(1) * (zp / zpus(1))**(GASCON*zlus(2)/GA)
         elseif (zp <= zpus(2) .and. zp > zpus(3)) then
            ztrs(jlev) = ztus(2) * (zp / zpus(2))**(GASCON*zlus(3)/GA)
         elseif (zp <= zpus(3) .and. zp > zpus(4)) then
            ztrs(jlev) = ztus(3) * (zp / zpus(3))**(GASCON*zlus(4)/GA)
         elseif (zp <= zpus(4) .and. zp > zpus(5)) then
            ztrs(jlev) = ztus(4) * (zp / zpus(4))**(GASCON*zlus(5)/GA)
         elseif (zp <= zpus(5) .and. zp > zpus(6)) then
            ztrs(jlev) = ztus(5) * (zp / zpus(5))**(GASCON*zlus(6)/GA)
         elseif (zp <= zpus(6) .and. zp > zpus(7)) then
            ztrs(jlev) = ztus(6) * (zp / zpus(6))**(GASCON*zlus(7)/GA)
         else
            ztrs(jlev) = ztus(7)
         endif
      enddo

!     2. Symmetric equator-pole forcing mode (DTEP) and
!     3. Asymmetric Npole-Spole forcing mode (DTNS)

!     sid(nlat) is sine of latitude, taking into account the nonequally
!     spaced Gaussian latitudes.
      do jlat=1,nlat
         zdtep(jlat) = -dtep * (sid(jlat)**2 - 1./3.)
         zdtns(jlat) =  dtns * sid(jlat) / 2.
      enddo

!     4. Factor modulating the DTEP and DTNS modes (f)

      zsigtp = zpus(1)/zpus(0)   ! sigma at tropopause
      zff(:) = 0.
      do jlev=1,nlev
         if (sigma(jlev) > zsigtp) then
            zff(jlev) = sin(0.5*PI * (sigma(jlev) - zsigtp)             &
     &                                                  / (1. - zsigtp))
         endif
      enddo

!     5. Vertical profile of stratospheric polar vortex forcing

      ztpv(:) = 0.
      if (pmaxpv == 0.) then
         do jlev=1,nlev
            if (sigma(jlev) <= zsigtp) then
               ztpv(jlev) = ztus(1) * (sigma(jlev)*PSURF / zpus(1))     &
     &                                               **(GASCON*alrpv/GA)
            endif
         enddo
      elseif (pmaxpv > 0.) then
         do jlev=1,nlev
            if (sigma(jlev) <= pmaxpv/PSURF) then
               ztpv(jlev) = ztus(1) * (sigma(jlev)*PSURF / pmaxpv)      &
     &                                               **(GASCON*alrpv/GA)
            else
               ztpv(jlev) = ztus(1)
            endif
         enddo
      endif

!     6. Factor confining the stratosph. polar vortex to high latitudes

      zphi(:) = dla(:) * 180. / PI
      zfw(:,:) = 0.
      if (edgepv > 0.) then
         do jlev=1,nlev
            if (sigma(jlev) <= pmaxpv/PSURF) then
               do jlat=1,nlat
                  zfw(jlat,jlev) =                                      &
     &                  0.5 * (1. - tanh((radpv - zphi(jlat)) / edgepv))
               enddo
            endif
         enddo
      endif

!     7. Lower stratospheric forcing

      zfls(:,:) = 0.
      if (nfls == 1) then
         do jlev=1,nlev
            zp =sigma(jlev) * PSURF
            do jlat=1,nlat
               if (zp > flsp0-flsdp .and. zp < flsp0+flsdp) then
                  zfls(jlat,jlev) = cos(0.5 * PI * (zp - flsp0) / flsdp)&
     &            * (flsoff - flsamp * cos(2. * zphi(jlat) * PI / 180.))
               endif
            enddo
         enddo
      endif

!     construct restoration temperature field

      do jlev=1,nlev
         do jlat=1,nlat
            zgr1(:,jlat,jlev) = ((1. - zfw(jlat,jlev)) * ztrs(jlev)     &
     &         + zfw(jlat,jlev) * ztpv(jlev) + zff(jlev) * zdtep(jlat)  &
     &     + (1. - zfw(jlat,jlev)) * zfls(jlat,jlev) - t0k(jlev)) / CT
            zgr2(:,jlat,jlev) = (zff(jlev) * zdtns(jlat)) / CT
         enddo
      enddo

      do jlev = 1 , nlev
         gtc(:,:,jlev) = t0k(jlev) + CT * zgr1(:,:,jlev)
         gtv(:,:,jlev) =             CT * zgr2(:,:,jlev)
      enddo

!     ---------- test output to control T_r field ----------
      open(112,file='tr_test.srv',form='unformatted')
      do jlev=1,nlev
         ip = int(sigma(jlev) * PSURF * 1000.)
         write(112) 121, ip, 0000, 00, nlon, nlat, 0, 0
         write(112) t0k(jlev) + (zgr1(:,:,jlev) + zgr2(:,:,jlev)) * CT
      enddo
      close(112)
!     ---------- test output to control T_r field ----------

      print *,'**************************************************'
      print *,'* Restoration Temperature set up for aqua planet *'
      print *,'* including stratosphere and polar vortex        *'
      print *,'**************************************************'
      return
      end
