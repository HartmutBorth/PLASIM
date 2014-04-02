      module pumamod

!     ****************************************************************
!     *                    Shallow water model  SAM                  *
!     ****************************************************************
!     *                       Thomas Frisius                         *
!     *                   Meteorologisches Institut                  *
!     *                     KlimaCampus Hamburg                      *
!     * 10-September-2009        GERMANY                             *
!     ****************************************************************
!
!                              based on the                    
!
!     ****************************************************************
!     *        Portable University Model of the Atmosphere           *
!     ****************************************************************
!
!                                    by
!
!     ****************************************************************
!     *          Frank Lunkeit & Edilbert Kirk & Torben Kunz         *
!     *                   Meteorologisches Institut                  *
!     *                      Universitaet Hamburg                    *
!     *                        Bundesstrasse 55                      *
!     *                          20146 HAMBURG                       *
!     * 27-Feb-2007                 GERMANY                          *
!     ****************************************************************

      integer :: npro = 1                  ! Number of processes (CPUs)

!     ****************************************************************
!     * Set the vertical resolution, needed by GUI                   *
!     ****************************************************************

      parameter(NLEV = 1)                  ! Number of levels

!     ****************************************************************
!     * Don't touch the following parameter definitions !            *
!     ****************************************************************

      integer, parameter :: PUMA   = 0     ! Model ID
      integer, parameter :: SAM    = 1     ! Model ID
      integer, parameter :: PLASIM = 2     ! Model ID
      
      integer :: nlat = 32
      integer :: nhpp = 16
      integer :: NLON = 64                 ! Number of longitudes
      integer :: NTRU = 21                 ! Triangular truncation
      integer :: NLPP = 32                 ! Latitudes per process
      integer :: NHOR = 2048               ! Horizontal part
      integer :: NUGP = 2048
      integer :: NPGP = 1024               ! Number of packed words
      integer :: NTP1 = 22                 ! Truncation + 1
      integer :: NRSP = 506                ! No of real    global modes
      integer :: NCSP = 253                ! No of complex global modes
      integer :: NSPP = 506                ! Modes per process
      integer :: NESP = 506                ! Dim of spectral fields
      integer :: NZOM = 44                 ! Number of zonal wavenumbers

      integer :: nud  = 6                  ! I/O unit for diagnostic output
      parameter(NROOT = 0)                 ! Master node

      parameter(EZ     = 1.632993161855452D0) ! ez = 1 / sqrt(3/8)
      parameter(GA     = 9.81)                ! Gravity
!      parameter(GA     = 22.9)               ! Gravity Jupiter
      parameter(PI     = 3.141592653589793D0) ! Pi
      parameter(TWOPI  = PI + PI)             ! 2 Pi
      parameter(PLARAD = 6371000.0)           ! Planet radius
!       parameter(PLARAD = 70000000.0)        ! Planet radius Jupiter
!       parameter(PLARAD = 71500000.0)        ! Planet radius Jupiter2
!      parameter(PLARAD = 0.495)              ! Planet radius tank
      parameter(PNU    = 0.02)                ! Time filter
      parameter(PNU21  = 1.0 - 2.0*PNU)       ! Time filter 2
      parameter(WW     = 0.00007292)          ! Rotation speed [1/sec]
!     parameter(WW     = 0.000176)            ! Rotation speed Jupiter [1/sec]
!      parameter(WW     = 2.5)                ! Rotation speed tank Experiment A
!     parameter(WW     = 1.44)                ! Rotation speed tank Experiment B
      parameter(CV     = PLARAD * WW)         ! cv
      parameter(RHO0   = 1.)                  ! Density of the shallow water
!      parameter(RHO0   = 1000.)              ! Density of the shallow water (Tank)
!      parameter(HS     = 5960.)              ! Mean Height of shallow water
      parameter(HS     = 10600.)              ! Mean Height of shallow water
!      parameter(HS     = 5000.)              ! Mean Height of shallow water
!      parameter(HS     = 20000.)             ! Mean Height of shallow water (Jupiter)
!      parameter(HS     = 0.01)               !  Mean Height of shallow water (Tank Experiment A))
!      parameter(HS     = 0.04)               ! Mean Height of shallow water (Tank Experiment B))
      parameter(PSURF  = GA*RHO0*HS)          ! Mean Surface pressure [Pa]
      parameter(FR     = CV**2/GA/HS)         ! Froude-Number 
      parameter(F0     = 1.e-4)               ! Mean Coriolis parameter (lequiv=.true.)
      parameter(FREQ   = CV**2*F0**2/WW**2/GA/HS)      ! Froude-Number for equivalent barotropic model (lequiv=.true.)

!     *******************
!     * Global Logicals *
!     *******************

      logical :: lrestart =  .false. ! Existing "sam_restart" sets to .true.
      logical :: lselect  =  .false. ! true: disable some zonal waves
      logical :: lspecsel =  .false. ! true: disable some spectral modes
      logical :: llid     =  .false. ! true: A rigid lid is introduced
      logical :: lequiv   =  .false. ! true: Equivalent barotropic model
      logical :: lbal     =  .false. ! true: The initial state is balanced
      
!     **************************
!     * Global Integer Scalars *
!     **************************

      integer :: model    = SAM
 
      integer :: nexp     =  4  ! Index for initial state
      integer :: kick     =  1 ! kick > 0 initializes eddy generation
      integer :: nafter   = 12 ! write data interval
      integer :: ncoeff   =  0 ! number of modes to print
      integer :: ndel     =  6 ! ndel
      integer :: ndiag    = 12 ! write diagnostics interval
      integer :: ngui     =  0 ! activate Graphical User Interface
      integer :: nguidbg  =  0 ! activate GUI debug mode
      integer :: nkits    =  3 ! number of initial timesteps
      integer :: noutput  =  1 ! global switch for output on (1) or off (0)
      integer :: nwspini  =  1 ! write sp_init after initialization
      integer :: nrun     =  0 ! if (nstop == 0) nstop = nstep + nrun
      integer :: nseedlen =  0 ! length of random seed
      integer :: nstep    =  0 ! current timestep
      integer :: nstep1   =  0 ! 1st. timestep
      integer :: nstop    =  0 ! finishing timestep
      integer :: ntspd    =144 ! number of timesteps per day
      integer :: ncu      =  0 ! check unit (debug output)
      integer :: nwrioro  =  1 ! controls output of orography
      integer :: nruido   =  0 ! 0: no noise added, 1: temporal (global
      integer :: nf       =  0 ! Forcing wavenumber (Nozawa and Yoden 1997)
      integer :: nyears   =  0 ! use instead of nrun for simulation time
      integer :: mmax     = 256  ! Maximum zonal waven. for climat. forcing
      integer :: nmax     = 256 ! Maximum total waven. for climat. forcing
!
!     THIS IS A CHANGE BY TF
!
      integer  :: irec   = 1  ! Index for the GRADS-output file

!
!     END CHANGE BY TF
!
!
!     ***********************
!     * Global Real Scalars *
!     ***********************

      real :: delt             ! 2 pi / ntspd timestep interval
      real :: delt2            ! 2 * delt
      real :: plavor=    EZ    ! planetary vorticity
      real :: psmean= PSURF    ! Mean of Ps 
      real :: rotspd=     1.0  ! rotation speed 1.0 = normal Earth rotation
      real :: tdiss =     0.25 ! diffusion time scale [days]
      real :: ah    =     0.   ! horizontal momentum exchange coefficient
!       real :: ah    =     1.e-6  ! horizontal momentum exchange coefficient (Tank! Viscosity of water)
      real :: tac   =       0. ! length of annual cycle [days] (0 = no cycle)
      real :: pac   =       0. ! phase of the annual cycle [days]
      real :: disp  =     0.0  ! noise dispersion
      real :: fnorm0 =    0.0  ! RMS-Amplitude of spectr. stoch. Forcing
                                ! (Nozawa and Yoden 1997) 
      real :: vortamp =   0.0  !Vortex amplitude for tank experiment
      real :: vortrad =   0.0  !Vortex radius for tank experiment

!     real :: pldm(NLPP,NCSP)                  ! pldp * m

!     **************************
!     * Global Spectral Arrays *
!     **************************

      real, allocatable :: sd(:)        ! Spectral Divergence
      real, allocatable :: sz(:)        ! Spectral Vorticity
      real, allocatable :: sp(:)        ! Spectral Pressure 
      real, allocatable :: sdo(:)      ! Spectral Divergence old
      real, allocatable :: szo(:)      ! Spectral Vorticity  old
      real, allocatable :: spo(:)      ! Spectral Pressure   old
      real, allocatable :: sp2(:)       ! Spectral Pressure at t-1
!                                    for extended output only
      real, allocatable :: sp3(:)       ! Spectral Pressure at t-2
!                                    for extended output only
      real, allocatable :: so(:)        ! Spectral Orography

      real, allocatable :: sdp(:)       ! Spectral Divergence  Partial
      real, allocatable :: szp(:)       ! Spectral Vorticity   Partial
      real, allocatable :: spp(:)       ! Spectral Pressure    Partial
      real, allocatable :: sop(:)       ! Spectral Orography   Partial

      real, allocatable :: sdt(:)       ! Spectral Divergence  Tendency
      real, allocatable :: szt(:)       ! Spectral Vorticity   Tendency
      real, allocatable :: spt(:)       ! Spectral Pressure    Tendency

      real, allocatable :: sdm(:)       ! Spectral Divergence  Minus
      real, allocatable :: szm(:)       ! Spectral Vorticity   Minus
      real, allocatable :: spm(:)       ! Spectral Pressure    Minus

      real, allocatable :: su(:)            ! Spectral angular momentum
      real, allocatable :: svort(:)        ! Spectral absolute vorticity for case lequiv=.true.
      real, allocatable :: svortp(:)       ! Spectral absolute vorticity for case lequiv=.true.
      real, allocatable :: senst(:)         ! EnstrophySpectral angular momentum
      real, allocatable :: spvo(:)         ! Potential vorticity
      real, allocatable :: sak(:)          ! Hyper diffusion
      real, allocatable :: salpha(:)       ! Restoration rate for climatological forcing
      real, allocatable :: srcn(:)         ! 1.0 / (n * (n+1))
      real, allocatable :: span(:)         ! Pressure for diagnostics
      real, allocatable :: spnorm(:)       ! Factors for output normalization
      real, allocatable :: ske(:)          ! Kinetic energy
      real, allocatable :: spsi(:)         ! Streamfunction
      real, allocatable :: spsip(:)        ! Streamfunction partial
      real, allocatable :: schi(:)         ! Velocity potential 
      real, allocatable :: schip(:)        ! Velocity potential 
      real, allocatable :: szr(:)          ! Vorticity of the climatological forcing
      real, allocatable :: szrp(:)         ! Vorticity of the climatological forcing
      real, allocatable :: spsir(:)        ! Streamfunction of the climatological forcing
      real, allocatable :: spsirp(:)       ! Streamfunction of the climatological forcing
      real, allocatable :: szforce(:)      ! Random Vorticity forcing
      real, allocatable :: szforcenew(:)   ! Random Vorticity forcing
      real, allocatable :: szforcenewg(:)  ! Random Vorticity forcing
      real, allocatable :: sztau(:)        ! Curl of Windstress
      real, allocatable :: sdtau(:)        ! Divergence of Windstress

      real :: sr1(5,2)        ! Dummy declaration for GUI compatibility
      real :: sr2(5,2)        ! Dummy declaration for GUI compatibility
      integer, allocatable  :: nindex(:)   ! Holds wavenumber
      integer, allocatable  :: nscatsp(:)  ! Used for reduce_scatter op
      integer, allocatable  :: nselect(:)  ! Enable/disable selected zonal waves
      integer, allocatable  :: nspecsel(:) ! Enable/disable slected spectral modes
      integer, allocatable  :: iszforce(:) ! Enable/disable selected random Vorticity forcing
      real, allocatable     :: riszforce(:) ! Enable/disable selected random Vorticity forcing
      real, allocatable     :: riszforcep(:) ! Enable/disable selected random Vorticity forcing, partial

!     ***************************
!     * Global Gridpoint Arrays *
!     ***************************

      real, allocatable    :: gd(:)         ! Divergence
      real, allocatable    :: gz(:)         ! Vorticity
      real, allocatable    :: gu(:)         ! Zonal velocity
      real, allocatable    :: gv(:)         ! MeridiSpectral Divergence
      real, allocatable    :: gp(:)         ! Pressure 
      real, allocatable    :: gum(:)        ! Zonal velocity old
      real, allocatable    :: gvm(:)        ! Meridional velocity old
      real, allocatable    :: go(:)         ! Orography
      real, allocatable    :: gop(:)        ! Orography, partial
      real, allocatable    :: gland(:)      ! Integer land mask 
      real, allocatable    :: gzforce(:)    ! Vorticity forcing
      integer, allocatable :: iland(:)      ! Land mask
      real, allocatable    :: gpsi(:)       ! Streamfunction
      real, allocatable    :: gchi(:)       ! Velocity potential 
      real, allocatable    :: gpvo(:)       ! Potential vorticity 
      real, allocatable    :: genst(:)      ! Enstrophy
      real, allocatable    :: gtaux(:)      ! Zonal wind stress
      real, allocatable    :: gtauy(:)      ! Meridional wind stress
      real*4, allocatable  :: outgrid2(:)   !
      real, allocatable    :: outgrid(:)    !
      real, allocatable    :: outgridp(:)   !
      real, target,  allocatable  :: gflam(:)  ! Zonal vector component times cos(phi)
      real, target,  allocatable  :: gfphi(:)  ! Meridional vector component times cos(phi)
      real, target,  allocatable  :: gtn(:)    ! Dummy
      real, target,  allocatable  :: gfu(:)    ! Term Fu in Primitive Equations
      real, target,  allocatable  :: gfv(:)    ! Term Fv in Primitive Equations
      real, target,  allocatable  :: gpu(:)    ! Term u * p
      real, target,  allocatable  :: gpv(:)    ! Term v * p
      real, target,  allocatable  :: gke(:)    ! Kinetic energy u * u + v * v
      real, allocatable :: rcsq(:)          ! 1 / cos2(phi)
      real, allocatable :: ruidop(:)        ! noise partial
      real, allocatable :: ruido(:)         ! noise

!     *********************
!     * Diagnostic Arrays *
!     *********************

      integer :: ndl = 0

      real, allocatable :: csu(:)
      real, allocatable :: csv(:)

!     *******************
!     * Latitude Arrays *
!     *******************

      character (3), allocatable  :: chlat(:)
      real, allocatable           :: csq(:)   ! cos2(phi)
      real(kind=8), allocatable   :: sid(:)   ! sin(phi)
      real, allocatable           :: sia(:)   ! sin(phi) alternating
      real(kind=8), allocatable   :: gwd(:)   ! Gaussian weight (phi)
      real, allocatable           :: rcs(:)   ! 1/cos(phi)

!     ****************
!     * Level Arrays *
!     ****************

      real :: tfrc   =   0.0  ! tau F [days]
      real :: fric   =   0.0  ! 1.0 / (2 Pi * tfrc  )
      real :: restim =   0.*10.  ! timescale for climatological forcing
                              ! in days
      real :: sigma(2) = 0.0  ! dummy array for GUI compatibility

!     ******************
!     * Parallel Stuff *
!     ******************

      integer :: myworld = 0
      integer :: mpinfo  = 0
      integer :: mypid   = 0
      real    :: tmstart = 0.0
      real    :: tmstop  = 0.0
      character(80), allocatable :: ympname(:)

!     ******************************************
!     * GUI (Graphical User Interface for X11) *
!     ******************************************

      parameter (NPARCS = 10)         ! Number of GUI parameters
      integer :: nshutdown = 0        ! Flag for shutdown request
      integer :: ndatim(6) = -1       ! Date & time array
      real(kind=4) :: parc(NPARCS)    ! Values of GUI parameters
      real(kind=4) :: crap(NPARCS)    ! Backup of parc(NPARCS)
      logical :: ldtep   = .FALSE.    ! DTEP changed by GUI
      logical :: ldtns   = .FALSE.    ! DTNS changed by GUI
      character(len=32) :: yplanet = "Earth"
      integer, allocatable :: nr(:)

      ! ***************
      ! * Random seed *
      ! ***************
      
      integer, allocatable :: meed(:)     ! machine dependent seed

      end module pumamod

      module radmod       ! Dummy declaration for compatibility
      use pumamod         ! with PLASIM (needed in guimod)
      end module radmod

      program puma_main

!     ***********
!     * History *
!     ***********

!     1972 - W. Bourke:
!            An efficient one-level primitive equation spectral model
!            Mon. Weath. Rev., 100, pp. 683-689

!     1975 - B.J. Hoskins and A.J. Simmons: 
!            A multi-layer spectral model and the semi-implicit method
!            Qart. J. R. Met. Soc., 101, pp. 637-655

!     1993 - I.N. James and J.P. Dodd:
!            A Simplified Global Circulation Model
!            Users' Manual, Dept. of Meteorology, University of Reading

!     1998 - Klaus Fraedrich, Edilbert Kirk, Frank Lunkeit
!            Portable University Model of the Atmosphere
!            DKRZ Technical Report No. 16

!
!     2007 - Back to the roots: SAM - An efficient one-level primitive 
!                                     equation spectral model. 
!
!


!     ******************
!     * Recent Changes *
!     ******************

use pumamod
      call mpstart
      call opendiag
      call read_resolution
      call resolution
      call prolog
      call master
      call epilog
      call guistop
      call mpstop
      stop
      end program puma_main


! ***********************
! * SUBROUTINE OPENDIAG *
! ***********************

subroutine opendiag
use pumamod

if (mypid == NROOT) then 
   open(nud,file='sam_diag')
endif

return
end


! ***************************
! * SUBROUTINE CLEAR_ARRAYS *
! ***************************

subroutine clear_arrays
use pumamod

allocate(sd(nesp))   ; sd(:)    = 0.0 ! Spectral Divergence
allocate(sz(nesp))   ; sz(:)    = 0.0 ! Spectral Vorticity
allocate(sp(nesp))   ; sp(:)    = 0.0 ! Spectral Pressure (ln Ps)
allocate(sdo(nesp))  ; sdo(:)   = 0.0 ! Spectral Divergence old
allocate(szo(nesp))  ; szo(:)   = 0.0 ! Spectral Vorticity old
allocate(spo(nesp))  ; spo(:)   = 0.0 ! Spectral Pressure old
allocate(sp2(nesp))  ; sp2(:)   = 0.0 ! Spectral Pressure (ln Ps) at t-1
allocate(sp3(nesp))  ; sp3(:)   = 0.0 ! Spectral Pressure (ln Ps) at t-2
allocate(so(nesp))   ; so(:)    = 0.0 ! Spectral Orography
allocate(sdp(nspp))  ; sdp(:)   = 0.0 ! Spectral Divergence  Partial
allocate(szp(nspp))  ; szp(:)   = 0.0 ! Spectral Vorticity   Partial
allocate(spp(nspp))  ; spp(:)   = 0.0 ! Spectral Pressure    Partial
allocate(sop(nspp))  ; sop(:)   = 0.0 ! Spectral Orography   Partial
allocate(sdt(nspp))  ; sdt(:)   = 0.0 ! Spectral Divergence  Tendency
allocate(szt(nspp))  ; szt(:)   = 0.0 ! Spectral Vorticity   Tendency
allocate(spt(nspp))  ; spt(:)   = 0.0 ! Spectral Pressure    Tendency
allocate(sdm(nspp))  ; sdm(:)   = 0.0 ! Spectral Divergence  Minus
allocate(szm(nspp))  ; szm(:)   = 0.0 ! Spectral Vorticity   Minus
allocate(spm(nspp))  ; spm(:)   = 0.0 ! Spectral Pressure    Minus


allocate(su(NESP))         ; su(:)      = 0.0  ! Spectral angular momentum
allocate(svort(NESP))      ; svort(:)   = 0.0  ! Spectral absolute vorticity for case lequiv=.true.
allocate(svortp(NSPP))     ; svortp(:)  = 0.0  ! Spectral absolute vorticity for case lequiv=.true.
allocate(senst(NESP))      ; senst(:)   = 0.0   ! EnstrophySpectral angular momentum
allocate(spvo(NESP))       ; spvo(:)    = 0.0 ! Potential vorticity
allocate(sak(NESP))        ; sak(:)     = 0.0 ! Hyper diffusion
allocate(salpha(NESP))     ; salpha(:)  = 0.0 ! Restoration rate for climatological forcing
allocate(srcn(NESP))       ; srcn(:)    = 0.0 ! 1.0 / (n * (n+1))
allocate(span(NESP))       ; span(:)    = 0.0 ! Pressure for diagnostics
allocate(spnorm(NESP))     ; spnorm(:)  = 0.0 ! Factors for output normalization
allocate(ske(NESP))        ; ske(:)     = 0.0 ! Kinetic energy
allocate(spsi(NESP))       ; spsi(:)    = 0.0 ! Streamfunction
allocate(spsip(NSPP))      ; spsip(:)   = 0.0 ! Streamfunction partial
allocate(schi(NESP))       ; schi(:)    = 0.0 ! Velocity potential 
allocate(schip(NSPP))      ; schip(:)   = 0.0 ! Velocity potential 
allocate(szr(NESP))        ; szr(:)     = 0.0 ! Vorticity of the climatological forcing
allocate(szrp(NSPP))       ; szrp(:)    = 0.0 ! Vorticity of the climatological forcing
allocate(spsir(NESP))      ; spsir(:)    = 0.0 ! Streamfunction of climatological forcing
allocate(spsirp(NSPP))     ; spsirp(:)    = 0.0 ! Streamfunction of climatological forcing
allocate(szforce(NSPP))    ; szforce(:) = 0.0 ! Random Vorticity forcing
allocate(szforcenew(NSPP)) ; szforcenew(:) = 0.0  ! Random Vorticity forcing
allocate(szforcenewg(NESP)); szforcenewg(:) = 0.0  ! Random Vorticity forcing
allocate(sztau(NSPP))      ; sztau(:)      = 0.0 ! Curl of Windstress
allocate(sdtau(NSPP))      ; sdtau(:)      = 0.0 ! Divergence of Windstress

allocate(nindex(nesp))    ; nindex(:)    = ntru ! Holds wavenumber
allocate(nscatsp(npro))   ; nscatsp(:)   = nspp ! Used for reduce_scatter op
allocate(nselect(0:ntru)) ; nselect(:)   =    1 ! Enable selected zonal waves
allocate(nspecsel(NCSP))  ; nspecsel(:)  =    1 ! Enable selected zonal waves
allocate(iszforce(NESP))  ; iszforce(:)  =    0 !Enable/disable selected random Vorticity forcing
allocate(riszforce(NESP)) ; riszforce(:) =   0. !Enable/disable selected random Vorticity forcing
allocate(riszforcep(NSPP)); riszforcep(:) =  0. !Enable/disable selected random Vorticity forcing
allocate(iland(NHOR))     ; iland(:)     = 0    ! Land mask

allocate(gd(nhor))    ; gd(:)    = 0.0 ! Divergence
allocate(gz(nhor))    ; gz(:)    = 0.0 ! Vorticity
allocate(gu(nhor))    ; gu(:)    = 0.0 ! u * cos(phi)
allocate(gv(nhor))    ; gv(:)    = 0.0 ! v * sin(phi)
allocate(gum(nhor))   ; gum(:)   = 0.0 ! u * cos(phi)
allocate(gvm(nhor))   ; gvm(:)   = 0.0 ! v * sin(phi)
allocate(gp(nhor))    ; gp(:)    = 0.0 ! Ln(Ps)
allocate(go(nugp))    ; go(:)    = 0.0 ! Orography
allocate(gop(nhor))   ; gop(:)   = 0.0 ! Orography, partial
allocate(gland(nhor)) ; gland(:) = 0.0 ! Land-Mask
allocate(gpsi(nhor))  ; gpsi(:)  = 0.0 ! Streamfunction
allocate(gchi(nhor))  ; gchi(:)  = 0.0 ! Velocity potential
allocate(gpvo(nhor))  ; gpvo(:)  = 0.0 ! Potential vorticity
allocate(genst(nhor)) ; genst(:) = 0.0 ! Enstrophy
allocate(gtaux(nhor)) ; gtaux(:) = 0.0 ! Zonal wind stress
allocate(gtauy(nhor)) ; gtauy(:) = 0.0 ! Merdional wind stress
allocate(gzforce(nhor)) ; gzforce(:) = 0.0 ! Vorticity forcing
allocate(gflam(nhor))   ; gflam(:)  = 0.0  ! Term Fu in Primitive Equations
allocate(gfphi(nhor))   ; gfphi(:)  = 0.0  ! Term Fv in Primitive Equations
allocate(gtn(nhor))   ; gtn(:)  = 0.0  ! Dummy
allocate(gfu(nhor))   ; gfu(:)  = 0.0  ! Term Fu in Primitive Equations
allocate(gfv(nhor))   ; gfv(:)  = 0.0  ! Term Fv in Primitive Equations
allocate(gpu(nhor))   ; gpu(:)  = 0.0  ! Term u * p
allocate(gpv(nhor))   ; gpv(:)  = 0.0  ! Term v * p
allocate(gke(nhor))   ; gke(:) = 0.0 ! Kinetic energy u * u + v * v
allocate(rcsq(nhor))  ; rcsq(:)  = 0.0 ! 1 / cos2(phi)
allocate(ruido(nugp)) ; ruido(:) = 0.0   ! noise
allocate(ruidop(nhor)); ruidop(:)= 0.0  ! noise partial

allocate(outgrid2(nugp))  ; outgrid2(:) = 0.0 !
allocate(outgrid(nugp))   ; outgrid(:) = 0.0 !
allocate(outgridp(nhor)) ; outgridp(:) = 0.0 !
!

allocate(csu(nlat))   ; csu(:) = 0.0
allocate(csv(nlat))   ; csv(:) = 0.0

allocate(chlat(nlat)) ; chlat(:) = '   '
allocate(sid(nlat))   ; sid(:)   = 0.0  ! sin(phi)
allocate(sia(nlat))   ; sia(:)   = 0.0  ! sin(phi) alternating
allocate(gwd(nlat))   ; gwd(:)   = 0.0  ! Gaussian weight (phi)
allocate(csq(nlat))   ; csq(:)   = 0.0  ! cos2(phi)
allocate(rcs(nlat))   ; rcs(:)   = 0.0  ! 1/cos(phi)

allocate(nr(ncsp))    ; nr(:)   = 0  ! 1/cos(phi)

return
end

!     =====================
!     SUBROUTINE INITRANDOM
!     =====================

      subroutine initrandom
      use pumamod
      integer :: i, clock

!     Set random number generator seed

      call random_seed(size=nseedlen)
      allocate(meed(nseedlen))
      call system_clock(count=clock)
      meed(:) = clock + 37 * (/(i,i=1,nseedlen)/)
      call random_seed(put=meed)
      return
      end


!     =================
!     SUBROUTINE PROLOG
!     =================

      subroutine prolog
      use pumamod

      character(8) :: cpuma
      character(80) :: pumaversion
      data cpuma/'PUMA-II '/   ! format for pumaburner version >= 2

      pumaversion = '1.1  (19-Dec-2011)'
 
      call clear_arrays

      if (mypid == NROOT) then
         call cpu_time(tmstart)

         write(nud,'(/," ****************************************")')
         write(nud,'(" * SAM  ",a31," *")') trim(pumaversion)
         write(nud,'(" ****************************************")')
         write(nud,'(" * NTRU =",i4,"  NLON = ",i4,"   NLAT =",i4," *")') &
            NTRU,NLON,NLAT
         write(nud,'(" ****************************************")')
         if (NPRO > 1) then
           write(nud,'(/," ****************************************************")')
           do jpro = 1 , NPRO
              write (nud,'(" * CPU",i4,1x,a40," *")') jpro-1,ympname(jpro)
           enddo
           write(nud,'(" ****************************************************")')
         endif
         call restart_ini(lrestart,'sam_restart')
         call inigau(NLAT,sid,gwd)
         call inilat
         call legpri
         call readnl
         call initpm
         call altlat(csq,NLAT)
         sia(:) = sid(:)
         call altlat(sia,NLAT)
         call initrandom
         if (ngui > 0) call guistart
         if (nruido /= 0) call initruido
         if (nstop  > 0) nrun = nstop-nstep
         if (nyears > 0) nrun = nyears * ntspd * 360
         if (noutput > 0) then
            open(40,file='sam_output',form='unformatted')
            write(40) CPUMA
            write(40) NTRU
            write(40) NLAT
            write(40) 1
            write(40) 0.5
         endif ! (noutput > 0)
         OPEN(36,FILE='EN.DAT',FORM='FORMATTED')
         open(35,FILE='SPECTRUM.DAT',FORM='formatted')
      endif ! (mypid == NROOT)
                              ! 2 = Polvani & Kushner
!
!     ***********************
!     * broadcast & scatter *
!     ***********************

      call mpscdn(sid,NHPP)  ! kind 8 - reg - partial
      call mpscdn(gwd,NHPP)  ! kind 8 - reg - partial
      call mpscrn(csq,NLPP)  ! real   - alt - partial
      call mpscrn(sia,NLPP)  ! real   - alt - partial

      do jlat = 1 , NLPP
         rcsq(1+(jlat-1)*NLON:jlat*NLON) = 1.0 / csq(jlat)
      enddo

      call mpbci(kick    ) ! add noise for kick > 0
      call mpbci(nafter  ) ! write data interval
      call mpbci(ncoeff  ) ! number of modes to print
      call mpbci(ndel    ) ! ndel
      call mpbci(noutput ) ! global output switch
      call mpbci(ndiag   ) ! write diagnostics interval
      call mpbci(ngui    ) ! GUI on (1) or off (0)
      call mpbci(nkits   ) ! number of initial timesteps
      call mpbci(nrun    ) ! if (nstop == 0) nstop = nstep + nrun
      call mpbci(nstep   ) ! current timestep
      call mpbci(nstop   ) ! finishing timestep
      call mpbci(ntspd   ) ! number of timesteps per day
      call mpbci(nyears  ) ! simulation time
      call mpbci(nruido )  ! no noise, temp. noise, spatio-temporal noise
      call mpbci(nf )      ! Forcing wavenumber
      call mpbci(mmax )    !  Maximum zonal waven. for climat. forcing
      call mpbci(nmax )    !  Maximum meridional waven. for climat. forcing
      call mpbci(nexp )    !  Experiment number

      call mpbcl(lrestart) ! true: read restart file, false: initial run
      call mpbcl(lselect ) ! true: disable some zonal waves
      call mpbcl(lspecsel) ! true: disable some spectral modes
      call mpbcl(lbal)     ! true: balancing the initial states
      call mpbcl(llid)     ! true: introduction of a rigid lid
      call mpbcl(lequiv)   ! true: equivalent barotropic model

      call mpbcr(tdiss   )
      call mpbcr(tac     )
      call mpbcr(pac     )
      call mpbcr(plavor  )
      call mpbcr(rotspd  )
      call mpbcr(disp    )
      call mpbcr(ah)
      call mpbcr(restim)

      call mpbcin(ndl  ,1)

      call mpbcrn(fric  ,1)
      call mpbcrn(tfrc  ,1)


      call mpscin(nindex,NSPP)
      call mpscin(nselect,NTP1)
      call mpscin(nspecsel,NCSP)
      call mpscrn(srcn  ,NSPP)
      call mpscrn(sak   ,NSPP)

      call legini(nlat,nlpp,nesp,nlev,plavor,sid,gwd)

      if (lrestart) then
         call read_atmos_restart
         call mpbcrn(so,NESP)
         call sp2fc(so,gop)
         call fc2gp(gop,NLON,NLPP)
         call mpgallsp(spo,spm,1)
         call mpgallsp(sdo,sdm,1)
         call mpgallsp(szo,szm,1)
         if (mypid == NROOT) then
!            if (kick > 10) call noise(kick-10)
         endif
      else
         call initfd
      endif
      call mpscrn(riszforce,NESP)
      call mpscsp(riszforce,riszforcep,1)
      call mpscrn(salpha ,NSPP)
      call mpbcrn(sp,NESP)
      call mpbcrn(sp3,NESP)
      call mpbcrn(sp2,NESP)
      call mpbcrn(sd,NESP)
      call mpbcrn(sz,NESP)
      call mpscsp(sd,sdp,1)
      call mpscsp(sz,szp,1)
      call mpscsp(sp,spp,1)
      call mpscsp(so,sop,1)
      call mpscsp(szr,szrp,1)
      call mpscgp(ruido,ruidop,1)

      return
      end

      subroutine specheck(csp)
      use pumamod
      real :: csp(NESP)
      real :: zspa(NESP)
      real :: zspb(NESP)
      real :: zg(NHOR)

      zspa = csp

      call sp2fc(zspa,zg)
      call fc2gp(zg,NLON,NLPP)
      call gp2fc(zg,NLON,NLPP)
      call fc2sp(zg,zspb)

      zdiff = zspa(1) - zspb(1)
      write (88,'(i5,3e15.6)') nstep,zspa(1),zspb(1),zdiff
      return
      end
      

!     =================
!     SUBROUTINE MASTER
!     =================

      subroutine master
      use pumamod

      real :: zgp(NHOR)

!     ***************************
!     * short initial timesteps *
!     ***************************

      ikits = nkits
      do 1000 jkits=1,ikits
         delt  = (TWOPI/ntspd) / (2**nkits)
         delt2 = delt + delt
         call gridpoint
         call spectral
         nkits = nkits - 1
 1000 continue
      delt  = TWOPI/ntspd
      delt2 = delt + delt

      nstep1 = nstep ! remember 1.st timestep
      do 2000 jstep=1,nrun
         nstep = nstep + 1
         call ntomin(nstep,ndatim(5),ndatim(4),ndatim(3),ndatim(2),ndatim(1))
         ndatim(2) = ndatim(2) + 1
         ndatim(3) = ndatim(3) + 1

!        ************************************************************
!        * calculation of non-linear quantities in grid point space *
!        ************************************************************

         call gridpoint

!
!        Preliminaries for energy
!
        if (mod(nstep,nafter) == 0 .and. noutput > 0) then
           gke   = gu * gu +  gv * gv
           gke   = gke * rcsq * 0.5
           gpvo  = gz / (1.+gp)
           zgp   = gu
           genst = gz * gz
           call gp2fc(gke  ,NLON,NLPP)
           call gp2fc(gpvo ,NLON,NLPP)
           call gp2fc(zgp  ,NLON,NLPP)
           call gp2fc(genst,NLON,NLPP)
           call fc2sp(gke  ,ske )
           call fc2sp(gpvo ,spvo)
           call fc2sp(zgp  ,su  )
           call fc2sp(genst,senst)
           call mpsumbcr(ske  ,NESP)
           call mpsumbcr(spvo ,NESP)
           call mpsumbcr(su   ,NESP)
           call mpsumbcr(senst,NESP)
        endif

         if (mypid == NROOT) then
            if (mod(nstep,nafter) == 0 .and. noutput > 0) call outsp
            if (mod(nstep,ndiag ) == 0 .or.  ngui    > 0) call diag
            if (ncu > 0) call checkunit
         endif
!        if (mod(nstep,nafter) == 0 .and. noutput > 0) call outgrads
         if (ngui > 0) call guistep_puma


!        ******************************
!        * adiabatic part of timestep *
!        ******************************


         call spectral
!        call specheck(sp) ! check accuray of spectral transform

         if (nshutdown > 0) return
 2000 continue
      return
      end

!     =================
!     SUBROUTINE EPILOG
!     =================

      subroutine epilog
      use pumamod

      real    (kind=8) :: zut,zst
      integer (kind=8) :: imem,ipr,ipf,isw,idr,idw

      if (mypid == NROOT) close(40) ! close output file

!     write restart file

      if (mypid == NROOT) then
         call restart_prepare('sam_status')
!        sp(1) = psmean ! save psmean
         call put_restart_integer('nstep'   ,nstep   )
         call put_restart_integer('nlat'    ,NLAT    )
         call put_restart_integer('nlon'    ,NLON    )
         call put_restart_integer('nrsp'    ,NRSP    )

!        Save current random number generator seed

         call random_seed(get=meed)
!        call put_restart_array('seed',meed,nseedlen,nseedlen,1)

         call put_restart_array('sz'  ,sz  ,NRSP,NESP,1)
         call put_restart_array('sd'  ,sd  ,NRSP,NESP,1)
         call put_restart_array('sp'  ,sp  ,NRSP,NESP,1)
         call put_restart_array('sp2' ,sp2 ,NRSP,NESP,1)
         call put_restart_array('sp3' ,sp3 ,NRSP,NESP,1)
         call put_restart_array('so'  ,so  ,NRSP,NESP,1)
         call put_restart_array('szr' ,szr ,NRSP,NESP,1)
         call put_restart_array('riszforce',riszforce,NRSP,NESP,1)
      endif

      call mpputgp('gland',gland,NHOR,1)
      call mpputsp('szm',szm,NSPP,1)
      call mpputsp('sdm',sdm,NSPP,1)
      call mpputsp('spm',spm,NSPP,1)

      if (mypid == NROOT) then
         call restart_stop
!        Get resource stats from function resources in file pumax.c
         ires = nresources(zut,zst,imem,ipr,ipf,isw,idr,idw)
         call cpu_time(tmstop)
         tmrun = tmstop - tmstart
         if (nstep > nstep1) then
            zspy = tmrun * 360.0 * real(ntspd) / (nstep - nstep1) ! sec / siy
            zypd = (24.0 * 3600.0 / zspy)                         ! siy / day
            write(nud,'(/,"****************************************")')
            if (zut > 0.0) &
            write(nud,  '("* User   time         : ", f10.3," sec *")') zut
            if (zst > 0.0) &
            write(nud,  '("* System time         : ", f10.3," sec *")') zst
            if (zut + zst > 0.0) tmrun = zut + zst
            write(nud,  '("* Total CPU time      : ", f10.3," sec *")') tmrun
            if (imem > 0) &
            write(nud,  '("* Memory usage        : ", f10.3," MB  *")') imem * 0.000001
            if (ipr > 0) &
            write(nud,  '("* Page reclaims       : ", i6," pages   *")') ipr
            if (ipf > 0) &
            write(nud,  '("* Page faults         : ", i6," pages   *")') ipf
            if (isw > 0) &
            write(nud,  '("* Page swaps          : ", i6," pages   *")') isw
            if (idr > 0) &
            write(nud,  '("* Disk read           : ", i6," blocks  *")') idr
            if (idw > 0) &
            write(nud,  '("* Disk write          : ", i6," blocks  *")') idw
            write(nud,'("****************************************")')
            if (zspy < 600.0) then
               write(nud,'("* Seconds per sim year: ",i6,9x,"*")') nint(zspy)
            else if (zspy < 900000.0) then
               write(nud,'("* Minutes per sim year  ",i6,9x,"*")') nint(zspy/60.0)
            else
               write(nud,'("* Days per sim year:    ",i6,5x,"*")') nint(zspy/86400.0)
            endif
            write(nud,'("* Sim years per day   :",i7,9x,"*")') nint(zypd)
            write(nud,'("****************************************")')
         endif
      endif

      return
      end

!     =============================
!     SUBROUTINE READ_ATMOS_RESTART
!     =============================

      subroutine read_atmos_restart
      use pumamod

!     read scalars and full spectral arrays

      if (mypid == NROOT) then
         call get_restart_integer('nstep',nstep)
!        call get_restart_array('seed',meed,nseedlen,1,1)
         call get_restart_array('sz' ,sz ,NRSP,NESP,1)
         call get_restart_array('sd' ,sd ,NRSP,NESP,1)
         call get_restart_array('sp' ,sp ,NRSP,NESP,1)
         call get_restart_array('sp2',sp2,NRSP,NESP,1)
         call get_restart_array('sp3',sp3,NRSP,NESP,1)
         call get_restart_array('so' ,so ,NRSP,NESP,1)
         call get_restart_array('szr',szr,NRSP,NESP,1)
         call get_restart_array('riszforce',riszforce,NRSP,NESP,1)
!        psmean = sp(1)
!        sp(1)  = 0.0
         call random_seed(put=meed)
      endif

      call mpbci(nstep)     ! broadcast current timestep
      call mpbcr(psmean)    ! broadcast mean surface pressure

!     read and scatter spectral arrays

      call mpgetsp('szm',szm,NSPP,1)
      call mpgetsp('sdm',sdm,NSPP,1)
      call mpgetsp('spm',spm,NSPP,1)
      call mpgetgp('gland',gland,NHOR,1)

      return
      end subroutine read_atmos_restart

!     =================
!     SUBROUTINE INITFD
!     =================

      subroutine initfd
      use pumamod
      real ampl(0:NTRU)
      real ampln(0:NTRU)
!     real*8 szinit(NESP)
      real*4 szinit(NESP)
      real zgp(NHOR)

      if (nkits < 1) nkits = 1

!     Look for start data and read them if there

      call read_grid('topo',so,    1,iread1)
      call read_grid('pres',sp,    1,iread2)

      if (mypid == NROOT) then
         if (iread1==0 .or. iread2==0) then
         else
            psmean = PSURF * (1.+ spnorm(1) * sp(1)) 
            sp(1)  = 0.0
            so(:) = so(:) * ga / (cv * cv) ! descale from [gpm]
            write(nud,'(a,f8.2,a)') ' Mean of Ps = ',0.01 * psmean, '[hPa]'
         endif
      endif

!
!     Set Vorticity Restoration Field of atmosphere at Rest
!
      szr(3)=szr(3)+plavor
!     
!     Test1: Global steady state nonlinear zonal geostrophic flow (Williamson et al. 1991)
!
      if (nexp.eq.1) then
        if (mypid == NROOT) then
          u0=1./12.
          sz(3)= plavor + 2.*u0/SQRT(3./2.) 
          sp(5)= - u0*(1.+u0/2.)/SQRT(5./8.)/3.*FR
        endif
        call mpscsp(sp,spm,1)
        call mpscsp(sz,szm,1)
      endif
!     
!     Test2: Spherical Gravity-Eigenoscillation without planet rotation
!
      if (nexp.eq.2) then
        if (mypid == NROOT) then
          sz(3)=plavor
          jr =-1
          jw = 0
          do jm=0,NTRU
            do jn=jm,NTRU
              jr=jr+2
              ji=jr+1
              if ((jm.eq.8).and.(jn.eq.jm+1)) then 
                sp(jr)=0.0001
              endif
            enddo 
          enddo 
        endif
        call mpscsp(sp,spm,1)
        call mpscsp(sz,szm,1)
      endif
!
!     Test3: Rossby Haurwitz Wave (Williamson et al. 1991)
!
      if (nexp.eq.3) then
        omega=0.7848e-5/WW 
        fk=omega*1.
        do jlat=1,NLPP
          do jlon=1,NLON
            coslasq=(1.-sia(jlat)**2)
            flam=2.*PI*real(jlon)/real(NLON)      
            gz(jlon+(jlat-1)*NLON)= 2.*omega*sia(jlat)                                   &
      &                              -fk*30.*sia(jlat)*coslasq**2*cos(4.*flam)  
            gp(jlon+(jlat-1)*NLON)=FR*omega/2.*(2.+omega)*coslasq                        &
      &                             +FR*fk**2/4.*coslasq**4*(5.*coslasq+26.-32./coslasq) &
      &                             +FR*(1.+omega)/15.*fk*coslasq**2*(26.-25.*coslasq)   &
      &                                                                    *cos(4.*flam) &
      &                             +FR/4.*fk**2*coslasq**4*(5.*coslasq-6.)*cos(8.*flam)  
          enddo
        enddo
        call gp2fc(gz,NLON,NLPP)
        call fc2sp(gz,sz)
        call mpsum(sz,1)
        call mpbcrn(sz,NESP)
        call gp2fc(gp,NLON,NLPP)
        call fc2sp(gp,sp)
        call mpsum(sp,1)
        if (mypid == NROOT) then
          sp(1)=0.
          sz(3)=sz(3)+plavor
        endif
        call mpbcrn(sp,NESP)
        call mpscsp(sp,spm,1)
        call mpscsp(sz,szm,1)
      endif
!
! 
!     Test4: Zonal flow over an isolated mountain (Williamson et al. 1991)
!
      if (nexp.eq.4) then
        u0=20./CV
        sz(3)= plavor + 2.*u0/SQRT(3./2.) 
        gop(:) = 0.0
        jhor = 0
        do jlat=1,NLPP
          do jlon=1,NLON
            jhor = jhor + 1
            distm=PI/9.
            dist=sqrt((2.*PI*real(jlon-3*NLON/4)/real(NLON))**2 &
       &                             +(asin(sia(jlat))-PI/6.)**2)
            if(dist.lt.distm) then
              gop(jhor)=2000./(CV**2/GA)*(1.-dist/distm)
            endif
            gp(jhor)=-FR*u0*(1.+u0/2.)*sia(jlat)**2 -FR*gop(jhor)
!             1.*0.05*exp(-0.1*(REAL(jlon-NLON/2)**2+REAL(jlat-NLAT/4)**2))
          enddo
        enddo

        call gp2fc(gop,NLON,NLPP)
        call gp2fc(gp ,NLON,NLPP)
        call fc2sp(gop,so)
        call fc2sp(gp ,sp)
        call mpsumbcr(so,NESP)
        call mpsumbcr(sp,NESP)
        call sp2fc(so,gop)
        call sp2fc(sp,gp )
        call fc2gp(gop,NLON,NLPP)
        call fc2gp(gp ,NLON,NLPP)
!       call gridpoint
!       sp=-srcn(1:NSPP)*sdt*FR-so*FR
        call mpscsp(sp,spm,1)
        call mpscsp(sz,szm,1)
      endif
!
!     Test5: Advection of cosine bell over the pole (Williamson et al. 1991)
!
      if (nexp.eq.5) then
        llid=.true.
        u0=1./12.
        sak(:)=0.
        if (mypid == NROOT) then
          sz(3)= plavor + 2.*u0/SQRT(3./2.)
        endif
        gp(:) = 0.0
        jhor = 0
        do jlat=1,NLPP
          do jlon=1,NLON
            jhor = jhor + 1
            alpha=PI/2.-0.05
            zsl  = sia(jlat)
            zlon = nlon
            gpsi(jhor)=-u0*(zsl*cos(alpha) &
       &        -sqrt(1.-zsl**2)*sin(alpha)*cos(2.*PI*real(jlon)/zlon))
            distm=PLARAD/3.
            dist=PLARAD*sqrt((2.*PI*(jlon-.75*zlon)/zlon)**2+(asin(zsl))**2)
            if(dist.lt.distm) then
              gp(jhor)=2000./(CV**2/GA)*(1.+cos(PI*dist/distm))
            endif
          enddo
        enddo   
        call gp2fc(gp  ,NLON,NLPP)
        call gp2fc(gpsi,NLON,NLPP)
        call fc2sp(gp  ,sp  )
        call fc2sp(gpsi,spsi)
        call mpsumbcr(sp  ,NESP)
        call mpsumbcr(spsi,NESP)
        call mpscsp(spsi,spsip,1)
        szm(:) = -nindex(1:NSPP) * (nindex(1:NSPP) + 1) * spsip(:)
        call mpgallsp(sz,szm,   1)
        if (mypid == NROOT) then
          sz(3)= sz(3) + plavor
        endif
        call mpscsp(sp,spm,1)
        call mpscsp(sz,szm,1)
!        itspd = 8 * 24 * nlat / 32 ! setup requires very short timestep
!        if (ntspd < itspd) ntspd = itspd
     endif
!
!     Test6: Drift of a tropical low
!
      if (nexp.eq.6) then
        llid=.true.
        do jlat=1,NLPP
          do jlon=1,NLON
            alpha=PI/2.-0.05
            distm=PLARAD/8.
            dist=PLARAD*sqrt((2.*PI*real(jlon-3*NLON/4)/real(NLON))**2 &
       &                             +(asin(sia(jlat))-PI/9.)**2)
            gz(jlon +(jlat-1)*NLON)=0.1*2.*(distm**2-dist**2)/distm**2*exp(-dist**2/distm**2)
          enddo
        enddo   
        call gp2fc(gz,NLON,NLPP)
        call fc2sp(gz,sz)
        call mpsum(sz,1)
        call mpbcrn(sz,NESP)
        if (mypid == NROOT) then
          sz(3)=sz(3)+plavor
        endif
!             
!          Balancing the initial state   (not parallel until now)
!
!         call gridpoint
!         spm=-srcn(1:NSPP)*sdt*FR
 !        sz=0.
 !        sz(3)=sz(3)+plavor
        call mpscsp(sz,szm,1)
!        call balance(1)
!        call mpgallsp(sp,spm,   1)
      endif
!
!     Test7: Jupiter Experiment by Nozawa and Yoden
!
      if (nexp.eq.7) then
!
!       Initial state: Atmosphere at rest
!
!        CALL RANDOM_SEED()    
        llid=.true.	  
!
        if (mypid == NROOT) then
          sz(3)=plavor
        endif
!
!       Forcing wavenumber
!
        nf=40
!
!       Momentum exchange coefficient
!
!
!       Forcing amplitude (dimension: 1/s^2)
!
!        fnorm0=7.85e-12
       fnorm0=2.18e-11*1.
!       fnorm0=2.18e-10*0.25
!
!     Spectral stochastic forcing after Nozawa and Yoden (1997)
!
        if (mypid == NROOT) then
          jr =-1
          jw = 0
          do jm=0,NTRU
            do jn=jm,NTRU
              jr=jr+2
              ji=jr+1
              if((jn.le.nf+2).and.(jn.ge.nf-2).and.(jm.ne.0)) then
!              if((jn.le.nf+4).and.(jn.ge.nf).and.(jm.eq.nf)) then
                iszforce(jr)=1
                iszforce(ji)=1
              else
                iszforce(jr)=0
                iszforce(ji)=0
              endif
            enddo
          enddo
          riszforce=real(iszforce)
        endif
        call mpscsp(sz,szm,1)
      endif
!
!     Test8: Forced unstable jet with mountain
!
      if (nexp.eq.8) then
         lequiv=.false.
!        llid=.true.
        restim=10.
        salpha(:)=1./restim/86400./WW
        if (mypid == NROOT) then
          salpha(3)=0.
        endif
!
!       Streamfunction of the restoration field
!
        jhor=0
        do jlat=1,NLPP
          do jlon=1,NLON
            jhor=jhor+1
            gpsi(jhor)=0.03*tanh(20./PI*(-asin(sia(jlat))+PI/4.))         &
     &                +0.03*tanh(20./PI*(-asin(sia(jlat))-PI/4.))         &
     &                +0.*0.005*sin(4.*2.*PI*real(jlon)/real(NLON))           &
     &                   *exp(-(-asin(sia(jlat))+PI/4.)**2*100./PI**2)
          enddo
        enddo
!
!      Inclusion of a mountain
!
        jhor=0
        do jlat=1,NLPP
          do jlon=1,NLON
            jhor=jhor+1
            distm=PI/9.
            dist=sqrt((2.*PI*real(jlon-3*NLON/4)/real(NLON))**2 &
       &                             +(asin(sia(jlat))-PI/6.)**2)
            if(dist.lt.distm) then
              gop(jhor)=2000./(CV**2/GA)*(1.-dist/distm)*1.
            else
              gop(jhor)=0.
            endif
            if(lequiv) then
            else
              gp(jhor)=-FR*gop(jhor)
            endif              
          enddo
        enddo   
        call gp2fc(gop,NLON,NLPP)
        call fc2sp(gop,so)
        call mpsumbcr(so,NESP)
        call sp2fc(so,gop)
        call fc2gp(gop,NLON,NLPP)
        call gp2fc(gp,NLON,NLPP)
        call fc2sp(gp,sp)
        call mpsumbcr(sp,NESP)
        call gp2fc(gpsi,NLON,NLPP)
        call fc2sp(gpsi,spsi)        
        call mpsumbcr(spsi,NESP)
        if (mypid == NROOT) then
          if(lequiv) then
            szr=-(real(nindex)*real(nindex+1)+FREQ)*spsi
          else 
            szr=-real(nindex)*real(nindex+1)*spsi 
          endif
          szr(1)=0.

          sz=szr 
          sz=0.
          sz(3)=sz(3)+plavor
          szr(3)=szr(3)+plavor
        endif
        call mpsumbcr(sz,NESP)
        call mpsumbcr(szr,NESP)
        call mpscsp(sz,szm,1)
        call mpscsp(sp,spm,1)
!        call balance(1)
      endif
!
!     Test9: 2d-Turbulence Experiment by Cho and Polvani 1996 (Tank-Experiment)
!
      if (nexp.eq.9) then
!        llid=.true.
!
!       Initial state: Turbulent flow with a specific energy spectrum
!
        if (mypid == NROOT) then
        ampl=0.
        ampln=0.
        sz=0.
!        CALL RANDOM_SEED()    
        jr =-1
        jw = 0
        do jm=0,NTRU
          do jn=jm,NTRU
            jr=jr+2
            ji=jr+1
            if(jm.eq.0) then
              CALL RANDOM_NUMBER(zr)
              sz(jr)=(zr-0.5)*2.
            else
              CALL RANDOM_NUMBER(zr)
              sz(jr)=(zr-0.5)
              CALL RANDOM_NUMBER(zr)
              sz(ji)=(zr-0.5)
            endif
            if(jm.eq.0) then		 	 
              ampl(jn)=ampl(jn)+sz(jr)**2
              ampl(jn)=ampl(jn)+sz(ji)**2
            else
              ampl(jn)=ampl(jn)+2.*sz(jr)**2
              ampl(jn)=ampl(jn)+2.*sz(ji)**2
            endif
          enddo
        enddo
        gamma=18
        jn0=14
        do jn=1,NTRU
          ampln(jn)=2.*dble(jn)**(gamma/2+1)/(dble(jn+jn0)**gamma)*dble(jn+1)/PLARAD**2
          write(nud,*) jn,ampln(jn)
        enddo
        jr =-1
        jw = 0
        do jm=0,NTRU
          do jn=jm,NTRU
            jr=jr+2
            ji=jr+1
            sz(jr)=sz(jr)*sqrt(ampln(jn)/ampl(jn))
            sz(ji)=sz(ji)*sqrt(ampln(jn)/ampl(jn))
          enddo
        enddo
        sz(1)=0.
        spsi(3:NESP)=srcn(3:NESP)*sz(3:NESP)
        zke   = dot_product(spsi(1:NZOM),sz(1:NZOM)) * 0.5       &
     &        + dot_product(spsi(NZOM+1:NRSP),sz(NZOM+1:NRSP))  
        zke=0.5*PLARAD**2*WW**2*zke
!       sz=sz*65./sqrt(2.*zke)
         sz=sz*0.011/sqrt(2.*zke)                 ! Tank! Experiment A
!       sz=sz*0.002/sqrt(2.*zke)                 ! Tank! Experiment B
!        sz=sz/sqrt(2.*zke)
!       sz=sz*200./sqrt(2.*zke)
        spsi(3:NESP)=srcn(3:NESP)*sz(3:NESP)
        zke   = dot_product(spsi(1:NZOM),sz(1:NZOM)) * 0.5       &
     &        + dot_product(spsi(NZOM+1:NRSP),sz(NZOM+1:NRSP)) 
        zke=0.5*zke*WW**2*PLARAD**2 
        write(nud,*) zke
        write(nud,*) 'U: ',sqrt(2.*zke)
        write(nud,*) tdiss
        write(nud,*) PLARAD/sqrt(2.*zke)*WW/TWOPI
        write(nud,*) NLAT
!
!       Store initial sz-field
!
        open(90,file='szinitsn.gra',FORM='unformatted',ACCESS='DIRECT',RECL=4*NESP)
        write(90,rec=1) (sz(jsp),jsp=1,NESP)
        endif
        call mpbcrn(sz,NESP)
        call sp2fc(sz,gz)
        call fc2gp(gz,NLON,NLPP)
!
!       Land-Sea mask (continent for phi<45N)
!
        jhor=0
        do jlat=1,NLPP
          do jlon=1,NLON
            jhor=jhor+1
            gland(jhor)=1.*0.5*(1.+tanh(160./PI*(-asin(sia(jlat))+1.*4./3.*PI/4.)))
            gz(jlon +(jlat-1)*NLON)=gz(jlon +(jlat-1)*NLON)*(1.-gland(jlon +(jlat-1)*NLON))   &
      &                             *(1.-cos(2.*(asin(sid(jlat))-PI/4.))**22)
          enddo
        enddo
        where (gland.lt.0.01)
          gland=0.
        endwhere
         call gp2fc(gz,NLON,NLPP)
        call fc2sp(gz,sz)
        call mpsum(sz,1)
        call mpbcrn(sz,NESP)
         call dv2uv(sd,sz,gu,gv)
         call fc2gp(gu  ,NLON,NLPP)
         call fc2gp(gv  ,NLON,NLPP)
         sd=0.
        call mpsum(sz,1)
        call mpbcrn(sz,NESP)
        if (mypid == NROOT) then
          sz(3)=sz(3)+plavor
        endif
        call mpscsp(sz,szm,1)
!        gland=0.   ! Full Sphere
         write(nud,*) mypid,FR,HS,PI
      endif
!
!     Test10: Decaying Rossby-Haurwitz Wave
!
      if (nexp.eq.10) then
        llid=.true.
        if (mypid == NROOT) then
          jr =-1
          jw = 0
          do jm=0,NTRU
            do jn=jm,NTRU
              jr=jr+2
              ji=jr+1
              if ((jm.eq.8).and.(jn.eq.jm+1)) then 
                sz(jr)=0.25
              endif
              if ((jm.eq.0).and.(jn.eq.3)) then 
                sz(jr)=-0.04*1.
              endif
!            if ((jm.eq.8).and.(jn.eq.jm+3)) then 
!              sz(ji)=-0.05
!            endif
            enddo
          enddo
          sz(3)=plavor
        endif
        call mpscsp(sz,szm,1)
      endif
!
!     Test11: Gravity wave propagation
!
      if (nexp.eq.11) then
        if (mypid == NROOT) then
          sz(3)= plavor
        endif
        flat0=10.*PI/180.
        flam0=180.*PI/180.
        jhor=0
        do jlat=1,NLPP
          do jlon=1,NLON
            jhor = jhor + 1
            flam=real(jlon-1)*2.*PI/real(NLON)-flam0
            gp(jhor)=-0.01/HS*exp(-(flam**2+1.*(asin(sia(jlat))-flat0)**2)/(2.*PI/24.)**2)
          enddo
        enddo
        call gp2fc(gp,NLON,NLPP)
        call fc2sp(gp,sp)
        call mpsum(sp,1)
        call mpbcrn(sp,NESP)
        call mpscsp(sp,spm,1)
        call mpscsp(sz,szm,1)
      endif
!
!     Test12: Idealized Gulf Stream simulations
!
      if (nexp.eq.12) then
        if (mypid == NROOT) then
          sz(3)= plavor
        endif
        ah=5.e5/WW/PLARAD**2
!        llid=.true.
        jhor=0
        do jlat=1,NLPP
          phi=asin(sia(jlat))
          if(phi.gt.PI/4.) phi=PI/4.
          do jlon=1,NLON/2
             jhor=jhor+1
             gland(jhor)=0.5*(1.+tanh((2.*PI*real(jlon-1)/real(NLON)-PI/2.)*50.))         &
                        +0.5*(1.-tanh(phi*50.))
          enddo
          do jlon=NLON/2+1,NLON
             jhor=jhor+1
             gland(jhor)=0.5*(1.-tanh((2.*PI*real(jlon-1)/real(NLON)-phi-3./2.*PI)*50.))  &
                        +0.5*(1.-tanh(phi*50.))
          enddo
        enddo
        where(gland.gt.1.) gland=1.
        jhor=0
!        gop=gland/HS
        do jlat=1,NLPP
          do jlon=1,NLON
            jhor=jhor+1
            gtaux(jhor)=-0.1*cos(3.*asin(sia(jlat)))
            gtauy(jhor)=0.
          enddo
        enddo
        gtaux=gtaux/(WW**2*RHO0*PLARAD)
        gtauy=gtauy/(WW**2*RHO0*PLARAD)
        call mpscsp(sz,szm,1)
      endif
      if(lequiv) then
!
!       Weighting orography with sin of latitude
!
        jhor=0
        do jlat=1,NLPP
          do jlon=1,NLON
            jhor=jhor+1
            gop(jhor)=gop(jhor)*2.*sia(jlat)/F0*WW
          enddo
        enddo
        call gp2fc(gop,NLON,NLPP)
        call fc2sp(gop,so)
        call fc2gp(gop,NLON,NLPP)
      endif
      if (mypid == NROOT) then
        if((llid).or.(lequiv)) then
          sd=0.
        else
          if(lbal) sd=0.
        endif
      endif
      call mpscsp(sd,sdm,1)
      if((llid).or.(lequiv)) then
      else
        if(lbal) call balance(1)
      endif
!      write(nud,*) 'mypid ',mypid,'after'
      call mpbcrn(sd,NESP)
      call mpbcrn(sz,NESP)
      call mpbcrn(sp,NESP)

      call dv2uv(sd,sz,gu,gv)
      call fc2gp(gu  ,NLON,NLPP)
      call fc2gp(gv  ,NLON,NLPP)
!     Add initial noise if wanted

      if (mypid == NROOT) then
         if (kick > 10) then
            call noise(kick-10)
         else
            call noise(kick)
         endif
      endif ! (mypid == NROOT)

      call mpscsp(sp,spm,1)
!
!        Preliminaries for energy
!
      gke   =  gu * gu +  gv * gv
      gke   =  gke*rcsq*0.5
      gpvo  =  gz / (1.+gp)
      genst = gz * gz
      call gp2fc(gke ,NLON,NLPP)
      call gp2fc(gpvo ,NLON,NLPP)
      call fc2sp(gke,ske)
      call mpsum(ske,1)
      call mpbcrn(ske,NESP)
      call fc2sp(gpvo,spvo)
      call mpsum(spvo,1)
      call mpbcrn(spvo,NESP)
      call gp2fc(gu,NLON,NLPP)
      call fc2sp(gu,su)
      call mpsum(su,1)
      call mpbcrn(su,NESP)
      call fc2gp(gu,NLON,NLPP)
      call gp2fc(genst ,NLON,NLPP)
      call fc2sp(genst,senst)
      call mpsum(senst,1)
      call mpbcrn(senst,NESP)
      if (mypid == NROOT) then
        CALL ENERGY
      endif
      return
      end

!     ==========================
!     SUBROUTINE READ_RESOLUTION
!     ==========================

      subroutine read_resolution
      use pumamod

      character (80) :: ylat 

      if (mypid == NROOT) then 
         call get_command_argument(1,ylat)
         read(ylat,*) nlat 
      endif

      call mpbci(nlat)
      return
      end  


!     =====================
!     SUBROUTINE RESOLUTION
!     =====================

      subroutine resolution
      use pumamod

      nlon = nlat + nlat ! Longitudes
      nlah = nlat / 2
      nlpp = nlat / npro
      nhpp = nlah / npro
      nhor = nlon * nlpp
      nugp = nlon * nlat
      npgp = nugp / 2

      ntru = (nlon - 1) / 3
      ntp1 = ntru + 1
      nzom = ntp1 + ntp1
      nrsp = (ntru + 1) * (ntru + 2)
      ncsp = nrsp / 2
      nspp = (nrsp + npro - 1) / npro
      nesp = nspp * npro

      return
      end


!     =================
!     SUBROUTINE READNL
!     =================

      subroutine readnl
      use pumamod

      namelist /sam_nl/ kick , nafter  , ncoeff  , ndel    , ndiag   &
                             , nkits   , nrun    , nwspini , noutput &
                   , nstep   , nstop   , ntspd   , ncu               &
                   , ngui    , nguidbg , nyears  , nmonths           &
                   , tdiss   , nruido  , disp    , ah      , restim  &
                   , ndl     , tac     , pac     , rotspd            &
                   , tfrc    , llid    , lequiv  , lbal    , nexp    &
                   , mmax    , nmax
      open(13,file='sam_namelist')
!       open(13,file='./somnamelisttank')
!      open(13,file='/home/thomas/SOM/namelistsomtsun')
      read (13,sam_nl)
      close(13)
      if (ntspd == 0) ntspd = (24 * nlat) / 32 ! automatic
      write(nud,sam_nl)
      return
      end


!     =============================
!     SUBROUTINE SELECT_ZONAL_WAVES
!     =============================

      subroutine select_zonal_waves
      use pumamod

      if (sum(nselect(:)) /= NTP1) then ! some wavenumbers disabled
         lselect = .true.
      endif
      return
      end

!     ================================
!     SUBROUTINE SELECT_SPECTRAL_MODES
!     ================================

      subroutine select_spectral_modes
      use pumamod

      if (sum(nspecsel(:)) /= NCSP) then ! some modes disabled
         lspecsel = .true.
      endif
      return
      end

!     =================
!     SUBROUTINE INITPM
!     =================

      subroutine initpm
      use pumamod

!       real           :: radea,zakk,zzakk
      real (kind=8) :: radea,zakk,zzakk

      radea  = PLARAD        ! Planet radius in high precision
      plavor = EZ * rotspd   ! Planetary vorticity

!     *************************************************************
!     * carries out all initialisation of model prior to running. *
!     * major sections identified with comments.                  *
!     * this s/r sets the model parameters and all resolution     *
!     * dependent quantities.                                     *
!     *************************************************************

      if (lrestart) nkits=0

!     ****************************************************
!     * Check for enabling / disabling zonal wavenumbers *
!     ****************************************************

      call select_zonal_waves
      call select_spectral_modes

      if (tfrc > 0.0) then
        fric = 1.0 / (TWOPI * tfrc)
      endif

!     annual cycle period and phase in timesteps

      if (tac > 0.0) tac = TWOPI / (ntspd * tac)
      pac = pac * ntspd

!     compute internal diffusion parameter

      jdelh = ndel/2
      if (tdiss > 0.0) then
         zakk = WW*(radea**ndel)/(TWOPI*tdiss*((NTRU*(NTRU+1.))**jdelh))
      else
         zakk = 0.0
      endif
      zzakk = zakk / (WW*(radea**ndel))
      ah=ah/WW/PLARAD**2
!
!    
!      if (nexp.eq.11) then
!        zzakk = 3.*1.e-17*(PI/2./7.)**ndel
!        zzakk = 3.*1.e-17
!      endif   
!
!
!      write(nud,*) 'NTRU',NTRU,(NTRU*(NTRU+1.))**jdelh

!     set coefficients which depend on wavenumber

      zrsq2 = 1.0 / sqrt(2.0)

      jr =-1
      jw = 0
      do jm=0,NTRU
         do jn=jm,NTRU
            jr=jr+2
            ji=jr+1
            jw=jw+1
            nindex(jr)=jn
            nindex(ji)=jn
            spnorm(jr)=zrsq2
            spnorm(ji)=zrsq2
            zsq = jn * (jn+1)
            if (jn > 0) then
               srcn(jr) = 1.0 / zsq
               srcn(ji) = srcn(jr)
            endif
            sak(jr) = -zzakk * zsq**jdelh
            sak(ji) = sak(jr)

!     Determine wavenumber dependent restoration rate of climatological forcing
!     (Lee and Anderson 1996)

            if (jm.le.mmax) then
              if (jn.le.nmax) then
                if(restim.gt.0.) then
!                  salpha(jr)=srcn(jr)/restim/86400./WW
!                  salpha(ji)=srcn(ji)/restim/86400./WW
                  salpha(jr)=1./restim/86400./WW
                  salpha(ji)=1./restim/86400./WW
                else
                  salpha(jr)=0.
                  salpha(ji)=0.
                endif                  
              endif
            endif

         enddo
         zrsq2=-zrsq2
      enddo
!      write(nud,*) sak
!      write(nud,*) 1./restim/86400./PLARAD**2

!     print out

      write(nud,8120)
      write(nud,8000)
      write(nud,8020) NTRU
      write(nud,8030) NLAT
      write(nud,8040) NLON
      if (zakk == 0.0) then
         write(nud,8060)
      else
         write(nud,8070) ndel
         write(nud,8080)
         write(nud,8090) zakk,ndel
         write(nud,8100) tdiss
      endif
      write(nud,8110) PNU
      write(nud,8130) HS
      write(nud,8140) PSURF/100.
      write(nud,8150) FR
      write(nud,8160) sqrt(GA*HS)
      write(nud,8000)
      write(nud,8120)
      return
 8000 format(' *****************************************************')
 8020 format(' * NTRU = ',i6,'   Triangular truncation             *')
 8030 format(' * NLAT = ',i6,'   Number of latitudes               *')
 8040 format(' * NLON = ',i6,'   Number of longitues               *')
 8060 format(' *                 No lateral dissipation            *')
 8070 format(' * ndel = ',i6,'   Lateral dissipation               *')
 8080 format(' * on vorticity and divergence                       *')
 8090 format(' * with diffusion coefficient = ',e13.4,' m**',i1,'/s *')
 8100 format(' * e-folding time for smallest scale is ',f7.3,' days *')
 8110 format(' * Robert time filter with parameter PNU =',f8.3,'   *')
 8130 format(' * Mean height of shallow water      = ',f9.2,' m   *')
 8140 format(' * Mean surface pressure             = ',f9.2,' hPa *')
 8150 format(' * Froude number                     = ',f9.2,'     *')
 8160 format(' * Surface gravity wave phase speed  = ',f9.2,' m/s *')
 8120 format(/)
      end


!     ====================
!     SUBROUTINE INITRUIDO
!     ====================

      subroutine initruido
      use pumamod

      integer (kind=4) :: iflag
      integer :: itime(8)
      real    :: zr

      call date_and_time(values=itime)
      iflag = -1*itime(5)*itime(8)*3.5+1

      if (nruido == 1) then
         zr = disp*gasdev(iflag)
         do jhor=1,NUGP
            ruido(jhor) = zr
         enddo
      elseif (nruido == 2) then
         do jhor=1,NUGP
            ruido(jhor) = disp*gasdev(iflag)
         enddo
      endif
      return
      end

!     ====================
!     SUBROUTINE STEPRUIDO
!     ====================

      subroutine stepruido
      use pumamod
     
      integer (kind=4) :: iflag
      real :: zr

      iflag = 1
      if (mypid == NROOT) then
         if (nruido == 1) then
            zr = disp*gasdev(iflag)
            do jhor=1,NUGP
              ruido(jhor) = zr
            enddo
         elseif (nruido == 2) then
            do jhor=1,NUGP
              ruido(jhor) = disp*gasdev(iflag)
            enddo
         endif
      endif
      call mpscgp(ruido,ruidop,1)
      return
      end
!
!     ====================
!     SUBROUTINE SPECFORCE
!     ====================
   
      subroutine specforce
      use pumamod
      real zrand
!
!     Spectral stochastic forcing after Nozawa and Yoden (1997)
!   
      do jspec=1,NSPP
        if(riszforcep(jspec).ge.0.5) then
         call random_number(zrand)
!          read(20,'(1X,F27.24)') zrand
         szforcenew(jspec)=(zrand-0.5)
        else
          szforcenew(jspec)=0.   
        endif
      enddo
      call mpgallsp(szforcenewg,szforcenew,   1)
      if (mypid == NROOT) then
!
!     Norm of forcing
! 
        fnorm=0.5*(dot_product(szforcenewg(1:NZOM),szforcenewg(1:NZOM)) &
     &   +dot_product(szforcenewg(NZOM+1:NRSP),szforcenewg(NZOM+1:NRSP)))
        fnorm=sqrt(fnorm)
        fnorm=fnorm0/WW**2/fnorm
      endif
      call mpbcr(fnorm)
      szforcenew=szforcenew*fnorm
!      fnorm=0.5*(dot_product(szforcenew(1:NZOM),szforcenew(1:NZOM)) &
!     &   +dot_product(szforcenew(NZOM+1:NRSP),szforcenew(NZOM+1:NRSP)))
!      fnorm=sqrt(fnorm)
!      WRITE(*,*) fnorm*WW**2
      fac=sqrt(1.-0.9604)
      szforce=szforce+(-0.02*szforce+fac*szforcenew)/1.
      return 
      end      

!     =============================
!     SUBROUTINE FILTER_ZONAL_WAVES
!     =============================

      subroutine filter_zonal_waves(pfc)
      use pumamod
      dimension pfc(2,NLON/2,NLPP)

      do jlat = 1 , NLPP
         pfc(1,1:NTP1,jlat) = pfc(1,1:NTP1,jlat) * nselect(:)
         pfc(2,1:NTP1,jlat) = pfc(2,1:NTP1,jlat) * nselect(:)
      enddo

      return
      end
      

!     ================================
!     SUBROUTINE FILTER_SPECTRAL_MODES
!     ================================

      subroutine filter_spectral_modes
      use pumamod

      j =  0
      k = -1
      do m = 0 , NTRU
         do n = m , NTRU
            k = k + 2
            j = j + 1
            if (nspecsel(j) == 0) then
               spp(k:k+1) = 0.0
               sdp(k:k+1) = 0.0
               spt(k:k+1) = 0.0
               sdt(k:k+1) = 0.0
               spm(k:k+1) = 0.0
               sdm(k:k+1) = 0.0
               if (n < NTRU) then
                  szp(k+2:k+3) = 0.0
                  szt(k+2:k+3) = 0.0
                  szm(k+2:k+3) = 0.0
               endif
            endif
         enddo
      enddo

      return
      end
      

!     ================
!     SUBROUTINE NOISE
!     ================

      subroutine noise(kickval)
      use pumamod

!     kickval = -1 : read ps from puma_sp_init
!     kickval =  0 : model runs zonally symmetric with no eddies
!     kickval =  1 : add white noise to Ps asymmetric hemispheres
!     kickval =  2 : add white noise to Ps symmetric to the equator
!     kickval =  3 : force mode(1,2) of Ps allowing reproducable runs
!     kickval =  4 : add white noise to symmetric zonal wavenumbers 7 of Ps

      integer :: kickval
      integer :: jsp, jsp1, jn, jm
      integer :: jr, ji, in
      real :: zr, zi, zscale, zrand

      zscale = 0.0001            ! amplitude of noise
      zr     = 0.001             ! kickval=3 value for mode(1,2) real
      zi     = 0.0005            ! kickval=3 value for mode(1,2) imag

      select case (kickval)
      case (-1)
         open(71, file='puma_sp_init',form='unformatted',iostat=iostat)
         if (iostat /= 0) then
            write(nud,*) ' *** kick=-1: needs file <puma_sp_init> ***'
            stop
         endif
         read(71,iostat=iostat) sp(:)
         if (iostat /= 0) then
            write(nud,*)' *** error reading file <puma_sp_init> ***'
            stop
         endif
         close(71)
         write(nud,*)'initial ps field read from puma_sp_init'
         return
      case (0)                  ! do nothing
      case (1)
         jsp1=2*NTP1+1
         do jsp=jsp1,NRSP
            call random_number(zrand)
            sp(jsp)=sp(jsp)+zscale*(zrand-0.5)
         enddo
         write(nud,*)'white noise added'
      case (2)
         jr=2*NTP1-1
         do jm=1,NTRU
            do jn=jm,NTRU
               jr=jr+2
               ji=jr+1
               if (mod(jn+jm,2) == 0) then
                  call random_number(zrand)
                  sp(jr)=sp(jr)+zscale*(zrand-0.5)
                  sp(ji)=sp(ji)+zscale*(zrand-0.5)
               endif
            enddo
         enddo
         write(nud,*)'symmetric white noise added'
      case (3)
         sp(2*NTP1+3) = sp(2*NTP1+3) + zr
         sp(2*NTP1+4) = sp(2*NTP1+4) + zi
         write(nud,*)'mode(1,2) of Ps set to (',sp(2*NTP1+3),',',sp(2*NTP1+4),')'
      case (4)
         jr=2*NTP1-1
         do jm=1,NTRU
            do jn=jm,NTRU
               jr=jr+2
               ji=jr+1
               if (mod(jn+jm,2) == 0 .and. jm == 7) then
                  call random_number(zrand)
                  sp(jr)=sp(jr)+zscale*(zrand-0.5)
                  sp(ji)=sp(ji)+zscale*(zrand-0.5)
               endif
            enddo
         enddo
         write(nud,*)'symmetric zonal wavenumbers 7 of Ps perturbed',   &
     &        'with white noise.'
      case default
         write(nud,*)'Value ',kickval  ,' for kickval not implemented.'
         stop
      end select

      if (nwspini == 1) then
         open(71, file='puma_sp_init', form='unformatted')
         write(71) sp(:)
         close(71)
      endif

      return
      end

!     ====================
!     SUBROUTINE READ_GRID
!     ====================

      subroutine read_grid(cname,psp,klev,kread)
      use pumamod

      logical :: lexist
      integer :: ihead(8)
      character(len=  *) :: cname
      character(len= 20) :: ynumber
      character(len=256) :: yfilename
      real :: psp(NESP,klev)
      real :: zgp(NUGP,klev)
      real :: zpp(NHOR,klev)

      kread = 0
      if (mypid == NROOT) then
         write(ynumber,'(I10)') NTRU
         yfilename = "puma_" // cname // "_t" // trim(adjustl(ynumber)) // ".txt"
         inquire(file=yfilename,exist=lexist)
      endif
      call mpbcl(lexist)
      if (.not. lexist) return

      if (mypid == NROOT) then
         open(65,file=yfilename,form='formatted')
         write(nud,*) 'Reading file <',trim(yfilename),'>'
         do jlev = 1 , klev
            read (65,'(8I10)') ihead(:)
            read (65,*) zgp(:,jlev)
         enddo
         close(65)
         if (cname == "pres") then
            write(nud,*) "Converting Ps to nomdimensional Ps-Anomaly"
            zscale   = log(100.0) - log(PSURF) ! Input [hPa] / PSURF [Pa]
            zgp(:,:) = (zgp(:,:)*zscale/PSURF-1.)
         endif
      endif ! (mypid == NROOT)
      
      call mpscgp(zgp,zpp,klev)
      call gp2fc(zpp,NLON,NLPP*klev)
      do jlev = 1 , klev
         call fc2sp(zpp(1,jlev),psp(1,jlev))
      enddo
      call mpsum(psp,klev)
      kread = 1
      return
      end subroutine read_grid


!     ===============
!     SUBROUTINE DIAG
!     ===============

      subroutine diag
      use pumamod
      if (noutput > 0 .and. mod(nstep,ndiag) == 0) then
         if (ncoeff > 0) call prisp
         call xsect
      endif
      call energy
      return
      end

!     ================
!     SUBROUTINE PRISP
!     ================

      subroutine prisp
      use pumamod

      character(30) :: title

      scale = 100.0
      title = 'Vorticity [10-2]'
      if (ndl.ne.0) call wrspam(sz,0,title,scale)

      title = 'Divergence [10-2]'
      if (ndl.ne.0) call wrspam(sd,0,title,scale)

      title = 'Pressure [10-3]'
      call wrspam(sp,0,title,scale)

      return
      end

!     ====================
!     SUBROUTINE POWERSPEC
!     ====================

      subroutine powerspec(pf,pspec)
      use pumamod
      real :: pf(2,NCSP)
      real :: pspec(NTP1)

      do j = 1 , NTP1
         pspec(j) = 0.5 * (pf(1,j) * pf(1,j) + pf(2,j) * pf(2,j))
      enddo

      j = NTP1 + 1
      do m = 2 , NTP1
         do l = m , NTP1
            pspec(l) = pspec(l) + pf(1,j) * pf(1,j) + pf(2,j) * pf(2,j)
            j = j + 1
         enddo
      enddo
      return
      end

!     =====================
!     SUBROUTINE POWERPRINT
!     =====================

      subroutine powerprint(text,pspec)
      use pumamod
      character(3) :: text
      real :: pspec(NTP1)

      zmax = maxval(pspec(:))
      if (zmax <= 1.0e-20) return
      zsca = 10 ** (4 - int(log10(zmax)))
      write(nud,1000)text,(int(pspec(j)*zsca),j=2,13)
      return
 1000 format(' * Power(',a3,') ',i8,11i5,' *')
      end




!     ==============
!     FUNCTION RMSSP
!     ==============

      function rmssp(pf)
      use pumamod
      real pf(NESP)

      zsum = 0.0
      zsum = zsum + &
     &          (dot_product(pf(1:NZOM),pf(1:NZOM)) * 0.5&
     &        +  dot_product(pf(NZOM+1:NRSP),pf(NZOM+1:NRSP)))
      rmssp = zsum
      return
      end

!     =================
!     SUBROUTINE ENERGY
!     =================

      subroutine energy
      use pumamod

      parameter (idim=5) ! Number of scalars for GUI timeseries

!     calculates various global diagnostic quantities
!     remove planetary vorticity so sz contains relative vorticity

      real :: spec(NTP1)
      real (kind=4) ziso(idim)


!    ***********************************************
!     calculate means - zpsitot rms vorticity
!                       zchitot rms divergence
!                       ztotp  pe    potential energy
!                       zke    ke    kinetic   energy
!                       ztot   te    total energy
!                       zenst  pens  potential enstrophy
!                       zamsp mean surface pressure
!     ***********************************************

      zsqrt2 = sqrt(2.0)
      zamsp  = 1.0 + sp(1) / zsqrt2
      ztout1 = FR
      ztout2 = 1. /2.
      zskew  = 0.0
      tke    = 0.0

      ztotp = ztout1*dot_product(sp(1:NZOM),so(1:NZOM)) * 0.5&
     &      + ztout1*dot_product(sp(NZOM+1:NRSP),so(NZOM+1:NRSP))&
     &      + ztout2*dot_product(sp(1:NZOM),sp(1:NZOM)) * 0.5&
     &      + ztout2*dot_product(sp(NZOM+1:NRSP),sp(NZOM+1:NRSP))
      zke   = ztout1*dot_product(sp(1:NZOM),ske(1:NZOM)) * 0.5&
     &      + ztout1*dot_product(sp(NZOM+1:NRSP),ske(NZOM+1:NRSP))&
     &      + ztout1*ske(1) * 0.5*sqrt(2.)
      zenst = 0.5*dot_product(sz(1:NZOM),spvo(1:NZOM)) * 0.5&
     &      + 0.5*dot_product(sz(NZOM+1:NRSP),spvo(NZOM+1:NRSP))
      ztot=zke+ztotp
      zke=zke*PSURF/GA*CV**2*PLARAD**2
      ztotp=ztotp*PSURF/GA*CV**2*PLARAD**2
      ztot=ztot*PSURF/GA*CV**2*PLARAD**2
      sz(3) = sz(3) - plavor

      if((llid).or.(lequiv)) then
        spsi(3:NESP)=srcn(3:NESP)*sz(3:NESP)
        zke   = dot_product(spsi(1:NZOM),sz(1:NZOM)) * 0.5
        tke   = dot_product(spsi(1:NZOM),sz(1:NZOM)) * 0.5      & 
     &        + dot_product(spsi(NZOM+1:NRSP),sz(NZOM+1:NRSP))  
        zenst  = dot_product(sz(1:NZOM),sz(1:NZOM)) * 0.5       &
     &        + dot_product(sz(NZOM+1:NRSP),sz(NZOM+1:NRSP))
        zskew  = dot_product(sz(1:NZOM),senst(1:NZOM)) * 0.5       &
     &         + dot_product(sz(NZOM+1:NRSP),senst(NZOM+1:NRSP))
        zke=0.5*PLARAD**2*WW**2*zke
        tke=0.5*PLARAD**2*WW**2*tke
        zskew=zskew/zenst/sqrt(zenst)
        zenst=0.5*WW**2*zenst
      endif    
    
      zpsitot = sqrt(rmssp(sz))
      zchitot = sqrt(rmssp(sd))

      ziso(1) = 0. ! T [C]
      ziso(2) = WW * zchitot * 1.0e6
      ziso(3) = 0.
      ziso(4) = 0. ! ztotp
      ziso(5) = sz(3)
      call guiput("SCALAR" // char(0) ,ziso,idim,1,1)

!     restore sz to absolute vorticity


      if (mod(nstep,ndiag) /= 0) return ! was called for GUI only
      write(nud,9001)
      write(nud,9002) nstep,zpsitot,zchitot,zke,ztotp,zamsp
      write(36,9002) nstep,zpsitot,zchitot,tke,zke,ztot,zenst,zskew,su(1),zamsp*PSURF
      write(nud,9002)
      write(nud,9011) (j,j=1,12)
      write(nud,9012)
      call powerspec(sp,spec)
      call powerprint('Pre',spec)
      call powerspec(sz,spec)
      call powerprint('Vor',spec)
      sz(3) = sz(3) + plavor
!
!     CHANGE TF
!
!     Write of energy-spectrum
!      
      do jspec=2,NTP1
        write(35,'(3E15.5)') real(nstep),real(jspec-1)            &
        ,spec(jspec)/real(jspec*(jspec-1))
      enddo
      write(35,'(2E15.5)')
!
!     END CHANG TF
!
      call powerspec(sd,spec)
      call powerprint('Div',spec)
      return
 9001 format(/,'     nstep     rms z       rms d       ke       &
     & pe+ie       msp')
 9002 format(i10,4x,9g13.5,g15.8)
!9009 format(' *',75(' '),' *')
!9010 format(' * Power(',a,') ',7e9.2,' *')
 9011 format(' * Wavenumber ',i8,11i5,' *')
 9012 format(' ',78('*'))
      end

!     =================
!     SUBROUTINE NTOMIN
!     =================

      subroutine ntomin(kstep,imin,ihou,iday,imon,iyea)
      use pumamod
      istep = kstep                          ! day [0-29] month [0-11]
      if (istep .lt. 0) istep = 0            ! min [0-59] hour  [0-23]
      imin = mod(istep,ntspd) * 1440 / ntspd ! minutes of current day
      ihou = imin / 60                       ! hours   of current day
      imin = imin - ihou * 60                ! minutes of current hour
      iday = istep / ntspd                   ! days    in this run
      imon = iday / 30                       ! months  in this run
      iday = iday - imon * 30                ! days    of current month
      iyea = imon / 12                       ! years   in this run
      imon = imon - iyea * 12                ! month   of current year
      return
      end

!     =================
!     SUBROUTINE NTODAT
!     =================

      subroutine ntodat(istep,datch)
      character(18) :: datch
      character(3) :: mona(12)
      data mona /'Jan','Feb','Mar','Apr','May','Jun',&
     &           'Jul','Aug','Sep','Oct','Nov','Dec'/
      call ntomin(istep,imin,ihou,iday,imon,iyea)
      write (datch,20030) iday+1,mona(imon+1),iyea+1,ihou,imin
20030 format(i2,'-',a3,'-',i4.4,2x,i2,':',i2.2)
      end


!     =================
!     SUBROUTINE WRSPAM
!     =================

      subroutine wrspam(ps,klev,title,scale)
      use pumamod
!
      dimension ps(NRSP)
      character(30) :: title
      character(18) :: datch

!     cab(i)=real(scale*sqrt(ps(i+i-1)*ps(i+i-1)+ps(i+i)*ps(i+i)))

      call ntodat(nstep,datch)
      write(nud,'(1x)')
      write(nud,20000)
      write(nud,20030) datch,title,klev
      write(nud,20000)
      write(nud,20020) (i,i=0,9)
      write(nud,20000)
      write(nud,20100) (cab(i),i=1,10)
      write(nud,20200) (cab(i),i=NTRU+2,NTRU+10)
      write(nud,20300) (cab(i),i=2*NTRU+2,2*NTRU+9)
      write(nud,20400) (cab(i),i=3*NTRU+1,3*NTRU+7)
      write(nud,20000)
      write(nud,'(1x)')

20000 format(1x,78('*'))
20020 format(' * n * ',10i7,' *')
20030 format(' *   * ',a18,2x,a30,'  Level ',i2,11x,'*')
20100 format(' * 0 *',f8.2,9f7.2,' *')
20200 format(' * 1 *',8x,9f7.2,' *')
20300 format(' * 2 *',15x,8f7.2,' *')
20400 format(' * 3 *',22x,7f7.2,' *')
      contains
      function cab(i)
         cab = scale * sqrt(ps(i+i-1)*ps(i+i-1)+ps(i+i)*ps(i+i))
      end function cab
      end

!     ===============
!     SUBROUTINE WRZS
!     ===============

      subroutine wrzs(zs,title,scale)
      use pumamod
!
      dimension zs(NLAT)
      character(30) :: title
      character(18) :: datch

      ip = NLAT / 16
      ia = ip/2
      ib = ia + 7 * ip
      id = NLAT + 1 - ia
      ic = id - 7 * ip

      call ntodat(nstep,datch)
      write(nud,'(1x)')
      write(nud,20000)
      write(nud,20030) datch,title
      write(nud,20000)
      write(nud,20020) (chlat(i),i=ia,ib,ip),(chlat(j),j=ic,id,ip)
      write(nud,20000)
      write(nud,20100) (int(zs(i)*scale),i=ia,ib,ip),&
     &               (int(zs(j)*scale),j=ic,id,ip)
  200 continue
      write(nud,20000)
      write(nud,'(1x)')

20000 format(1x,78('*'))
20020 format(' *    * ',16(1x,a3),' *    *')
20030 format(' *    * ',a18,2x,a30,20x,'*')
20100 format(' *    * ',16i4,' *    *')
      end

!     ================
!     SUBROUTINE XSECT
!     ================

      subroutine xsect
      use pumamod
      character(30) :: title

      scale = 10.0
      title = 'Zonal Wind [0.1 m/s]'
      call wrzs(csu,title,scale)
      title = 'Meridional Wind [0.1 m/s]'
      call wrzs(csv,title,scale)
      return
      end

!     =================
!     SUBROUTINE PACKGP
!     =================

      subroutine packgp(pu,pmin,psca,kp)
      use pumamod
      real         :: pu(2,NPGP)   ! in:  array to pack
      real(kind=4) :: pmin         ! out: minimum
      real(kind=4) :: psca         ! out: scaling factor
      integer(kind=4) :: kp(NPGP)  ! out: array for packed data
      integer(kind=4) :: ir,ii
      
      zmax = maxval(pu(:,:))
      pmin = minval(pu(:,:))
      zran = zmax - pmin           ! range of values
      if (zran > 1.0e-25) then
         psca = 64000.0 / zran
      else
         psca = 1.0
      endif

      do j = 1 , NPGP
         ir = ishft(nint((pu(1,j) - pmin) * psca),  16)
         ii = ibits(nint((pu(2,j) - pmin) * psca),0,16)
         kp(j) = ior(ir,ii)
      enddo

      return
      end

!     =================
!     SUBROUTINE PACKSP
!     =================

      subroutine packsp(pu,pp,kp)
      use pumamod
      real(kind=4)    :: pu(2,NCSP)       ! in:  array to pack
      real(kind=4)    :: pp(NTP1+1)       ! out: modes for m=0
      integer(kind=4) :: kp(NTP1+1:NCSP)  ! out: array for packed data
      integer(kind=4) :: ir,ii

      pp(1:NTP1) = pu(1,1:NTP1) ! store first modes unpacked

      zmax = maxval(pu(:,NTP1+1:NCSP))
      zmin = minval(pu(:,NTP1+1:NCSP))
      zabs = max(abs(zmax),abs(zmin))
      if (zabs > 1.0e-25) then
         zsca = 32000.0 / zabs
      else
         zsca = 1.0
      endif
      pp(NTP1+1) = zsca

      do j = NTP1+1 , NCSP
         ir = ishft(32768+nint(pu(1,j) * zsca),  16)
         ii = ibits(32768+nint(pu(2,j) * zsca),0,16)
         kp(j) = ior(ir,ii)
      enddo

      return
      end


!     ===================
!     SUBROUTINE UNPACKSP
!     ===================

      subroutine unpacksp(pu,pp,kp)
      use pumamod
      real             :: pu(2,NCSP)
      real    (kind=4) :: pp(NTP1+1)
      integer (kind=4) :: kp(NTP1+1:NCSP)
      integer (kind=4) :: ir,ii
      integer (kind=4),parameter :: ipata = 65535
      integer (kind=4),parameter :: ipatb = 32768

      pu(1,1:NTP1) = pp(1:NTP1) ! first modes at full precision
      zsca = pp(NTP1+1)

      do j = NTP1+1 , NCSP
         ir = iand(ishft(kp(j),-16),ipata) - ipatb
         ii = iand(kp(j),ipata) - ipatb
         pu(1,j) = ir / zsca
         pu(2,j) = ii / zsca
      enddo

      return
      end

!     ==================
!     SUBROUTINE WRITESP
!     ==================

      subroutine writesp(kunit,pf,kcode,klev,pscale,poff)
      use pumamod
      real :: pf(NRSP)
!     The pumaburner reads IEEE 32-bit words only
      real    (kind=4) :: zf(NRSP)
      integer (kind=4) :: ihead(8)
      integer (kind=4) :: la(NTP1+1:NCSP)
      real    (kind=4) :: za(NTP1+1)

      istep = nstep - nafter
      call ntomin(istep,nmin,nhour,nday,nmonth,nyear)
      nyear  = nyear  + 1
      nmonth = nmonth + 1
      nday   = nday   + 1

      ihead(1) = kcode
      ihead(2) = klev
      ihead(3) = nday + 100 * nmonth + 10000 * nyear
      ihead(4) = nmin + 100 * nhour
      ihead(5) = NRSP
      ihead(6) = 1
      ihead(7) = 1
      ihead(8) = 0

!     normalize ECHAM compatible and scale to physical dimensions

      zf(:) = pf(:) * spnorm(1:NRSP) * pscale
      zf(1) = zf(1) + poff ! Add offset if necessary
      write (kunit) ihead  ! IEEE 32-bit integer
      call packsp(zf,za,la)
      write (kunit) za     ! IEEE 32-bit real
      write (kunit) la     ! IEEE 32-bit integer

      return
      end

!     ================
!     SUBROUTINE OUTSP
!     ================

      subroutine outsp
      use pumamod
      real zsr(NESP)

      if (nwrioro == 1) then
         call writesp(40,so,129,0,CV*CV,0.0)
         nwrioro = 0
      endif

!     ************
!     * pressure *
!     ************

      call writesp(40,log(sp+1.),152,0,1.0,log(psmean))

!     **************
!     * divergence *
!     **************

      call writesp(40,sd,155,0,WW,0.0)

!     *************
!     * vorticity *
!     *************

!
!     Equivalent barotropic atmosphere
!
      if (lequiv) then
        sz(3)=sz(3)-plavor
        svort(3:nesp)=(real(nindex(3:nesp))*real(nindex(3:nesp)+1))                                  &
     &      /(real(nindex(3:nesp))*real(nindex(3:nesp)+1)+FREQ)*sz(3:nesp)
        svort(3)=svort(3)+plavor
        sz(3)=sz(3)+plavor
      else
        svort=sz
      endif

      zsave = svort(3)
      svort(3) = svort(3) - plavor
      call writesp(40,svort,138,0,WW,0.0)
      svort(3) = zsave

      return
      end
!
!     ====================
!     SUBROUTINE CHECKUNIT
!     ====================

      subroutine checkunit
      use pumamod

      write (ncu,1000) nstep,'sp(  1)',sp(1),sp(1)*spnorm(1)+log(psmean)
      write (ncu,1000) nstep,'sd(  1)',sd(1),sd(1)*spnorm(1)*WW
      write (ncu,1000) nstep,'sz(  1)',sz(1),sz(1)*spnorm(1)*WW

      write (ncu,1000) nstep,'sd(  1)',sd(1),sd(1)*spnorm(1)*WW
      write (ncu,1000) nstep,'sz(  1)',sz(1),sz(1)*spnorm(1)*WW

      if (100 < NRSP) then
      write (ncu,1000) nstep,'sp(100  )',sp(100),sp(100)*spnorm(100)
      write (ncu,1000) nstep,'sd(100)',sd(100),sd(100)*spnorm(100)*WW
      write (ncu,1000) nstep,'sz(100)',sz(100),sz(100)*spnorm(100)*WW
      endif

      return
 1000 format(i5,1x,a,1x,2f14.7)
      end

!     =====================
!     * SUBROUTINE LEGPRI *
!     =====================

      subroutine legpri
      use pumamod

      write(nud,231)
      write(nud,232)
      write(nud,233)
      write(nud,232)
      do 14 jlat = 1 , NLAT
         zalat = asin(sid(jlat))*180.0/PI
         write(nud,234) jlat,zalat,csq(jlat),gwd(jlat)
   14 continue
      write(nud,232)
      write(nud,231)
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
      real (kind=8) :: zcsq

      do jlat = 1 , NLAT
         zcsq       = 1.0 - sid(jlat) * sid(jlat)
         csq(jlat)  = zcsq
         rcs(jlat)  = 1.0 / sqrt(zcsq)
      enddo
      do jlat = 1 , NLAT/2
         ideg = nint(180.0/PI * asin(sid(jlat)))
         write(chlat(jlat),'(i2,a1)') ideg,'N'
         write(chlat(NLAT+1-jlat),'(i2,a1)') ideg,'S'
      enddo
      return
      end


!     ====================
!     SUBROUTINE GRIDPOINT
!     ====================

      subroutine gridpoint
      use pumamod

!*    Calculation of tendencies due to nonlinear products
!     by applying the transform method 

      real sdf(NESP)
      real szf(NESP)
      real spf(NESP)
      real zgp(NLON,NLAT)
      real (kind=4) :: zcs(NLAT)
      real (kind=4) :: zsp(NRSP)
!
!     Equivalent barotropic atmosphere
!
      if (lequiv) then
        if (mypid == NROOT) then
          sz(3)=sz(3)-plavor
          svort(3:nesp)=(real(nindex(3:nesp))*real(nindex(3:nesp)+1))                                  &
     &      /(real(nindex(3:nesp))*real(nindex(3:nesp)+1)+FREQ)*sz(3:nesp)
          svort(3)=svort(3)+plavor
          sz(3)=sz(3)+plavor
        endif
      else
        if (mypid == NROOT) then
          svort=sz
        endif
      endif
      call mpbcrn(svort,NESP)
!      
!     1. Transformation from spectral to Fourier space
!        and solving Helmholtz's equations

      call sp2fc(sd,gd)
      call sp2fc(svort,gz)
      call sp2fc(sp,gp)

      call dv2uv(sd,svort,gu,gv)

      if (lselect) then
         call filter_zonal_waves(gp)
         call filter_zonal_waves(gu)
         call filter_zonal_waves(gv)
         call filter_zonal_waves(gd)
         call filter_zonal_waves(gz)
      endif

      if (ngui > 0 .or. mod(nstep,ndiag) == 0) then
        do jlat = 1 , NLPP
          sec = CV / sqrt(csq(jlat))
          csu(jlat) = gu(1+(jlat-1)*NLON) * sec
          csv(jlat) = gv(1+(jlat-1)*NLON) * sec
        enddo
      endif

!     2. Transformation from  Fourier to gridpoint space

      call fc2gp(gu  ,NLON,NLPP)
      call fc2gp(gv  ,NLON,NLPP)
!     call landmask
      call fc2gp(gd  ,NLON,NLPP)
      call fc2gp(gz  ,NLON,NLPP)
      call fc2gp(gp  ,NLON,NLPP)

!     3. Calculation of specific kinetic energy 

      gke = gu * gu + gv * gv

!     4. Calculation of vorticity fluxes and orographic forcing (lequiv=.true.)

      if(lequiv) then
        gfv = -gu * (gz + FR*F0/WW*gop)
        gfu =  gv * (gz + FR*F0/WW*gop)
      else
        gfv = -gu * gz
        gfu =  gv * gz
      endif

!     5. Calculation of mass fluxes

      gpu =   -gu * gp
      gpv =   -gv * gp

!
!     6. Friction by obstacles 
!
!
!     Equivalent barotropic atmosphere
!
      if (lequiv) then
        if (mypid == NROOT) then
          szo(3)=szo(3)-plavor
          svort(3:nesp)=(real(nindex(3:nesp))*real(nindex(3:nesp)+1))                                  &
     &      /(real(nindex(3:nesp))*real(nindex(3:nesp)+1)+FREQ)*szo(3:nesp)
          svort(3)=svort(3)+plavor
          szo(3)=szo(3)+plavor
        endif
      else
        if (mypid == NROOT) then
          svort=szo
        endif
      endif
      call mpbcrn(svort,NESP)
      call dv2uv(sdo,svort,gum,gvm)

      call fc2gp(gum  ,NLON,NLPP)
      call fc2gp(gvm  ,NLON,NLPP)
      where(gland.gt.0.01)
!        gfu =gfu  -1.* gu/TWOPI*ntspd/4.
!        gfv =gfv  -1.* gv/TWOPI*ntspd/4.
       gfu = gfu  -1.* gum/delt/4.*gland**2 +gtaux/sqrt(rcsq)/(HS-gop*HS)
       gfv = gfv  -1.* gvm/delt/4.*gland**2 +gtauy/sqrt(rcsq)/(HS-gop*HS)
      elsewhere
       gfu = gfu  + gtaux/sqrt(rcsq)/(HS-gop*HS)
       gfv = gfv  + gtauy/sqrt(rcsq)/(HS-gop*HS)
      endwhere   

!     7. Transformation from  gridpoint to Fourier space

      call gp2fc(gpu ,NLON,NLPP)
      call gp2fc(gpv ,NLON,NLPP)
      call gp2fc(gfv ,NLON,NLPP)
      call gp2fc(gfu ,NLON,NLPP)
      call gp2fc(gke ,NLON,NLPP)

      if (lselect) then
         call filter_zonal_waves(gpu)
         call filter_zonal_waves(gpv)
         call filter_zonal_waves(gfv)
         call filter_zonal_waves(gfu)
         call filter_zonal_waves(gke)
      endif

!     8. Finalizing tendencies due to nonlinear products

      call mktend(sdf,spf,szf,gtn,gfu,gfv,gke,gpu,gpv)
      if (nruido /= 0) call stepruido
      call mpsumsc(spf,spt,1)
      call mpsumsc(sdf,sdt,1)
      call mpsumsc(szf,szt,1)

!     9. Gui calculations

      if (ngui > 0 .or. mod(nstep,ndiag) == 0) then
         call sp2fc(spsi,gpsi)
         call fc2gp(gpsi,NLON,NLPP)
         call sp2fc(schi,gchi)
         call fc2gp(gchi,NLON,NLPP)

         call guigv("GU"    // char(0),gu*(1.-gland))
         call guigv("GV"    // char(0),gv*(1.-gland))
         call guihor("GZ"   // char(0),gz*(1.-gland),1,WW*1.e4,0.0)
         call guihor("GPSI" // char(0),gpsi*(1.-gland),1,PLARAD*PLARAD*WW*1.e-9,0.0)
         call guihor("GP"   // char(0),(gp+FR*gop)*(1.-gland),1,HS,0.0)
         if ((llid).or.(lequiv)) then
         else
           call guihor("GD"   // char(0),gd*(1.-gland),1,WW*1.e6,0.0)
           call guihor("GCHI" // char(0),gchi*(1.-gland),1,PLARAD*PLARAD*WW*1.e-7,0.0)
         endif
         call guihor("GO" // char(0),gop,1,CV*CV/GA,-10.0)
         call gp2fc(gp,NLON,NLPP)
         call fc2sp(gp,span)
         call mpsum(span,1)               ! span = Ps spectral
         call mpgacs(csu)
         call mpgacs(csv)
         if (mypid == NROOT) then 
            call altcs(csu)
            call altcs(csv)
            if (ngui > 0 .and. mypid == NROOT) then 
               zsp(:) = span(1:NRSP)
               call guiput("SPAN" // char(0) ,zsp ,NCSP,-NTP1,1)
            endif
         endif
      endif ! (ngui > 0 .or. mod(nstep,ndiag) == 0)
      return
      end


!     ===================
!     SUBROUTINE SPECTRAL
!     ===================

      subroutine spectral
      use pumamod

!*    Add adiabatic and diabatic tendencies - perform leapfrog

!     The adiabatic tendencies are added using a semi implicit scheme like that 
!     of Hoskins and Simmons 1975 (Q.J.R.Meteorol.Soc.,101,637-655) (HS75)

!     Name rule for global arrays <abc>:
!     a : representation (s=spectral, g=grid, z=local)
!     b : variable (p=pressure, d=divergence, z=vorticity, t=temperature)
!     c : modifier (m=previous timestep, p=present timestep, t=tendency)

!     global arrays variable in time
!     ------------------------------

!     spt - pressure    tendency HS75 (10)
!     sdt - divergence  tendency HS75 ( 8)
!     szt - vorticity   tendency

!     spm - pressure    at previous timestep
!     sdm - divergence  at previous timestep
!     szm - vorticity   at previous timestep

!     spp - pressure    at present timestep
!     sdp - divergence  at present timestep
!     szp - vorticity   at present timestep

!     global arrays constant in time
!     ------------------------------

!     sak(NSPP)      -     = hyper diffusion
!     sop(NSPP)      - g*  = orography as geopotential
!     nindex(NSPP)   - n   = total wavenumber n for spectral modes
!     srcn(NSPP)   - 1/Cn  = 1.0 / (n * (n+1))
!     fric        1/tau F  = time constant for Rayleigh friction
!     ah                   = kinematic viscosity

      real zpm(NSPP)      ! new spm
      real zdm(NSPP)      ! new sdm
      real zzm(NSPP)      ! new szm

      real zgt(NSPP)      ! work array

!     0. Special code for experiments with mode filtering

      if (lspecsel) call filter_spectral_modes

!     1. Initialize local arrays

      zpm = spp
      zdm = sdp
      zzm = szp


!     2. Calculate divergence on timelevel t (sdt) HS75 (17)
!        which will replace the divergence tendency sdt
!        (semi implicit scheme)

!         z0 : Reverse of Froude number squared
!         zq : 1.0 / Cn
!         zt :  - script P / Fr  (pressure tendency due to nonlinear terms)
!         zm :  + Ps(t-dt) / Fr 

!         (note that geopotential is needed in HS75 (17) and, therefore,
!         the surface geopotential phi* [sop] is added

      z0 = 1./FR
      do jsp=1,NSPP
        zq = srcn(jsp)                   ! 1.0 / (n * (n + 1))
        zt = - z0 * spt(jsp)
        zm = + z0 * spm(jsp)
        za = sdt(jsp) * zq
        zb = sdm(jsp) * zq
        zgt(jsp) = zb + delt * (za + zm + sop(jsp) + zt * delt)
      enddo

      if (mypid == NROOT) then
        jsmin=3
      else
        jsmin=1
      endif
      do jsp = jsmin , NSPP
        jn = nindex(jsp)            ! total wavenumber n
        sdt(jsp) = zgt(jsp)/(delt * delt / FR + 1.0 / real(jn*(jn+1)))
      enddo

!     3. Calculate surface pressure tendency -ps HS75 (15)


      if ((llid).or.(lequiv)) then
      else
        spt = spt + sdt
      endif

!     4. Add tendencies

      spp = spm - delt2 * spt              ! spt = -ps tendency
      sdp =   2.0 * sdt - sdm              ! sdt = sdm + delt * tend.
      szp = delt2 * szt + szm              ! vorticity

!
!     Equivalent barotropic atmosphere
!
      if (lequiv) then
        if (mypid == NROOT) then
          szp(3)=szp(3)-plavor
        endif
        svortp(1:NSPP)=real(nindex(1:NSPP))*real(nindex(1:NSPP)+1)                                  &
     &          /(real(nindex(1:NSPP))*real(nindex(1:NSPP)+1)+FREQ)*szp(1:NSPP)
        spsip(3:NSPP)=-1./(real(nindex(3:NSPP))*real(nindex(3:NSPP)+1)+FREQ)*szp(3:NSPP)
        if (mypid == NROOT) then
          szp(3)=szp(3)+plavor
          svortp(3)=svort(3)+plavor
        endif
      else
        svortp=szp
        if (mypid == NROOT) then
          szp(3)=szp(3)-plavor
        endif
        spsip(3:NSPP)=-1./(real(nindex(3:NSPP))*real(nindex(3:NSPP)+1))*szp(3:NSPP)
        schip(3:NSPP)=-1./(real(nindex(3:NSPP))*real(nindex(3:NSPP)+1))*sdp(3:NSPP)
       if (mypid == NROOT) then
          szp(3)=szp(3)+plavor
        endif
      endif

!     5. Calculate  friction and biharmonic diffusion
!        sak  = diffusion
!        fric = friction

!      Annual cycle (currently deactivated)
!      zampl = cos((real(nstep)-pac)*tac)

      sdt=0.
      szt=0.
      
      sdt = sdt + sdp * (sak(1:NSPP) - fric - ah*(real(nindex(1:NSPP))*real(nindex(1:NSPP)+1)))
      szt = szt + svortp * (sak(1:NSPP) - fric - ah*(real(nindex(1:NSPP))*real(nindex(1:NSPP)+1)-2.))

!     6. Add climatological forcing
!        salpha=restoration rate (wavenumber dependent)
      if (lequiv) then
!        szt=szt + salpha(1:NSPP)*(szrp(1:NSPP) + FREQ * (spsip - spsirp))
        szt = szt + salpha(1:NSPP)*(szrp(1:NSPP) - szp)
      else
        szt = szt + salpha(1:NSPP)*(szrp(1:NSPP) - szp)
      endif
!
!     7. Add stochastic forcing

      nforce=1
      if ((nf.gt.0).and.(mod(nstep,nforce) == 0)) then
!        szforceold=szforce
        call specforce
!        szt = szt + 0.5*(szforce+szforceold)
!      else
      endif 
      szt = szt + szforce 

!        Conserve ps by forcing mode(0,0) to zero (not necessary)
!        Correct vorticity by canceling the friction and diffusion
!        applied to planetary vorticity
!        Only root node processes the first NSPP modes

      if (mypid == NROOT) then
!         spp(1) = 0.0                   ! No artificial mass conservation
         spp(2) = 0.0
         szt(3) = szt(3) + plavor * (fric - sak(3)+ ah*(real(nindex(3))*real(nindex(3)+1)-2.))
         sdt(1) = 0.
         sdt(2) = 0.
!         szt(3) = szt(3) + plavor * (fric - sak(3))
      endif

!     8. Add friction and diffusion tendencies

      sdp = sdp + delt2 * sdt
      szp = szp + delt2 * szt


!     9. Apply Robert Asselin time filter
!       (not to be done on short initial timesteps)

!                             t           t+1   t-1
      if (nkits == 0) then ! ---          ---   ---
         spm =       PNU21 * zpm + PNU * (spp + spm)
         sdm =       PNU21 * zdm + PNU * (sdp + sdm)
         szm =       PNU21 * zzm + PNU * (szp + szm)
      endif
      if((llid).or.(lequiv)) then
        sdm=0.
        sdp=0.
      endif

!     10. Gather spectral modes from all processes

      call mpgallsp(sp,spp,   1)
      call mpgallsp(sd,sdp,   1)
      call mpgallsp(sz,szp,   1)

!     Gather spectral modes from old spectral fields

      call mpgallsp(spo,spm,   1)
      call mpgallsp(sdo,sdm,   1)
      call mpgallsp(szo,szm,   1)

!     Streamfunction (for diagnostics)

      call mpgallsp(spsi,spsip,   1)
      call mpgallsp(schi,schip,   1)


      return
      end

!     =================
!     SUBROUTINE GASDEV
!     =================

!     Gaussian noise generator with zero mean and unit variance.

      function gasdev(idum) 
      implicit none
      real :: gasdev, fac, v1, v2, r, ran1
      real , save :: gset
      integer (kind=4) :: idum
      integer, save :: iset=0

      if (iset == 0) then
        r=2.0
        do while (r >= 1.0)
          v1=2.*ran1(idum)-1.
          v2=2.*ran1(idum)-1.
          r=v1**2+v2**2
        enddo
        fac=sqrt(-2.*log(r)/r)
        gset=v1*fac
        gasdev=v2*fac
        iset=1
      else
        gasdev=gset
        iset=0
      endif

      return
      end

!     =============
!     FUNCTION RAN1
!     =============

!     New random number generator written completely in Fortran 90. 

      function ran1(idum)
      implicit none
      integer, parameter :: k4b=selected_int_kind(9)
      integer(k4b), intent(inout) :: idum
      real :: ran1
      integer(k4b), parameter :: ia=16807, im=2147483647, iq=127773, ir=2836
      real, save :: am
      integer(k4b), save :: ix=-1,iy=-1,k
      integer(k4b) :: absidum
      integer(k4b),parameter :: ipata = 888889999
      integer(k4b),parameter :: ipatb = 777755555
      integer(k4b),parameter :: ione  = 1
      if (idum <= 0 .or. iy<0) then
         absidum=abs(idum)
         am=nearest(1.0,-1.0)/im
         iy=ior(ieor(ipata,absidum),ione)
         ix=ieor(ipatb,absidum)
         idum=absidum+1
      endif
      ix=ieor(ix,ishft(ix,13))
      ix=ieor(ix,ishft(ix,-17))
      ix=ieor(ix,ishft(ix,5))
      k=iy/iq
      iy=ia*(iy-k*iq)-ir*k
      if (iy < 0) iy=iy+im
      ran1=am*ior(iand(im,ieor(ix,iy)),ione)
      end function ran1

!     ==================
!     SUBROUTINE BALANCE
!     ==================

!
!     This subroutine generates a balanced initial state
!
      subroutine balance(itype_bal)
      use pumamod
      real zgt(NSPP)      ! work array
      real zpm(NSPP)      ! new spm
      real zdm(NSPP)      ! new sdm
      real zzm(NSPP)      ! new szm
      real gzt(NHOR)
      integer itype_bal    ! itype_bal=1 pressure is adapted to the vorticity field
                           ! itype_bal=2 vorticity is adapted to the pressure field     
!          
!     Balancing the initial state
!
      if(itype_bal.eq.1) then
        delt  = TWOPI/ntspd
        sdm=0.
        call gridpoint
        spm=-srcn(1:NSPP)*sdt*FR
      endif
!        write(nud,*) sdt

      call mpgallsp(sp,spm,   1)
      return
      end

      subroutine landmask
      use pumamod
!
!     Set velocity to zero over land
!
     gu=(1.-gland)*gu
     gv=(1.-gland)*gv
      where(gland.ge.0.99)
        gu=0.
        gv=0.
      endwhere
      gflam=gu
      gfphi=gv
      call gp2fc(gflam ,NLON,NLPP)
      call gp2fc(gfphi ,NLON,NLPP)
      call uv2dv(gflam,gfphi,sd,sz)
!      sz(3)=sz(3)+plavor
!      call sp2fc(sz,gz)
      return
      end
