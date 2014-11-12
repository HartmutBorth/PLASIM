      module lsgvar
      implicit none

!     PD configuration & Eurafrican Mediterranean Sea with depth=1500 m

      integer, parameter :: ien    =   72
      integer, parameter :: jen    =   76
      integer, parameter :: ken    =   22
!!      integer, parameter :: matrx  = 10080
      integer, parameter :: ken1   = ken+1
      integer, parameter :: kenm1  = ken-1
      integer, parameter :: ienjen = ien*jen
      integer, parameter :: nstoen =   10
      integer, parameter :: myear  =   12
      integer, parameter :: kb     = 4*ien+4
      integer, parameter :: km     = kb+1
      integer, parameter :: kbm    = kb+km
!!      integer, parameter :: matot  = matrx+km
      integer :: matrx
      integer :: matot 

!     common /lsgadv/
!     -----------------------------------------------------------------
!     kiterm    number of advection-diffusion sub-cycles per time step
!     a_tv      depth-dependent vertical diffusion coefficient
!     a_th      depth-dependent horizontal diffusion coefficient
!     quick_x   coefficients for zonal advection
!     curv_xp   coefficients for zonal advection
!     curv_xn   coefficients for zonal advection
!     quick_y   coefficients for meridional advection
!     curv_yp   coefficients for meridional advection
!     curv_yn   coefficients for meridional advection
!     quick_z   coefficients for vertical advection
!     curv_zp   coefficients for vertical advection
!     curv_zn   coefficients for vertical advection

      integer :: kiterm = 4
      real (kind=8) :: a_tv(ken)
      real (kind=8) :: a_th(ken)
      real (kind=8) :: quick_x(ien,2)
      real (kind=8) :: curv_xp(ien,3)
      real (kind=8) :: curv_xn(ien,3)
      real (kind=8) :: quick_y(jen,2)
      real (kind=8) :: curv_yp(jen,3)
      real (kind=8) :: curv_yn(jen,3)
      real (kind=8) :: quick_z(ien,jen,ken,2)
      real (kind=8) :: curv_zp(ien,jen,ken,3)
      real (kind=8) :: curv_zn(ien,jen,ken,3)


!     common /lsgadv2/
!     -----------------------------------------------------------------

      logical :: first_call
      logical :: debug_advect
      logical :: no_melting
!PHCONS
      logical, parameter :: debug_conservation = .false.
!PHCONS

!     common /lsgbbl/
!     -----------------------------------------------------------------

      integer :: ibbl(ien,jen)
      integer :: ibo(ien,jen)
      integer :: ibo2(ien,jen)
      real (kind=8) :: rhobo(ien,jen)
      real (kind=8) :: sbott(ien,jen)
      real (kind=8) :: tbott(ien,jen)
      real (kind=8) :: bblthick       = 300.0

!     for BBL model
!     ibbl      index indicating direction of flow in z = ibo

!     common /lsgbou/

      real (kind=8) :: temin(ien,jen)
      real (kind=8) :: salin(ien,jen)
      real (kind=8) :: hflin(ien,jen)
      real (kind=8) :: wflin(ien,jen)
      real (kind=8) :: taxin(ien,jen)
      real (kind=8) :: tayin(ien,jen)
      real (kind=8) :: sicein(ien,jen)

!     temin    input temperature at surface (Kelvin).
!     salin    input salinity at surface (0/00).
!     hflin    input heat flux at surface, positive heat release of the
!              ocean (W/m**2).
!     wflin    input net fresh water flux at surface ,positive
!              values for net precipitation (m/s).
!     taxin    input zonal wind stress at surface (Pa).
!     tayin    input meridional wind stress at surface (Pa).
!     sicein   input sea ice obtained from surface model

!     common /lsgcon/

      real (kind=8) :: g
      real (kind=8) :: rhonul
      real (kind=8) :: roref(ken)
      real (kind=8) :: sippo(ken)
      real (kind=8) :: sippow(ken)
      real (kind=8) :: tkelvin
      real (kind=8) :: tfreez
      real (kind=8) :: cp
      real (kind=8) :: entmel
      real (kind=8) :: day         = 86400.0
      real (kind=8) :: frih
      real (kind=8) :: salicc
      real (kind=8) :: erdof

!     g         gravitation, g = 9.80665 m/s.
!     rhonul    density        = 1020 kg/m**3
!     roref     reference-density (tref,sref) as function of
!               pressure (g/m**3)
!     tkelvin   temperature in Kelvin equivalent to 0 Celsius
!                 = 273.15 K (note: original value = 273.16 K)
!     sippo     difference between in-situ temperature and pot.
!               temperature in u-points (deg C).
!     sippow    same as sippo, but in w-points (deg C).
!     day       length of day in seconds.
!     frih      coefficient of horizontal friction (m**2/s)
!     cp        specific heat of sea water = 4000 Ws/(deg kg)
!     entmel    melt enthalpy of sea ice = 80*4186 Ws/kg
!     salicc    salinity of sea ice
!
!     lsgddr.h: Headers for data sets, see manual.

!     common /lsgddr/

      integer (kind=8) :: nddr(512) ! integer header
      real    (kind=8) :: oddr(512) ! real    header

!     common /lsgdia/

      integer :: nrca(ken)
      integer :: mindi
      integer :: mindj
      real (kind=8) :: psimax
      real (kind=8) :: trup(ken)
      real (kind=8) :: trdo(ken)
      real (kind=8) :: tm(ken)
      real (kind=8) :: sm(ken)

!     psimax  maximum value of barotropic stream function
!     mindi   index i of psimax
!     mindj   index j of psimax
!     trup    sum of upward transports through a certain layer
!     trdo    sum of downward transports through a certain layer
!     sm      layer-averaged salinity
!     tm      layer-averaged potential temperature
!     nrca    number of convective adjustment events in the layer

!     new common /lsgdiv/ for dividing do_lsg_step and do_lsg_init:

      integer :: monit(ien,jen)
      integer :: moni2(ien,jen)
      integer :: moni3(ien,jen)
      integer :: moni4(ien,jen)
      integer :: mons
      real (kind=8) :: psimer(jen,ken+2,4)
      real (kind=8) :: rtime

!     common /lsgdum/ Dummies for storage of auxiliary fields

      real (kind=8) :: vrdb(ienjen,ken)
      real (kind=8) :: weto(ien,jen)
      real (kind=8) :: wetu(ien,jen)
      real (kind=8) :: s2(ien,jen)
      real (kind=8) :: t2(ien,jen)
      real (kind=8) :: th1(ien,jen),th2(ien,jen)
      real (kind=8) :: sh1(ien,jen),sh2(ien,jen)
      real (kind=8) :: rh1(ien,jen),rh2(ien,jen)

!     common /lsgfie/
!     -----------------------------------------------------------------
!     utot      zonal      component of the total velocities (m/s)
!     vtot      meridional     "     "   "    "        "       "   
!     w         vertical       "     "   "    "        "       "   
!     ub        zonal      component of the barotropic velocities (m/s)
!     vb        meridonal      "     "   "      "          "        "  
!     psi       barotropic stream function (m**2/s).
!     taux      zonal      component of wind stress divided by *ronul*
!     tauy      meridional     "     "    "     "      "    "    "    
!     t         potential temperature (deg C)
!     tbound    boundary value of t (deg C)
!     s         salinity (0/00)
!     sbound    boundary value of s (0/00)
!     rh        normalized density (actual density/ronul)
!     rhdif     difference in potential density between 2-layers,
!               related to the intermediate w-level divided by *ronul*
!     convad    contains information about appliance of convective adjust-
!               ment

      real (kind=8) :: utot(ien,jen,ken) ! 3-dim
      real (kind=8) :: utota(ienjen,ken) ! 2-dim
      equivalence(utot,utota)

      real (kind=8) :: vtot(ien,jen,ken) ! 3-dim
      real (kind=8) :: vtota(ienjen,ken) ! 2-dim
      equivalence(vtot,vtota)

      real (kind=8) :: vta(ienjen,ken)
      real (kind=8) :: w(ien,jen,ken)
      real (kind=8) :: ub(ien,jen)
      real (kind=8) :: vb(ien,jen)
      real (kind=8) :: psi(ien,jen)
      real (kind=8) :: taux(ien,jen),tauxa(ienjen)
      real (kind=8) :: tauy(ien,jen),tauya(ienjen)
      real (kind=8) :: tbound(ien,jen),tba(ienjen)
      real (kind=8) :: sbound(ien,jen)
      real (kind=8) :: convad(ien,jen,ken),convaa(ienjen,ken)
      real (kind=8) :: convef(ien,jen)
      real (kind=8) :: t(ien,jen,ken),ta(ienjen,ken)
      real (kind=8) :: s(ien,jen,ken),sa(ienjen,ken)
      real (kind=8) :: rh(ien,jen,ken)
      real (kind=8) :: rhdif(ien,jen,kenm1),rhdifa(ienjen,kenm1)
      real (kind=8) :: sbound2(ien,jen)

      equivalence(t,ta)
      equivalence(s,sa)
      equivalence(vtot,vta)
      equivalence(taux,tauxa)
      equivalence(tauy,tauya)
      equivalence(tbound,tba)
      equivalence(convad,convaa)
      equivalence(rhdif,rhdifa)


!     common /lsgfil/
!     -----------------------------------------------------------------

      character(len= 7) :: filbac1 = "kleiin1"
      character(len= 7) :: filbac2 = "kleiin2"
      character(len= 7) :: filauf  = "kleiauf"
      character(len= 7) :: filpost = "datpost"
      character(len= 7) :: filswit = "kleiswi"
      character(len=11) :: filpnew

!     common /lsggri/
!     -----------------------------------------------------------------
!     ntyear    number of time steps / year
!     dt        timestep [seconds]
!     dtdays    timestep [days]
!     dphi      meridional gridsize (m)
!     dpin      1./dphi
!     dl        zonal gridsize as a function of latitude (m)
!     dli       1./dl
!     du        depth of the u-points (m)
!     dw        depth of the w-points (m)
!     ddu       layer thickness between 2 u-points (m)
!     ddw       layer thickness between 2 w-points (m)

      character(len= 2) :: exper = "xy"
      integer :: ntyear = 36
      real (kind=8) :: dt        = 86400.0
      real (kind=8) :: dtdays    =    10.0
      real (kind=8) :: dphi
      real (kind=8) :: dpin
      real (kind=8) :: dl(jen)
      real (kind=8) :: dli(jen)
      real (kind=8) :: ff(jen)
      real (kind=8) :: du(ken1)
      real (kind=8) :: dw(ken1)
      real (kind=8) :: ddu(ken1)
      real (kind=8) :: ddw(ken1)


!     common /lsgh1/
!     -----------------------------------------------------------------
!     u         zonal      component of the baroclinic modes (m**2/s)
!     v         meridional     "     "   "      "        "   (m**2/s)

      real (kind=8) :: u(ien,jen,ken) ! 3-dim
      real (kind=8) :: uc(ienjen,ken) ! 2-dim
      real (kind=8) :: ua(ienjen*ken) ! 1-dim
      equivalence(u,uc,ua)

      real (kind=8) :: v(ien,jen,ken) ! 3-dim
      real (kind=8) :: vc(ienjen,ken) ! 2-dim
      real (kind=8) :: va(ienjen*ken) ! 1-dim
      equivalence(v,vc,va)


!     common /lsgh2/
!     -----------------------------------------------------------------
!     t1        dummies for storage of auxiliary fields.
!     s1        

      real (kind=8) :: s1(ien,jen,ken) ! 3-dim
      real (kind=8) :: s1a(ienjen,ken) ! 2-dim
      equivalence(s1,s1a)

      real (kind=8) :: t1(ien,jen,ken) ! 3-dim
      real (kind=8) :: t1a(ienjen,ken) ! 2-dim
      equivalence(t1,t1a)


!     common /lsgmat/
!     -----------------------------------------------------------------
!     This common block contains variables only active when *nsve*=2.
!     These variables are used to solve the equation system for the
!     barotropic velocities. They are computed in subroutine *matrix*.

!!      real (kind=8) :: elim(kb,matrx)   ! elimination factors
!!      real (kind=8) :: trisys(km,matrx) ! triangularised band matrix
!!      real (kind=8) :: skal(matrx)      !  scaling factors
      real(kind=8),allocatable :: elim(:,:)   ! elimination factors
      real(kind=8),allocatable :: trisys(:,:) ! triangularised band matrix
      real(kind=8),allocatable :: skal(:)     ! scaling factors


!     common /lsgphi/
!     -----------------------------------------------------------------
      real (kind=8) :: phimax(4)
      real (kind=8) :: phimin(4)
      real (kind=8) :: phimin1(4)
      real (kind=8) :: phimax1(4)
      real (kind=8) :: phimax2(4)
      real (kind=8) :: phimin2(4)
      real (kind=8) :: phimax3(4)
      real (kind=8) :: phimin3(4)

!     common /lsgpre/
!     -----------------------------------------------------------------
!     p         normalized pressure , defined as actual pressure
!               normalized by division through rhonul

      real (kind=8) :: p(ien,jen,ken)
      real (kind=8) :: pa(ienjen,ken)
      equivalence(p,pa)


!     common /lsgref/
!     -----------------------------------------------------------------

      real (kind=8) :: sref,tref

!     tref      reference temperature (=2.5 deg S)
!     sref      reference salinity (=35. 0/00)

!     common /lsgriv/

      integer :: iriv(ien,jen)
      integer :: irbek(ien,jen)
      integer :: irbek2(ien,jen)
      real (kind=8) :: stor(ien,jen)
      real (kind=8) :: riv(ien,jen)
      real (kind=8) :: stomax          = 0.3
      real (kind=8) :: rglacc(ien,jen)

!     for runoff model

!     stor      storage on land points (in m)
!     riv       runoff in river
!     stomax    bucket depth
!     iriv      index indicating direction of flow

!     rglacc    for glac mask (in m)

!     common /lsgsur/
!     -----------------------------------------------------------------
!     zeta      surface elevation (m).
!     zetado    time derivative of zeta (m/s).
!     sice      ice thickness (m).
!     fluhea    heat flux (W/m**2).
!     flukhea   heat flux due to Newtonian cooling (W/m**2).
!     fluwat    net freshwater flux (m/s).
!     flukwat   fresh water flux due to newtonian cooling (m/s).

      real (kind=8) :: zeta(ien,jen),zetaa(ienjen)
      real (kind=8) :: zetado(ien,jen)
      real (kind=8) :: sice(ien,jen),sicea(ienjen)
      real (kind=8) :: fluhea(ien,jen)
      real (kind=8) :: fluwat(ien,jen)
      real (kind=8) :: flukhea(ien,jen)
      real (kind=8) :: flukwat(ien,jen)
      real (kind=8) :: brine(ien,jen)

      equivalence(zeta,zetaa)
      equivalence(sice,sicea)


!     common /lsgtap/

      integer :: itop
      integer :: idat
      integer :: itapin
      integer :: ivel
      integer :: ipsi
      integer :: ist
      integer :: idin
      integer :: inorm
      integer :: ibac1
      integer :: ibac2
      integer :: iswit
      integer :: no71
      integer :: no73
      integer :: nowat
      integer :: nohea
      integer :: no6
      integer :: no76
      integer :: no75
      integer :: nopost
      integer :: no37
      integer :: no70

!     Input files.
!     itop      number of the tapeunit with topographic data
!     idat        "    "   "     "      "   start data
!     itapin      "    "   "     "      "   parameters
!
!     Output files.
!     ivel      number of the tapeunit with velocity data
!     ipsi        "    "   "      "     "   stream functions
!     ist         "    "   "      "     "   temperature and
!               salinity data
!     idin      number of the tapeunit last written backup on
!     ibac1     number of tapeunit of backupfile1
!     ibac2       "    "     "     "  backupfile2
!     iswit       "    "     "     written on which backupfile
!                was last in use.

!     common /lsgtim/ runt

      real (kind=8) :: runt
      real (kind=8), parameter :: startt = 87000.0
      real (kind=8), parameter :: leadt  =     0.0
      real (kind=8), parameter :: diagt  =   100.0
      real (kind=8), parameter :: durat  =   500.0

!     variables and parameters for transient simulations
!     runt:      model run time
!     startt:  start time in years
!     leadt:   preliminary lead time to (transient) experiments
!     diagt:   period of diagnosing pc14
!     durat:   duration of transient experiments    

!     common /lsgtio/
!     -----------------------------------------------------------------
!     newst     number of time steps to be computed in this run.
!     ntcont    number of time steps between control outputs.
!     ntnout    number of the time step to write next output.
!     ntout     number of time steps between data outputs.
!     ntaver    number of time steps between mean data outputs.

      integer :: newst
      integer :: ntcont =   1
      integer :: ntnout =   1
      integer :: ntout
      integer :: ntaver =   1
      integer :: ntback = 120
      integer :: ntsurf =   1


!     common /lsgtop/
!     -----------------------------------------------------------------
!     depth     depth of the bottom in u-points (m).
!     depp      depth of the bottom in p-points (m).
!               maximum of surrounding u-points.
!     wet       3-d array, =1., when cell is wet, else =0.
!     wetvec    3-d array, =1., when wet vectorpoint,else=0.
!     numx      numbers of the grid-points in the ocean.
!               land points: =0
!     ddz       individual layer thicknes in p-points (m).
!     cofric    coefficient of horizontal friction *dx*dy (1/s).
!     fricne    0.5 * (cofric(ilo,jla) + cofric(ilo,jm1))    
!     fricse    0.5 * (cofric(ilo,jla) + cofric(ip1,jp1))
!     fricsw    0.5 * (cofric(ilo,jla) + cofric(ilo,jp1))
!     fricnw    0.5 * (cofric(ilo,jla) + cofric(im1,jm1))

      real (kind=8) :: wet(ien,jen,ken)
      real (kind=8) :: weta(ienjen,ken)
      real (kind=8) :: wetb(ienjen*ken)
      equivalence(wet,weta,wetb)

      logical :: lsmv(ien,jen) ! land-sea-mask vector points (sea = .true.)
      integer :: nweten
      integer :: num(ien,jen)
      integer :: numx(ien,jen)
      real (kind=8) :: depth(ien,jen),deptha(ienjen)
      real (kind=8) :: depp(ien,jen),deppa(ienjen)
      real (kind=8) :: ddz(ien,jen,ken),ddza(ienjen,ken)
      real (kind=8) :: cofric(ien,jen),cofrica(ienjen)
      real (kind=8) :: wetvec(ien,jen,ken)
      real (kind=8) :: fricne(ien,jen)
      real (kind=8) :: fricse(ien,jen)
      real (kind=8) :: fricsw(ien,jen)
      real (kind=8) :: fricnw(ien,jen)

      equivalence(cofric,cofrica)
      equivalence(depth,deptha)
      equivalence(depp,deppa)
      equivalence(ddz,ddza)


!     common /lsgtrd/
!     -----------------------------------------------------------------
!     delta     individual layer thickness in u-points (m).
!     hred      reduced depth (m).

      real (kind=8) :: delta(ien,jen,ken),deltaa(ienjen,ken)
      real (kind=8) :: hred(ien,jen,ken),hreda(ienjen,ken)

      equivalence(delta,deltaa)
      equivalence(hred,hreda)


!     common /lsgvec/
!     -----------------------------------------------------------------
!     all variables in this commonblock are introduced for the
!     sake of vectorization.
!     dlh       zonal grid step (2-dim).
!     dlih      1./dlh
!     depthl    depth in the next u-point to the left (m).

      real (kind=8) :: dlih(ien,jen),dliha(ienjen)
      real (kind=8) :: dlh(ien,jen),dlha(ienjen)
      real (kind=8) :: depthl(ien,jen),depthla(ienjen)

      equivalence(dlh,dlha)
      equivalence(dlih,dliha)
      equivalence(depthl,depthla)


!     common /lsgver/ nsve,nsmix,nt0,nt,nsflu,nscoup,grad
!     -----------------------------------------------------------------
!     nsve      gives the method for computation of barotropic
!               velocities.
!               =2    computation of the triangularised band matrix
!                     and direct computation of the barotropic
!                     velocities.
!     nscoup    switch for coupling to agcm.
!               =0    uncoupled,boundary fields are read from files.
!               =1    coupled,boundary fields are taken from common
!                     block lsgbou.
!     nt0       start time step.
!     nt        actual time step.
!     grad      grid step in degrees.

      integer :: nsve    = 2
      integer :: nsmix   = 3
      integer :: nt0
      integer :: nt
      integer :: nsflu   = 2
      integer :: nscoup  = 0
      integer :: ndbox,ndboy   ! Gridpoint for debug output
      integer :: mdatim(7) = 0 ! date and time copied from PlaSim
      integer :: n1styear  = 0 ! PlaSim start year
      integer :: ndcoup    = 0 ! timestep in days from PlaSim
      integer :: ncoupling = 0 ! coupling method from Plasim
      integer :: naqua     = 0 ! switch for aqua planet (etc.) setup 
                               ! i.e. no aditional friction + diffusion
      real (kind=8) :: grad
      character (len=3) :: ymon(12) = (/'Jan','Feb','Mar','Apr','May',  &
                           'Jun','Jul','Aug','Sep','Oct','Nov','Dec'/)
!
!     new arrays for flux correction diagnostics (see readbou_fl.f)
!
      real(kind=8) :: fldtaux(ien,jen),fldtauy(ien,jen) 
      real(kind=8) :: fldsst(ien,jen)
      real(kind=8) :: fldpme(ien,jen)
      real(kind=8) :: fldice(ien,jen)
!
!     arrays to enable perfect restart
!
      real(kind=8) :: utotold(ien,jen,ken),vtotold(ien,jen,ken)
      real(kind=8) :: wold(ien,jen,ken)
!
!     additional diagnostics
!
      integer :: nprint  = 0               ! switch for extra printout
      integer :: npri    = 1               ! i index for extra diag
      integer :: nprj    = 1               ! j index for extra diag
      integer :: nprk    = 1               ! k index for extra diag
      integer :: ndiagfl = 0               ! switch for diag output
      real (kind=8) :: dtdt(ien,jen,ken,4) ! temperature tendencies
!!FL
      real (kind=8) :: tote,toted,toteu     ! sum(cp*rho*t) (diagnostics)
!
      end module lsgvar
      subroutine actdate(kdate,ktime)
      use lsgvar
      implicit none
!     ------------------------------------------------------------------
!
!**** *actdate* gives actual date and time
!
!     by U. Mikolajewiz 12/87.
!     using standard fortran90 f90/f95 - S. Lorenz 01/2005
!
!     Output.
!     -------
!     kdate     actual date format yymmdd (integer).
!     ktime     actual time format hhmmss (integer).
!
!     Interface.
!     ----------
!     *call* *actdate*(*kdate*,*ktime*)
!
      integer, intent(out) :: ktime
      integer, intent(out) :: kdate
      integer :: idtim(8)

      call date_and_time (values=idtim)

      ktime = idtim(5)*10000 + idtim(6)*100 + idtim(7)
      kdate = idtim(1)*10000 + idtim(2)*100 + idtim(3)

      write (6,*) 'Current date:',kdate,' time:',ktime
      return
      end subroutine actdate
!     ==================================================================

      subroutine adv_quick
      use lsgvar
      implicit none
!     ------------------------------------------------------------------
!
!**** *adv_quick*.
!
!     Based on code by E. Maier-Reimer and U. Mikolajewicz,
!     modified by A. Paul 11/2000.
!     Adapted for tracers by Martin Butzin 01/2002.
!     - added a couple of corrections with respect to masking out
!       sub-bathymetry (also in equation of state), and taking into
!       account the variable thickness of the first layer due to
!       variable SSH in the calculation of advective and diffusive
!       tracer fluxes. 
!     - ensure use of absolute temperatures (K) instead of degC when
!       calculating adiabatic heat content changes due to time
!       evolution of SSH.
!     Peter Herrmann, 08/2010. Watch out for "!PHCONS" for the respective 
!     code changes.
!
!FL: modified to include T*dh/dt into the sub cycles to enhance conservation
!  
!     Purpose.
!     --------
!
!     *adv_quick* solves the advection-diffusion equation for
!     temperature *t* and salinity *s* using a numerical scheme
!     based on the 1 scheme by Leonard (1979). The actual
!     formulation is due to Farrow and Stevens (1995), e.g., it is
!     explicit in time.
!
!     Called from *step*.
!
!**   Input.
!     ------
!     t        old potential temperature.
!     s        old salinity.
!     utot     zonal      component of the total velocity.
!     vtot     meridional     "     "   "    "      "    .
!     w        vertical velocity.
!     tracer   old tracers.
!
!     Output.
!     -------
!     t        new potential temperature.
!     s        new salinity.
!     tracer   new tracers.
!
!     all via  common/ lsgfie /.
!
!
!     Interface.
!     ----------
!     *call* *adv_quick*.
!
!     ------------------------------------------------------------------
!     Declaration of local variables.
!     -------------------------------
      integer :: i,j,k,l,i1,j1,kiter,im1,im2,ip1,ip2,jm1,jm2,jm4
      integer :: jp1,jp2,kp1,kp2,ko,ku,km1,il,ir,ia

      real (kind=8) :: salsua,temsua,sumica,thicko,eps,dts,dth
      real (kind=8) :: quart,half,one,two,zero,pi
      real (kind=8) :: volum,volol,volt2,vols2,dus,dvs,dws,volua
      real (kind=8) :: tsum,tvol,ssum,svol,tsum0,ssum0
      real (kind=8) :: upos,uneg,vpos,vneg,wpos,wneg
      real (kind=8) :: totvel,difh,salsun,temsun,sumicn
      real (kind=8) :: temerr,salerr,told,sold,thick,alph
      real (kind=8) :: difzet,volz2,gewnw,gewne,gewsw,gewse,gewmm
      real (kind=8) :: difsmo,ttranw,ttrane,ttrasw,ttrase
      real (kind=8) :: stranw,strane,strasw,strase

      real (kind=8) :: ze2(ien,jen)
      real (kind=8) :: area(ken), dzt(ken)
      real (kind=8) :: utot1(ien,jen,ken)
      real (kind=8) :: utot2(ien,jen,ken)
      real (kind=8) :: vtot1(ien,jen,ken),vtot2(ien,jen,ken)
      real (kind=8) :: w1(ien,jen,ken),w2(ien,jen,ken)
      save utot1,utot2,vtot1,vtot2,w1,w2

!     adv_fe   advective transport across the eastern face of a cell
!     adv_fn   advective transport across the northern face of a cell
!     adv_fb   advective transport across the bottom face of a cell
!
      real (kind=8) :: adv_te(ien,ken,jen,2)
      real (kind=8) :: adv_tn(ien,ken,jen,2)
      real (kind=8) :: adv_tb(ien,ken,jen,2)
!
      integer       :: ipm1,ipp1,jpm1,jpp1,ipm2,ipp2,jpm2,jpp2,kpm1
      real (kind=8) :: dzetadt(ien,jen)
      real (kind=8) :: zeta1(ien,jen)
!
!     diagnostics
!
      real (kind=8) :: ztold(ien,jen,ken)
!
!*    1.       Initialisation.
!     ---------------
!
      ztold(:,:,:)=t(:,:,:)
      sice(:,:)=sice(:,:)*wet(:,:,1)
      if(nsmix < 2) sice(:,:)=0.
!
!PHCONS: just to be safe: a few additional inititalisations,
!        to ensure zero data on sub-bathymetry points.
!
      adv_te(:,:,:,:) = 0.0
      adv_tn(:,:,:,:) = 0.0
      adv_tb(:,:,:,:) = 0.0
      t1(:,:,:)       = 0.0
      s1(:,:,:)       = 0.0
!
!*    1.1      Separation of ice and sea water
!     ----------------------------------------    
!
      do j=1,jen
        do i=1,ien
          if (wet(i,j,1)<0.5) cycle
!
!         Protection to avoid negative layer thicknesses.
!
          if (nsmix>=2) sice(i,j)=dmin1(sice(i,j),dw(1)-2.)
          zeta(i,j)=dmax1(zeta(i,j),sice(i,j)-dw(1)+1.)
!
!         subtract sea ice (can not be advected) 
!
          zeta(i,j)=zeta(i,j)-sice(i,j)
        end do
      end do
!
!*    1.2     Set initial values and constants.
!     -----------------------------------------
!
!     eps avoids division 0./0.
!
      eps=1.e-20
      dts=dt/kiterm
      dth=0.5*dts
      quart=0.25
      half=0.5
      one=1.
      two=2.
      zero=0.
      pi=4.0*atan(1.0)
!
      l=(ien*(jen-4))-1
!
!*    3.        Computation of new pot temperature and salinity.
!     ------------------------------------------------
!
!-----------------------------------------------------------------------
!     Begin of loop over advection-diffusion sub-cycles.
!-----------------------------------------------------------------------
!
      do kiter=1,kiterm
!
!*    3.1 compute interpolated velocities
!     note: in the original lsg they are computet at t+1/2 and t+1
!     while in Farro & Stevens its t an t+1/2 (see comments) 

        do k=1,ken
          do j=1,jen
            do i=1,ien
!             velocity increments for sub-cycle
              dus=(utot(i,j,k)-utotold(i,j,k))/kiterm
              dvs=(vtot(i,j,k)-vtotold(i,j,k))/kiterm
              dws=(w(i,j,k)-wold(i,j,k))/kiterm
!
!FL: original LSG
!
!             velocities at center of sub-cycle
              utot1(i,j,k)=utotold(i,j,k)+(2.0*kiter-1.0)/2.0*dus
              vtot1(i,j,k)=vtotold(i,j,k)+(2.0*kiter-1.0)/2.0*dvs
              w1(i,j,k)   =wold(i,j,k)   +(2.0*kiter-1.0)/2.0*dws
!             velocities at end of sub-cycle
              utot2(i,j,k)=utotold(i,j,k)+kiter*dus
              vtot2(i,j,k)=vtotold(i,j,k)+kiter*dvs
              w2(i,j,k)   =wold(i,j,k)+kiter*dws
!
!FL: like Farrow and Stevens
!
!             velocities at center of sub-cycle
!              utot2(i,j,k)=utotold(i,j,k)+(2.0*kiter-1.0)/2.0*dus
!              vtot2(i,j,k)=vtotold(i,j,k)+(2.0*kiter-1.0)/2.0*dvs
!              w2(i,j,k)   =wold(i,j,k)   +(2.0*kiter-1.0)/2.0*dws
!             velocities at begin sub-cycle
!              utot1(i,j,k)=utotold(i,j,k)+(kiter-1)*dus
!              vtot1(i,j,k)=vtotold(i,j,k)+(kiter-1)*dvs
!              w1(i,j,k)   =wold(i,j,k)+(kiter-1)*dws
            end do
          end do
        end do
!
!*      3.2.      Computation of intermediate temperature and salinity.
!       -----------------------------------------------------
!
!-----------------------------------------------------------------------
!       This is the predictor step (second-order in space,
!       explicit in time; advection only).
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!       Calculate 2*advective transport across eastern face of "t" cells.
!-----------------------------------------------------------------------
!
        do j=1,jen
          do k=1,ken
            do i=1,ien
              if (i>1) then
                im1=i-1
              else
                im1=ien
              end if
              if (i<ien) then
                ip1=i+1
             else
                ip1=1
              end if
              if (i<(ien-1)) then
                ip2=i+2
              else
                ip2=i-ien+2
              end if
!
!PHCONS: totvel = 0.5*subcycling_timestep*centered_zonal_velocity*
!                 delta_z_at_u*delta_phi
!    Note that all horizontal quick_coefficients are one for this 
!    particular global E-grid. Masking of lateral boundaries done in 
!    calculation of the fluxes themselves. Note that now zetao is
!    taken into account in the first level.
!
              totvel=dth*utot1(i,j,k)*dphi*delta(i,j,k)
!
!PHCONS: advective_flux=totvel*(tracer_at_i + tracer_at_i+1) , only if both
!        are wet ocean cells. 
!
              adv_te(i,k,j,1)=totvel*(quick_x(i,1)*(t(i,j,k)+tkelvin)   &
     &                               +quick_x(i,2)*(t(ip1,j,k)+tkelvin))
              adv_te(i,k,j,2)=totvel*(quick_x(i,1)*s(i,j,k)+quick_x(i,2)&
     &                        *s(ip1,j,k))*wet(i,j,k)*wet(ip1,j,k)

            end do
          end do
        end do
!
!-----------------------------------------------------------------------
!       Calculate 2*advective transport across northern face of "t" cells.
!-----------------------------------------------------------------------
!
        do j=1,jen
          jm4=max(j-4,1)
          jm2=max(j-2,1)
          jm1=max(j-1,1)
          jp2=min(j+2,jen)
          do k=1,ken
            do i=1,ien
              if (i>2) then
                im2=i-2
              else
                im2=i+ien-2
              end if
              if (i>1) then
                im1=i-1
              else
                im1=ien
              end if
              if (i<ien) then
                ip1=i+1
              else
                ip1=1
              end if
              if (i<(ien-1)) then
                ip2=i+2
              else
                ip2=i-ien+2
              end if
!
!PHCONS: analogous to the zonal flux, we use the velocity at the northern
!    corner of the tracer cell (i,j), which is v(im1,jm1); tracers are 
!    again interpolated between (i,j) and the cell excactly north of it,
!    which is (im1,jm2). zetao taken into account for k=1 again.
!
              totvel=dth*vtot1(im1,jm1,k)*dlh(im1,jm1)                  &
     &              *delta(im1,jm1,k)
              adv_tn(i,k,j,1)=totvel*(quick_y(j,1)*(t(i,j,k)+tkelvin)   &
     &                               +quick_y(j,2)*(t(im1,jm2,k)+tkelvin))
              adv_tn(i,k,j,2)=totvel*(quick_y(j,1)*s(i,j,k)+quick_y(j,2)&
     &                        *s(im1,jm2,k))*wet(i,j,k)*wet(im1,jm2,k)
            end do
          end do
        end do
!
!-----------------------------------------------------------------------
!       Calculate 2*advective transport across bottom face of "t" cells.
!-----------------------------------------------------------------------
!
        do j=1,jen
          do k=1,ken-1
            kp1=min(k+1,ken)
            do i=1,ien
!
!PHCONS: totvel as vertical velocity times twice the area of the individual
!    gridcel (times half the subcycling timestep)
!
              totvel=w1(i,j,k)*dphi*dlh(i,j)*dth
!
              adv_tb(i,k,j,1)=totvel*(quick_z(i,j,k,1)*(t(i,j,k)+tkelvin) &
     &                               +quick_z(i,j,k,2)*t(i,j,k+1)+tkelvin)
              adv_tb(i,k,j,2)=totvel*(quick_z(i,j,k,1)*s(i,j,k)           &
     &                               +quick_z(i,j,k,2)*s(i,j,k+1))
            end do
          end do
        end do
!
!FL     compute dh/dt (dzetadt) and new (preliminary) zeta (zeta1)  
!       from subcycle (interpolated) u,v,w   
!
        call mkdzeta(utot1,vtot1,w1,zeta,dth,dzetadt)
        zeta1(:,:)=zeta(:,:)+dth*dzetadt(:,:)
!
        do k=1,ken
          ko=max(k-1,1)
          ku=min(k+1,ken)
          do j=1,jen
            do i=1,ien
              weto(i,j)=wet(i,j,ko)
              wetu(i,j)=wet(i,j,ku)
            end do
          end do
          if (k==1) then
            do j=1,jen
              do i=1,ien
                weto(i,j)=zero
              end do
            end do
          else if (k==ken) then
            do j=1,jen
              do i=1,ien
                wetu(i,j)=zero
              end do
            end do
          end if
!
          km1=max(k-1,1)
          do j=3,jen-2
            jm1=j-1
            jp2=j+2
            do i=1,ien
              if (wet(i,j,k)<0.5) cycle
              if (i<ien) then
                ip1=i+1
              else
                ip1=1
              end if
              if (i>1) then
                im1=i-1
              else
                im1=ien
              end if
              if (k==1) then
                volum=(zeta1(i,j)+ddz(i,j,k))*dlh(i,j)*dphi
              else
                volum=ddz(i,j,k)*dlh(i,j)*dphi
              end if
!
!PHCONS: Now taken into account SSH in the first layer in both the present
!    state and in the lateral fluxes, we safely calculate the new box-
!    averaged tracers.
!
              volt2=2.*(t(i,j,k)+tkelvin)*volum                         &
     &              -(adv_te(i,k,j,1)-adv_te(im1,k,j,1)                 &
     &              )-(adv_tn(i,k,j,1)-adv_tn(ip1,k,jp2,1))             &
     &              -(weto(i,j)*adv_tb(i,km1,j,1)-wetu(i,j)             &
     &              *adv_tb(i,k,j,1))
!
!FL: correct for dh/dt  note: use zn on left hand side and zo on the right
!                      (d(ht)=adv => hdt=adv-tdh)
!               
              if(k == 1) then
               volt2=volt2-2.*dth*(t(i,j,k)+tkelvin)*dzetadt(i,j)       &
     &                    *dlh(i,j)*dphi

              endif
!
!PHCONS: ensure zero on dry points
!
              t1(i,j,k)=0.5*volt2*wet(i,j,k)/MAX(volum,eps)
!
              vols2=2.*s(i,j,k)*volum-(adv_te(i,k,j,2)-adv_te(im1,k,j,2)&
     &              )-(adv_tn(i,k,j,2)-adv_tn(ip1,k,jp2,2))             &
     &              -(weto(i,j)*adv_tb(i,km1,j,2)-wetu(i,j)             &
     &              *adv_tb(i,k,j,2))
              if(k == 1) then
               volt2=volt2-2.*dth*s(i,j,k)*dzetadt(i,j)                 &
     &                    *dlh(i,j)*dphi
              endif
!
              s1(i,j,k)=0.5*vols2*wet(i,j,k)/MAX(volum,eps)
!
            end do
          end do
        end do
!
!*      3.3.      Computation of final temperature and salinity.
!       ----------------------------------------------
!
!-----------------------------------------------------------------------
!       This is the corrector step (third-order in space,
!       explicit in time; advection and diffusion).
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!       In computing advective transports, use temperature and salinity
!       at center of time step (i.e., from predictor step).
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!       Calculate 2*advective transport across eastern face of "t" cells.
!-----------------------------------------------------------------------
!
        do j=1,jen
          do k=1,ken
            do i=1,ien
              if (i>1) then
                im1=i-1
              else
                im1=ien
              end if
              if (i<ien) then
                ip1=i+1
              else
                ip1=1
              end if
              if (i<(ien-1)) then
                ip2=i+2
              else
                ip2=i-ien+2
              end if
              totvel=dts*utot2(i,j,k)*dphi*delta(i,j,k)
!
!PHCONS: mask out sub-bathymetric points as in predictor timestep
!
              totvel=totvel*wet(i,j,k)*wet(ip1,j,k)
!
!PHCONS: if u>0 then upos=totvel if adjacent cells wet, otherwise 0
!
              upos=0.5*(totvel+abs(totvel))*wet(im1,j,k)*wet(i,j,k)     &
     &             *wet(ip1,j,k)
!
!PHCONS: if u<0 then uneg=totvel if adjacent cells wet, otherwise 0
!
              uneg=0.5*(totvel-abs(totvel))*wet(ip2,j,k)*wet(ip1,j,k)   &
     &             *wet(i,j,k)
!
!PHCONS: normal centered diff. fluxes as above (at end of subcycle), plus
!    upstream correction flux of delta_T = .25*Tm - .5*T + .25*Tp
!
              adv_te(i,k,j,1)=totvel*(quick_x(i,1)*t1(i,j,k)            &
     &                        +quick_x(i,2)*t1(ip1,j,k))                &
     &                        -upos*(curv_xp(i,1)*t1(ip1,j,k)           &
     &                        +curv_xp(i,2)*t1(i,j,k)+curv_xp(i,3)      &
     &                        *t1(im1,j,k))                             &
     &                        -uneg*(curv_xn(i,1)*t1(ip2,j,k)           &
     &                        +curv_xn(i,2)*t1(ip1,j,k)+curv_xn(i,3)    &
     &                        *t1(i,j,k))
              adv_te(i,k,j,2)=totvel*(quick_x(i,1)*s1(i,j,k)            &
     &                        +quick_x(i,2)*s1(ip1,j,k))                &
     &                        -upos*(curv_xp(i,1)*s1(ip1,j,k)           &
     &                        +curv_xp(i,2)*s1(i,j,k)+curv_xp(i,3)      &
     &                        *s1(im1,j,k))                             &
     &                        -uneg*(curv_xn(i,1)*s1(ip2,j,k)           &
     &                        +curv_xn(i,2)*s1(ip1,j,k)+curv_xn(i,3)    &
     &                        *s1(i,j,k))
            end do
          end do
        end do
!
!-----------------------------------------------------------------------
!       Calculate 2*advective transport across northern face of "t" cells.
!-----------------------------------------------------------------------
!
        do j=1,jen
          jm4=max(j-4,1)
          jm2=max(j-2,1)
          jm1=max(j-1,1)
          jp2=min(j+2,jen)
          do k=1,ken
            do i=1,ien
              if (i>2) then
                im2=i-2
              else
                im2=i+ien-2
              end if
              if (i>1) then
                im1=i-1
              else
                im1=ien
              end if
              if (i<ien) then
                ip1=i+1
              else
                ip1=1
              end if
              if (i<(ien-1)) then
                ip2=i+2
              else
                ip2=i-ien+2
              end if
!
              totvel=dts*vtot2(im1,jm1,k)*dlh(im1,jm1)                  &
     &              *delta(im1,jm1,k)
!
!PHCONS: mask out sub-bathymetric points as in predictor step
!
              totvel=totvel*wet(i,j,k)*wet(im1,jm2,k)
              vpos=0.5*(totvel+abs(totvel))*wet(ip1,jp2,k)*wet(i,j,k)   &
     &             *wet(im1,jm2,k)
              vneg=0.5*(totvel-abs(totvel))*wet(im2,jm4,k)              &
     &             *wet(im1,jm2,k)*wet(i,j,k)
              adv_tn(i,k,j,1)=totvel*(quick_y(j,1)*t1(i,j,k)            &
     &                        +quick_y(j,2)*t1(im1,jm2,k))              &
     &                        -vpos*(curv_yp(j,1)*t1(im1,jm2,k)         &
     &                        +curv_yp(j,2)*t1(i,j,k)+curv_yp(j,3)      &
     &                        *t1(ip1,jp2,k))                           &
     &                        -vneg*(curv_yn(j,1)*t1(im2,jm4,k)         &
     &                        +curv_yn(j,2)*t1(im1,jm2,k)+curv_yn(j,3)  &
     &                        *t1(i,j,k))
              adv_tn(i,k,j,2)=totvel*(quick_y(j,1)*s1(i,j,k)            &
     &                        +quick_y(j,2)*s1(im1,jm2,k))              &
     &                        -vpos*(curv_yp(j,1)*s1(im1,jm2,k)         &
     &                        +curv_yp(j,2)*s1(i,j,k)+curv_yp(j,3)      &
     &                        *s1(ip1,jp2,k))                           &
     &                        -vneg*(curv_yn(j,1)*s1(im2,jm4,k)         &
     &                        +curv_yn(j,2)*s1(im1,jm2,k)+curv_yn(j,3)  &
     &                        *s1(i,j,k))
            end do
          end do
        end do
!
!-----------------------------------------------------------------------
!       Calculate 2*advective transport across bottom face of "t" cells.
!-----------------------------------------------------------------------
!
        do j=1,jen
          k=1
          do i=1,ien
            totvel=w2(i,j,k)*dphi*dlh(i,j)*dts
            wpos=0.5*(totvel+abs(totvel))*wet(i,j,k+2)*wet(i,j,k+1)     &
     &           *wet(i,j,k)
!
!PHCONS: Couldn't we also assume here that at the surface dtracer/dz = 0?
!
            adv_tb(i,k,j,1)=totvel*(quick_z(i,j,k,1)*t1(i,j,k)          &
     &                      +quick_z(i,j,k,2)*t1(i,j,k+1))              &
     &                      -wpos*(curv_zn(i,j,k,1)*t1(i,j,k+2)         &
     &                      +curv_zn(i,j,k,2)*t1(i,j,k+1)               &
     &                      +curv_zn(i,j,k,3)*t1(i,j,k))
            adv_tb(i,k,j,2)=totvel*(quick_z(i,j,k,1)*s1(i,j,k)          &
     &                      +quick_z(i,j,k,2)*s1(i,j,k+1))              &
     &                      -wpos*(curv_zn(i,j,k,1)*s1(i,j,k+2)         &
     &                      +curv_zn(i,j,k,2)*s1(i,j,k+1)               &
     &                      +curv_zn(i,j,k,3)*s1(i,j,k))
          end do

          do k=2,ken-2
            kp1=min(k+1,ken)
            do i=1,ien
              totvel=w2(i,j,k)*dphi*dlh(i,j)*dts
              wpos=0.5*(totvel+abs(totvel))*wet(i,j,k+2)*wet(i,j,k+1)   &
     &             *wet(i,j,k)
              wneg=0.5*(totvel-abs(totvel))*wet(i,j,k-1)*wet(i,j,k)     &
     &             *wet(i,j,k+1)
              adv_tb(i,k,j,1)=totvel*(quick_z(i,j,k,1)*t1(i,j,k)        &
     &                        +quick_z(i,j,k,2)*t1(i,j,k+1))            &
     &                        -wneg*(curv_zp(i,j,k,1)*t1(i,j,k+1)       &
     &                        +curv_zp(i,j,k,2)*t1(i,j,k)               &
     &                        +curv_zp(i,j,k,3)*t1(i,j,k-1))            &
     &                        -wpos*(curv_zn(i,j,k,1)*t1(i,j,k+2)       &
     &                        +curv_zn(i,j,k,2)*t1(i,j,k+1)             &
     &                        +curv_zn(i,j,k,3)*t1(i,j,k))
              adv_tb(i,k,j,2)=totvel*(quick_z(i,j,k,1)*s1(i,j,k)        &
     &                        +quick_z(i,j,k,2)*s1(i,j,k+1))            &
     &                        -wneg*(curv_zp(i,j,k,1)*s1(i,j,k+1)       &
     &                        +curv_zp(i,j,k,2)*s1(i,j,k)               &
     &                        +curv_zp(i,j,k,3)*s1(i,j,k-1))            &
     &                        -wpos*(curv_zn(i,j,k,1)*s1(i,j,k+2)       &
     &                        +curv_zn(i,j,k,2)*s1(i,j,k+1)             &
     &                        +curv_zn(i,j,k,3)*s1(i,j,k))
            end do
          end do
          k=ken-1
          do i=1,ien
            totvel=w2(i,j,k)*dphi*dlh(i,j)*dts
            wneg=0.5*(totvel-abs(totvel))*wet(i,j,k-1)*wet(i,j,k)       &
     &           *wet(i,j,k+1)
!
!PHCONS: same question here: couldn't we assume that dtracer/dz vanishes
!    at the bottom (as we presumably do in calculating vertical
!    diffusion)?
!
            adv_tb(i,k,j,1)=totvel*(quick_z(i,j,k,1)*t1(i,j,k)          &
     &                      +quick_z(i,j,k,2)*t1(i,j,k+1))              &
     &                      -wneg*(curv_zp(i,j,k,1)*t1(i,j,k+1)         &
     &                      +curv_zp(i,j,k,2)*t1(i,j,k)+curv_zp(i,j,k,3)&
     &                      *t1(i,j,k-1))
            adv_tb(i,k,j,2)=totvel*(quick_z(i,j,k,1)*s1(i,j,k)          &
     &                      +quick_z(i,j,k,2)*s1(i,j,k+1))              &
     &                      -wneg*(curv_zp(i,j,k,1)*s1(i,j,k+1)         &
     &                      +curv_zp(i,j,k,2)*s1(i,j,k)+curv_zp(i,j,k,3)&
     &                      *s1(i,j,k-1))

          end do
        end do
!
!       Copy old temperature and salinity:
!       Compute diffusion for whole time step dts
!
        do k=1,ken
          do j=1,jen
            do i=1,ien
              t1(i,j,k)=t(i,j,k)
              s1(i,j,k)=s(i,j,k)
            end do
          end do
        end do
!
!FL     compute dz/dt and new zeta 
!       from subcycle (interpolated) u,v,w
!
        call mkdzeta(utot1,vtot1,w1,zeta,dts,dzetadt)
        zeta(:,:)=zeta(:,:)+dts*dzetadt(:,:)*wet(:,:,1)
!
        do k=1,ken
!
!       Horizontal diffusion *difh*.
!
          difh=a_th(k)*dts
          ko=max(k-1,1)
          ku=min(k+1,ken)
          do j=1,jen
            do i=1,ien
              weto(i,j)=wet(i,j,ko)
              wetu(i,j)=wet(i,j,ku)
            end do
          end do
          if (k==1) then
            do j=1,jen
              do i=1,ien
                weto(i,j)=zero
              end do
            end do
          else if (k==ken) then
            do j=1,jen
              do i=1,ien
                wetu(i,j)=zero
              end do
            end do
          end if
!
          km1=max(k-1,1)
          do j=3,jen-2
            jm1=j-1
            jp1=j+1
            jp2=j+2
            do i=1,ien
              if (wet(i,j,k)<0.5) cycle
              if (i<ien) then
                ip1=i+1
              else
                ip1=1
              end if
              if (i>1) then
                im1=i-1
              else
                im1=ien
              end if
              if (k==1) then
                volum=(zeta(i,j)+ddz(i,j,k))*dlh(i,j)*dphi
              else
                volum=ddz(i,j,k)*dlh(i,j)*dphi
              end if
!
!FL add effect of dh/dt (dzetadt)
!
              if (k==1) then
!               nw - ne - se - sw
                volt2=2.*(t1(i,j,k)+tkelvin)*volum                      &
     &               -(adv_te(i,k,j,1)-adv_te(im1,k,j,1))               &
     &                -(adv_tn(i,k,j,1)-adv_tn(ip1,k,jp2,1))            &
     &                -(weto(i,j)*adv_tb(i,km1,j,1)-wetu(i,j)           &
     &                *adv_tb(i,k,j,1))                                 &
     &               -2.*dts*(t1(i,j,k)+tkelvin)*dzetadt(i,j)           &
     &                *dlh(i,j)*dphi                                    &
     &               +4.*difh*wet(im1,jm1,k)                            &
     &                *(t1(im1,jm1,k)-t1(i,j,k))                        &
     &                *(delta(im1,jm1,k)+delta(im1,j,k)+                &
     &                  zeta(im1,jm1)+zeta(i,j) )*half                  &
     &               +4.*difh*wet(i,jm1,k)*(t1(i,jm1,k)-t1(i,j,k))      &
     &                *(delta(im1,jm1,k)+delta(i,j,k)+                  &
     &                  zeta(i,jm1)+zeta(i,j) )*half                    &
     &               +4.*difh*wet(ip1,jp1,k)                            &
     &                *(t1(ip1,jp1,k)-t1(i,j,k))                        &
     &                *(delta(i,jp1,k)+delta(i,j,k)+                    &
     &                  zeta(ip1,jp1)+zeta(i,j) )*half                  &
     &               +4.*difh*wet(i,jp1,k)*(t1(i,jp1,k)-t1(i,j,k))      &
     &                *(delta(i,jp1,k)+delta(im1,j,k)+                  &
     &                  zeta(i,jp1)+zeta(i,j) )*half                    &
     &               +2.*a_tv(ko)*weto(i,j)*dphi*dlh(i,j)*dts           &
     &                *(t1(i,j,ko)-t1(i,j,k))*2./(ddu(ko)+ddu(k))       &
     &               +2.*a_tv(k)*wetu(i,j)*dphi*dlh(i,j)*dts            &
     &                *(t1(i,j,ku)-t1(i,j,k))*2./(ddu(k)+ddu(ku))
              else
!               nw - ne - se - sw
                volt2=2.*(t1(i,j,k)+tkelvin)                            &
     &                *volum-(adv_te(i,k,j,1)-adv_te(im1,k,j,1))        &
     &                -(adv_tn(i,k,j,1)-adv_tn(ip1,k,jp2,1))            &
     &                -(weto(i,j)*adv_tb(i,km1,j,1)-wetu(i,j)           &
     &                *adv_tb(i,k,j,1))                                 &
     &                +4.*difh*wet(im1,jm1,k)                           &
     &                *(t1(im1,jm1,k)-t1(i,j,k))                        &
     &                *(delta(im1,jm1,k)+delta(im1,j,k))*half           &
     &                +4.*difh*wet(i,jm1,k)*(t1(i,jm1,k)-t1(i,j,k))     &
     &                *(delta(im1,jm1,k)+delta(i,j,k))*half             &
     &                +4.*difh*wet(ip1,jp1,k)                           &
     &                *(t1(ip1,jp1,k)-t1(i,j,k))                        &
     &                *(delta(i,jp1,k)+delta(i,j,k))*half               &
     &                +4.*difh*wet(i,jp1,k)*(t1(i,jp1,k)-t1(i,j,k))     &
     &                *(delta(i,jp1,k)+delta(im1,j,k))*half             &
     &                +2.*a_tv(ko)*weto(i,j)*dphi*dlh(i,j)*dts          &
     &                *(t1(i,j,ko)-t1(i,j,k))*2./(ddu(ko)+ddu(k))       &
     &                +2.*a_tv(k)*wetu(i,j)*dphi*dlh(i,j)*dts           &
     &                *(t1(i,j,ku)-t1(i,j,k))*2./(ddu(k)+ddu(ku))
              endif
!
!PHCONS: first 4 lines represent the divergence of the advective fluxes of
!    the encapsulating lon-lat box, whereas the following 10 lines
!    represent the divergence of the diffusive fluxes of the original
!    inclined faces of the tracer cell. Why is advction not calculated
!    in that (more accurate) way? The last lines represent vertical
!    diffusion.
!
              t(i,j,k)=0.5*volt2*wet(i,j,k)/MAX(volum,eps)
!
              if (k==1) then
!               nw - ne - se - sw
                vols2=2.*s1(i,j,k)*volum                                &
     &               -(adv_te(i,k,j,2)-adv_te(im1,k,j,2))               &
     &                -(adv_tn(i,k,j,2)-adv_tn(ip1,k,jp2,2))            &
     &                -(weto(i,j)*adv_tb(i,km1,j,2)                     &
     &                 -wetu(i,j)*adv_tb(i,k,j,2))                      &
     &               -2.*dts*s1(i,j,k)*dzetadt(i,j)                     &
     &                *dlh(i,j)*dphi                                    &
     &               +4.*difh*wet(im1,jm1,k)                            &
     &                *(s1(im1,jm1,k)-s1(i,j,k))                        &
     &                *(delta(im1,jm1,k)+delta(im1,j,k)+                &
     &                  zeta(im1,jm1)+zeta(i,j) )*half                &
     &               +4.*difh*wet(i,jm1,k)*(s1(i,jm1,k)-s1(i,j,k))      &
     &                *(delta(im1,jm1,k)+delta(i,j,k)+                  &
     &                  zeta(i,jm1)+zeta(i,j) )*half                  &
     &               +4.*difh*wet(ip1,jp1,k)                            &
     &                *(s1(ip1,jp1,k)-s1(i,j,k))                        &
     &                *(delta(i,jp1,k)+delta(i,j,k)+                    &
     &                  zeta(ip1,jp1)+zeta(i,j) )*half                &
     &               +4.*difh*wet(i,jp1,k)*(s1(i,jp1,k)-s1(i,j,k))      &
     &                *(delta(i,jp1,k)+delta(im1,j,k)+                  &
     &                  zeta(i,jp1)+zeta(i,j) )*half                  &
     &               +2.*a_tv(ko)*weto(i,j)*dphi*dlh(i,j)*dts           &
     &                *(s1(i,j,ko)-s1(i,j,k))*2./(ddu(ko)+ddu(k))       &
     &               +2.*a_tv(k)*wetu(i,j)*dphi*dlh(i,j)*dts            &
     &                *(s1(i,j,ku)-s1(i,j,k))*2./(ddu(k)+ddu(ku))
              else
!               nw - ne - se - sw
                vols2=2.*s1(i,j,k)                                      &
     &                *volum-(adv_te(i,k,j,2)-adv_te(im1,k,j,2))        &
     &                -(adv_tn(i,k,j,2)-adv_tn(ip1,k,jp2,2))            &
     &                -(weto(i,j)*adv_tb(i,km1,j,2)-wetu(i,j)           &
     &                *adv_tb(i,k,j,2))                                 &
     &                +4.*difh*wet(im1,jm1,k)                           &
     &                *(s1(im1,jm1,k)-s1(i,j,k))                        &
     &                *(delta(im1,jm1,k)+delta(im1,j,k))*half           &
     &                +4.*difh*wet(i,jm1,k)*(s1(i,jm1,k)-s1(i,j,k))     &
     &                *(delta(im1,jm1,k)+delta(i,j,k))*half             &
     &                +4.*difh*wet(ip1,jp1,k)                           &
     &                *(s1(ip1,jp1,k)-s1(i,j,k))                        &
     &                *(delta(i,jp1,k)+delta(i,j,k))*half               &
     &                +4.*difh*wet(i,jp1,k)*(s1(i,jp1,k)-s1(i,j,k))     &
     &                *(delta(i,jp1,k)+delta(im1,j,k))*half             &
     &                +2.*a_tv(ko)*weto(i,j)*dphi*dlh(i,j)*dts          &
     &                *(s1(i,j,ko)-s1(i,j,k))*2./(ddu(ko)+ddu(k))       &
     &                +2.*a_tv(k)*wetu(i,j)*dphi*dlh(i,j)*dts           &
     &                *(s1(i,j,ku)-s1(i,j,k))*2./(ddu(k)+ddu(ku))
              endif
              s(i,j,k)=0.5*vols2*wet(i,j,k)/MAX(volum,eps)
            end do
          end do
        end do
!
!       convert t to C again
! 
        t(:,:,:)=(t(:,:,:)-tkelvin)*wet(:,:,:)
!
!-----------------------------------------------------------------------
!       End of loop over advection-diffusion sub-cycles.
!-----------------------------------------------------------------------
!
      end do
!
!     Recombination of sea ice and water.
!
      zeta(:,:)=zeta(:,:)+sice(:,:)
!
!     Save velocity.
!
      do k=1,ken
        do j=1,jen
          do i=1,ien
            utotold(i,j,k)=utot(i,j,k)
            vtotold(i,j,k)=vtot(i,j,k)
            wold(i,j,k)=w(i,j,k)
          end do
        end do
      end do
!
!*    4.1      Diffusion of *zeta* for fast subgrid gravity waves.
!     ---------------------------------------------------
!
!     difzet horizontal diffusion for zeta due to fast, subgrid
!     gravity waves (m**/s)
!
      difzet=1000.
      difh=difzet*dt
      do j=1,jen
        do i=1,ien
          ze2(i,j)=zeta(i,j)
        end do
      end do
!
      do j=3,jen-2
        do i=1,ien
          if (wet(i,j,1)<0.5) cycle
          il=i-1
          if (il<1) il=ien
          ir=i+1
          if (ir>ien) ir=1
          volz2=2.*ze2(i,j)*dl(j)*dphi+2.*difh*wet(il,j-1,1)            &
     &          *(ze2(il,j-1)-ze2(i,j))+2.*difh*wet(i,j-1,1)            &
     &          *(ze2(i,j-1)-ze2(i,j))+2.*difh*wet(ir,j+1,1)            &
     &          *(ze2(ir,j+1)-ze2(i,j))+2.*difh*wet(i,j+1,1)            &
     &          *(ze2(i,j+1)-ze2(i,j))
          zeta(i,j)=.5*volz2/(dl(j)*dphi)
        end do
      end do
!
!     Effect of zeta exchange on T and S at surface.
!
      do j=1,jen
        do i=1,ien
          t2(i,j)=t(i,j,1)
          s2(i,j)=s(i,j,1)
        end do
      end do
!
      do j=3,jen-2
        do i=1,ien
          if (wet(i,j,1)<0.5) cycle
          il=i-1
          if (il<1) il=ien
          ir=i+1
          if (ir>ien) ir=1
          volum=dl(j)*dphi*(ddz(i,j,1)+zeta(i,j))
          volol=dl(j)*dphi*(ddz(i,j,1)+ze2(i,j))
          gewnw=0.5*(abs(ze2(il,j-1)-ze2(i,j))+(ze2(il,j-1)-ze2(i,j)))  &
     &          *wet(il,j-1,1)*difh
          gewne=0.5*(abs(ze2(i,j-1)-ze2(i,j))+(ze2(i,j-1)-ze2(i,j)))    &
     &          *wet(i,j-1,1)*difh
          gewsw=0.5*(abs(ze2(i,j+1)-ze2(i,j))+(ze2(i,j+1)-ze2(i,j)))    &
     &          *wet(i,j+1,1)*difh
          gewse=0.5*(abs(ze2(ir,j+1)-ze2(i,j))+(ze2(ir,j+1)-ze2(i,j)))  &
     &          *wet(ir,j+1,1)*difh
          gewmm=0.5*(abs(ze2(il,j-1)-ze2(i,j))-(ze2(il,j-1)-ze2(i,j)))  &
     &          *wet(il,j-1,1)                                          &
     &          *difh+0.5*(abs(ze2(i,j-1)-ze2(i,j))-(ze2(i,j-1)-ze2(i,j)&
     &          ))*wet(i,j-1,1)                                         &
     &          *difh+0.5*(abs(ze2(i,j+1)-ze2(i,j))-(ze2(i,j+1)-ze2(i,j)&
     &          ))*wet(i,j+1,1)                                         &
     &          *difh+0.5*(abs(ze2(ir,j+1)-ze2(i,j))-(ze2(ir,j+1)       &
     &          -ze2(i,j)))*wet(ir,j+1,1)*difh
          s2(i,j)=(volol*s(i,j,1)+gewse*s(ir,j+1,1)+gewnw*s(il,j-1,1)   &
     &            +gewne*s(i,j-1,1)+gewsw*s(i,j+1,1)-gewmm*s(i,j,1))    &
     &            /(volum)
          t2(i,j)=(volol*t(i,j,1)+gewse*t(ir,j+1,1)+gewnw*t(il,j-1,1)   &
     &            +gewne*t(i,j-1,1)+gewsw*t(i,j+1,1)-gewmm*t(i,j,1))    &
     &            /(volum)
        end do
      end do
      do j=1,jen
        do i=1,ien
          s(i,j,1)=s2(i,j)
          t(i,j,1)=t2(i,j)
        end do
      end do
!     call check_nan("adv_quick.f(979):t",t,ien,jen,ken)
!
!     Smoothe in isolated points
!     with diffusion coefficient *difsmo* (m**2/s).
!
      difsmo=2000.
!
!     for aqua planet (etc.)
!
      if(naqua == 0) then
       do ia=1,5
        if (ia==1) then
          i=15
          j=25
        end if
        if (ia==2) then
          i=13
          j=25
        end if
        if (ia==3) then
          i=11
          j=25
        end if
        if (ia==4) then
          i=67
          j=29
        end if
        if (ia==5) then
          i=44
          j=12
        end if
        il=i-1
        if (il<1) il=ien
        ir=i+1
        if (ir>ien) ir=1
        do k=1,ken
          difh=difsmo*dt/ddw(k)
          if (wet(i,j,k)<0.5) cycle
          ttranw=difh*dw(k)*wet(il,j-1,k)*(t(il,j-1,k)-t(i,j,k))
          ttrane=difh*dw(k)*wet(i,j-1,k)*(t(i,j-1,k)-t(i,j,k))
          ttrasw=difh*dw(k)*wet(i,j+1,k)*(t(i,j+1,k)-t(i,j,k))
          ttrase=difh*dw(k)*wet(ir,j+1,k)*(t(ir,j+1,k)-t(i,j,k))
          stranw=difh*dw(k)*wet(il,j-1,k)*(s(il,j-1,k)-s(i,j,k))
          strane=difh*dw(k)*wet(i,j-1,k)*(s(i,j-1,k)-s(i,j,k))
          strasw=difh*dw(k)*wet(i,j+1,k)*(s(i,j+1,k)-s(i,j,k))
          strase=difh*dw(k)*wet(ir,j+1,k)*(s(ir,j+1,k)-s(i,j,k))
          if (k==1) then
            volum=(ddz(i,j,1)+zeta(i,j))*dl(j)*dphi
          else
            volum=ddz(i,j,k)*dl(j)*dphi
          end if
!
          s(i,j,k)=(volum*s(i,j,k)+stranw+strane+strasw+strase)/volum
          t(i,j,k)=(volum*t(i,j,k)+ttranw+ttrane+ttrasw+ttrase)/volum
          if (wet(il,j-1,k)>0.5) then
            if (k==1) then
              volum=(ddz(il,j-1,1)+zeta(il,j-1))*dl(j-1)*dphi
            else
              volum=ddz(il,j-1,k)*dl(j-1)*dphi
            end if
            s(il,j-1,k)=(s(il,j-1,k)*volum-stranw)/volum
            t(il,j-1,k)=(t(il,j-1,k)*volum-ttranw)/volum
          end if
          if (wet(i,j-1,k)>0.5) then
            if (k==1) then
              volum=(ddz(i,j-1,1)+zeta(i,j-1))*dl(j-1)*dphi
            else
              volum=ddz(i,j-1,k)*dl(j-1)*dphi
            end if
            s(i,j-1,k)=(s(i,j-1,k)*volum-strane)/volum
            t(i,j-1,k)=(t(i,j-1,k)*volum-ttrane)/volum
          end if
          if (wet(i,j+1,k)>0.5) then
            if (k==1) then
              volum=(ddz(i,j+1,1)+zeta(i,j+1))*dl(j+1)*dphi
            else
              volum=ddz(i,j+1,k)*dl(j+1)*dphi
            end if
            s(i,j+1,k)=(s(i,j+1,k)*volum-strasw)/volum
            t(i,j+1,k)=(t(i,j+1,k)*volum-ttrasw)/volum
          end if
          if (wet(ir,j+1,k)>0.5) then
            if (k==1) then
              volum=(ddz(ir,j+1,1)+zeta(ir,j+1))*dl(j+1)*dphi
            else
              volum=ddz(ir,j+1,k)*dl(j+1)*dphi
            end if
            s(ir,j+1,k)=(s(ir,j+1,k)*volum-strase)/volum
            t(ir,j+1,k)=(t(ir,j+1,k)*volum-ttrase)/volum
          end if
!
        end do
       end do
      endif
!
!     diagnostics
!
      dtdt(:,:,:,3)=(t(:,:,:)-ztold(:,:,:))/dt
!
!     debug diagnostic by FL
!
      if(nprint == 1) then
       write(6,*) 'End of adv_quick'
       call prfldiag
       write(6,*) 'Max dt/dt: ',maxval(dtdt(:,:,:,3))                       &
     &       ,' at ',maxloc(dtdt(:,:,:,3))
       write(6,*) 'Min dt/dt: ',minval(dtdt(:,:,:,3))                       &
     &       ,' at ',minloc(dtdt(:,:,:,3))
       call inventory
      elseif(nprint == 2) then
       write(6,*) 'At flux transfere into ml'
       write(6,*) 's: ',s(npri,nprj,nprk)
       write(6,*) 't: ',t(npri,nprj,nprk)
      endif

      return
      end subroutine adv_quick

      subroutine cont2
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
!
!     Interface.
!     ----------
!     *call* *cont1*.
!
!     ------------------------------------------------------------
!
!     Declaration of local variables.
!     -------------------------------
      integer :: i,j,k,m
      real (kind=8) :: hordiv
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
      do k=ken,2,-1
        do j=3,jen-2
          do i=1,ien
            m = mod(i+ien-2,ien) + 1
            hordiv = utot(m,j  ,k) * delta(m,j  ,k) * dphi              &
     &             - utot(i,j  ,k) * delta(i,j  ,k) * dphi              &
     &             + vtot(i,j+1,k) * delta(i,j+1,k) * dl(j+1)           &
     &             - vtot(m,j-1,k) * delta(m,j-1,k) * dl(j-1)
            w(i,j,k-1) = w(i,j,k) + hordiv / (dl(j) * dphi) * wet(i,j,k)
          end do
        end do
      end do
      return
      end subroutine cont2
!     ===================================
      subroutine dbgoprt (var,nidim,njdim,ni,nj,string)
!     ===================================

!     debug routine for ocean arrays

      implicit none 
      character(len=*), intent(in) :: string        ! info string
      real (kind=8), intent(in) :: var(nidim,njdim) ! printout variable
      integer, intent(in) :: nidim, njdim           ! dimension of var
      integer, intent(in) :: ni, nj                 ! index of printout
      integer             :: ina                    ! start index 

      if (ni > nidim) then 
         ina = nidim - 3
      else 
         ina = ni 
      endif

      write(6,'("|",a18,"(",i3,",",i3,")=",4f12.4," |")')               &
     &      trim(string),ina,nj,var(ina:ina+3,nj)
      return
      end subroutine dbgoprt
!     ==================================================================
!     ------------------------------------------------------------------
      subroutine dens
      use lsgvar
      implicit none
!
!     ------------------------------------------------------------------
!**** *dens*
!
!     by E. Maier-Reimer.
!     Last modified by U. Mikolajewicz   9/87.
!
!     Purpose.
!     --------
!     *dens* computes the actual density-field *rh* and the pot
!     density differences *rhdif* from potential temperature *t*
!     and salinity *s*. It adjusts the stratification due to
!     convective adjustment.
!
!**   Input.
!     ------
!     t       potential temperature field (in degree celsius).
!     s       salinity (in 0/00).
!        both via *common* .
!
!     sippow  difference between in-situ and pot temperature.
!        via *common* .
!
!     Output.
!     -------
!
!     t      adjusted potential temperature.
!     s      adjusted salinity.
!     rh     normalized density (non dimensional).
!     rhdif  normalized density-difference between 2 layers.
!     convad information about the appliance of convective adjustment.
!       all via *common* .
!     nrca   number of convective adjustment event per layer.
!       via *common*
!
!     Interface.
!     ----------
!     *call* *dens*.
!
!*    Externals.
!     ----------
!     *subroutine* *rhof1*
!     *function* *rho*.
!
!
!     Declaration of local variables.
!     -------------------------------

      interface
         real (kind=8) function rho(s,t,p)
            real (kind=8), intent(in) :: s,t,p
         end function rho
      end interface
 
      integer :: i,j,k,klo,il,ir,ii,jj,kk
!!FL
      integer :: jlat
      integer :: iifl(ien,jen),jjfl(ien,jen)
      real (kind=8) :: zt(ien,jen,ken),zs(ien,jen,ken)
!!FL

      real (kind=8) :: thousi,zero,tenth,hal,one,two,tenm6
      real (kind=8) :: thicko,z,delrh,tupper,tlower,supper,slower
      real (kind=8) :: thicku,tm1,tm2,zminn,thicku1,thicku2
      real (kind=8) :: tun,sun,thicko1,thicko2,tob,sob
!!FL
      real (kind=8) :: zzs1,zzs2,zsinn,zsino,zsin
!!FL
!
!     diagnostics
!
      real (kind=8) :: ztold(ien,jen,ken)
!
!*    1.     Set constants and default values.
!     ---------------------------------


      thousi=0.001
      zero=0.
      tenth=0.1
      hal=0.5
      one=1.
      two=2.
      tenm6=0.000001
      do k=1,ken
        do i=1,ien*jen
          convaa(i,k)=zero
        end do
      end do
      do k=1,ken
        nrca(k)=0
      end do
!
      ztold(:,:,:)=t(:,:,:)
!
!     Separation of ice and sea water.
!
      do j=1,jen
        do i=1,ien
          if (wet(i,j,1)<0.5) cycle
!         Protection to avoid negative layer thicknesses.
          sice(i,j)=dmin1(sice(i,j),dw(1)-2.)
          zeta(i,j)=dmax1(zeta(i,j),sice(i,j)-dw(1)+1.)
!PHCONS: regard salt as in liquid part of column only, don't do rescaling
!        to dw(1)
!          thicko=ddz(i,j,1)+zeta(i,j)
!          s(i,j,1)=(s(i,j,1)*thicko-sice(i,j)*salicc)/(thicko-sice(i,j))
          zeta(i,j)=zeta(i,j)-sice(i,j)
        end do
      end do
!
!*    2.     Computation of the actual density field.
!     ----------------------------------------
!
!
!*    2.1.   Computation for all layers except the lowest one.
!
      do k=2,ken
        z=dw(k-1)*rhonul/g*1.e-3
!
!       Calculation of in-situ temperatures.
!
        th1(:,:) = 0.0
        th2(:,:) = 0.0
        sh1(:,:) = 0.0
        sh2(:,:) = 0.0
        do j=1,jen
          do i=1,ien
            if (wet(i,j,k-1)<=0.5) cycle
            th1(i,j)=t(i,j,k-1)+sippow(k-1)
            sh1(i,j)=s(i,j,k-1)
          end do
        end do
        do j=1,jen
          do i=1,ien
            if (wet(i,j,k)<=0.5) cycle
            th2(i,j)=t(i,j,k)+sippow(k-1)
            sh2(i,j)=s(i,j,k)
          end do
        end do
!
!       Search for instabilities and appliance of the convective
!       adjustment if necessary.
!
        call rhof1(th1,sh1,z,rh1,ienjen)
        call rhof1(th2,sh2,z,rh2,ienjen)
        do j=1,jen
          do i=1,ien
            if (wet(i,j,k)<=0.5) cycle
            if (rh2(i,j)>=rh1(i,j)) cycle
            delrh=rh1(i,j)-rh2(i,j)
            tupper=t(i,j,k-1)
            tlower=t(i,j,k)
            supper=s(i,j,k-1)
            slower=s(i,j,k)
            thicku=ddz(i,j,k-1)
            if (k==2) thicku=thicku+zeta(i,j)
            if (ddz(i,j,k)>thicku) then
              t(i,j,k-1)=tlower
              s(i,j,k-1)=slower
              t(i,j,k)=tlower+(tupper-tlower)*thicku/ddz(i,j,k)
              s(i,j,k)=slower+(supper-slower)*thicku/ddz(i,j,k)
            else
              t(i,j,k)=tupper
              s(i,j,k)=supper
              t(i,j,k-1)=tupper+(tlower-tupper)*ddz(i,j,k)/thicku
              s(i,j,k-1)=supper+(slower-supper)*ddz(i,j,k)/thicku
            end if
            tm1=t(i,j,k-1)+sippow(k-1)
            rh1(i,j)=rho(s(i,j,k-1),tm1,z)
            tm2=t(i,j,k)+sippow(k-1)
            rh2(i,j)=rho(s(i,j,k),tm2,z)
            convad(i,j,k)=1.
            nrca(k-1)=nrca(k-1)+1
          end do
        end do

!
      end do
!
!     ------------------------------------------------------------------
!
!     Search for instabilities and appliance of the BBL
!     adjustment if necessary.
!
      tbott(:,:) = 0.0
      sbott(:,:) = 0.0
      do j=3,jen-2
        do i=1,ien
!     avoid t(i,j,0) for array boundary check (S. Lorenz 2005/08)
          k=ibo(i,j)
!          if (k==0) sbott(i,j)=35.  !PHCONS: zero outside of wet points
          if (k==0) cycle
          sbott(i,j)=s(i,j,k)
          tbott(i,j)=t(i,j,k)

        end do
      end do
!
      zminn=bblthick*dl(jen/2)
!
!  bblthick is a constant in namelist param
!     write(6,*) "XXX bblthick= ",bblthick
!
      z=0.
!
      call rhof1(tbott,sbott,z,rhobo,ienjen)
!
!!FL
      iifl(:,:)=0
      jjfl(:,:)=0
      zt(:,:,:)=t(:,:,:)
      zs(:,:,:)=s(:,:,:)
!!

      do j=3,jen-2
       do i=1,ien
        k=ibo(i,j)
        if (k < 2) cycle
        if (ibbl(i,j)==0) cycle
        if (wet(i,j,1)<0.5) cycle
        il=i-1
        if (il<1) il=ien
        ir=i+1
        if (ir>ien) ir=1
!
! southern hemisphere
!
        if(j > jen/2) then
! to the east and to the pol
         if (ibbl(i,j)==1.and.(rhobo(i,j)>rhobo(ir,j+1))) then
          iifl(i,j)=ir
          jjfl(i,j)=j+1
! to the east and to the equator
         else if (ibbl(i,j)==2.and.(rhobo(i,j)>rhobo(i,j-1))) then
          iifl(i,j)=i
          jjfl(i,j)=j-1
! to the west and to the equator
         else if (ibbl(i,j)==3.and.(rhobo(i,j)>rhobo(il,j-1))) then
          iifl(i,j)=il
          jjfl(i,j)=j-1
! to the west and to the pol
         else if (ibbl(i,j)==4.and.(rhobo(i,j)>rhobo(i,j+1))) then
          iifl(i,j)=i
          jjfl(i,j)=j+1
         end if
        else
!
! northern hemisphere
!
! to the east and to the pol
         if (ibbl(i,j)==1.and.(rhobo(i,j)>rhobo(i,j-1))) then
          iifl(i,j)=i
          jjfl(i,j)=j-1
! to the east and to the equator
         else if (ibbl(i,j)==2.and.(rhobo(i,j)>rhobo(ir,j+1))) then
          iifl(i,j)=ir
          jjfl(i,j)=j+1
! to the west and to the equator
         else if (ibbl(i,j)==3.and.(rhobo(i,j)>rhobo(i,j+1))) then
          iifl(i,j)=i
          jjfl(i,j)=j+1
! to the west and to the pol
         else if (ibbl(i,j)==4.and.(rhobo(i,j)>rhobo(il,j-1))) then
          iifl(i,j)=il
          jjfl(i,j)=j-1
         end if
        endif
       enddo
      enddo
!
      do jlat=3,jen/2
       j=jlat
       do i=1,ien
        k=ibo(i,j)
        ii=iifl(i,j)
        jj=jjfl(i,j)
        if(k < 2) cycle
        if(ibbl(i,j)==0) cycle
        if(wet(i,j,1)<0.5) cycle
        if(ii==0 .or. jj==0) cycle
        kk=ibo(ii,jj)
        if (wet(ii,jj,kk)<0.5) stop "false in bblm"
        thicku=ddz(ii,jj,kk)*dl(jj)       ! unten (german) = down
        thicku2=dmin1(zminn,thicku)       ! bbl is max 200 m thick
        thicku1=thicku-thicku2
        tun=t(ii,jj,kk)
        sun=s(ii,jj,kk)
        thicko=ddz(i,j,k)*dl(j)           ! oben (german) = above
        thicko2=dmin1(zminn,thicko)       ! BBL is max 200 m thick
        thicko1=thicko-thicko2
        tob=t(i,j,k)
        sob=s(i,j,k)
        if (thicko2>thicku2) then
         t(ii,jj,kk)=tob
         s(ii,jj,kk)=sob
         t(i,j,k)=tob+(tun-tob)*thicku2/thicko2
         s(i,j,k)=sob+(sun-sob)*thicku2/thicko2
        else
         t(i,j,k)=tun
         s(i,j,k)=sun
         t(ii,jj,kk)=tun+(tob-tun)*thicko2/thicku2
         s(ii,jj,kk)=sun+(sob-sun)*thicko2/thicku2
        end if
        t(ii,jj,kk)=(tun*thicku1+thicku2*t(ii,jj,kk))/thicku
        s(ii,jj,kk)=(sun*thicku1+thicku2*s(ii,jj,kk))/thicku
        t(i,j,k)=(tob*thicko1+thicko2*t(i,j,k))/thicko
        s(i,j,k)=(sob*thicko1+thicko2*s(i,j,k))/thicko
       end do
       j=jen+1-jlat
       do i=1,ien
        k=ibo(i,j)
        ii=iifl(i,j)
        jj=jjfl(i,j)
        if(k < 2) cycle
        if(ibbl(i,j)==0) cycle
        if(wet(i,j,1)<0.5) cycle
        if(ii==0 .or. jj==0) cycle
        kk=ibo(ii,jj)
        if (wet(ii,jj,kk)<0.5) stop "false in bblm"
        thicku=ddz(ii,jj,kk)*dl(jj)       ! unten (german) = down
        thicku2=dmin1(zminn,thicku)       ! bbl is max 200 m thick
        thicku1=thicku-thicku2
        tun=t(ii,jj,kk)
        sun=s(ii,jj,kk)
        thicko=ddz(i,j,k)*dl(j)           ! oben (german) = above
        thicko2=dmin1(zminn,thicko)       ! BBL is max 200 m thick
        thicko1=thicko-thicko2
        tob=t(i,j,k)
        sob=s(i,j,k)
        if (thicko2>thicku2) then
         t(ii,jj,kk)=tob
         s(ii,jj,kk)=sob
         t(i,j,k)=tob+(tun-tob)*thicku2/thicko2
         s(i,j,k)=sob+(sun-sob)*thicku2/thicko2
        else
         t(i,j,k)=tun
         s(i,j,k)=sun
         t(ii,jj,kk)=tun+(tob-tun)*thicko2/thicku2
         s(ii,jj,kk)=sun+(sob-sun)*thicko2/thicku2
        end if
        t(ii,jj,kk)=(tun*thicku1+thicku2*t(ii,jj,kk))/thicku
        s(ii,jj,kk)=(sun*thicku1+thicku2*s(ii,jj,kk))/thicku
        t(i,j,k)=(tob*thicko1+thicko2*t(i,j,k))/thicko
        s(i,j,k)=(sob*thicko1+thicko2*s(i,j,k))/thicko
       end do
      end do
!
!     ------------------------------------------------------------------
!
      rhdif(:,:,:) = 0.0
      do k=2,ken
        z=dw(k-1)*rhonul/g*1.e-3
!
!       Calculation of in-situ temperatures.
!
        th1(:,:) = 0.0
        th2(:,:) = 0.0
        sh1(:,:) = 0.0
        sh2(:,:) = 0.0
        do j=1,jen
          do i=1,ien
            if (wet(i,j,k-1)<0.5) cycle
            th1(i,j)=t(i,j,k-1)+sippow(k-1)
            sh1(i,j)=s(i,j,k-1)
          end do
        end do
        do j=1,jen
          do i=1,ien
            if (wet(i,j,k)<0.5) cycle
            th2(i,j)=t(i,j,k)+sippow(k-1)
            sh2(i,j)=s(i,j,k)
          end do
        end do
        call rhof1(th1,sh1,z,rh1,ienjen)
        call rhof1(th2,sh2,z,rh2,ienjen)
!
!       Computation of pot density-differences.
! 
        do j=1,jen
          do i=1,ien
            if (wet(i,j,k-1)<0.5) cycle
            rhdif(i,j,k-1)=(rh2(i,j)-rh1(i,j))/rhonul
            rhdif(i,j,k-1)=hal*(rhdif(i,j,k-1)+tenm6+abs(rhdif(i,j,k-1)-&
     &                     tenm6))
          end do
        end do
        if (k==2) then
!
!         Recombination of sea ice and water.
!
          do j=1,jen
            do i=1,ien
              if (wet(i,j,1)<0.5) cycle
              zeta(i,j)=zeta(i,j)+sice(i,j)
!PHCONS: don't do scaling with dw(1)
!              s(i,j,1)=(s(i,j,1)*(dw(1)+zeta(i,j)-sice(i,j))+sice(i,j)  &
!     &                 *salicc)/(dw(1)+zeta(i,j))
            end do
          end do
        end if
!
!
!       Computation of the density.
!
        th1(:,:) = 0.0
        sh1(:,:) = 0.0
        z=du(k-1)*rhonul/g*1.e-3
        do j=1,jen
          do i=1,ien
            if (wet(i,j,k-1)<0.5) cycle
            th1(i,j)=t(i,j,k-1)+sippo(k-1)
            sh1(i,j)=s(i,j,k-1)
          end do
        end do
        call rhof1(th1,sh1,z,rh1,ienjen)
        do j=1,jen
          do i=1,ien
            rh(i,j,k-1)=rh1(i,j)/rhonul
          end do
        end do
!
      end do
!
!*    2.2.    Computation of density in the lowest layer.
!
      th1(:,:) = 0.0
      sh1(:,:) = 0.0
      z=du(ken)*rhonul/g*1.e-3
      do j=1,jen
        do i=1,ien
          if (wet(i,j,ken)<0.5) cycle
          th1(i,j)=t(i,j,ken)+sippo(ken)
          sh1(i,j)=s(i,j,ken)
        end do
      end do
!
      call rhof1(th1,sh1,z,rh1,ienjen)
!
      do j=1,jen
        do i=1,ien
          rh(i,j,ken)=rh1(i,j)/rhonul
        end do
      end do
!
!     diagnostics
!
      dtdt(:,:,:,4)=(t(:,:,:)-ztold(:,:,:))/dt
      if(nprint == 1) then
       write(6,*) 'At end of dens:'
       call prfldiag
       call inventory
      endif
!
      end subroutine dens
!     ==================================================================
!     ------------------------------------------------------------------
      subroutine densin
      use lsgvar
      implicit none
!
!     ------------------------------------------------------------------
!**** *densin*
!
!     by E. Maier-Reimer.
!     Last modified by U. Mikolajewicz   9/87.
!
!     Purpose.
!     --------
!     *densin* computes the actual density-field *rh* and the
!     pot density differences *rhdif* from potential temperature *t*
!     and salinity *s*.
!
!**   Input.
!     ------
!
!     t       potential temperature field (in degree celsius).
!     s       salinity (in 0/00).
!       both via *common* .
!
!     sippow  difference between pot and in-situ temperature.
!       via common /lsgcon/.
!
!     Output.
!     -------
!
!     *rh*     normalized density(nondimensional).
!     *rhdif*  normalized density-difference between 2 layers.
!
!     all via *common* .
!
!*    Interface.
!     ----------
!
!     *call* *densin*.
!
!
!     uses the *subroutine* *rhof1*
!
!     Declaration of local variables.
!     -------------------------------
      integer :: i,j,k

      real (kind=8) :: thousi,tenm6,z
!
!*    1.     Set constants and default values.
!     ---------------------------------
      thousi=0.001
      tenm6=1.e-30
!
!*    2.     Computation of the actual density field.
!     ----------------------------------------
!
!*    2.1.   Computation for all layers except the lowest one.
!
      rhdif(:,:,:) = 0.0
      do k=2,ken
        z=dw(k-1)*rhonul/g*1.e-3
!
!       Calculation of in-situ temperatures.
!
        th1(:,:) = 0.0
        th2(:,:) = 0.0
        sh1(:,:) = 0.0
        sh2(:,:) = 0.0
        do j=1,jen
          do i=1,ien
            if (wet(i,j,k-1)<0.5) cycle
            th1(i,j)=t(i,j,k-1)+sippow(k-1)
            sh1(i,j)=s(i,j,k-1)
          end do
        end do
        do j=1,jen
          do i=1,ien
            if (wet(i,j,k)<0.5) cycle
            th2(i,j)=t(i,j,k)+sippow(k-1)
            sh2(i,j)=s(i,j,k)
          end do
        end do
!
        call rhof1(th1,sh1,z,rh1,ienjen)
        call rhof1(th2,sh2,z,rh2,ienjen)
!
!       Computation of pot density-differences  *rhdif*.
!
        
        do j=1,jen
          do i=1,ien
            if (wet(i,j,k-1)<0.5) cycle
            rhdif(i,j,k-1)=(rh2(i,j)-rh1(i,j))/rhonul
            rhdif(i,j,k-1)=dmax1(rhdif(i,j,k-1),tenm6)
          end do
        end do
!
!       Computation of the density in cgs-units.
!
        z=du(k-1)*rhonul/g*1.e-3
        do j=1,jen
          do i=1,ien
            if (wet(i,j,k-1)<0.5) cycle
            th1(i,j)=t(i,j,k-1)+sippo(k-1)
            sh1(i,j)=s(i,j,k-1)
          end do
        end do
        call rhof1(th1,sh1,z,rh1,ienjen)
        do j=1,jen
          do i=1,ien
            if (wet(i,j,k-1)<0.5) cycle
            rh(i,j,k-1)=rh1(i,j)/rhonul
          end do
        end do
      end do
!
!*    2.2.    Computation of density in the lowest layer.
!
      z=du(ken)*rhonul/g*1.e-3
      do j=1,jen
        do i=1,ien
          if (wet(i,j,ken)<0.5) cycle
          th1(i,j)=t(i,j,ken)+sippo(ken)
          sh1(i,j)=s(i,j,ken)
        end do
      end do
!
      call rhof1(th1,sh1,z,rh1,ienjen)
!
      do j=1,jen
        do i=1,ien
          if (wet(i,j,ken)<0.5) cycle
          rh(i,j,ken)=rh1(i,j)/rhonul
        end do
      end do
!
      end subroutine densin
!     ==================================================================
!     ------------------------------------------------------------------
!
      subroutine diva2
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
!                    *uvtrop2*   computes barotropic velocities *ub*,*vb*.
!                    *zetapsi*   computes surface elevation *zeta*.
!                    *geost2*    computes baroclinic velocities.
!                    *cont2*     computes vertical velocities *w*.
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
!     trup        sum of vertical transports up (m**3/s).
!     trdo         "  "     "          "     down   "   .

!
!    -------------------------------------------------------------------
!
!     Declaration of local constants/variables.
!     -----------------------------------------
!
      integer :: i,j,k
!
      real (kind=8) :: zu,zv,hordiv
!
!
!*    1.        Computation of normalized pressure.
!     -----------------------------------
!
      call press
!
!*    2.        Computation of velocities.
!     --------------------------
!
!     Compute barotropic velocities.
!
      call uvtrop2
!
!     Compute baroclinic velocities.
!
!     call geost3  ! replacement for geost2 under construction (EK)
      call geost2
!
!     predictor-corrector for sea ice/sea level
!
      call uvtrop2
      call zetapsi
!
      do j=3,jen-2
        do i=1,ien
          if (lsmv(i,j)) then
            zu =ub(i,j)-dot_product(utot(i,j,:),delta(i,j,:))/depth(i,j)
            zv =vb(i,j)-dot_product(vtot(i,j,:),delta(i,j,:))/depth(i,j)
            utot(i,j,:) = (utot(i,j,:) + zu) * wetvec(i,j,:)
            vtot(i,j,:) = (vtot(i,j,:) + zv) * wetvec(i,j,:)
          endif
        enddo
      enddo
!
!
!*    3.        Compute vertical velocities.
!     ----------------------------
!
      call cont2
!
!*    4.        Computation of control parameters.
!     ----------------------------------
!     *trup* and *trdo* describe the sum of the transports up and down.
!
      do k = 1 , ken
        trup(k) = sum(w(:,:,k)*dlh(:,:),mask=w(:,:,k)>0.0) * dphi
        trdo(k) = sum(w(:,:,k)*dlh(:,:),mask=w(:,:,k)<0.0) * dphi
      enddo
      return
      end subroutine diva2
!=======================================================================
!     ------------------------------------------------------------------
!
      subroutine do_lsg_visual
      use lsgvar
      implicit none
!
!     ------------------------------------------------------------------
!
!**** *do_lsg_visual* plots LSG fields in GUI
!
!     by JvH  01/2007.
!
!
!     Interface:
!     -----------------
!     *call* *do_lsg_visual* called by interface routine *lsgstep* of PlaSim
!
!     Type declaration.
!     -----------------
!
!      Force data transfer to GUI to 32 bit kind

       real (kind=4) :: zm,za
       real (kind=4) :: zf(ien,jen,ken)
       real (kind=4) :: zw(ien,jen,ken)

       zm = 1.0
       za = 0.0
       zw(:,:,:) = wet(:,:,:)
       zf(:,:,:) = t(:,:,:)
       call  guihorlsg("OCET" // char(0),zf,zw,ien,jen,ken,zm,za)
       zf(:,:,:) = s(:,:,:)
       call  guihorlsg("OCES" // char(0),zf,zw,ien,jen,ken,zm,za)
       return
       end subroutine do_lsg_visual
!     ==================================================================
!     ------------------------------------------------------------------
!
      subroutine geost2
      use lsgvar
      implicit none
!
!     ------------------------------------------------------------------
!
!**** *geostr*.
!
!     by E. Maier-Reimer.
!     Last modified by U. Mikolajewicz 11/93.
!
!     Purpose.
!     --------
!     computes baroclinic velocities by decomposition of
!     the field of motion into different baroclinic modes.
!     the system is solved by iteration.
!     the iteration is done *itera* times.
!
!**   Input.
!     ------
!     p           norm pressure          (common /lsgpre/).
!     taux,tauy   windstress             (common /lsgfie/).
!     rhdif       pot density difference (common /lsgfie/).
!
!     Output.
!     -------
!     utot        ) baroclinic horizontal velocity.
!     vtot        )
!        via common /lsgfie/.
!
!
!     Interface.
!     ----------
!     *call* *geost2*.
!
!     ------------------------------------------------------------------
!
!     Declaration of local constants/variables.
!     -----------------------------------------
!
      integer, save :: itera = 8

      integer :: i,j,k,l,jn1,jn2,kk,k1,kit,l1,iadv

      real (kind=8) :: dti,gewmin,dphiqq,zero,half,one,four
      real (kind=8) :: udi,vdi,udif,vdif,tko,wfac,eps

      real (kind=8) :: gew(ien,jen,ken),gewr(ien,jen),prd(ien,jen)
      real (kind=8) :: pd(ien,jen)
      real (kind=8) :: u1(ien,jen),u1r(ien,jen),u1l(ien,jen),v1(ien,jen)
      real (kind=8) :: v1r(ien,jen),v1l(ien,jen),a11(ien,jen)
      real (kind=8) :: a22(ien,jen)
      real (kind=8) :: rex(ien,jen),rey(ien,jen),utes(ien,jen)
      real (kind=8) :: vtes(ien,jen)
      real (kind=8) :: sumu(ien,jen),sumv(ien,jen),tratr(ien,jen)
      real (kind=8) :: ho(ien,jen)
      real (kind=8) :: ffh(ien,jen)
      real (kind=8) :: walt(ien,jen,ken),rhalt(ien,jen,ken)
      real (kind=8) :: gewa(ien*jen,ken)
      real (kind=8) :: u1la(ien*jen)
      real (kind=8) :: v1la(ien*jen)
      real (kind=8) :: u1ra(ien*jen),v1ra(ien*jen),gewra(ien*jen)
      real (kind=8) :: prda(ien*jen)
      real (kind=8) :: pda(ien*jen),utesa(ien*jen)
      real (kind=8) :: vtesa(ien*jen),u1a(ien*jen),v1a(ien*jen)
      real (kind=8) :: rexa(ien*jen),zx(ien*jen)
      real (kind=8) :: reya(ien*jen)
      real (kind=8) :: a11a(ien*jen),a22a(ien*jen),ffha(ien*jen)
      real (kind=8) :: hoa(ien*jen)
      real (kind=8) :: sumua(ien*jen),sumva(ien*jen)
      real (kind=8) :: tratra(ien*jen)
      real (kind=8) :: cofricl(ien,jen),cofricr(ien,jen)
      real (kind=8) :: cofricla(ienjen)
      real (kind=8) :: cofricra(ienjen)

      equivalence (cofricl(1,1),cofricla(1))
      equivalence (cofricr(1,1),cofricra(1))
!
      equivalence (gewa,gew)
      equivalence (u1la,u1l)
      equivalence (v1la,v1l)
      equivalence (u1ra,u1r)
      equivalence (v1ra,v1r)
      equivalence (gewra,gewr)
      equivalence (prda,prd)
      equivalence (pda,pd)
      equivalence (utesa,utes)
      equivalence (vtesa,vtes)
      equivalence (u1a,u1)
      equivalence (v1a,v1)
      equivalence (rexa,rex)
      equivalence (reya,rey)
      equivalence (a11a,a11)
      equivalence (a22a,a22)
      equivalence (ffha,ffh)
      equivalence (hoa,ho)
      equivalence (sumu,sumua)
      equivalence (sumva,sumv)
      equivalence (tratra,tratr)
!
!     Number of iterations .
!
      data itera/8/
!
!*    1.        Set initial constants and values.
!     ---------------------------------
!
!*    1.1.      Set initial constants.
!     ----------------------
      zero=0.
      half=0.5
      one=1.
      four=4.
      jn1=jen-1
      jn2=jen-2
      dti=1.0/dt
      gewmin=0.05
      dphiqq=dphi*dphi
      rhalt(:,:,:) = rh(:,:,:)
!
!*    1.2.      Set initial values on arrays.
!     -----------------------------
!
      t1(:,:,:) = u(:,:,:)
      s1(:,:,:) = v(:,:,:)

      u(:,:,1) = ub(:,:)
      v(:,:,1) = vb(:,:)
!
!     *ffh* Coriolis parameter.
!
      do i=1,ien
         ffh(i,:) = ff(:)
      end do

      gew(:,:,:) = 0.0
!
!
!*    2.        Computation of the weight-factors.
!     ----------------------------------
      l=jen*ien
      do k=1,kenm1
        do kk=1,k
          do i=1,l
            if (weta(i,k+1)>half) gewa(i,k)=gewa(i,k)+g*rhdifa(i,kk)    &
     &          *(dw(kk)/dw(k))
          end do
        end do
        k1=k+1
        if (k1>kenm1) cycle
        do kk=k1,kenm1
          do i=1,l
            if (dw(kk)<deppa(i)) gewa(i,k)=gewa(i,k)+g*rhdifa(i,kk)     &
     &          *weta(i,kk+1)*(deppa(i)-dw(kk))/(deppa(i)-dw(k))
          end do
        end do
      end do
      do k=1,kenm1
        do i=1,ien*jen
          gewa(i,k)=dmax1(gewa(i,k),gewmin)
        end do
      end do
!
!*    3.        Iterative solution.
!     -------------------
!
!     Main loop.
!
      do kit=1,itera
!       do 3700 kit=1,30
        udi=zero
        vdi=zero
!
!       Loop for the different depths.
!
        do k=2,ken
          if (k/=2) cycle
!
          udif=0.
          vdif=0.
          tko=zero
          if (k==2) tko=one
!
!         3.1.      Give values to the arrays describing
!         neighbour-cells.
!         -----------------------------------------------------
!
          l=ien*jen-1
          do i=2,l+1
            u1la(i)=uc(i-1,k)
            v1la(i)=vc(i-1,k)
            u1ra(i-1)=uc(i,k)
            v1ra(i-1)=vc(i,k)
            gewra(i-1)=gewa(i,k-1)
            prda(i-1)=(pa(i,k)-pa(i,k-1))
            cofricla(i)=cofrica(i-1)
            cofricra(i-1)=cofrica(i)
          end do
!
          l=l+1
          do i=1,l
            u1a(i)=uc(i,k)
            v1a(i)=vc(i,k)
            pda(i)=pa(i,k)-pa(i,k-1)
          end do
!
          do j=1,jen
            u1l(1,j)=u(ien,j,k)
            v1l(1,j)=v(ien,j,k)
            u1r(ien,j)=u(1,j,k)
            v1r(ien,j)=v(1,j,k)
            gewr(ien,j)=gew(1,j,k-1)
            prd(ien,j)=p(1,j,k)-p(1,j,k-1)
            cofricl(1,j)=cofric(ien,j)
            cofricr(ien,j)=cofric(1,j)
          end do
!
          l=ien*jen
          do i=1,l
            utesa(i)=uc(i,k)
            vtesa(i)=vc(i,k)
          end do
!
!*        3.2.      Compute the baroclinic modes.
!         -----------------------------
!
          zx(:) = 0.0
          l=(jen-4)*ien
          do i=2*ien+1,2*ien+l
            if (hreda(i,k)>half) then
              rexa(i)=dliha(i)*hreda(i,k)*(pda(i)-prda(i))-tko*tauxa(i) &
     &                +dt*dliha(i)*hreda(i,k)                           &
     &                *(gewra(i)*dliha(i)*(u1ra(i)                      &
     &                +dpin*(dlha(i-ien)*v1a(i-ien)-dlha(i+ien)         &
     &                *v1ra(i+ien)))+gewa(i,k-1)*dliha(i)               &
     &                *(u1la(i)+dpin*(dlha(i+ien)*v1a(i+ien)-dlha(i-ien)&
     &                *v1la(i-ien))))                                   &
     &                +0.5*(u1la(i-ien)*(cofrica(i)+cofricla(i-ien))    &
     &                +u1a(i-ien)*(cofrica(i)+cofrica(i-ien))+u1a(i+ien)&
     &                *(cofrica(i)+cofrica(i+ien))+u1ra(i+ien)          &
     &                *(cofrica(i)+cofricra(i+ien)))+dti*t1a(i,k)
              reya(i)=hreda(i,k)*dpin*(prda(i+ien)-pda(i-ien))          &
     &                -tko*tauya(i)+dt*dpin*hreda(i,k)                  &
     &                *(gewa(i-ien,k-1)*dliha(i-ien)                    &
     &                *(u1a(i-ien)-u1la(i-ien)                          &
     &                +dpin*(dlha(i-2*ien)*v1la(i-2*ien)))+gewra(i+ien) &
     &                *dliha(i+ien)                                     &
     &                *(u1a(i+ien)-u1ra(i+ien)+dpin*(dlha(i+2*ien)      &
     &                *v1ra(i+2*ien))))                                 &
     &                +0.5*(v1la(i-ien)*(cofrica(i)+cofricla(i-ien))    &
     &                +v1a(i-ien)*(cofrica(i)+cofrica(i-ien))+v1a(i+ien)&
     &                *(cofrica(i)+cofrica(i+ien))+v1ra(i+ien)          &
     &                *(cofrica(i)+cofricra(i+ien)))+dti*s1a(i,k)
              a11a(i)=0.5*(four*cofrica(i)+cofrica(i-ien)+cofrica(i+ien)&
     &                +cofricla(i-ien)+cofricra(i+ien))+dt*dliha(i)     &
     &                *dliha(i)*hreda(i,k)*(gewra(i)+gewa(i,k-1))+dti
              a22a(i)=0.5*(four*cofrica(i)+cofrica(i-ien)+cofrica(i+ien)&
     &                +cofricla(i-ien)+cofricra(i+ien))                 &
     &                +dt*dpin**2*dlha(i)                               &
     &                *(dliha(i-ien)*gewa(i-ien,k-1)+dliha(i+ien)       &
     &                *gewra(i+ien))*hreda(i,k)+dti
              uc(i,k)=(a22a(i)*rexa(i)+ffha(i)*reya(i))                 &
     &                /(a11a(i)*a22a(i)+ffha(i)*ffha(i))
              vc(i,k)=(a11a(i)*reya(i)-ffha(i)*rexa(i))                 &
     &                /(a11a(i)*a22a(i)+ffha(i)*ffha(i))
              vc(i,k)=half*(vc(i,k)+v1a(i))
              uc(i,k)=half*(uc(i,k)+u1a(i))
            end if
          end do
!
!         Compute control parameters.
!
          do j=1,jen
            do i=1,ien
              utes(i,j)=abs(u(i,j,k)-utes(i,j))
              vtes(i,j)=abs(v(i,j,k)-vtes(i,j))
            end do
          end do
!
!         Computation of control parameters describing the convergence.
!
          l1=ien*jen
          do i=1,l1
            udif=udif+utesa(i)
            vdif=vdif+vtesa(i)
          end do
          udi=udi+udif*float(k)
          vdi=vdi+vdif*float(k)
        end do
      end do
!
!*    4.        Computation of the baroclinic velocities.
!     -----------------------------------------
!
!*    4.1.      Set initial values.
!     -------------------
      utot(:,:,:) = 0.0
      vtot(:,:,:) = 0.0
!
      l=ien*jen
      do i=1,l
        hoa(i)=zero
      end do
!
!*    4.2.      Composition of the modes.
!     -------------------------
!
      do k=2,ken
        do i=1,l
          if (deltaa(i,k)>half) hoa(i)=hoa(i)+deltaa(i,k-1)
        end do
!
!       Contribution from the layers below.
!
        k1=k-1
        do kk=1,k1
          do i=1,l
            if (deltaa(i,kk)>half) then
              utota(i,kk)=utota(i,kk)-uc(i,k)/hoa(i)
              vtota(i,kk)=vtota(i,kk)-vc(i,k)/hoa(i)
            end if
          end do
        end do
!
!       Contribution from the layers above.
!
        do kk=k,ken
          do i=1,l
            if (deltaa(i,kk)>half) then
              utota(i,kk)=utota(i,kk)+uc(i,k)/(deptha(i)-hoa(i))
              vtota(i,kk)=vtota(i,kk)+vc(i,k)/(deptha(i)-hoa(i))
            end if
          end do
        end do
!
      end do
!
!*    4.3.      Correction for possible cutoff-errors.
!     --------------------------------------
!
!     Guarantees that the integral over the velocities with depth
!     is zero. The corrections are usually very small.
!
      l=ien*jen
      do i=1,l
        tratra(i)=zero
        sumua(i)=zero
        sumva(i)=zero
      end do
!
      do k=1,ken
        do i=1,l
          sumua(i)=sumua(i)+utota(i,k)*deltaa(i,k)
          sumva(i)=sumva(i)+vtota(i,k)*deltaa(i,k)
          tratra(i)=tratra(i)+deltaa(i,k)
        end do
      end do
!
      do i=1,l
        if (tratra(i)>half) then
          sumua(i)=sumua(i)/tratra(i)
          sumva(i)=sumva(i)/tratra(i)
        end if
      end do
!
      do k=1,ken
        do i=1,l
          utota(i,k)=utota(i,k)-sumua(i)
          vtota(i,k)=vtota(i,k)-sumva(i)
        end do
!
!       Set velocities on land-cells equal to zero.
!
        do i=1,l
          if (deltaa(i,k)<half) then
            utota(i,k)=zero
            vtota(i,k)=zero
          else
            utota(i,k)=utota(i,k)+ua(i)
            vtota(i,k)=vtota(i,k)+va(i)
          end if
        end do
      end do
!
!     Loop over predictor corrector scheme.
!
      wfac=1.0
      do iadv=1,1
!
!       Compute vertical velocities
!
        call cont2

        do k=1,ken
          do j=1,jen
            do i=1,ien
              if (iadv>1) w(i,j,k)=0.5*w(i,j,k)+0.5*walt(i,j,k)
              walt(i,j,k)=w(i,j,k)
              rh(i,j,k)=rhalt(i,j,k)
            end do
          end do
        end do
!
!       Advection.
!
        eps=1.e-20
        do k=2,ken-1
          do j=3,jen-2
            do i=1,ien
              if (wet(i,j,k)<0.5) cycle
              if (wet(i,j,k+1)<0.5.or.k==ken) rhdif(i,j,k)=0.
              rh(i,j,k)=rh(i,j,k)+wfac*0.5*(1./(ddz(i,j,k)+eps))        &
     &                  *(rhdif(i,j,k)                                  &
     &                  *dmin1(ddz(i,j,k),dt*(w(i,j,k)+abs(w(i,j,k))))  &
     &                  +rhdif(i,j,k-1)                                 &
     &                  *dmax1(-ddz(i,j,k),dt*(w(i,j,k-1)-              &
     &                  abs(w(i,j,k-1)))))
            end do
          end do
        end do
!
        do k=1,1
          do j=3,jen-2
            do i=1,ien
              rh(i,j,k)=rh(i,j,k)+wfac*0.5*(1./(ddz(i,j,k)+eps))        &
     &                  *(rhdif(i,j,k)                                  &
     &                  *dmin1(ddz(i,j,k),dt*(w(i,j,k)+abs(w(i,j,k)))))
            end do
          end do
        end do
        do k=ken,ken
          do j=3,jen-2
            do i=1,ien
              if (wet(i,j,k)<0.5) cycle
              rh(i,j,k)=rh(i,j,k)+wfac*0.5*(1./(ddz(i,j,k)+eps))        &
     &                  *rhdif(i,j,k-1)                                 &
     &                  *dmax1(-ddz(i,j,k),dt*(w(i,j,k-1)-              &
     &                  abs(w(i,j,k-1))))
            end do
          end do
        end do
!
!
!       Computation of new pressure field.
!
        call press
!
!       Iteration of the baroclinic modes with new density field.
!
!
!       Main loop.
!
        do kit=1,itera
          udi=zero
          vdi=zero
!
!         Loop for the different depths.
!
          do k=2,ken
!
            udif=0.
            vdif=0.
            tko=zero
            if (k==2) tko=one
!
!           3.1.      Give values to the arrays describing
!           neighbour-cells.
!           -----------------------------------------------------
!
            l=ien*jen-1
            do i=2,l+1
              u1la(i)=uc(i-1,k)
              v1la(i)=vc(i-1,k)
              u1ra(i-1)=uc(i,k)
              v1ra(i-1)=vc(i,k)
              gewra(i-1)=gewa(i,k-1)
              prda(i-1)=(pa(i,k)-pa(i,k-1))
            end do
!
            l=l+1
            do i=1,l
              u1a(i)=uc(i,k)
              v1a(i)=vc(i,k)
              pda(i)=pa(i,k)-pa(i,k-1)
            end do
!
            do j=1,jen
              u1l(1,j)=u(ien,j,k)
              v1l(1,j)=v(ien,j,k)
              u1r(ien,j)=u(1,j,k)
              v1r(ien,j)=v(1,j,k)
              gewr(ien,j)=gew(1,j,k-1)
              prd(ien,j)=p(1,j,k)-p(1,j,k-1)
            end do
!
            l=ien*jen
            do i=1,l
              utesa(i)=uc(i,k)
              vtesa(i)=vc(i,k)
            end do
!
!*          3.2.      Compute the baroclinic modes.
!           -----------------------------
!
            l=(jen-4)*ien
            do i=2*ien+1,2*ien+l
              if (hreda(i,k)>half) then
                rexa(i)=dliha(i)*hreda(i,k)*(pda(i)-prda(i))            &
     &                  -tko*tauxa(i)+dt*dliha(i)*hreda(i,k)            &
     &                  *(gewra(i)*dliha(i)                             &
     &                  *(u1ra(i)+dpin*(dlha(i-ien)*v1a(i-ien)          &
     &                  -dlha(i+ien)*v1ra(i+ien)))+gewa(i,k-1)*dliha(i) &
     &                  *(u1la(i)                                       &
     &                  +dpin*(dlha(i+ien)*v1a(i+ien)-dlha(i-ien)       &
     &                  *v1la(i-ien))))                                 &
     &                  +0.5*(u1la(i-ien)*(cofrica(i)+cofricla(i-ien))  &
     &                  +u1a(i-ien)*(cofrica(i)+cofrica(i-ien))         &
     &                  +u1a(i+ien)*(cofrica(i)+cofrica(i+ien))         &
     &                  +u1ra(i+ien)*(cofrica(i)+cofricra(i+ien)))      &
     &                  +dti*t1a(i,k)
                reya(i)=hreda(i,k)*dpin*(prda(i+ien)-pda(i-ien))        &
     &                  -tko*tauya(i)+dt*dpin*hreda(i,k)                &
     &                  *(gewa(i-ien,k-1)*dliha(i-ien)                  &
     &                  *(u1a(i-ien)-u1la(i-ien)                        &
     &                  +dpin*(dlha(i-2*ien)*v1la(i-2*ien)))            &
     &                  +gewra(i+ien)*dliha(i+ien)                      &
     &                  *(u1a(i+ien)-u1ra(i+ien)                        &
     &                  +dpin*(dlha(i+2*ien)*v1ra(i+2*ien))))           &
     &                  +0.5*(v1la(i-ien)*(cofrica(i)+cofricla(i-ien))  &
     &                  +v1a(i-ien)*(cofrica(i)+cofrica(i-ien))         &
     &                  +v1a(i+ien)*(cofrica(i)+cofrica(i+ien))         &
     &                  +v1ra(i+ien)*(cofrica(i)+cofricra(i+ien)))      &
     &                  +dti*s1a(i,k)
                a11a(i)=0.5*(four*cofrica(i)+cofrica(i-ien)+cofrica(i+  &
     &                  ien)+cofricla(i-ien)+cofricra(i+ien))           &
     &                  +dt*dliha(i)*dliha(i)*hreda(i,k)                &
     &                  *(gewra(i)+gewa(i,k-1))+dti
                a22a(i)=0.5*(four*cofrica(i)+cofrica(i-ien)+cofrica(i+  &
     &                  ien)+cofricla(i-ien)+cofricra(i+ien))           &
     &                  +dt*dpin**2*dlha(i)                             &
     &                  *(dliha(i-ien)*gewa(i-ien,k-1)+dliha(i+ien)     &
     &                  *gewra(i+ien))*hreda(i,k)+dti
                uc(i,k)=(a22a(i)*rexa(i)+ffha(i)*reya(i))               &
     &                  /(a11a(i)*a22a(i)+ffha(i)*ffha(i))
                vc(i,k)=(a11a(i)*reya(i)-ffha(i)*rexa(i))               &
     &                  /(a11a(i)*a22a(i)+ffha(i)*ffha(i))
                vc(i,k)=half*(vc(i,k)+v1a(i))
                uc(i,k)=half*(uc(i,k)+u1a(i))
              end if
            end do
!
!           Compute control parameters.
!
            do j=1,jen
              do i=1,ien
                utes(i,j)=abs(u(i,j,k)-utes(i,j))
                vtes(i,j)=abs(v(i,j,k)-vtes(i,j))
              end do
            end do
!
!           Computation of control parameters describing the
!           convergence.
            l1=ien*jen
            do i=1,l1
              udif=udif+utesa(i)
              vdif=vdif+vtesa(i)
            end do
            udi=udi+udif*float(k)
            vdi=vdi+vdif*float(k)
          end do
        end do
!
!*      4.        Computation of the baroclinic velocities.
!       -----------------------------------------
!
!*      4.1.      Set initial values.
!       -------------------
        utot(:,:,:) = 0.0
        vtot(:,:,:) = 0.0
!
        l=ien*jen
        do i=1,l
          hoa(i)=zero
        end do
!
!*      4.2.      Composition of the modes.
!       -------------------------
!
        do k=2,ken
          do i=1,l
            if (deltaa(i,k)>half) hoa(i)=hoa(i)+deltaa(i,k-1)
          end do
!
!         Contribution from the layers below.
!
          k1=k-1
          do kk=1,k1
            do i=1,l
              if (deltaa(i,kk)>half) then
                utota(i,kk)=utota(i,kk)-uc(i,k)/hoa(i)
                vtota(i,kk)=vtota(i,kk)-vc(i,k)/hoa(i)
              end if
            end do
          end do
!
!         Contribution from the layers above.
!
          do kk=k,ken
            do i=1,l
              if (deltaa(i,kk)>half) then
                utota(i,kk)=utota(i,kk)+uc(i,k)/(deptha(i)-hoa(i))
                vtota(i,kk)=vtota(i,kk)+vc(i,k)/(deptha(i)-hoa(i))
              end if
            end do
          end do
!
        end do
!
!*      4.3.      Correction for possible cutoff-errors.
!       --------------------------------------
!
!       Guarantees, that the integral over the velocities with depth
!       is zero. The corrections are usually very small.
!
        l=ien*jen
        do i=1,l
          tratra(i)=zero
          sumua(i)=zero
          sumva(i)=zero
        end do
!
        do k=1,ken
          do i=1,l
            sumua(i)=sumua(i)+utota(i,k)*deltaa(i,k)
            sumva(i)=sumva(i)+vtota(i,k)*deltaa(i,k)
            tratra(i)=tratra(i)+deltaa(i,k)
          end do
        end do
!
        do i=1,l
          if (tratra(i)>half) then
            sumua(i)=sumua(i)/tratra(i)
            sumva(i)=sumva(i)/tratra(i)
          end if
        end do
!
        do k=1,ken
          do i=1,l
            utota(i,k)=utota(i,k)-sumua(i)
            vtota(i,k)=vtota(i,k)-sumva(i)
          end do
!
!         Set velocities on land-cells equal to zero.
!
          do i=1,l
            if (deltaa(i,k)<half) then
              utota(i,k)=zero
              vtota(i,k)=zero
            else
              utota(i,k)=utota(i,k)+ua(i)
              vtota(i,k)=vtota(i,k)+va(i)
            end if
          end do
        end do
!
      end do

      end subroutine geost2
!     ==================================================================
!     ------------------------------------------------------------------
!
      subroutine inicoz
      use lsgvar
      implicit none
!
!     ------------------------------------------------------------------
!
!**** *inicoz*.
!
!     by E. Maier-Reimer.
!     Last modified by U. Mikolajewicz 9/87.
!     modified by Larry 08/08
!
!     Purpose.
!     --------
!     *inicoz* initializes the commonblock /lsgcon/.
!     all constants are defined by datacards.
!
!**   Input.
!     ------
!     not required.
!
!     Output.
!     -------
!     commonblock /lsgcon/.
!
!     Interface.
!     ----------
!     *call* *inicoz*.
!
!     function *rho * is used.
!
      interface
         real (kind=8) function rho(s,t,p)
            real (kind=8), intent(in) :: s,t,p
         end function rho 
      end interface

      integer :: i,j,k
      real (kind=8) :: sippoc,z
!
      g=9.80665
      rhonul=1020.
      sippoc=0.00
      tkelvin=273.15
      tref=2.
      sref=35.
      tfreez=-1.79
      cp=4000.
      frih=4.e5
      salicc=0.
      entmel=80.*4186.
!
!     Values changed according to PlaSim:
!
      tkelvin=273.16
      tfreez=-1.91
      cp=4180.
      rhonul=1030.
      entmel=80.*cp
!
!     Computation of reference-densities.
!
      do k=1,ken
        z=du(k)*rhonul/g*1.e-3
        roref(k)=rho(sref,tref,z)
      end do
!
!     Computation of the difference between true and pot temperature.
!
      do k=1,ken
        sippo(k)=du(k)/1000.*sippoc
        sippow(k)=dw(k)/1000.*sippoc
      end do
      erdof=0.
      do j=3,jen-2
        do i=1,ien
          erdof=erdof+0.5*dl(j)*dphi
        end do
      end do
!
      return
      end subroutine inicoz
!     ==================================================================
!     ------------------------------------------------------------------
!
      subroutine inigr
      use lsgvar
      implicit none
!
!     ------------------------------------------------------------------
!
!**** *inigr*.
!
!     by E. Maier-Reimer.
!     last modified by U. Mikolajewicz 9/87.
!     modified by Larry 08/08
!
!     Purpose.
!     --------
!     *inigr* initializes commonblocks /lsggri/ and /lsgvec/,
!     both containing variables only dependent on the grid and
!     the Coriolis parameter.
!
!**   Input.
!     ------
!     meridional gridstep in degrees (*grad*),
!     number of gridpoints in meridional (*jen*)and zonal (*ien*)
!     direction, number of layers (*ken*) and depth
!     of the layers (*du*) from namelist *param*.
!
!     Output.
!     -------
!     commonblocks  and .
!
!     Interface.
!     ----------
!     *call* *inigr*.
!
!     ------------------------------------------------------------------
      integer :: j,k

      real (kind=8) :: erdrad,dphi1,pi,four,erdrot,phi,pifac
!
!     Set initial values and constants.
!    
      four=4.
      pi=2.*asin(1.)
      erdrad=6371000.00
      dphi1=2.*pi*erdrad/360.
      erdrot=four*pi/86164.
!
!     Computation of mer. gridstep(*dphi*) and the inv(*dpin*).
!   
      dphi=dphi1*grad
      dpin=1./dphi
!
!     Coriolis param *ff*, zonal gridstep *dl* and inv *dli*.
!  
      pifac = pi / 360.0 * grad
      do j=1,jen
        phi = pifac * (jen/2-j+0.5)
        ff(j)=erdrot*sin(phi)
        dl(j)=dphi*cos(phi)
        dli(j)=1./dl(j)
!       write (6,*) 'INIGR: ',j,phi,sin(phi),cos(phi)
      end do
!
!     Computation of parameters, describing the vertical spacing.
! 
      do k=1,ken
        ddu(k)=du(k+1)-du(k)
        dw(k)=du(k)+ddu(k)/2.
      end do
      ddu(ken+1)=ddu(ken)
      dw(ken)=du(ken)+ddu(ken)*0.5
!
      ddw(1)=dw(1)
      do k=2,ken
        ddw(k)=dw(k)-dw(k-1)
      end do
!
!     *dl* and *dli* as 2d-array (for vectorization).
!     
      do j=1,jen
         dlih(:,j)=dli(j)
         dlh(:,j)=dl(j)
      end do
!
      return
      end subroutine inigr
!     ==================================================================
!     ------------------------------------------------------------------
!
      subroutine inipar
      use lsgvar
      implicit none
!
!     ------------------------------------------------------------------
!
!     by E. Maier-Reimer.
!     Last modified by U. Mikolajewicz 9/87.
!
!     Purpose.
!     --------
!     *inipar* reads important parameters from tapeunit *itapin*.
!
!**   Input.
!     ------
!     itapin   number of the tape to read from (common /lsgtap/).
!              type integer.
!     on local file
!     grad1    gridsize in degree.
!     phinor   latitude of the first point (northernmost).
!     dt       timestep , positive in days, negative in seconds.
!     du       array containing the depths of the various layers.
!              (in m).
!     iyear    year to stop the execution at (integer).
!     idate    date to stop the execution at (integer, format:-mmdd).
!     newst    number of timesteps to be executed in this run (integer).
!     nsmix    select the upper forcing (integer).
!     nsve     select, which method to use computation of barotropic
!              velocities (integer).
!     nscoup   select coupled or uncoupled mode.
!              (0 uncoupled , 1 coupled)
!     filauf   name of old restart file (character(len=7)).
!     filbac1  name of new restart file (character(len=7)).
!     filbac2  name of backup restart file (character(len=7)).
!     filswi   file name, control which of the files to use(character(len=7)).
!
!     Output.
!     -------
!     variables named in input plus:
!     grad     gridsize in degree.
!     dt       timestep in seconds.
!     ntyear   number of time steps /year.
!     ntcont   number of timesteps between control outputs (integer).
!     newst    number of timesteps to be executed in this run (integer).
!     nsmix    select the upper forcing (integer).
!     nscoup   select coupled or uncoupled mode.
!              (0 uncoupled , 1 coupled) (integer).
!     nsve     select, which method to use computation of barotropic
!              velocities (integer).
!     nsflu    decide, whether fluxes are used or not.
!     filpost  local file name for output data (character(len=7)).
!
!     the data are transferred via commonblocks /lsggri/,/lsgtio/
!     /lsgver/ and /lsgfil/.
!
!     Interface.
!     ----------
!     *call* *inipar*
!
!     Input parameters via namelist param.
!     ------------------------------------
      integer :: k,iyear,idate

      real (kind=8) :: grad1  =  5.0
      real (kind=8) :: phinor = 93.75
      real (kind=8) :: rmaty

      namelist /param / nsmix,nsve,nsflu,newst,ntcont,iyear,idate,      &
     &  nscoup,ntout,ntaver,ntback,grad1,phinor,dt,du,kiterm,bblthick,  &
     &  exper,ntsurf,ndiagfl,naqua,nprint,npri,nprj,nprk

!     some initial values, overwritten by input

      nscoup=1
      newst=120
      iyear=10
      idate=-1230

      if (ken == 11) then
        du(01)=25.
        du(02)=75.
        du(03)=150.
        du(04)=250.
        du(05)=450.
        du(06)=700.
        du(07)=1000.
        du(08)=2000.
        du(09)=3000.
        du(10)=4000.
        du(11)=5000.
        du(12)=7000.
      elseif (ken == 22) then
        du(01)=  25.
        du(02)=  75.
        du(03)= 125.
        du(04)= 175.
        du(05)= 225.
        du(06)= 275.
        du(07)= 350.
        du(08)= 450.
        du(09)= 550.
        du(10)= 650.
        du(11)= 750.
        du(12)= 850.
        du(13)= 950.
        du(14)=1100.
        du(15)=1300.
        du(16)=1500.
        du(17)=1800.
        du(18)=2250.
        du(19)=2750.
        du(20)=3500.
        du(21)=4500.
        du(22)=5500.
        du(23)=6500.
      endif
!
!     Read the parameters from tape itapin.
!     -------------------------------------
      read (itapin,nml=param)
      write(no6,*) ' INIPAR: Read namelist input :'
      write(no6,param)
!
!     Length of timestep taken from PlaSim [sec]
!     ------------------------------------------
      dtdays = ndcoup       ! Copy from PlaSim
      dt     = ndcoup * day ! LSG timestep [seconds]
      ntyear = 360 / ndcoup ! number of timesteps/year with 360 days
      write(no6,9131) dt,ndcoup
 9131 format(' INIPAR: DT set by PlaSim:   DT =',f11.0,                 &
     &       ' seconds   ( =',i7,' days )')
!
!  1 LSG timestep per coupling fixed
!
      grad=grad1
!
!     Write out time interval for writing matrix in years:
!     ----------------------------------------------------
!       - ntback is used for frequency of call of matrix
      rmaty=real(ntback)/real(ntyear)
      write (no6,9135) ntback,rmaty
 9135 format(' INIPAR: matrix called every NTBACK=',i7,                 &
     &       '  timesteps ( =',f7.1,' years)')
!
!     Write parameters for control purpose on file output.
!     ----------------------------------------------------
      write (no6,9102) dt
 9102 format (" ","timestep (seconds):      ",f10.2)
      write (no6,9152) 'nsmix=',nsmix,' ,nsve= ',nsve,                  &
     &              ' ,nsflu=',nsflu,' nscoup= ',nscoup
 9152 format(2x,'Selections:'/2x,4(a,i2))
      write (no6,9103) du(1)
 9103 format (" ","levels of depth (m):     ",f10.2)
      do k=2,ken1
        write (no6,9104) du(k)
 9104   format (" ",25x,f10.2)
      end do
!
!     Check switches.
!     ---------------
      if (nsmix==2.and.nsflu==1) write (no6,*)                          &
     &               " nsmix=2 and nsflu=1 is notallowed at the moment."
      if (nsmix<1.or.nsmix>3) then
        write (no6,*) " wrong selection of nsmix"
        stop "select"
      end if
      if (nsflu<0.or.nsflu>2) then
        write (no6,*) " wrong selection of nsflu."
        stop "select"
      end if
!
!     Control if *nsve* is selected correct.
!
      if (nsve==2) return
      write (no6,9876) nsve
 9876 format (1x,"nsve may have the value 2, you demanded",i6)
      stop "select"
      end subroutine inipar
!     ==================================================================
!     -----------------------------------------------------------------
!
      subroutine initop
      use lsgvar
      implicit none
!
!     -----------------------------------------------------------------
!
!**** *initop*
!
!     by E. Maier-Reimer.
!     Last modified by U. Mikolajewicz 5/88.
!     modified by Larry
!
!     Purpose.
!     --------
!     Reads the topography from a file and computes several
!     topography-dependent variables.
!     initializes commonblocks/ lsgtop/ and /lsgtrd/.
!
!**   Input.
!     ------
!     itop    number of tapeunit with topographic data. (integer).
!     - several parameters computed in subroutine inigr.
!
!     Output.
!     -------
!     /lsgtop/
!      depth      depth in u-points.
!      depp       depth in p-points.
!      wet        =0 then dry, =1 then wet, defined in scalar points.
!      wetvec     as wet,but in vector-points.
!      numx       number of the grid-points, only water-points have
!                 numbers different from zero. (integer).
!      ddz        individual layerthickness in p-points.
!      cofric     coefficient of horizontal friction /(dx*dy) (1/s).
!     /lsgtrd/
!      delta      individual layerthickness in u-points.
!      hred       reduced depth.
!     /lsgvec/.
!      depthl     depth in the next u-point to the left.
!
!     Interface.
!     ----------
!     *call* *initop*.
!
!
!     ------------------------------------------------------------------
!
!     Declaration of local variables.
!
      integer :: ilo,jla,im1,ip1,jm1,jp1
      integer :: i,j,k,l,jen1,nomsk,len,jl,il,iejeke,iter,iiter
      integer :: ifehl,ir,iffst,izz,ii
      real (kind=8) :: zero,one,half,onehalf,tmin,zlen,tmp,ztm,tp
      real (kind=8) :: zellsu,tmm,hmin
      real (kind=8) :: zhead(20)
      character(len=80) form
      character(len=1) ch*1
      character(len=1) textopo(ien,jen)

!      real :: ddepth22   ! only necessary if ken==22
!
!     *fricin* field for increased friction in straits.
!     (0 normal,>1 increased friction by factor).
!
      real (kind=8) :: fricin(ien,jen),hilf(ien,jen),hilf2(ien,jen)
      real (kind=8) :: hilf3(ien,jen)
!
!*    1. Set initial constants and values.
!   
      zero=0.
      one=1.
      half=0.5
      onehalf=1.5
      jen1=jen-1
      tmin=dw(2)+onehalf
!
!*    1.2 Open file for output of control masks:
!    
!!FL      nomsk=91
!!FL      open (nomsk,file="contmask",form="formatted")
!
!*    2. Read topography from file.
!    
      read (itop,"(a80)") form
      read (itop,form) zlen
      len=nint(zlen)
      if (len>20) then
        write (no6,*) " topography header longer than 20!"
        stop "setup"
      end if
      write (no6,'(a,f10.1)') '  zlen = ',zlen
      read (itop,"(a80)") form
      read (itop,form) (zhead(jl),jl=1,len)
      read (itop,"(a80)") form
      read (itop,form) depth
      write (no6,*) '  depth in u-points read from file'
!
!     Make sure that topography does fit to grid.
!
      hilf3(:,:)=0.
      do j=1,jen
        do i=1,ien
          if (j<=3.or.j>=jen-2) depth(i,j)=zero
          if (depth(i,j)<onehalf) cycle
          if (depth(i,j)>du(ken+1)-onehalf) depth(i,j)=du(ken+1)-onehalf
          if (depth(i,j)<dw(2)+onehalf) depth(i,j)=dw(2)+onehalf
          do k=3,ken
            if (abs(depth(i,j)-du(k))<=one) depth(i,j)=du(k)+onehalf
          end do
        end do
      end do
!
!*    2.2 Read and check *fricin* and *depp*.
!
!     read *fricin*.
!
      read (itop,form,err=2900,end=2900) fricin
!
      if(naqua == 0) then
       if (ien==72.and.ken==22) then
!
!       Bering Strait
!
        fricin(44,13)=fricin(44,13)*3.3
        fricin(43,12)=fricin(44,13)
        fricin(44,14)=fricin(44,13)/6.
!       Banda
        fricin(44,40)=dmax1(100.,fricin(44,40))
        fricin(45,40)=dmax1(100.,fricin(45,40))
        fricin(44,41)=100.
        fricin(45,41)=100.
        fricin(46,41)=100.
       end if
!
!      Mediterranean Sea
!
       do j=23,27
        do i=9,11
         fricin(i,j)=1.
        end do
       end do
      endif
!
!     Check *fricin*.
!
      do j=1,jen
        do i=1,ien
          if (fricin(i,j)<one) fricin(i,j)=one
        end do
      end do
!
!     Read *depp*.
!
      read (itop,form,end=2950,err=2950) depp
!     write (no6,*) " read depp"
      write (no6,*) '  depp in t-points read from file'
!
!     Check *depp*.
!
      do i=1,ien
        il=i-1
        if (il<1) il=ien
        do j=3,jen-2
          tmp=dmax1(depth(i,j),depth(il,j),depth(il,j-1),depth(i,j+1))
          if (tmp<=onehalf) then
            depp(i,j)=dmin1(depp(i,j),-10.)
          else
            depp(i,j)=dmax1(tmin,tmp)
          end if
        end do
        depp(i,1)=zero
        depp(i,2)=zero
        depp(i,jen-1)=zero
        depp(i,jen)=zero
      end do
!
!
      goto 3000
!
!*    2.9       Assign values to *fricin* and *depp* if not given.
!     --------------------------------------------------
 2900 continue
      do j=1,jen
        do i=1,ien
          fricin(i,j)=one
        end do
      end do
!
      if(naqua == 0) then
       fricin(44,13)=30.0
       fricin(45,40)=300.0
       fricin(45,41)=300.0
       fricin(46,41)=300.0
       write (6,*) "xxx fricin gesetzt! xxx"
       write (no6,*) '  fricin changed manually in some grid points'
      endif
      if (ken==22) then
        fricin=fricin*1.1
        write (no6,*) '   ***ken=22 -> fricin=fricin*1.1***'
      end if
!
 2950 continue
!     write (6,*) "depp nicht gelesen, sondern berechnet "
      write (no6,*) '  depp calculated from depth'
      do i=1,ien
        depp(i,1)=zero
        depp(i,jen)=zero
      end do
!
      do j=3,jen-2
        do i=1,ien
          if (i>1) il=i-1
          if (i==1) il=ien
!         Vanya--------------------------------------
          if ((depth(i,j)>0.).or.(depth(il,j)>0.).or.(depth(il,j-1)>0.) &
     &        .or.(depth(i,j+1)>0)) then
            tmp=depth(i,j)
            if (tmp<depth(il,j)) tmp=depth(il,j)
            if (tmp<depth(il,j-1)) tmp=depth(il,j-1)
            if (tmp<depth(i,j+1)) tmp=depth(i,j+1)
          else
            tmp=depth(i,j)
            if (tmp>=depth(il,j)) tmp=depth(il,j)
            if (tmp>=depth(il,j-1)) tmp=depth(il,j-1)
            if (tmp>=depth(i,j+1)) tmp=depth(i,j+1)
          end if
!         Vanya---------------------------------------
!
          depp(i,j)=tmp
          hilf(i,j)=depp(i,j)
          if (depp(i,j)<zero) depp(i,j)=zero
          if (depth(i,j)<zero) depth(i,j)=zero
        end do
      end do
!
!*    3.        Compute values for *wet* and *ddz*.
!     -----------------------------------
 3000 continue
      do k=1,ken
        do j=1,jen
          do i=1,ien
            wet(i,j,k)=zero
            tmp=depp(i,j)
            if (tmp>du(k)) wet(i,j,k)=one
            ddz(i,j,k)=ddw(k)
            if (k>1) then
              if (tmp<du(k+1)) ddz(i,j,k)=tmp-dw(k-1)
            else
              if (tmp<du(k+1)) ddz(i,j,k)=tmp
            end if
            if (tmp<du(k)) ddz(i,j,k)=zero
          end do
        end do
      end do
!
!*    4.        Compute reduced depth *hred*.
!     -----------------------------
      do k=2,ken
        do j=1,jen
          do i=1,ien
            ztm=depth(i,j)
            hred(i,j,1)=ztm
            hred(i,j,k)=zero
            tp=du(k)
            if (ztm<tp) cycle
            hred(i,j,k)=(ztm-dw(k-1))*dw(k-1)/ztm
          end do
        end do
      end do
!
!*    5.        Computation of layer thickness *delta*.
!     ---------------------------------------
      do k=1,ken
        do j=1,jen
          do i=1,ien
            ztm=depth(i,j)
            delta(i,j,k)=zero
            tp=du(k)
            if (ztm<tp) cycle
            delta(i,j,k)=ddw(k)
            if (ztm<du(k+1)) delta(i,j,k)=ztm-dw(k-1)
          end do
        end do
      end do
!!FL        open(unit=23,file='delta.dat')
      do k=1,ken
        do j=1,jen
          do i=1,ien
!!FL            write(23,*) delta(i,j,k) 
            if (delta(i,j,k)>0.1) then
              wetvec(i,j,k)=1.
            else
              wetvec(i,j,k)=0.
            end if
          end do
        end do
      end do
!!FL        close(23)
!
!*    6.        Count and number the cells.
!     ---------------------------
!
!*    6.1.      Count the wet cells (3-d).
!     --------------------------
      iejeke=ien*jen*ken
      zellsu=0.
      do i=1,iejeke
        zellsu=zellsu+wetb(i)
      end do
      write (no6,9811) zellsu
!
 9811 format ("  sum of wet ",f12.0)
!
!*    6.2.      Number the sea-points (2-d).
!     ----------------------------
      numx(:,:) = 0
      num(:,:)  = 0
      lsmv(:,:)  = .false.

      l=0
      do j=3,jen-2
        do i=1,ien
          if (depth(i,j)<=zero) cycle
          l=l+1
          numx(i,j)=l
          lsmv(i,j) = .true.
        end do
      end do

      l=0
      do j=3,jen-2
        do i=1,ien
          if (depp(i,j)<=zero) cycle
          l=l+1
          num(i,j)=l
        end do
      end do
      nweten=l
!
!
!*    7.        Computation of *dephtl*.
!     ------------------------
      l=ien*jen-1
      do i=2,l+1
        depthla(i)=deptha(i-1)
      end do
      do j=1,jen
        depthl(1,j)=depth(ien,j)
      end do
!
!*    8.        Computation of frction *cofric*.
!     --------------------------------
      do j=1,jen
        do i=1,ien
          cofric(i,j)=frih*dpin*dlih(i,j)*fricin(i,j)
        end do
      end do
!
      do jla = 2 , jen - 1
        jm1 = jla - 1
        jp1 = jla + 1
        do ilo = 1 , ien
          im1 = mod(ilo+ien-2,ien) + 1
          ip1 = mod(ilo      ,ien) + 1
          fricne(ilo,jla) = 0.5 * (cofric(ilo,jla) + cofric(ilo,jm1))
          fricse(ilo,jla) = 0.5 * (cofric(ilo,jla) + cofric(ip1,jp1))
          fricsw(ilo,jla) = 0.5 * (cofric(ilo,jla) + cofric(ilo,jp1))
          fricnw(ilo,jla) = 0.5 * (cofric(ilo,jla) + cofric(im1,jm1))
        enddo
      enddo
!     
      close (itop)
!
!*    9.         Compute parameters for runoff-model
!     -----------------------------------
      write (no6,*) '  compute parameters for runoff-model '
      do j=1,jen
        do i=1,ien
!         new for runoff
          hilf(i,j)=-hilf(i,j)
!
          if (j<=2.or.j>=jen-1) hilf(i,j)=1.e6
          iriv(i,j)=0                   ! iriv = 0  over the oceans
        end do
      end do
!
      iiter=0
 9019 continue
      iiter=iiter+1
      ifehl=0
!
      do j=3,jen-2
        do i=1,ien
          hilf2(i,j)=hilf(i,j)
          if (hilf(i,j)<0.) cycle
          il=i-1
          if (il<1) il=ien
          ir=i+1
          if (ir>ien) ir=1
          tmp=dmin1(hilf(il,j-1),hilf(i,j-1),hilf(i,j+1),hilf(ir,j+1))
          tmm=dmax1(hilf(il,j-1),hilf(i,j-1),hilf(i,j+1),hilf(ir,j+1))
          if (tmp>=hilf(i,j)-0.001) then
            ifehl=ifehl+1
            hilf2(i,j)=0.5*(tmp+hilf(i,j)+0.02*tmm)+1.000001
          end if
        end do
      end do
!
      do j=3,jen-2
        do i=1,ien
          hilf(i,j)=hilf2(i,j)
        end do
      end do
!             write (6,*) " ifehl: ",ifehl
              if (iiter==1) iffst=ifehl
              if (ifehl/=0) goto 9019           ! loop until there are no sinks
      if (iiter>400) write(6,*) "iiter too high "
      write (6,'(a,i5,a,i5)') "   ifehl=",iffst,"    at iiter=    1"
      write (6,'(a,i5,a,i5)') "   ifehl=",ifehl," after iiter=",iiter
!
      do j=3,jen-2
        do i=1,ien
          hilf3(i,j)=hilf(i,j)
          if (hilf(i,j)<0.) cycle       ! for ocean points
          il=i-1
          if (il<1) il=ien
          ir=i+1
          if (ir>ien) ir=1
          hmin=1.e9
          if (hilf(i,j-1)<hmin) then    ! to N
            iriv(i,j)=1
            hmin=hilf(i,j-1)
          end if
          if (hilf(ir,j+1)<hmin) then   ! to S
            iriv(i,j)=2
            hmin=hilf(ir,j+1)
          end if
          if (hilf(i,j+1)<hmin) then    ! to W
            iriv(i,j)=3
            hmin=hilf(i,j+1)
          end if
          if (hilf(il,j-1)<hmin) iriv(i,j)=4 ! to E
!
!
        end do
      end do
!!FL      write (nomsk,*) " IRIV"
!!FL      write (nomsk,"(1x,72i1)") iriv
!     ------------------------------------------------------------------

!
!*    10.         Compute parameters for BBL-model
!     -----------------------------------
      write(6,*) "  "
      write(6,*) " compute parameters for bbl-model "
!
      do j=1,jen
        do i=1,ien
          ibo(i,j)=0
          ibo2(i,j)=0
          if (wet(i,j,1)<0.5) goto 10012
          do k=1,ken
            if (wet(i,j,k)>0.5) ibo(i,j)=ibo(i,j)+1
          end do
          ibo2(i,j)=ibo(i,j)-2
10012     continue
          hilf(i,j)=-depp(i,j)          ! this is the ocean
          if (j<=2.or.j>=jen-1) hilf(i,j)=1.e6*1.0
          if (depp(i,j)<=25.) hilf(i,j)=1.e6*1.0
          ibbl(i,j)=0                   ! land
        end do
      end do
!
      izz=ken
      iiter=0
!
10019 continue
      iiter=iiter+1
      izz=izz-nint(real(iiter)/250.)
      ifehl=0
!
      do j=3,jen-2
        do i=1,ien
          hilf2(i,j)=hilf(i,j)
          if (hilf(i,j)>=1.e5) cycle
          il=i-1
          if (il<1) il=ien
          ir=i+1
          if (ir>ien) ir=1
          tmp=dmin1(hilf(il,j-1),hilf(i,j-1),hilf(i,j+1),hilf(ir,j+1))
          tmm=dmax1(hilf(il,j-1),hilf(i,j-1),hilf(i,j+1),hilf(ir,j+1))
          if (ibo(i,j)<izz) then
            if (tmp>=hilf(i,j)-0.001) then
              ifehl=ifehl+1
              if (0.5*(hilf(i,j)+tmp)+real(ibo(i,j))+110.<0.) hilf2(i,j)&
     &            =0.5*(hilf(i,j)+tmp)+real(ibo(i,j))+100.
              if (tmm>0.) then
                if (0.5*hilf(i,j)+120.<0.) hilf2(i,j)=0.5*hilf(i,j)+120.
              end if
            end if
          end if
        end do
      end do
!
      do j=3,jen-2
        do i=1,ien
          hilf(i,j)=hilf2(i,j)
        end do
      end do

      if (iiter<500.and.ifehl/=0) goto 10019    ! loop until there are no sinks
!
      do j=3,jen-2
        do i=1,ien
          if (hilf(i,j)>=1.e5) then
            ibbl(i,j)=iriv(i,j)
            ibbl(i,j)=0
            cycle       ! for land points
          end if
          il=i-1
          if (il<1) il=ien
          ir=i+1
          if (ir>ien) ir=1
          hmin=1.e9
          hmin=1500.
!
! southern hemisphere
!
          if(j>jen/2) then
           if (hilf(i,j-1)<hmin) then   ! to E and to Eq.
            ibbl(i,j)=2
            hmin=hilf(i,j-1)
           endif
           if (hilf(il,j-1)<hmin) then  ! to W and to Eq.
            ibbl(i,j)=3
            hmin=hilf(il,j-1)
           end if
           if (hilf(i,j+1)<hmin) then   ! to W and to Pol
            ibbl(i,j)=4
            hmin=hilf(i,j+1)
           end if
           if (hilf(ir,j+1)<hmin) then  ! to E and to Pol
            ibbl(i,j)=1
            hmin=hilf(ir,j+1)
           end if
          else
!
! southern hemisphere
!
           if (hilf(ir,j+1)<hmin) then   ! to E and to Eq.
            ibbl(i,j)=2
            hmin=hilf(ir,j+1)
           endif
           if (hilf(i,j+1)<hmin) then    ! to W and to Eq.
            ibbl(i,j)=3
            hmin=hilf(i,j+1)
           end if
           if (hilf(il,j-1)<hmin) then   ! to W and to Pol
            ibbl(i,j)=4
            hmin=hilf(il,j-1)
           end if
           if (hilf(i,j-1)<hmin) then    ! to E and to Pol
            ibbl(i,j)=1
            hmin=hilf(i,j-1)
           end if
          endif
          if (hilf(i,j)<hmin) ibbl(i,j)=5
        end do
      end do
!
!!FL      write (nomsk,*) " IBO2"
!!FL      write (nomsk,"(1x,72i1)") ibo2
!!FL      write (nomsk,*) " IBBL"
!!FL      write (nomsk,"(1x,72i1)") ibbl
!
!!FL      do j=1,jen
!!        do i=1,ien
!!          if (depp(i,j)<=0.) then
!!            textopo(i,j)="#"
!!          else
!!            ii=int(depp(i,j)/1000.)
!!            write (ch(1:1),"(i1)") ii
!!            textopo(i,j)=ch
!!          end if
!!        end do
!!FL      end do
!
!!FL      write (nomsk,*) "  "
!!FL      write (nomsk,*) " ocean topography: "
!!FL      write (nomsk,54642) (mod(i,10),i=1,ien)
54642 format (13x,72i1)
!
!!FL      write (nomsk,*) " "
!!FL      do j=3,jen-4
!!FL        write (nomsk,54596) j-2,j,(textopo(i,j),i=1,ien)
!!FL      end do
!!FL      write (nomsk,*) " "
!!FL      write (nomsk,54641) (mod(i,10),i=1,ien)
!
!!FL      write (nomsk,*) "  "
!
!!FL      do j=1,jen
!!        do i=1,ien
!!          if (hilf3(i,j)<=0) then
!!           textopo(i,j)="#"
!!          else
!!            ii=int(hilf3(i,j)/1000.)
!!            write (ch(1:1),"(i1)") ii
!!            textopo(i,j)=ch
!!          end if
!!        end do
!!FL      end do
!
!!FL      write (nomsk,*) "  "
!!FL      write (nomsk,*) " atmospheric orography: "
!!FL      write (nomsk,54641) (mod(i,10),i=1,ien)
!!FL      write (nomsk,*) " "
!
!!FL      do j=3,jen-4
!!FL        write (nomsk,54596) j-2,j,(textopo(i,j),i=1,ien)
!!FL      end do
!
!!FL      write (nomsk,*) " "
!!FL      write (nomsk,54641) (mod(i,10),i=1,ien)
!!FL      write (nomsk,*) " "
!
!     ------------------------------------------------------------------
!
!!FL      write (nomsk,*)
!!FL      write (nomsk,*) 'WET:'
!!FL      write (nomsk,'(72i1)') int(wet(:,:,1)+0.1)
!!FL      write (nomsk,*) 'WETVEC:'
!!FL      write (nomsk,'(72i1)') int(wetvec(:,:,1)+0.1)

!!FL      write (nomsk,*)
!!FL      write (nomsk,*) '   K      DU     DDU      DW     DDW'
!!FL      write (nomsk,*) '------------------------------------'
!!FL      write (nomsk,'(i5,4f8.1)') (k,du(k),ddu(k),dw(k),ddw(k),k=1,ken)
!!FL      write (nomsk,*)

!     insl=nomsk

      return

54641 format (13x,72i1)
54596 format (1x,i4,2x,i4,2x,72a1)
!
      end subroutine initop
!     ==================================================================
!     ------------------------------------------------------------------
!
      subroutine inival2(noread,yfile)
      use lsgvar
      implicit none
!
!     ------------------------------------------------------------------
!
!**** *inival2*.
!
!     by U. Mikolajewicz 12/87.
!
!     Purpose.
!     --------
!     *inival2* gives initial values to variables.
!     reads initial values for *nt*,*u*,*v*,*zeta*,*ub*,*vb*,
!     *t* and *s* from tapeunuit *noread*.
!     The subroutine *outback* creates an file, that can be used
!     as input file for the next run.
!
!**   Input.
!     ------
!     noread   number of the tape to read from.  (integer)
!     yfile    name of the file to read from (created by *outback*).
!
!     Output.
!     -------
!     /lsgfie/.
!      t        pot temperature.
!      s        salinity.
!      ub,vb    barotropic velocities.   )
!     /lsgh1/.                           )  all equal zero.
!      u,v      baroclinic modes.        )  if abnormal
!     /lsgsur/.                           )  startup.
!      zeta     surface elevation.       )
!     /lsgver/.
!      nt0      number of timesteps computed in former runs.
!                (integer).
!     Interface.
!     ----------
!     *call* *inival2*(*noread*,*yfile*).
!
!     ------------------------------------------------------------------
!
!     Declaration of local variables.
!
      integer :: i,j,k,noread,kcdat,kctim,ksdat,kstim,iyy,imm,idd
      integer :: kk,k1
      real (kind=8) :: zero,half,zmax
      character(len=7) yfile
      real (kind=8) :: zhelp(ien,jen)
      real (kind=8) :: zprel,zlen,zcode,zlev,zdum
      real (kind=8) :: ho(ien,jen)
!
!*    1.        Set initial value and open file.
!     --------------------------------
!               Data are all kind=8 (64bit)
       write(6,*) 'noread=',noread
      zero=0.0
      open (noread,file=yfile,access="sequential",form="unformatted")
!
!*    2.        Normal input data.
!     ------------------
!
!     2.1       Read *nddr* and *oddr*.
!
      rewind noread
      read (noread) nddr

!  machine independent: do not use integers for character output:

      nddr(300:399) = 0
      read (noread) oddr
!
!     Writes actual switches on nddr.
!
      nddr(400)=nsve
      nddr(401)=nsmix
      nddr(402)=nsflu
      nddr(403)=ntyear
!
!     Computes actual date and time.
!
      call actdate(kcdat,kctim)
      nddr(8)=-1    ! kcdat
      nddr(7)=-1    ! kctim
!
!     Store simulation year and date (fetched from PlaSim)
!
      ksdat = mdatim(1) *1000 + mdatim(2) * 100 + mdatim(3)
      kstim = 0
      nddr( 9) = -1  ! nddr(11)
      nddr(10) = -1  ! nddr(12)
      nddr(11) = ksdat
      nddr(12) = kstim
!
!     Compute current time step from PlaSim date
!
      iyy = 360 * (mdatim(1) - n1styear) 
      imm =  30 * (mdatim(2) -        1)
      idd =        mdatim(3) -        1
      nt  = (iyy + imm + idd) / dt
      nt0 = nt

      write (no6,*) 'INIVAL2: Read input from <',trim(yfile),'>'
      write (no6,*) 'INIVAL2: Current date:',kcdat,kctim
      write (no6,*) 'INIVAL2: Sim     date:',ksdat,kstim
      write (no6,*) 'INIVAL2: Timestep    :',nt
!
!*    prepare for write in kleiauf: update nddr
!
!     Check depth distribution.
!
      half=.5
      do k=1,ken+1
        if (abs(du(k)-oddr(9+k))>half) then
          write (no6,*) " mismatch in depth distribution  "
          write (no6,9811) (kk,du(kk),oddr(kk+9),kk=1,ken+1)
 9811     format (1x,i3,".layer, namelist:",g10.3," header of data:",   &
     &            g10.3)
          stop "data 1"
        end if
      end do
!
!     Check dimensions.
!
      if (nddr(18)/=ien*jen) then
        write (no6,*) " mismatch in dimensions, ien*jen=",ien*jen,",",  &
     &                nddr(18)
        stop "data 2"
      end if
      if (nddr(19)/=ien) then
        write (no6,*) " mismatch in dimensions, ien=",ien,",",nddr(19)
        stop "data 3"
      end if
      if (nddr(20)/=jen) then
        write (no6,*) " mismatch in dimensions, jen=",jen,",",nddr(20)
        stop "data 4"
      end if
      if (nddr(21)/=ken) then
        write (no6,*) " mismatch in dimensions, ken=",ken,",",nddr(21)
        stop "data 5"
      end if
!
!     Read temperature field.
!
      write (no6,*) " read temperature"
      do k=1,ken
        read (noread) zprel,zlen,zcode,zlev,zdum,zdum
        if (abs(zcode+62.)>0.01.and.abs(zcode+2)>0.01) then
          write (no6,*) " wrong codenumber, -62 expected"
          stop "data 6"
        end if
        read (noread) zhelp
        t(:,:,k) = zhelp(:,:)
!!
!!        do j=1,jen
!!          do i=1,ien
!!            if (t(i,j,k)>100.) t(i,j,k)=t(i,j,k)-tkelvin
!!
!  no manipulation of t(i,j,k)
!           if (t(i,j,k)<tfreez) t(i,j,k)=tfreez
!!
!!          end do
!!        end do
      end do
!
!     Read salinity field.
!
      write (no6,*) " read salinity"
      do k=1,ken
        read (noread) zprel,zlen,zcode,zlev,zdum,zdum
        if (abs(zcode+5.)>0.01) then
          write (no6,*) " wrong codenumber, -5  expected"
          stop "data 7"
        end if
        read (noread) zhelp
        s(:,:,k) = zhelp(:,:)
!
! Attention - no more salinity correction here - warnings in STEP (2005/06)
!           if (s(i,j,k)<5.) s(i,j,k)=34.5
!           if (s(i,j,k)<20.) s(i,j,k)=34.5
      end do
!
!     Initialize velocities and baroclinic modes.
!
      do k=1,ken
        do j=1,jen
          do i=1,ien
            u(i,j,k)=0.
            v(i,j,k)=0.
            utot(i,j,k)=0.
            vtot(i,j,k)=0.
            w(i,j,k)=0.
          end do
        end do
      end do
!
!     Read baroclinic modes.
!
      write (no6,*) " read baroclinic modes"
      do k=1,ken
        read (noread,err=3000,end=3000) zprel,zlen,zcode,zlev,zdum,zdum
        read (noread) zhelp
        if (abs(zcode+63)>0.01) goto 2149
        u(:,:,k) = zhelp(:,:)
      end do
      do k=1,ken
        read (noread,err=3000,end=3000) zprel,zlen,zcode,zlev,zdum,zdum
        if (abs(zcode+64)>0.01) goto 3000
        read (noread) zhelp
        v(:,:,k) = zhelp(:,:)
      end do
      do k=1,ken
        do j=1,jen
          do i=1,ien
            if (delta(i,j,k)<=0.5) then
              u(i,j,k)=zero
              v(i,j,k)=zero
            end if
          end do
        end do
      end do
!
!     Read barotropic velocities.
!
 2149 continue
      write (no6,*) " read barotropic velocities"
      read (noread,err=3000,end=3000) zprel,zlen,zcode,zlev,zdum,zdum
      read (noread) zhelp
      if (abs(zcode+37)>0.01) then
        write (no6,*) "code-number -37 expected, found ",zcode
        goto 2149
      end if
      ub(:,:) = zhelp(:,:)
      read (noread,err=3000,end=3000) zprel,zlen,zcode,zlev,zdum,zdum
      if (abs(zcode+38)>0.01) goto 3000
      read (noread) zhelp
      vb(:,:) = zhelp(:,:)
      do j=1,jen
        do i=1,ien
          if (depth(i,j)>0.5) cycle
          ub(i,j)=zero
          vb(i,j)=zero
        end do
      end do
!
!     Read surface elevation.
!
      write (no6,*) " read surface elevation"
      read (noread,err=3000,end=3000) zprel,zlen,zcode,zlev,zdum,zdum
      if (abs(zcode+1)>0.01) goto 3000
      read (noread) zhelp
      zmax=10.
      do j=1,jen
        do i=1,ien
          zeta(i,j)=zhelp(i,j)
!!          if (abs(zeta(i,j))>zmax) then
!            write (6,*)                                                 &
!     &                " sea level field is extreme in startup file!"
!            write (6,*) " sea level: i=",i," j=",j," zeta=",zeta(i,j)
!            zeta(i,j)=zeta(i,j)*0.90
!            zeta(i,j)=dmax1(zeta(i,j),-du(1))
!            zeta(i,j)=dmax1(zeta(i,j),-zmax)
!            zeta(i,j)=dmin1(zeta(i,j),zmax)
!!          end if
        end do
      end do
!
!     Read ice thickness.
!
      write (no6,*) " read ice thickness"
      read (noread,err=3000,end=3000) zprel,zlen,zcode,zlev,zdum,zdum
      if (abs(zcode+13.)>0.01) goto 3000
      read (noread) zhelp
      sice(:,:) = zhelp(:,:)
!
!
!     Read data for runoff model.
!
      write (no6,*) " read soil storage (runoff model)"
      read (noread,err=3330,end=3000) zprel,zlen,zcode,zlev,zdum,zdum
      if (abs(zcode+88.)>0.01) goto 3000
      read (noread) zhelp
      stor(:,:) = zhelp(:,:)
!
      write (no6,*) " read runoff (riv)"
      read (noread,err=3000,end=3000) zprel,zlen,zcode,zlev,zdum,zdum
      if (abs(zcode+89.)>0.01) goto 3000
      read (noread) zhelp
      riv(:,:) = zhelp(:,:)
!
!     added for perfect restart:
!
!     a) dummy (unused depth) 
!
      read (noread,err=3000,end=3000) zprel,zlen,zcode,zlev,zdum,zdum
      write(no6,*) 'read dummy, code= ',zcode
      if (abs(zcode+99)>0.01) then
       write(no6,*) 'old restart file found, restart completed'
       goto 4000    
      endif
      read (noread) zhelp
!
!     b) Read 'old' velocities for adv_quick
!
      do k=1,ken
        read (noread,err=3330,end=3330) zprel,zlen,zcode,zlev,zdum,zdum
        read (noread) zhelp
        if (abs(zcode+3)>0.01) then
         write(no6,*) 'old restart found, code= ',zcode,' continue'
         goto 4000
        endif
        utotold(:,:,k) = zhelp(:,:)
      end do
      do k=1,ken
        read (noread,err=3330,end=3330) zprel,zlen,zcode,zlev,zdum,zdum
        if (abs(zcode+4)>0.01) then
         write(no6,*) 'ERROR: found code ',zcode,' searched for 4'
         stop
        endif
        read (noread) zhelp
        vtotold(:,:,k) = zhelp(:,:)
      end do
      do k=1,ken
        read (noread,err=3330,end=3330) zprel,zlen,zcode,zlev,zdum,zdum
        if (abs(zcode+7)>0.01) then
         write(no6,*) 'ERROR: found code ',zcode,' searched for 7'
         stop
        endif
        read (noread) zhelp
        wold(:,:,k) = zhelp(:,:)
      end do
      write(no6,*) 'old velocities read'
 3330 continue
      goto 4000         ! start without error
!
!
!*    3.        Startup with errors in reading startup file.
!     --------------------------------------------
 3000 continue
      write (no6,*) " error in reading tapeunit ",noread,               &
     &              "abnormal start"
      inorm=1

      do k=1,ken
        do j=1,jen
          do i=1,ien
            u(i,j,k)=zero
            v(i,j,k)=zero
          end do
        end do
      end do
      do j=1,jen
        do i=1,ien
          ub(i,j)=zero
          vb(i,j)=zero
          zeta(i,j)=zero
          sice(i,j)=zero
          stor(i,j)=stomax
          riv(i,j)=0.
        end do
      end do
!
!
!*    4.        Initialisation of fluxes.
!     -------------------------
 4000 continue
      do k=1,ken
        do j=1,jen
          do i=1,ien
            utot(i,j,k)=0.
            vtot(i,j,k)=0.
            w(i,j,k)=0.
            fluhea(i,j)=zero
            fluwat(i,j)=zero
          end do
        end do
      end do
      do j=1,jen
        do i=1,ien
          ho(i,j)=zero
        end do
      end do
!
!*    4.2.      Composition of the modes.
!     -------------------------
!
      do k=2,ken
        do j=1,jen
          do i=1,ien
            if (delta(i,j,k)>half) ho(i,j)=ho(i,j)+delta(i,j,k-1)
          end do
        end do
!
!       Contribution from the layers below.
!
        k1=k-1
        do kk=1,k1
          do j=1,jen
            do i=1,ien
              if (delta(i,j,kk)>half) then
                utot(i,j,kk)=utot(i,j,kk)-u(i,j,k)/ho(i,j)
                vtot(i,j,kk)=vtot(i,j,kk)-v(i,j,k)/ho(i,j)
              end if
            end do
          end do
        end do
!
!       Contribution from the layers above.
!
        do kk=k,ken
          do j=1,jen
            do i=1,ien
              if (delta(i,j,kk)>half) then
                utot(i,j,kk)=utot(i,j,kk)+u(i,j,k)/(depth(i,j)-ho(i,j))
                vtot(i,j,kk)=vtot(i,j,kk)+v(i,j,k)/(depth(i,j)-ho(i,j))
              end if
            end do
          end do
        end do
!
      end do
!
!*    5.        Close *yfile*.
!     --------------
      close (noread)
!
      end subroutine inival2
!     ==================================================================
!
      subroutine inp
!
!     ------------------------------------------------------------------
!
!**** *inp*.
!
!     by E. Maier-Reimer.
!     Last modified by U. Mikolajewicz 9/87.
!
!     Purpose.
!     --------
!     *inp* controls initial input.
!     Calls subroutines *inicoz*,*inigr*,*inipar* and *initop*.
!
!**   Input.
!     ------
!     itop   number of tapeunit with topographic data. )  integer
!     itapin                         inputparameters.  )
!
!     Output.
!     -------
!     initializes global variables in lsgvar
!
!     Interface.
!     ----------
!     *call* *inp*.
!
!     -----------------------------------------------------------------
!
!
!*    1.       Initialisation of parameters.
!     -----------------------------
      call inipar
!
!*    2.       Initialisation of grid-dependent variables.
!     -------------------------------------------
      call inigr
!
!*    3.       Initialisation of constants.
!     ----------------------------
      call inicoz
!
!*    4.       Input of topography.
!     --------------------
      call initop
!
!*    5.       Input of advection-diffusion.
!     ----------------------------
      call set_quick
!

      end subroutine inp
!     =================
!     subroutine LSGINI
!     =================
!
      subroutine lsgini(kdcoup,kcpl,kdbox,kdboy)
      use lsgvar
!
      real :: zlsms(ien,jen) = 0. ! land sea mask (0/1) for scalars 
      real :: zlsmv(ien,jen) = 0. ! land sea mask (0/1) for vectors
!
!*    Print Version Number
!
      write(6,*)
      write(6,*) '**************************************************'
      write(6,*) '* Large Scale Geostrophic (LSG) Ocean Model      *'
      write(6,*) '* ---------------------------------------------- *'
      write(6,*) '* Author        : Ernst Maier-Reimer     (MPIfM) *'
      write(6,*) '* Modified by   : Uwe  Mikolajewicz      (MPIfM) *'
      write(6,*) '* Mixg     by   : Gerrit Lohmann      (AWI-Brhv) *'
      write(6,*) '* QuickAdv by   : Andre Paul             (UniHB) *'
      write(6,*) '* A-O coupling  : Stephan Lorenz         (MPIfM) *'
      write(6,*) '*    and        : Jost von Hardenberg (ISAC-CNR) *'
      write(6,*) '* PlaSim version: Frank Lunkeit       (MI UniHH) *'
      write(6,*) '*    and        : Edilbert Kirk       (MI UniHH) *'
      write(6,*) '**************************************************'
      write(6,*)
!
!     copy input from cpl
!
      ndcoup = kdcoup   ! LSG timestep [days]
      ndbox  = kdbox    ! dbug 
      ndboy  = kdboy    ! dbug
      ncoupling = kcpl  ! coupling method
!
!     Open and rewind tapeunits.
!
      call opfile
!
!     Initialize program.
!    
      call start
!
      if (nddr(404)==0) then
       nddr(404)=-9
       write (6,*) " Attention, no random number: nddr(404)=",nddr(404)
      end if
      write (6,*) " nddr(404) ",nddr(404)
!
      nt=nt0
      rewind 79
      do j=1,jen
        do i=1,ien
          monit(i,j)=0
        end do
      end do
      rewind 79
      if (ien==72) then
        read (79,7900) ((monit(i,j),i=1,ien),j=3,jen-2)
      else if (ien==64) then
        read (79,"(1x,64i1)") ((monit(i,j),i=1,ien),j=3,jen-2)
      else
        read (79,"(1x,128i1)") ((monit(i,j),i=1,ien),j=3,jen-2)
      end if
 7900 format (1x,72i1)
!
!     monit has to be changed
!
      write(6,*) "Modified beklas."
      do j=1,jen
        do i=1,ien
          moni3(i,j)=monit(i,j)
          if (depp(i,j)<20.) monit(i,j)=0
          if (depp(i,j)>20..and.monit(i,j)==0.and.j>9) monit(i,j)=1
          if (depp(i,j)>20..and.monit(i,j)==0.and.j<=9) monit(i,j)=4
          irbek2(i,j)=monit(i,j)
          moni4(i,j)=monit(i,j)
          if (j>39.and.monit(i,j)==1) moni4(i,j)=5
          if (j>39.and.monit(i,j)==2) moni4(i,j)=6
        end do
      end do
!
      monn=(jen-4)/3+2
      mons=56   ! Southern Ocean is now defined as southward of 45S
      do j=1,jen
        do i=1,ien
          moni2(i,j)=monit(i,j)
          if (moni2(i,j)==4) moni2(i,j)=1
          if (moni2(i,j)==0) goto 79771
          if (j>monn.and.j<mons) moni2(i,j)=3
          if (j>=mons) moni2(i,j)=4
79771     continue
          if (moni3(i,j)==4) moni3(i,j)=1
          if (moni3(i,j)==0) goto 79772
          if (j>monn.and.j<mons) moni3(i,j)=3
          if (j>=mons) moni3(i,j)=4
          irbek(i,j)=moni3(i,j)         ! for runoff routine Gerrit
79772     continue
        end do
      end do
!
!     open larrys output
!
      open(93,file='lsg_output',form='unformatted')
!
!     put land sea mask to cpl 
!
      zlsms(:,:)=wet(:,:,1)
      zlsmv(:,:)=wetvec(:,:,1)
      call cputlsmo(zlsms,zlsmv,ien,jen)
!
      write (no6,*) ' '
      write (no6,*) ' LSG_INIT:  *** INITIALISATION FINISHED *** '
      write (no6,*) ' '
!
      end subroutine lsgini
!
!     ==================
!     subroutine LSGSTEP
!     ==================
!
      subroutine lsgstep
      use lsgvar
      implicit none
!
!     ------------------------------------------------------------------
!
!**** *lsgstep* controles main loop of LSG when coupled to PlaSim
!
!     by S. Lorenz  11/2005. (old do_lsg_step)
!     modified by larry 06/2008
!
!
!     Type declaration.
!     -----------------
!
      integer :: kdatim(7)
      integer :: i,j,k,ibek,ioc,mom,ke1,ke2,je1,je2,je3,je4
      integer :: idra,jdra,kcdat,kctim,ksynt

      real (kind=8) :: cc
      real (kind=8) :: geshfl(6),gesice(6),geswfl(6),gescen(6)
!
      real :: flxdeep(ien,jen) 
      real :: zsst(ien,jen),zmld(ien,jen)
!
!
!*    1. Executing the program.
!     -------------------------
!
!     increase timestep
!
      nt = nt + 1
!
!*    get date/time array
!
      call cgettime(kdatim)
      mdatim(:) = kdatim(:)
!
!*    3.1 Execute one timestep.
!     -------------------------
!
      call step
!
!*    3.2 coupling with PlaSim atmosphere
!     -----------------------------------
!
!     compute heatflux going into the PlaSim mixed layer (ML)
!     by accounting for the following facts and assumptions:
!     - ML depth is set dw(1) (=const) in PlaSim  
!     - see ice has an energy of cp*rho*tfreeze*sice+l*sice 
!       and this was not changed by lsg
!     - only the temperature below ice was changed during the lsg step
!       therefor only the energy in this part of the collumn (i.e. dw(1)-sice) 
!       needs to be considered
!     - not only advection but also mixing in mixg needs to be considered
!
      cc=(rhonul*cp)/dt
!
      do j=1,jen
       do i=1,ien
        if(zeta(i,j) >= 0) then
         fluhea(i,j)=(t(i,j,1)+tkelvin)*(ddz(i,j,1)-sice(i,j))  
        else
         fluhea(i,j)=(t(i,j,1)+tkelvin)*(ddz(i,j,1)+zeta(i,j)-sice(i,j))&
     &              -(t(i,j,2)+tkelvin)*zeta(i,j)*wet(i,j,2)
        endif
        fluhea(i,j)=cc*(fluhea(i,j)                                     &
     &       -(tbound(i,j)+tkelvin)*(ddz(i,j,1)-sice(i,j))*wet(i,j,1))
       enddo
      enddo

      fluhea(:,1)=0.
      fluhea(:,2)=0.
      fluhea(:,jen-1)=0.
      fluhea(:,jen)=0.
!
!     Convert all values to default kind
!
!     add flux correction and put hflx
!     to coupling interface
!
      if(ncoupling < 2) then
       flxdeep(:,:)=fluhea(:,:)
       call cflukob(flxdeep)
       call cputhflo(flxdeep)
      else
       zsst(:,:)=t(:,:,1)
       call cputssto(zsst)
      endif
!
!     debug diagnostic by FL
!
      if(nprint == 1) then
       write(6,*) 'At flux transfere into ml'
       call prfldiag
       write(6,*) 'Max flux into ml: ',maxval(flxdeep),' at ',maxloc(flxdeep)
       write(6,*) 'Min flux into ml: ',minval(flxdeep),' at ',minloc(flxdeep)
       write(6,*) 'Max temp (k=1): ',maxval(t(:,:,1)),' at ',maxloc(t(:,:,1))
       write(6,*) 'Min temp (k=1): ',minval(t(:,:,1)),' at ',minloc(t(:,:,1))
       call inventory
      elseif(nprint == 2) then
       write(6,*) 'At flux transfere into ml'
       write(6,*) 'flux into ml: ',flxdeep(npri,nprj)
       write(6,*) 't (k=1): ',t(npri,nprj,1)
      endif
!
!-----------------------------------------------------------------------
!

!
! Error!     do ibek=1,6 - phimax(4) only
        do ibek=1,6
          gesice(ibek)=0.
          geswfl(ibek)=0.
          gescen(ibek)=0.
          geshfl(ibek)=0.
        enddo
        do ibek=1,4
          phimax1(ibek)=0.
          phimin1(ibek)=0.
          phimax2(ibek)=0.
          phimin2(ibek)=0.
          phimax3(ibek)=0.
          phimin3(ibek)=0.
        end do
!
        do j=1,jen
          do i=1,ien
            if (moni4(i,j)<1) cycle
            geshfl(moni4(i,j))=geshfl(moni4(i,j))+fluhea(i,j)*dl(j)     &
     &                         *dphi*0.5
            geswfl(moni4(i,j))=geswfl(moni4(i,j))+fluwat(i,j)*dl(j)     &
     &                         *dphi*0.5
            gescen(moni4(i,j))=gescen(moni4(i,j))+convad(i,j,1)*dl(j)   &
     &                         *dphi*.5
            gesice(moni4(i,j))=gesice(moni4(i,j))+sice(i,j)*dl(j)       &
     &                         *dphi*0.5
          end do
        end do
!
        do k=1,ken+2
          do ioc=1,3
            do j=1,jen
              psimer(j,k,ioc)=0.
            end do
          end do
        end do
!
        do k=ken,1,-1
          do j=3,jen-3
            do i=1,ien
              mom=monit(i,j)
              if (mom<1) cycle
              psimer(j,k,mom)=psimer(j,k,mom)+vtot(i,j,k)*delta(i,j,k)  &
     &                        *dl(j)
            end do
            do ioc=1,3
              psimer(j,k,ioc)=psimer(j,k,ioc)+psimer(j,k+1,ioc)
            end do
          end do
        end do
!

! @@@     > Achtung - pac/ind ist 0.
!         > neu machen: phiout(50) -> nacheinander und in outdiag ausgeben
! @@@     > UND: auch als Jahresmittel ausgeben (auch outdiagm)!
!
!       For diagnostics calculate max and min overturning:
!         > psimer<0 is southward bottom, northward surface !
!         > maximum prod.  of NADW 46-66N is phimin1(1) - position  4
!         > maximum prod.  of NADW 16-44N is phimin2(1) - position 10
!         > maximum export of NADW at 30N is phimin3(1) - position 16
!         > maximum inflow of AABW 16-44N is phimax2(1) - position 11
!
        if (ken==11) then
          ke1=6    ! 700m
          ke2=5    ! 450m
        else if (ken==22) then
          ke1=11   ! 750m
          ke2=8    ! 450m
        end if

        je1=12     ! 66.25 N
        je2=20     ! 46.25 N
        je3=32     ! 16.25 N
        je4=51     ! 31.75 S (southernmost boundary with land)
!
!       1) NADW/AABW at 46-66N max/min below 700 m
!
        do ioc=1,3
          do k=ke1,ken         ! below 750 m
            do j=je1,je2       ! 46-66N
              phimax1(ioc)=max(phimax1(ioc),0.5*(psimer(j,k,ioc)+       &
     &                     psimer(j+1,k,ioc)))
              phimin1(ioc)=min(phimin1(ioc),0.5*(psimer(j,k,ioc)+       &
     &                     psimer(j+1,k,ioc)))
            end do
          end do
        end do
!
!       2) NADW/AABW at 16-44N max/min below 700 m
!
        do ioc=1,3
          do k=ke1,ken         ! below 750 m
            do j=je2+1,je3     ! 16-44N
              phimax2(ioc)=max(phimax2(ioc),0.5*(psimer(j,k,ioc)+       &
     &                     psimer(j+1,k,ioc)))
              phimin2(ioc)=min(phimin2(ioc),0.5*(psimer(j,k,ioc)+       &
     &                     psimer(j+1,k,ioc)))
            end do
          end do
        end do
!
!       3) Export to Southern Ocean
!
        do ioc=1,3
          do k=ke1,ken          ! below 750 m
            j=je4               ! 31.25 S
            phimax3(ioc)=max(phimax3(ioc),0.5*(psimer(j,k,ioc)+psimer(j+&
     &                  1,k,ioc)))
            phimin3(ioc)=min(phimin3(ioc),0.5*(psimer(j,k,ioc)+psimer(j+&
     &                  1,k,ioc)))
          end do
        end do

!  @@@ hier: statt tm/sm - Gebietsmittel der Becken ablegen!
!      - tropisch (A/P/I); subtropisch (3*s/n); highlats: nNA/sSA/sSP/sSI
!        > 13 Gebiete in 4 Tiefen: 25/250/2000/4000 = 52 Zahlen?
!      - alternativ: 3 Schnitte (jeweils jeden 2. Punkt, 5 Tiefen?)
!

! not used?
!  near Drake Passage stream function (near maximum)
!       T42: jdra=71
        if (ien==64) then
          idra=8
          jdra=40
        end if
!
        idra=20
        jdra=61
!
        rtime=real(nt)/real(ntyear)
!
!
! ascii control output for 11/22 levels
!
!!FL        open  (32,file="LSGtims_phi",form="formatted",position="append")
!!        write (32,9041) rtime,tm(1),sm(1),                              &
!!     &        (-phimin1(ioc)*1.0e-6,-phimax1(ioc)*1.0e-6,ioc=1,3),      &
!!     &        (-phimin2(ioc)*1.0e-6,-phimax2(ioc)*1.0e-6,ioc=1,3),      &
!!     &        (-phimin3(ioc)*1.0e-6,-phimax3(ioc)*1.0e-6,ioc=1,3),      &
!!     &        (psi(idra,jdra)+psi(idra,jdra+1))*0.5
!!        close (32)
!!        open  (32,file="LSGtims_ice",form="formatted",position="append")
!!        write (32,9042) rtime,gescen,geshfl,geswfl,gesice
!!        close (32)
!
!!        open  (32,file="LSGtims_tem",form="formatted",position="append")
!!        write (32,9043) rtime,tm
!!        close (32)
!
!!        open  (32,file="LSGtims_sal",form="formatted",position="append")
!!        write (32,9043) rtime,sm
!!        close (32)
!
 9041   format (f10.3,30f12.5)
 9042   format (f10.3,24e11.3)
 9043   format (f10.3,22f9.4)
!
! diagnose stepping:
!   - with ntcont and normal stop
!   - should move to outdiag !!
!
      if (mod(nt,ntcont).eq.0) then

        call actdate(kcdat,kctim)
        nddr(8)=-1 ! kcdat
        nddr(7)=-1 ! kctim
        ksynt=mod(nt-1,ntyear)+1
!
! write out where the model is at end of each step:
!
!        open  (37,file='run_status_LSG',status='replace')
!        rewind(37)
!        write (37,8090) ktyear,-ktdate,nt,ksynt,kcdat,kctim
! 8090   format (3x,i7,i5.4,3x,'nt=',i10,4x,'nt(year)=',i4,/             &
!     &          3x,'  YYYYY MMDD',                                      &
!     &          3x,'LSG MODEL RUN STATUS (simulation time)',/           &
!     &          18x,'executed at:  yyyymmdd  hhmmss',/                  &
!     &          32x,i8.8,2x,i6.6)
!        close (37)
      endif
!
!
!*      3.2       Write control diagnostics on file output.
!       -----------------------------------------
        if (mod(nt,ntcont)==0) call outdiaglsg
!
!

!       Calculate data average:
!       subroutine deactivated and moved to ../etc
!       ------------------------
!       calculate averages
!       ntaver=3
!       call aver(nt,ntaver)
!
!
!*      3.3       Store the output.
!       -----------------
        if (mod(nt,ntout)==0) then
         call outpostsrv
        endif
!
!*      3.4       Backup.
!       -------

!        > no runtime backup when coupled,
!        > at end of loop only:
!           -  using switch=0 -> inorm=0 for input data

!
!*      3.5.1     Write imatread=0 (compute matrit) in file mat76.sw
!       -------
!  preliminary using ntback
        if (mod(nt,ntback)==0) then
          open (43,file="mat76.sw",form="formatted")
          rewind 43
          write (43,1152)
          write (no6,1153) nt,ntback
 1153     format (' + MAIN: NT=',i8,' NTBACK=',i6,' mod(nt,ntback)=0')
          write (no6,1152)
          close (43)
        endif
 1152     format (' 0 = compute matrit and write to file "mat76"')
!
!       CAUTION: LSG is reproducable only when using same matrit file mat76 !
!                (Stephan, 14.09.2005)
!

!       No End time step loop jnt - controlled by do_lsg_stop
!     end do

      end subroutine lsgstep
!     ==================
!     subroutine LSGSTOP
!     ==================
!
      subroutine lsgstop
      use lsgvar
      implicit none
!
!     ------------------------------------------------------------------
!
!**** *do_lsg_stop* controles end of LSG when coupled to PlaSim
!
!     by S. Lorenz  11/2005.
!
!
!     Interface:
!     -----------------
!     *call* *do_lsg_init* called by interface routine *lsgini* of PlaSim
!     *call* *do_lsg_step* called by interface routine *lsgstep* of PlaSim
!     *call* *do_lsg_stop* called by interface routine *lsgstop* of PlaSim
!
!
!     Type declaration.
!     -----------------
!
      character(len=7) file
      integer :: i,j,k,ioc,mom,kstep
      real (kind=8) :: psisu
!
! new common for dividing do_lsg_step and do_lsg_init:
!
!
!*    4.        Final processing.
!     -----------------
!
!     Prepare file to start the next run with.
!
      write(6,*) "end of time step loop, make backup"
!-----------------------------------------------------------------------
!
      file=filbac1
      call outback(ibac1,iswit,0,file)

!
      rewind 79
      do j=1,jen
        do i=1,ien
          monit(i,j)=0
        end do
      end do
      rewind 79
      if (ien==72) then
        read (79,7900) ((monit(i,j),i=1,ien),j=3,jen-2)
      else if (ien==64) then
        read (79,"(1x,64i1)") ((monit(i,j),i=1,ien),j=3,jen-2)
      else
        read (79,"(1x,128i1)") ((monit(i,j),i=1,ien),j=3,jen-2)
      end if
      ioc=1
      do j=1,jen
        do k=1,ken+2
          psimer(j,k,ioc)=0.
        end do
      end do
!
      do k=ken,1,-1
        do j=1,jen
          psisu=0.
          do i=1,ien
            mom=monit(i,j)
            if (mom==4) mom=1
            if (mom/=1) cycle
            psisu=vtot(i,j,k)*1.e-6*delta(i,j,k)*dl(j)+psisu
          end do
          ioc=1
          psimer(j,k,ioc)=psimer(j,k+1,ioc)+psisu
        end do
      end do
!
! for 22 levels - every second only:
!  - should be changed to output of long-term mean only
!  - or: mean values instead of max/min due to time-mean biases (2006/03)
!
      kstep=1
!     if (ken>=19) kstep=2
      write (no6,*) " meridional circulation Atlantic "
!     write (no6,4711) 0,0,(nint(dw(k)),k=1,ken)
      write (no6,4711) 0,0,(nint(dw(k)),k=1,ken,kstep)
      do j=3,jen-2,2
        write (no6,4711) j,(nint((psimer(j,k,1)+psimer(j+1,k,1))/2.),   &
     &                   k=1,ken+1,kstep)
      end do
      do j=1,jen
        do k=1,ken+2
          psimer(j,k,ioc)=0.
        end do
      end do
!
      do k=ken,1,-1
        do j=1,jen
          psisu=0.
          do i=1,ien
            mom=monit(i,j)
            if (mom==4) mom=1
            if (mom/=2) cycle
            psisu=vtot(i,j,k)/1000000.*delta(i,j,k)*dl(j)+psisu
          end do
          ioc=1
          psimer(j,k,ioc)=psimer(j,k+1,ioc)+psisu
        end do
      end do
      write (no6,*) " meridional circulation Pacific     "
      write (no6,4711) 0,0,(nint(dw(k)),k=1,ken,kstep)
      do j=3,jen-2,2
        write (no6,4711) j,(nint((psimer(j,k,1)+psimer(j+1,k,1))/2.),   &
     &                   k=1,ken+1,kstep)
      end do
!!
!     close larrys output
!
      close(93)
!
!
!4711 format (1x,i2,12(i4))
!4711 format (1x,i2,23(i5))
 4711 format (1x,i2,23(i4))
!

!
!     write (no6,'(a,f10.3,a,i10)')
!    &    ' normal end of lsg run - time =',rtime,' - step =',nt
      write (no6,'(a,f10.3,a,i10)')                                     &
     &    '  %%% Normal end of LSG run - time =',rtime,                 &
     &    '    last timestep =',nt
 7900 format (1x,72i1)

      end subroutine lsgstop
!     ==================================================================
      subroutine matrix
      use lsgvar
      implicit none
!     ------------------------------------------------------------------
!     by E. Maier-Reimer.
!     Last modified by U. Mikolajewicz 9/87.
!     Optimised by E. Kirk 07-Mar-2008
!
!     Purpose.
!     --------
!     *matrit* is only active when nsve=2 !!!!!!!!
!
!     Computes the matrix of the left hand side of the equation
!     for the barotropic velocities. The matrix is triangularised
!     by Gauss" elimination. The elimination factors,
!     scaling factors and the triangularised matrix are
!     used in subroutine *uvtrop* to solve
!     the abovementioned equation.
!
!**   Input.
!     ------
!     rhdif     difference in pot density
!
!     Output.
!     -------
!     trisys    triangularised matrix in compressed form.
!     elim      eliminationfactors.
!     skal      scalingfactors.
!
!     Interface.
!     ----------
!     *call* *matrix*.
!     ------------------------------------------------------------------
!
!     Declaration of local constants/variables.
!     -----------------------------------------
      logical :: l77       ! existence of file <mat77>
      logical :: lnn,lne,lee,lse,lss,lsw,lww,lnw
      integer :: ilo       ! index for longitude
      integer :: jla       ! index for latitude
      integer :: k,l,m
      integer :: im1,ip1
      integer :: jm2,jm1,jp1,jp2
      integer :: jcc,jnn,jne,jee,jse,jss,jsw,jww,jnw
      integer :: ll
      integer :: itri,kvar
      integer :: ja,je,k1
      real (kind=8) :: ghnn,ghne,ghse,ghss,ghsw,ghnw
      real (kind=8) :: audtne,audtse,audtsw,audtnw
      real (kind=8) :: ztm,hm,wilkin,zsum,gt2
      real (kind=8) :: zf,dz,gre,grw,grn,grs,skalil
      real (kind=8) :: re,r1,zdlisq,zdpinsq
      real (kind=8) :: a(kbm,km)
      real (kind=8) :: zrow(-kb:kb)
      real          :: cpu1,cpu2

      call cpu_time(cpu1)
!
!*    1.        Set initial values and constants.
!     ---------------------------------
!
      zdpinsq = dpin * dpin ! (1.0 / delta phi) ** 2
      wilkin = 1.0 ! max. element of elimination-array
      itri=0
!
!
!*    1.1   Read matrices if <mat77> exists
!
      inquire(file="mat77",exist=l77)
      if (l77) then
         open (42,file="mat77",form="unformatted")
         read (42) elim
         read (42) trisys
         read (42) skal
         close(42)
         write(6,*) ">>> Read matrix data from file <mat77> <<<"
         return ! Success
      else
         write(6,*) ">>> No file <mat77> : computing matrix <<<"
      endif ! (l77)
!
!
!*    2.        Computation of the matrix.
!     --------------------------
!
!     l running index
!
      l=0
!
!     Main loops.
!
      do jla=3,jen-2
!
        jm2 = jla-2
        jm1 = jla-1
        jp1 = jla+1
        jp2 = jla+2
!
        zdlisq = dli(jla) * dli(jla) ! (1.0 / delta lambda) ** 2
!
        do ilo=1,ien
!
!         Computes only, when sea-point.
!
          if (numx(ilo,jla) == 0) cycle
!
!*        2.1.     Set initial values dependent on for ilo and jla.
!         --------------------------------------------
!
!         Definition of surrounding indices.
!
          im1 = ilo-1
          ip1 = ilo+1

          if (im1 <   1) im1 = ien
          if (ip1 > ien) ip1 =   1

          lnn = numx(im1,jm2) > 0
          lne = numx(ilo,jm1) > 0
          lee = numx(ip1,jla) > 0
          lse = numx(ip1,jp1) > 0
          lss = numx(ip1,jp2) > 0
          lsw = numx(ilo,jp1) > 0
          lww = numx(im1,jla) > 0
          lnw = numx(im1,jm1) > 0
!
!         Determination of the surrounding depths t and the inverses h .
!
          ztm = depth(ilo,jla)
          hm  = 1.0 / ztm
!
!         Determination of the numbers of the surrounding points.
!
          jcc = 2 * numx(ilo,jla)
          jnn = 2 * numx(im1,jm2) - jcc
          jne = 2 * numx(ilo,jm1) - jcc
          jee = 2 * numx(ip1,jla) - jcc
          jse = 2 * numx(ip1,jp1) - jcc
          jss = 2 * numx(ip1,jp2) - jcc
          jsw = 2 * numx(ilo,jp1) - jcc
          jww = 2 * numx(im1,jla) - jcc
          jnw = 2 * numx(im1,jm1) - jcc
!
!         Factors.
!
          gt2  = g * dt * dt
          ghnn = dpin * gt2 * depth(im1,jm2)
          ghne = dpin * gt2 * depth(ilo,jm1)
          ghse = dpin * gt2 * depth(ip1,jp1)
          ghss = dpin * gt2 * depth(ip1,jp2)
          ghsw = dpin * gt2 * depth(ilo,jp1)
          ghnw = dpin * gt2 * depth(im1,jm1)
!
!         Horizontal friction
!
          audtne = 0.5 * dt * (cofric(ilo,jla) + cofric(ilo,jm1))
          audtse = 0.5 * dt * (cofric(ilo,jla) + cofric(ip1,jp1))
          audtsw = 0.5 * dt * (cofric(ilo,jla) + cofric(ilo,jp1))
          audtnw = 0.5 * dt * (cofric(ilo,jla) + cofric(im1,jm1))
!
!         Main loop!!!!!!
!
!
!         *kvar*=1: computation for *ub*.
!         *kvar*=2:     "        "  *vb*.
!
          do kvar=1,2
!
!           *l* gives the number of lines of the matrix.
!
            l = l + 1
            m = l
            if (l > km) m = km
            zrow(:) = 0.0
!
!           Computation of *ub* (kvar == 1) and *vb* (kvar == 2).
!
            if (kvar == 1) then
!
!*          2.2.      computation of *ub*-terms.
!           --------------------------
!
!           Determination of zonal component of the terms containing
!           sea-surface variation, horizontal friction and Coriolis
!           force.
            zrow(0) = ztm*gt2*zdlisq*2.0+audtnw+audtne+audtsw+audtse
            zrow(1) = -ff(jla)*dt
            if (lne) zrow(jne)   = -audtne
            if (lse) zrow(jse)   = -audtse
            if (lsw) zrow(jsw)   = -audtsw
            if (lnw) zrow(jnw)   = -audtnw
            if (lee) zrow(jee)   = -gt2*depth(ip1,jla)*zdlisq
            if (lww) zrow(jww)   = -gt2*depth(im1,jla)*zdlisq
            if (lne) zrow(jne+1) = -ghne*zdlisq*dl(jm1)
            if (lse) zrow(jse+1) =  ghse*zdlisq*dl(jp1)
            if (lsw) zrow(jsw+1) = -ghsw*zdlisq*dl(jp1)
            if (lnw) zrow(jnw+1) =  ghnw*zdlisq*dl(jm1)
!
!           There was a bug in the original code for gridpoint
!           (2,5) the neighboure (1,4) was classified as land
!           instead of water ( numx > 1) instead of ( numx >= 1)
!
!           Terms including the density differences for the zonal
!           component.
!
            zf = g * hm * (dt * dt) * (dli(jla) * dli(jla))
!
            do k=2,ken
              dz = zf * delta(ilo,jla,k)
              do ll=2,k
                gre=rhdif(ip1,jla,ll-1)
                grw=rhdif(ilo,jla,ll-1)
                zrow(0)=zrow(0)+dz*delta(ilo,jla,ll)*(gre+grw)
                if (lee) zrow(jee)=zrow(jee)-dz*delta(ip1,jla,ll)*gre
                if (lww) zrow(jww)=zrow(jww)-dz*delta(im1,jla,ll)*grw
                if (lne) zrow(jne+1)=zrow(jne+1)                        &
     &              -dz*dpin*dl(jm1)*delta(ilo,jm1,ll)*gre
                if (lse) zrow(jse+1)=zrow(jse+1)                        &
     &              +dz*dpin*dl(jp1)*delta(ip1,jp1,ll)*gre
                if (lsw) zrow(jsw+1)=zrow(jsw+1)                        &
     &              -dz*dpin*dl(jp1)*delta(ilo,jp1,ll)*grw
                if (lnw) zrow(jnw+1)=zrow(jnw+1)                        &
     &              +dz*dpin*dl(jm1)*delta(im1,jm1,ll)*grw
              end do
            end do
            endif ! (kvar == 1)

            if (kvar == 2) then
!
!*          2.3.      Computation of *vb*-terms.
!           --------------------------
!
!           Computation of the meridional component of the terms
!           containing sea-surface variation and horizontal friction
!           and Coriolis force.
!
            zrow(0)=ztm*gt2*zdpinsq*dl(jla)*(dli(jm1)+dli(jp1))         &
     &               +audtnw+audtne+audtsw+audtse
            zrow(-1)=ff(jla)*dt
            if (lne) zrow(jne)  =-audtne
            if (lse) zrow(jse)  =-audtse
            if (lsw) zrow(jsw)  =-audtsw
            if (lnw) zrow(jnw)  =-audtnw
            if (lnn) zrow(jnn)  =-ghnn*dpin*dl(jm2)*dli(jm1)
            if (lss) zrow(jss)  =-ghss*dpin*dl(jp2)*dli(jp1)
            if (lne) zrow(jne-1)=-ghne*dli(jm1)
            if (lse) zrow(jse-1)= ghse*dli(jp1)
            if (lsw) zrow(jsw-1)=-ghsw*dli(jp1)
            if (lnw) zrow(jnw-1)= ghnw*dli(jm1)
!
!           Terms including the density differences for the merid
!           component.
!
            zf=g*hm*dt**2
            do k=2,ken
              dz = zf * delta(ilo,jla,k)
              do ll=2,k
                grn=rhdif(ilo,jm1,ll-1)
                grs=rhdif(ip1,jp1,ll-1)
                zrow(0)=zrow(0)+dz*delta(ilo,jla,ll)*zdpinsq*           &
     &                   dl(jla)*(dli(jm1)*grn+dli(jp1)*grs)
                if (lnn) zrow(jnn)=zrow(jnn)-dz*                        &
     &              delta(im1,jm2,ll)*zdpinsq*dl(jm2)*dli(jm1)*grn
                if (lss) zrow(jss)=zrow(jss)-dz*                        &
     &              delta(ip1,jp2,ll)*zdpinsq*dl(jp2)*dli(jp1)*grs
                if (lne) zrow(jne-1)=zrow(jne-1)-                       &
     &              dpin*dli(jm1)*grn*dz*delta(ilo,jm1,ll)
                if (lsw) zrow(jsw-1)=zrow(jsw-1)-                       &
     &              dpin*dli(jp1)*grs*dz*delta(ilo,jp1,ll)
                if (lnw) zrow(jnw-1)=zrow(jnw-1)+                       &
     &              dpin*dli(jm1)*grn*dz*delta(im1,jm1,ll)
                if (lse) zrow(jse-1)=zrow(jse-1)+                       &
     &              dpin*dli(jp1)*grs*dz*delta(ip1,jp1,ll)
              end do
            end do
            endif ! (kvar == 2)
!
!*          2.4.      Computation of *elim*, *trisys* and *skal*.
!           -------------------------------------------
!           Scaling.
!
            skalil = 1.0 / maxval(abs(zrow(:)))
            skal(l)=skalil
!
!           Write a new lowest row of the rhombus *a*.
!
            a(:,m) = zrow(:) * skalil
!
!           If the rhombus not complete, get a new row.
!
            if (l<=kb) cycle
!
!           If all lines are computed, handle the remaining rhombus.
!
            if (l < matrx) then !  Elimination of the 1st column
              re=a(km,1)
              itri=itri+1
              do k=1,kb
                ja=km-k
                je=ja+kb
                r1=a(ja,k+1)/re
                elim(k,itri)=r1
                wilkin=max(wilkin,abs(r1))
                a(ja:je,k+1) = a(ja:je,k+1) - r1 * a(ja+k:je+k,1)
              end do
              trisys(:,itri) =  a(km:kbm,1) !  Write a new row
              a(:,1:kb) = a(:,2:km)     ! Shift *a* one row up
              a(:,km) = 0.0  ! Set the last line equal to zero
!
            else !  Final processing on the rhombus
!
              do k1=1,kb
                itri=itri+1
                re=a(km,k1)
!
!               Elimination.
!
                do k=1,km-k1
                  ja=km-k
                  je=ja+kb
                  r1=a(ja,k+k1)/re
                  elim(k,itri)=r1
                  wilkin=max(wilkin,abs(r1))
                  a(ja:je,k+k1) = a(ja:je,k+k1) - r1 * a(ja+k:je+k,k1)
                end do ! k
                trisys(:,itri) = a(km:kbm,k1)
              end do ! k1
            endif ! (l < matrx)
!
!           End of the main loops.
!
          end do ! kvar
        end do ! ilo
      end do ! jla
!
!     Processing on the last row.
!
      trisys(1   ,matrx) = a(km,km)
      trisys(2:km,matrx) = 0.0

      call cpu_time(cpu2)
      write(6,'("CPU time for building mat77:",f10.4," *")') cpu2-cpu1
      write(6,*) 'rhdif min = ',minval(rhdif)
!
!*    2.2   write computed matrit data (lsgmat)
!
      open  (42,file="mat77",form="unformatted")
      write (42) elim
      write (42) trisys
      write (42) skal
      close (42)
!
!     Write info and check for correct dimension
!
      write(6,'("* Maximum elimination element:",g20.12)') wilkin
      write(6,'("* Matrix rows = ",i7," *")') l

      if (l /= matrx) then
         write(6,'("* matrx    = ",i7," *")') matrx
         write(6,*) ">>> Abort: rows must be equal to (matrx) <<<"
         stop "wrong matrx"
      endif
      return
      end subroutine matrix
!     ==================================================================
!     ------------------------------------------------------------------
!
      subroutine mixg
      use lsgvar
      implicit none
!
!     ------------------------------------------------------------------
!
!**** subroutine *mixg*
!
!     by G. Lohmann, MPI, nov. 1998
!     modified by larry 12.03.2008
!
!**   Purpose.
!     --------
!     *mixg*  connects the ocean surface
!     to atmospheric temperature tendencies and freshwater fluxes
!     These tendencies come from surf-routines or data
!
!     New boundary values for temperature and salinity are calculated
!
!     for nsmix*=1 : *mix* (computes changes the surface layer, without ice,
!                                must be driven by water temperatures.)
!     for nsmix*=2 : *mixatt* (computes changes of *t* and *s*, with ice,
!                                  must be driven by air temperatures)
!     for nsmix*=3 : *mixg* (computes changes of *t* and *s*, without ice,
!                                  driven by atmospheric model/data)
!
!     Input.
!     ------
!     t       old temperature.
!     s       old salinity.
!     tbound  boundary tempearture.
!     sbound  boundary salinity.
!             both in common block /lsgfie/.
!     fluwat  fresh water flux [m/s].
!             both in common block /lsgsur/.
!     zeta    old surface elevation.
!
!     Output.
!     -------
!     t       new pot temperature on the surface.
!     s       new salinity on the surface.
!     flukhea heat flux due to newtonian coupling.
!     flukwat fresh water flux due to newtonian coupling.
!
!     Interface.
!     ----------
!     *call* *mixg*.
!
!     ------------------------------------------------------------------

      integer :: i,j,k
      real (kind=8) :: eps,sold,dn,zdold,zdnew
!
      real (kind=8) :: told(ien,jen)

      real (kind=8) :: flumean, fluarea
      real (kind=8) :: ztbound
!
!     diagnostics
!
      real (kind=8) :: ztold(ien,jen,ken)
!
!*    Set initial values and constants.
!
      eps=1.e-10
!
      ztold(:,:,:)=t(:,:,:)
!
!*    Computation of new temperature and salinity.
!     --------------------------------------------
!
!     change surface parameter (uppermost level) due to boundary
!     fluxes and values
!
!PHCONS:
      if (debug_conservation) then
        flumean=0.0
        fluarea=0.0
        do j=3,jen-2
          do i=1,ien
            if (wet(i,j,1)<0.5) cycle
            flumean=flumean+fluwat(i,j)*dlh(i,j)*dphi
            fluarea=fluarea+dlh(i,j)*dphi
          end do
        end do
        flumean=flumean/fluarea
        do j=3,jen-2
          do i=1,ien
            if (wet(i,j,1)<0.5) cycle
            fluwat(i,j)=fluwat(i,j)-flumean
          end do
        end do
      endif
!
      flukhea(:,:)=0.
!
      do j=3,jen-2
       do i=1,ien
        if (wet(i,j,1)<0.5) cycle
        ztbound=tbound(i,j)
!
!     sea level change
!
        zdold=ddz(i,j,1)+zeta(i,j)-sice(i,j)
        zeta(i,j)=zeta(i,j)+dt*fluwat(i,j)
!
!     limit sea level (to omit crashes)
!
        if(abs(zeta(i,j)) > du(1)) then
         write(6,*) 'WARNING: abs(zeta) > du(1) detected at ',i,j
         zeta(i,j)=dmax1(zeta(i,j),-du(1))
         zeta(i,j)=dmin1(zeta(i,j),du(1))
        endif 
        zdnew=ddz(i,j,1)+zeta(i,j)-sice(i,j)
!
!     replace uppermost 50m (dw1) temperature by boundary condition
!     and diagnose heat flux
!     (if zeta > 0 mix the remaining water, 
!      if zeta < 0 mix new water into 2nd layer)
!
!PHCONS: no heat exchange with PlaSim when debug_conservation

        if (.not.debug_conservation) then 
         if(zeta(i,j) >= 0.0) then
          told(i,j)=t(i,j,1)
          t(i,j,1)=(zeta(i,j)*t(i,j,1)+ddz(i,j,1)*ztbound)              &
     &            /(ddz(i,j,1)+zeta(i,j))
          flukhea(i,j)=(t(i,j,1)-told(i,j))*(ddz(i,j,1)+zeta(i,j))      &
     &                *rhonul*cp/dt 
         else
          if(wet(i,j,2) < 0.5) then
           write(6,*) 'schiefgegangen'
           stop
          endif
          told(i,j)=t(i,j,2)
          t(i,j,2)=((ddz(i,j,2)+zeta(i,j))*t(i,j,2)                     &
     &              -zeta(i,j)*ztbound)/ddz(i,j,2)

          flukhea(i,j)=((ztbound-t(i,j,1))*(ddz(i,j,1)+zeta(i,j))       &
     &                 +(t(i,j,2)-told(i,j))*ddz(i,j,2))                &
     &                *rhonul*cp/dt
          t(i,j,1)=ztbound
         endif
        endif
!
!     salinity change due to sea level change
!     and diagnose fresh water flux
!
        s(i,j,1)=s(i,j,1)*zdold/zdnew
        flukwat(i,j)=fluwat(i,j)
       end do
      end do
!
!     diagnostics
!
      dtdt(:,:,:,2)=(t(:,:,:)-ztold(:,:,:))/dt
!
!     debug diagnostic by FL
!
      if(nprint == 1) then
       write(6,*) 'At end of mixg:'
       call prfldiag
       write(6,*) 'Max s (k=1): ',maxval(s(:,:,1),MASK=(wet(:,:,1) > 0.5))  &
     &       ,' at ',maxloc(s(:,:,1),MASK=(wet(:,:,1) > 0.5))
       write(6,*) 'Min s (k=1): ',minval(s(:,:,1),MASK=(wet(:,:,1) > 0.5))  &
     &       ,' at ',minloc(s(:,:,1),MASK=(wet(:,:,1) > 0.5))
       write(6,*) 'Max t (k=1): ',maxval(t(:,:,1),MASK=(wet(:,:,1) > 0.5))  &
     &,' at ',maxloc(t(:,:,1),MASK=(wet(:,:,1) > 0.5))
       write(6,*) 'Min t (k=1): ',minval(t(:,:,1),MASK=(wet(:,:,1) > 0.5))  &
     &       ,' at ',minloc(t(:,:,1),MASK=(wet(:,:,1) > 0.5))
       call inventory
      elseif(nprint == 2) then
       write(6,*) 'At end of mixg:'
       write(6,*) 's (k=1): ',s(npri,nprj,1)
       write(6,*) 't (k=1): ',t(npri,nprj,1)
       write(6,*) 'z: ',zeta(npri,nprj)
       write(6,*) 'ice: ',sice(npri,nprj)
      endif
!
      return
!
      end subroutine mixg
!     =================================================================
!     ------------------------------------------------------------------
      subroutine opfile
      use lsgvar
      implicit none
!
!     ------------------------------------------------------------------
!
!**** *opfile*- opens files and initialises common blocks /lsgtap.
!
!*    by U. Mikolajewicz 1/88.
!     modified by Larry 08/08
!
!*    Input.
!     ------
!     not required.
!
!     Output.
!     -------
!     common blocks *lsgtap* in module *lsgvar*
!
!     Interface.
!     ----------
!     *call* *opfile*.
!
!     ------------------------------------------------------------------
!
!*    1. Initialize tape numbers
!
      itapin=5    ! namelist
      no6=6       ! print out
      idat=52     ! old restart
      ibac1=84    ! new restart
      ibac2=85    ! backup restart
      iswit=86    ! tape unit for the last written backup
      nopost=53   ! postprocessor output
      itop=51     ! topography
      no71=72     ! sst climatology (switched to 72 because 71 is used)
      no73=73     ! sss climatology
      nowat=74    ! pme climatology
      no75=75     ! v-stress climatology
      no76=76     ! u-stress climatology
      nohea=77    ! heat flux climatology
      no37=54     !
      no70=70     ! sea ice climatology
!
!*    2. Open and rewind files.
!
!*    2.1 Files with input data.
!
!     namelist.
!
      open (itapin,file="input",access="sequential",form="formatted")
      rewind itapin
!
!     depth distribution.
!
      open (itop,file="topogr",access="sequential",form="formatted")
      rewind itop
!
!*    2.2 Files for output.
!    
      open (79,file="beklas",access="sequential",form="formatted")
!
!*    2.3 Boundary conditions
!
!     sst
!
      open (unit=no71,file="montem")
!
!     sss
!
      open (unit=no73,file="monsal")
!
!     sea ice
!
      open (unit=no70,file="seaice")
!
!     wind stress
!
      open (unit=no76,file="wistrx")
      open (unit=no75,file="wistry")
!
!     fresh water flux
!
      open (unit=nowat,file="salflu")
!
!     heat flux
!
      if (nsflu==1) open (unit=nohea,file="heaflu")
!
      return
      end subroutine opfile
!     ==================================================================
!     ------------------------------------------------------------------
!
      subroutine outback(ktape,ktswit,kswit,yfile)
      use lsgvar
      implicit none
      integer :: ktape,ktswit,kswit
      character(len=7) :: yfile
!
!     ------------------------------------------------------------------
!
!**** *outback*
!
!     by U. Mikolajewicz 12/87.
!
!     Purpose.
!     --------
!     *outback* writes a restart file on file *yfile*.
!
!     Input.
!     ------
!     ktape   number of the tapeunit to write the backup on.
!     ktswit  number of tapeunit to write switches on.
!     kswit   value of the switch.
!     yfile   file to write on.
!
!     ------------------------------------------------------------------
!
!     Declaration of local variables.
!

      integer :: jlev,i
      real (kind=8) :: zprel,zlen,zcode,zdum,zlev,z_du,zero
      real (kind=8) :: z_dp(ien,jen)
!
!*    1.        Actualize *nddr*.
!     -----------------
      open (ktape,file=yfile,access="sequential",form="unformatted")
      zero=0.
      zdum=zero
      nddr(11) = -(100 * mdatim(2) + mdatim(3))
      nddr(12) = mdatim(1)
!
!*    2.        Write backup on file.
!     ---------------------
!     write *nddr* and *oddr* on file.
!
      rewind ktape
      write (ktape) nddr
      write (ktape) oddr
      zlen=ien*jen
      zprel=6.
!
!     Write temperature.
!
      zcode=-62
      do jlev=1,ken
        z_du = du(jlev)
        z_dp(:,:) = t(:,:,jlev)
        write (ktape) zprel,zlen,zcode,z_du,zdum,zdum
        write (ktape) z_dp
      end do
!
!     Write salinity.
!
      zcode=-5
      do jlev=1,ken
        z_du = du(jlev)
        z_dp(:,:) = s(:,:,jlev)
        write (ktape) zprel,zlen,zcode,z_du,zdum,zdum
        write (ktape) z_dp
      end do
!
!     Write baroclinic modes
!
      zcode=-63
      do jlev=1,ken
        z_du = du(jlev)
        z_dp(:,:) = u(:,:,jlev)
        write (ktape) zprel,zlen,zcode,z_du,zdum,zdum
        write (ktape) z_dp
      end do
!
      zcode=-64
      do jlev=1,ken
        z_du = du(jlev)
        z_dp(:,:) = v(:,:,jlev)
        write (ktape) zprel,zlen,zcode,z_du,zdum,zdum
        write (ktape) z_dp
      end do
!
!     Write barotropic velocities.
!
      zlev=-100.
      zcode=-37.
      z_dp(:,:) = ub(:,:)
      write (ktape) zprel,zlen,zcode,zlev,zdum,zdum
      write (ktape) z_dp
!
      zcode=-38.
      z_dp(:,:) = vb(:,:)
      write (ktape) zprel,zlen,zcode,zlev,zdum,zdum
      write (ktape) z_dp
!
!     Write surface elevation and ice thickness.
!
      zcode=-1.
      z_dp(:,:) = zeta(:,:)
      write (ktape) zprel,zlen,zcode,zlev,zdum,zdum
      write (ktape) z_dp
!
      zcode=-13.
      z_dp(:,:) = sice(:,:)
      write (ktape) zprel,zlen,zcode,zlev,zdum,zdum
      write (ktape) z_dp
!
      zcode=-88.
      z_dp(:,:) = stor(:,:)
      write (ktape) zprel,zlen,zcode,zlev,zdum,zdum
      write (ktape) z_dp
      zcode=-89.
      z_dp(:,:) = riv(:,:)
      write (ktape) zprel,zlen,zcode,zlev,zdum,zdum
      write (ktape) z_dp
      zcode=-99.
      z_dp(:,:) = depth(:,:)
      write (ktape) zprel,zlen,zcode,zlev,zdum,zdum
      write (ktape) z_dp
!
!     added for perfect restart
!
!     Write 'old' velocities for adv_quick
!
      zcode=-3
      do jlev=1,ken
        z_du = du(jlev)
        z_dp(:,:) = utotold(:,:,jlev)
        write (ktape) zprel,zlen,zcode,z_du,zdum,zdum
        write (ktape) z_dp
      end do
!
      zcode=-4
      do jlev=1,ken
        z_du = du(jlev)
        z_dp(:,:) = vtotold(:,:,jlev)
        write (ktape) zprel,zlen,zcode,z_du,zdum,zdum
        write (ktape) z_dp
      end do
!
      zcode=-7
      do jlev=1,ken
        z_du = du(jlev)
        z_dp(:,:) = wold(:,:,jlev)
        write (ktape) zprel,zlen,zcode,z_du,zdum,zdum
        write (ktape) z_dp
      end do
!
      close (ktape)
!
!
!
!*    3.        Actualize file *filswit*.
!     -------------------------
      if (ktape==idat) goto 4000
 9860 format (4x,'OUTBACK: actualise: ',3a,i2,a,i2,a,i3)
      write(no6,9860) 'filswit=',filswit,' idat=',idat,' kswit=',kswit, &
     &                ' ktape=',ktape
      open (ktswit,file=filswit,form="formatted")
      rewind ktswit
      write (ktswit,*) ktape,kswit
      close (ktswit)
!
!*    4.        Send message to dayfile.
!     ------------------------
 4000 continue
      write (no6,9870) nddr(12),nddr(11),yfile
 9870 format (" +  OUTBACK  YYYYY -MMDD :",i5.2,i6.4,                   &
     &        " backup written on file ",a7,"  +")
!
      end subroutine outback
!     ==================================================================
      subroutine outdiaglsg
      use lsgvar
      implicit none
!     ------------------------------------------------------------------
!
!**** *outdiag* -writes diagnostics on file *output*.
!
!*    by U. Mikolajewicz 1/88.
!
!     Interface.
!     ----------
!     *call* *outdiag*.
!     ------------------------------------------------------------------
!
!     Declaration of local variables.
!
      interface
         real (kind=8) function rho(s,t,p)
            real (kind=8), intent(in) :: s,t,p
         end function rho
      end interface
            
      integer :: i,j,k,kend5,imin,imax,jmin,jmax,jn,js

      real (kind=8) :: xlat,xlon,vol,are,zsima,zsimi,zsime
      real (kind=8) :: sumekin,sumarea,hflsum,wflsum,eps
      real (kind=8) :: zetmin,zetmax,arn,ars,von,vos,salinv,volinv
      real (kind=8) :: phiout(20)
      real (kind=8) :: area(ken),ekin(ken)
      real (kind=8) :: field(ien,jen)
      real (kind=8) :: zetbar
      character(80) :: yf9975 = "(1x,a,i2,a,i2,a,g10.3,a,g10.3,a,f3.0)"
!
!*    1.        Set initial constants.
!     ----------------------
      kend5 = max(1,ken/5)
!
!*    2.        Write diagnostics.
!     ------------------
!
!*    2.1       Timestep and date
!
 9801 format("* LSG timestep",i7,"   date: ",i2,"-",a,"-",i4.4," *")
      write(no6,9801) nt,mdatim(3),ymon(mdatim(2)),mdatim(1)
!
!*    2.2       Maximum of barotropic stream function.
!     (zero at Antarctica, computed in *uvtrop*).
!
      xlat=oddr(52-1+mindj)
      xlon=oddr(52-1+jen+mindj)+float(mindi-1)*oddr(50)
      if (xlon>180.) xlon=xlon-360.
      if (xlon<-180.) xlon=xlon+360.
      write (no6,9808) psimax,xlon,xlat
 9808 format (1x,"Maximum of barotropic streamfunction in Sv",g10.3,    &
     &        " at lon=",f5.0," lat=",f5.0)
!
!*    2.3       Average potential temperature and salinity.
!     (computed in *advect*).
      write (no6,9802) (du(k),k=2,ken,kend5)
 9802 format (1x,"Layer depth in m      ",6f9.2)
      write (no6,9803) (tm(k),k=2,ken,kend5)
 9803 format (1x,"Av.pot.temp.(C)       ",6f9.4)
      write (no6,9813) (sm(k),k=2,ken,kend5)
 9813 format (1x,"Av.salinity in 0/00   ",6f9.4)
      write (no6,9713) ((rho(sm(k),tm(k),0.0_8)-1000.),k=2,ken,kend5)
 9713 format (1x,"Av.sigma-theta        ",6f9.4)
      write (no6,9714) ((rho(sm(k),tm(k),.1*du(k))-1000.),k=2,ken,kend5)
 9714 format (1x,"Av.sigma-in-situ      ",6f9.4)
!
!
!*    2.4       Average kinetic energy.
!
      do k=1,ken
        sumekin = 0.0
        sumarea = 0.0
        do j=1,jen
          do i=1,ien
            sumekin=sumekin+(utot(i,j,k)**2+vtot(i,j,k)**2)             &
     &             *wetvec(i,j,k)*rhonul*dphi*dlh(i,j)*0.5
            sumarea=sumarea+dphi*dlh(i,j)*wetvec(i,j,k)
          end do
        end do
        ekin(k) = sumekin
        area(k) = sumarea
        if (sumarea > 1.0) ekin(k) = sumekin / sumarea
      end do
      write (no6,9807) (ekin(k),k=2,ken,kend5)
 9807 format (1x,"Av. kin. energy J/m**3 ",8f9.4)
!
!*    2.5       Sum of upward and downward transports and conv. adj.
!     (computed in *cont1* and *dens*).
!
!     Depth in w-points.
!
      write (no6,9804) (dw(k),k=1,ken,kend5)
 9804 format (1x,"Depth in m (vector-P) ",6f9.2)
      write (no6,9805) (trup(k)*1.e-6,k=1,ken,kend5)
 9805 format (1x,"Upw.transports in Sv  ",6f9.2)
!     write (no6,9806) (abs(trdo(k))*1.e-6,k=1,ken-1,kend5)
 9806 format (1x,"Downw.transports in Sv",6f9.2)
      write (no6,9823) (nrca(k),k=1,ken,kend5)
 9823 format (1x,"Conv.adjustm. events: ",6i9)
!
!*    2.6       Ice-diagnostics.
!
      vol=0.0
      are=0.0
      do j=1,jen
        do i=1,ien
          if (sice(i,j)*wet(i,j,1)<=0.01) cycle
          vol=vol+dlh(i,j)*dphi*sice(i,j)*0.5
          are=are+dlh(i,j)*dphi*0.5
        end do
      end do
      zsime=0.0
      if (are>0.0) then
        zsime=vol/are
        write (no6,9810) vol,are,zsime
 9810   format (1x,"Icevol. m**3 ",e10.3," icecov.area m**2 ",e10.3,    &
     &          " Av.thickness in m ",f10.3)
      end if
      zsima=sicea(1)
      zsimi=sicea(1)
      do i=1,ienjen
        if (sicea(i)*wetb(i)>zsima) zsima=sicea(i)
        if (sicea(i)*wetb(i)<zsimi) zsimi=sicea(i)
      end do
      write (no6,9811) zsimi,zsima,zsime
 9811 format (1x,"Range icethickness       from",f10.4,"  to",f10.4,    &
     &        "   average ",f9.4," m")
!
!
!*    2.7       Temperature in the surface layer.
!     (computation of extremes).
      do j=1,jen
        do i=1,ien
          field(i,j)=t(i,j,1)
          if (wet(i,j,1)<0.5) field(i,j)=tm(1)
          if ((field(i,j)<-3.).or.(field(i,j)>50.)) then
            write(no6,yf9975) 'Extr.Temp. at i=',i,' j=',j,             &
     &      ' T=',field(i,j),' Tb=',tbound(i,j),' wet=',wet(i,j,1)
!           stop "STOP -- in computation of extremes"
          end if
        end do
      end do
      zsima=maxval(field)
      zsimi=minval(field)
      write (no6,9812) zsimi,zsima,tm(1)
 9812 format (1x,"Range temperat. at surf. from",f10.4,"  to",f10.4,    &
     &        "   average ",f9.4," C")
!
!*    2.8       Salinity in the surface layer.
!     (computation of extremes).
      do j=1,jen
        do i=1,ien
          field(i,j)=s(i,j,1)
          if (wet(i,j,1)<0.5) field(i,j)=sm(1)
        end do
      end do
      zsima=maxval(field)
      zsimi=minval(field)
      write (no6,9814) zsimi,zsima,sm(1)
 9814 format (1x,"Range salinity at surf. from ",f10.4,"  to",f10.4,    &
     &        "   average ",f9.4," 0/00")
!
!
!
!*    2.9       Flux diagnostics.
!     -----------------
!     Global averaged fluxes.
      hflsum=0.
      wflsum=0.
      are=0.
      do j=1,jen
        do i=1,ien
          if (wet(i,j,1)<0.5) cycle
          are=are+dlh(i,j)*dpin
          hflsum=hflsum+flukhea(i,j)*dlh(i,j)*dpin
          wflsum=wflsum+flukwat(i,j)*dlh(i,j)*dpin
        end do
      end do
!  for real*4: e-20 is enough!
!     eps=1.e-60
      eps=1.e-20
      hflsum=hflsum/(are+eps)
      wflsum=wflsum/(are+eps)
      write (no6,9874) hflsum,wflsum*30.*day*1000.
 9874 format ("Glob.av.heatflux ",1pe10.3,                              &
     &        "W/m**2, glob.av.net freshwaterflux ",1pe10.3,"mm/month")
!
!
      zetbar=0.0
      are=0.0
      zetmin=zeta(1,3)
      zetmax=zeta(1,3)
      imin=1
      imax=1
      jmin=3
      jmax=3
      do j=1,jen
        do i=1,ien
          if (wet(i,j,1)<0.5) cycle
          are=are+dlh(i,j)*dpin
          zetbar=zetbar+dlh(i,j)*dpin*zeta(i,j)
          if (zeta(i,j)<zetmin) then
            zetmin=zeta(i,j)
            imin=i
            jmin=j
          end if
          if (zeta(i,j)>zetmax) then
            zetmax=zeta(i,j)
            imax=i
            jmax=j
          end if
        end do
      end do
      zetbar=zetbar/MAX(are,eps)
      write (no6,4010) zetmin,imin,jmin,zetmax,imax,jmax
 4010 format (1x,"Min. sealevel: ",f8.2,2i4," max.: ",f8.2,2i4)
      write (no6,4011) zetbar
 4011 format (1x,"Average sealevel: ",f18.13)
!
!
!     Icearea
!
      arn=0.
      von=0.
      ars=0.
      vos=0.
      do j=1,jen/2
        jn=j
        js=jen+1-j
        do i=1,ien
          if (wet(i,jn,1)>0.5.and.sice(i,jn)>1.e-2) then
            arn=arn+dl(jn)*dphi*0.5
            von=von+dl(jn)*dphi*0.5*sice(i,jn)
          end if
          if (wet(i,js,1)>0.5.and.sice(i,js)>1.e-2) then
            ars=ars+dl(js)*dphi*0.5
            vos=vos+dl(js)*dphi*0.5*sice(i,js)
          end if
        end do
      end do
      write (6,5002) arn*1.e-12,ars*1.e-12,von*1.e-12,vos*1.e-12
 5002 format (" Iceareas: ",4g10.3)
!
!
!
!*    3.        Final processing.
!     -----------------
      salinv=0.
      volinv=0.
      do k=1,ken
        do j=3,jen-2
          do i=1,ien
!PHCONS: take into account salinity in wet column only
            if (k==1) then
              salinv=salinv+wet(i,j,k)*s(i,j,k)*dl(j)                   &
     &               *dphi*0.5*(ddz(i,j,k)+zeta(i,j)-sice(i,j))
!    &               *dphi*0.5*(ddz(i,j,k)+zeta(i,j))
              volinv=volinv+wet(i,j,k)*dl(j)                            &
     &               *dphi*0.5*(ddz(i,j,k)+zeta(i,j)-sice(i,j))
!    &               *dphi*0.5*(ddz(i,j,k)+zeta(i,j))
            else
              salinv=salinv+wet(i,j,k)*s(i,j,k)*dl(j)                   &
     &               *dphi*0.5*(ddz(i,j,k))
              volinv=volinv+wet(i,j,k)*dl(j)*dphi*0.5*(ddz(i,j,k))
            end if
          end do
        end do
      end do
!
      write (no6,*) "Salt inv: ",salinv
      write (no6,*) "Vol inv:  ",volinv
      write (no6,*) "Global mean salinity: ",salinv/volinv

!
!     Overturning diagnostics:
!
      phiout(1)=-phimax1(1)*1.0e-6
      phiout(2)=-phimax2(1)*1.0e-6
      phiout(3)=-phimax3(1)*1.0e-6
      phiout(4)=-phimin1(1)*1.0e-6
      phiout(5)=-phimin2(1)*1.0e-6
      phiout(6)=-phimin3(1)*1.0e-6
      phiout(7)=-phimax1(2)*1.0e-6
      phiout(8)=-phimax2(2)*1.0e-6
      phiout(9)=-phimax3(2)*1.0e-6
      phiout(10)=-phimin1(2)*1.0e-6
      phiout(11)=-phimin2(2)*1.0e-6
      phiout(12)=-phimin3(2)*1.0e-6
      write (no6,*)                                                     &
     &   "  MOC below 700m at         66-46N        16-44N          30S"
      write (no6,9711) "ATL max (NADW)   :",(phiout(i),i=4,6)
      write (no6,9711) "ATL min (AABW)   :",(phiout(i),i=1,3)
      write (no6,9711) "PAC max (outflow):",(phiout(i),i=10,12)
      write (no6,9711) "PAC min (inflow) :",(phiout(i),i=7,9)
 9711 format(5x,a,3f14.5)
      return
      end subroutine outdiaglsg
!     ------------------------------------------------------------------
!
      subroutine outpostsrv
      use lsgvar
      implicit none
!
!     ------------------------------------------------------------------
!
!**** *outpostsrv*
!     
!     by J. von Hardenberg
!     modified by larry 250208
!
!     Purpose.
!     --------
!     *outpostsrv writes output data on LSG_out.srv in srv format
!     ihead(8), yhelp(ien,jen)
!
!     Input.
!     ------
!     ktyear    actual year YYYYY (integer)
!     ktdate    actual date -MMDD (integer)
!
!     Output.
!     -------
!     File with output data.
!
!     Interface.
!     ----------
!     *call* *outpostsrv*(*ktyear*,*ktdate*)
!
!     ------------------------------------------------------------------
!
      integer(kind=4)  :: jcode,jlev,jdate,jtime,jlon,jlat,jdum,jc
      real(kind=4)     :: yhelp(ien,jen)
!
!*    write output to file.
!     ---------------------
!
!     Header for output
!
      jdate = mdatim(1) * 10000 + mdatim(2) * 100 + mdatim(3)
      jtime=0
      jlon=ien
      jlat=jen
      jdum=0

!
!     Write land sea mask for scalar points
!
      jcode=-40
      do jlev=1,ken
         write(93) jcode,jlev,jdate,jtime,jlon,jlat,jdum,jdum
         yhelp(:,:)=wet(:,:,jlev)
         write(93) yhelp
      end do
!
!     Write land sea mask for vector points
!
      jcode=-41
      do jlev=1,ken
         write(93) jcode,jlev,jdate,jtime,jlon,jlat,jdum,jdum
         yhelp(:,:)=wetvec(:,:,jlev)
         write(93) yhelp
      end do
!
!     Write potential  temperature.
!
      jcode=-2
      do jlev=1,ken
         write(93) jcode,jlev,jdate,jtime,jlon,jlat,jdum,jdum
         yhelp(:,:)=t(:,:,jlev)+tkelvin
         write(93) yhelp
      end do
!
!     Salinity.
!
      jcode=-5
      do jlev=1,ken
         write(93) jcode,jlev,jdate,jtime,jlon,jlat,jdum,jdum
         yhelp(:,:)=s(:,:,jlev)
         write(93) yhelp
      end do
!
!     Velocities.
!
      jcode=-3
      do jlev=1,ken
         write(93) jcode,jlev,jdate,jtime,jlon,jlat,jdum,jdum
         yhelp(:,:)=utot(:,:,jlev)
         write(93) yhelp
      end do
!
      jcode=-4
      do jlev=1,ken
         write(93) jcode,jlev,jdate,jtime,jlon,jlat,jdum,jdum
         yhelp(:,:)=vtot(:,:,jlev)
         write(93) yhelp
      end do
!
      jcode=-7
      do jlev=1,ken
         write(93) jcode,jlev,jdate,jtime,jlon,jlat,jdum,jdum
         yhelp(:,:)=w(:,:,jlev)
         write(93) yhelp
      end do
!
!     Barotropic velocities.
!
      jlev=1
      jcode=-37
      write(93) jcode,jlev,jdate,jtime,jlon,jlat,jdum,jdum
      yhelp(:,:)=ub(:,:)
      write(93) yhelp
!
      jlev=1
      jcode=-38
      write(93) jcode,jlev,jdate,jtime,jlon,jlat,jdum,jdum
      yhelp(:,:)=vb(:,:)
      write(93) yhelp
!
!     Surface elevation and ice thickness.
!
      jlev=1
      jcode=-1
      write(93) jcode,jlev,jdate,jtime,jlon,jlat,jdum,jdum
      yhelp(:,:)=zeta(:,:)
      write(93) yhelp
!
      jlev=1
      jcode=-13
      write(93) jcode,jlev,jdate,jtime,jlon,jlat,jdum,jdum
      yhelp(:,:)=sice(:,:)
      write(93) yhelp
!
!     Horizontal barotropic stream function
      jlev=1
      jcode=-27
      write(93) jcode,jlev,jdate,jtime,jlon,jlat,jdum,jdum
      yhelp(:,:)=psi(:,:)
      write(93) yhelp
!
!     Store topography in vector-points.
!
      jlev=1
      jcode=-99
      write(93) jcode,jlev,jdate,jtime,jlon,jlat,jdum,jdum
      yhelp(:,:)=depth(:,:)
      write(93) yhelp
!
!     Store topography in scalar-points.

      jlev=1
      jcode=-98
      write(93) jcode,jlev,jdate,jtime,jlon,jlat,jdum,jdum
      yhelp(:,:)=depp(:,:)
      write(93) yhelp
!
!     Store *fluwat*
!
      jlev=1
      jcode=-65
      write(93) jcode,jlev,jdate,jtime,jlon,jlat,jdum,jdum
      yhelp(:,:)=fluwat(:,:)
      write(93) yhelp
!
!     Store *flukwat*.
!
      jlev=1
      jcode=-67
      write(93) jcode,jlev,jdate,jtime,jlon,jlat,jdum,jdum
      yhelp(:,:)=flukwat(:,:)
      write(93) yhelp
!
!     Store *flukhea*.
!
      jlev=1
      jcode=-68
      write(93) jcode,jlev,jdate,jtime,jlon,jlat,jdum,jdum
      yhelp(:,:)=flukhea(:,:)
      write(93) yhelp
!
!     Store convective ice adjustments.
!
      jlev=1
      jcode=-66
      do jlev=1,ken
         write(93) jcode,jlev,jdate,jtime,jlon,jlat,jdum,jdum
         yhelp(:,:)=convad(:,:,jlev)
         write(93) yhelp
         jcode=-69
      end do
!
!     Store *taux*.
!
      jlev=1
      jcode=-52
      write(93) jcode,jlev,jdate,jtime,jlon,jlat,jdum,jdum
      yhelp(:,:)=taux(:,:)*rhonul
      write(93) yhelp
!
!     Store *tauy*.
!
      jlev=1
      jcode=-53
      write(93) jcode,jlev,jdate,jtime,jlon,jlat,jdum,jdum
      yhelp(:,:)=tauy(:,:)*rhonul
      write(93) yhelp
!
!     Store *tbound*.
!
      jlev=1
      jcode=-92
      write(93) jcode,jlev,jdate,jtime,jlon,jlat,jdum,jdum
      yhelp(:,:)=tbound(:,:)+tkelvin
      write(93) yhelp
!
!     Store *fluhea*.
!
      jlev=1
      jcode=-18
      write(93) jcode,jlev,jdate,jtime,jlon,jlat,jdum,jdum
      yhelp(:,:)=fluhea(:,:)
      write(93) yhelp
!
!     Store *sst-differences lsg-plasim*.
!
      jlev=1
      jcode=-80
      write(93) jcode,jlev,jdate,jtime,jlon,jlat,jdum,jdum
      yhelp(:,:)=fldsst(:,:)
      write(93) yhelp
!
!     Store *ice-differences lsg-plasim*.
!
      jlev=1
      jcode=-81
      write(93) jcode,jlev,jdate,jtime,jlon,jlat,jdum,jdum
      yhelp(:,:)=fldice(:,:)
      write(93) yhelp
!
!     Store *pme-differences lsg-plasim*.
!
      jlev=1
      jcode=-82
      write(93) jcode,jlev,jdate,jtime,jlon,jlat,jdum,jdum
      yhelp(:,:)=fldpme(:,:)
      write(93) yhelp
!
!     Store *taux-differences lsg-plasim*.
!
      jlev=1
      jcode=-83
      write(93) jcode,jlev,jdate,jtime,jlon,jlat,jdum,jdum
      yhelp(:,:)=fldtaux(:,:)
      write(93) yhelp
!
!     Store *taux-differences lsg-plasim*.
!
      jlev=1
      jcode=-84
      write(93) jcode,jlev,jdate,jtime,jlon,jlat,jdum,jdum
      yhelp(:,:)=fldtauy(:,:)
      write(93) yhelp
!
      if(ndiagfl==1) then
!
!     Write temperature tendencies
!
       do jc=1,4
        jcode=-90-jc
        do jlev=1,ken
         write(93) jcode,jlev,jdate,jtime,jlon,jlat,jdum,jdum
         yhelp(:,:)=dtdt(:,:,jlev,jc)
         write(93) yhelp
        end do
       end do
      end if
!
      end subroutine outpostsrv


!     ==================================================================
!     ------------------------------------------------------------------
!
      subroutine outsurf
      use lsgvar
      implicit none
!
!     ------------------------------------------------------------------
!
!**** *outsurf*
!
!     by S. Lorenz 10/2005.
!
!     Purpose.
!     --------
!     *outsurf writes surface output data on tapeunit *nopost*.
!     format is srv4 - service afterburner format:
!     ihead(8), yhelp(ien,jen)
!
!     Input.
!     ------
!     common blocks /lsgfie/ and /lsgsur/
!
!     Output.
!     -------
!     File with output data.
!
!     Interface.
!     ----------
!     *call* *outsurf*
!
!     ------------------------------------------------------------------
!
!
!     Dimension of local variables.
!     -----------------------------
      character(len=7)  filesrf
      integer :: jz1,jz2,jz

      integer(kind=4)  :: jcode,jlev,jdate,jtime,jlon,jlat,jdum
      real(kind=4)     :: yhelp(ien,jen)



!
!
!*    2.        Write output on file *filesrf*.
!     -------------------------------
!
      filesrf='LSG_srf'
      open (nopost,file=filesrf,form="unformatted",position="append")
!
!     Header for output
!
      jlev=oddr(10) ! sl17iw
!     jlev=0
      jdate=mdatim(1) * 10000 + mdatim(2) * 100 + mdatim(3)
      jtime=0
      jlon=ien
      jlat=jen
      jdum=0
!
!     write (no6,*) " OUTSURF: surface output on file ",filesrf         &
!    &, " Date: ",jdate
 9880 format (' OUTSURF: surface output ',a,' on file ',a7,             &
     &        ' YYYYYMMDD :',i10,'  +')

      write (no6,9880) '(real*4)',filesrf,jdate



!
!     Write temperature.
!
      jcode=-2
      write (nopost) jcode,jlev,jdate,jtime,jlon,jlat,jdum,jdum
      yhelp(:,:)=t(:,:,1)+tkelvin
      write (nopost) yhelp
!
!     Write two subsurface layers:
!       o  jz=3 and 5 - 125/275m/L22, 150/450m/L11
!       o  jz=2 and 4 - 75/175m/L22, 75/250m/L11
!
      jz1=2
      jz2=4
!
      jz=jz1
      jlev=oddr(10+jz-1)
      write (nopost) jcode,jlev,jdate,jtime,jlon,jlat,jdum,jdum
      yhelp(:,:)=t(:,:,jz)+tkelvin
      write (nopost) yhelp
!
      jz=jz2
      jlev=oddr(10+jz-1)
      write (nopost) jcode,jlev,jdate,jtime,jlon,jlat,jdum,jdum
      yhelp(:,:)=t(:,:,jz)+tkelvin
      write (nopost) yhelp
!
!     Salinity.
!
      jz=1
      jlev=oddr(10+jz-1)
      jcode=-5
      write (nopost) jcode,jlev,jdate,jtime,jlon,jlat,jdum,jdum
      yhelp(:,:)=s(:,:,1)
      write (nopost) yhelp
!
!     Write two subsurface layers
!
      jz=jz1
      jlev=oddr(10+jz-1)
      write (nopost) jcode,jlev,jdate,jtime,jlon,jlat,jdum,jdum
      yhelp(:,:)=s(:,:,jz)
      write (nopost) yhelp
!
      jz=jz2
      jlev=oddr(10+jz-1)
      write (nopost) jcode,jlev,jdate,jtime,jlon,jlat,jdum,jdum
      yhelp(:,:)=s(:,:,jz)
      write (nopost) yhelp
!
!     Surface elevation and ice thickness.
!
      jz=1
      jlev=oddr(10+jz-1)
      jcode=-1
      write (nopost) jcode,jlev,jdate,jtime,jlon,jlat,jdum,jdum
      yhelp(:,:)=zeta(:,:)
      write (nopost) yhelp
!
      jcode=-13
      write (nopost) jcode,jlev,jdate,jtime,jlon,jlat,jdum,jdum
      yhelp(:,:)=sice(:,:)
      write (nopost) yhelp
!
!     Horizontal barotropic stream function
      jcode=-27
      write (nopost) jcode,jlev,jdate,jtime,jlon,jlat,jdum,jdum
      yhelp(:,:)=psi(:,:)
      write (nopost) yhelp
!
!     Store topography in scalar-points.
!
      jcode=-99
      write (nopost) jcode,jlev,jdate,jtime,jlon,jlat,jdum,jdum
      yhelp(:,:)=wet(:,:,1)
!     yhelp(:,:)=depth(:,:)
      write (nopost) yhelp
!
!     Store topography in vector-points.
!
      jcode=-98
      write (nopost) jcode,jlev,jdate,jtime,jlon,jlat,jdum,jdum
      yhelp(:,:)=wetvec(:,:,1)
!     yhelp(:,:)=depp(:,:)
      write (nopost) yhelp
!
!     Store *fluwat*
!
      jcode=-65
      write (nopost) jcode,jlev,jdate,jtime,jlon,jlat,jdum,jdum
      yhelp(:,:)=fluwat(:,:)
      write (nopost) yhelp
!
!     Store *flukwat*.
!
      jcode=-67
      write (nopost) jcode,jlev,jdate,jtime,jlon,jlat,jdum,jdum
      yhelp(:,:)=flukwat(:,:)
      write (nopost) yhelp
!
!     Store *flukhea*.
!
      jcode=-68
      write (nopost) jcode,jlev,jdate,jtime,jlon,jlat,jdum,jdum
      yhelp(:,:)=flukhea(:,:)
      write (nopost) yhelp
!
!     Store convect ice adjustments.
!
      jcode=-66
      write (nopost) jcode,jlev,jdate,jtime,jlon,jlat,jdum,jdum
      yhelp(:,:)=convad(:,:,1)
      write (nopost) yhelp
!
!     Store *taux*.
!
      jcode=52
      write (nopost) jcode,jlev,jdate,jtime,jlon,jlat,jdum,jdum
      yhelp(:,:)=taux(:,:)*rhonul
      write (nopost) yhelp
!
!     Store *tauy*.
!
      jcode=53
      write (nopost) jcode,jlev,jdate,jtime,jlon,jlat,jdum,jdum
      yhelp(:,:)=tauy(:,:)*rhonul
      write (nopost) yhelp
!
!     Store *tbound*.
!
      jcode=-92
      write (nopost) jcode,jlev,jdate,jtime,jlon,jlat,jdum,jdum
      yhelp(:,:)=tbound(:,:)+tkelvin
      write (nopost) yhelp
!
!     Store *fluhea*.
!
      jcode=-18
      write (nopost) jcode,jlev,jdate,jtime,jlon,jlat,jdum,jdum
      yhelp(:,:)=fluhea(:,:)
      write (nopost) yhelp
!
!*    3.        Store output on file *filpnew*.
!     -------------------------------
      close (nopost)
      end subroutine outsurf
      subroutine press
      use lsgvar
      implicit none
!
!     ------------------------------------------------------------------
!
!**** *press*
!
!     by E. Maier-Reimer.
!     Last modified by U. Mikolajewicz 9/87.
!
!     Purpose.
!     --------
!     Computes the normalized pressure p* from the actual density rh*.
!     The normalized pressure defined as (actual pressure
!     divided by the density *rhonul*.
!
!**   Input.
!     ------
!     Normalized density rh*(3-d-array)  via common /lsgfie/.
!
!     Output.
!     -------
!     Normalized  pressure *p* via common /lsgpre/.
!
!     Interface.
!     ----------
!     *call* *press*.
!
!     Is called from *diva* and *geost2*
!
!     ------------------------------------------------------------------
      integer :: k
      real (kind=8) :: zfact
!
!*    1.      Computation of norm pressure the top layer.
!     ---------------------------------------------------
      zfact = g * du(1)
      p(:,:,1) = zfact * rh(:,:,1)
!
!*    2.      Computation of norm pressure the other layers.
!     ------------------------------------------------------
!
      do k = 2 , ken
        zfact = 0.5 * g * ddu(k-1)
        p(:,:,k) = p(:,:,k-1) + zfact * (rh(:,:,k) + rh(:,:,k-1))
      enddo
      return
      end subroutine press
!     =================
!     subroutine READBF
!     =================

      subroutine readbf(kpar,ktape,kmonth,pfield)
      use lsgvar
      implicit none
      integer :: kpar,ktape,kmonth
!
!     ------------------------------------------------------------------
!
!**** *readbf*.
!
!     by U. Mikolajewicz 1/88.
!     modified by larry 4/2008
!
!**   Purpose.
!     --------
!     *readbf* reads values for fields of boundary values from tapeunit
!     *kfield*.
!
!     Input.
!     ------
!     kpar      field code of the variable.
!     kmonth    relevant month.
!     ktape     number of tapeunit to read from.
!     ntyear    number of time steps/year (common *lsgver*).
!
!*    Output.
!     -------
!     pfield    field with boundary values.
!
!*    Interface.
!     ----------
!     *call* *readbf* (*kpar*,*ktape*,*kmonth*,*pfield*)
!
!
!     Declaration of parameter.
!     -------------------------
      real (kind=8) :: pfield(ien*jen)
!
!     Declaration of local variable.
!
      integer :: inum,jn,numrea,jread,mon,mon2

      real (kind=8) :: anum
      real (kind=8) :: zhead(20),zprel(6)
      integer :: ntcode(200)
      data ntcode/200*0/
      save ntcode
!
!*    1. Header handling.
!
 1000 continue
      if (ktape<=0.or.ktape>100) then
        write (no6,*) " illegal i/o :",ktape
        stop "b-files"
      end if
!
!*    1.1       Read header.
!
!     Number of used header variables.
!
      rewind ktape
      read (ktape,*,end=6666) anum
      inum=nint(anum)
      inum=min0(inum,20)
      inum=max0(inum,8)
!
!     Header.
!
      read (ktape,*) (zhead(jn),jn=1,inum)
!
!
!*    1.2 Check header.
!
      if (nint(zhead(1))/=kpar) then
        write (no6,9801) ktape,kpar,nint(zhead(1))
 9801   format (1x,"mismatch in header of tapeunit",i6,                 &
     &          " expected field code:",i6," found field code:",i6)
        stop "bfield"
      end if
!
      ntcode(ktape)=nint(zhead(3))
!
      if(ntcode(ktape) /= 1 .and. ntcode(ktape) /= 12) then
       write(no6,*) 'ERROR:',ntcode(ktape),' months in data tape ',ktape
       stop
      endif
!
      if (nint(zhead(4))/=ien.or.nint(zhead(5))/=jen) then
        write (no6,9804) kpar,ien,jen,nint(zhead(4)),nint(zhead(5))
 9804   format (1x,"mismatch in header of parameter",i6,                &
     &          " expected          dimensions:",i4," *",i4," found:",  &
     &          i4," *",i4)
        stop "bfield"
      end if
!
      if (nint(zhead(6))/=1) then
        write (no6,9806) kpar
 9806   format (1x,"mismatch in header of parameter",i6,                &
     &          " the input data    are not on an arakawa-e grid. ")
        stop "bfield"
      end if
!
!*    1.3 Determine,whether data must be read.
! 
 1300 continue
!
      if(ntcode(ktape) == 1) then
       numrea=1
      else
       numrea=kmonth
      endif
!
!     *numrea* number of data sets to be read in this call.
!
!      - ntcode(ktape)=zhead(3) is not the tape number but the number
!        of data sets in this file - see DKRZ-Report 2
!      - ntcode(ktape)=zhead(3)=12. for monthly arrays in an annual file,
!        or            zhead(3)= 1. for monthly array in a single file
!      - as before ascii-input is used (Gerrit)
!        (S. Lorenz 06/2005)
!
!*    2. Read data.
!     
      do jread=1,numrea
        read (ktape,*) zprel
        if (kpar==-65)                                                  &
     &  write(6,9715) jread,kpar,ktape,zprel(1),zprel(2)
 9715   format (' READBF: jread=',i3,' kpar=',i3,' ktape=',i3,          &
     &          ' code=',f6.1,' record no: ',f4.1)
!
        mon=nint(zprel(2)+0.0001)
        mon2=jread
        if (mon/=mon2) then
          write (no6,9810) kpar,numrea,mon
 9810     format (1x,"wrong data record for parameter no",i4,           &
     &            " did expect  record for timestep number",i4,         &
     &            " found record number ",i4)
          stop "boundaryfield"
        end if
!
        read (ktape,7100) pfield
 7100   format (4e20.10)
      end do
!
      return
 6666 continue
      write(no6,*) 'no data on tape ktape'
      pfield(:)=0.
      return
      end subroutine readbf
!     ==================
!     subroutine READBOU
!     ==================

      subroutine readbou
      use lsgvar
      implicit none
!
!     ------------------------------------------------------------------
!
!**** *readbou*
!
!     by U. Mikolajewicz 1/88.
!     strongly modified by larry 4/2008 
!
!     Purpose.
!     --------
!     *readbou* computes boundary values using plasim and file-input
!
!     Input.
!     ------
!     nsflu     switch describing wether fluxes shall be read.
!               numbers of tapeunits.
!
!     Output.
!    -------
!     tbound   boundary temperature in Celsius.
!     sbound   boundary salinity in 0/00.
!     taux     zonal wind stress.
!     tauy     meridional wind stress.
!     fluwat   net fresh water flux in m/s. 
!
!     Interface.
!     ----------
!     *call* *readbou*
!
!     Calls subroutine *readbf*,*runoff*,*mktin*,*cfluko*
!                      and *cpfac*
!     Is called by the subroutines start and step
!
!     ------------------------------------------------------------------
!
      integer :: i,j,k,jj
      integer :: kmonth
!
      real(kind=8) :: zztinb(ien,jen)   ! T    (file)
      real(kind=8) :: zziinb(ien,jen)   ! ice  (file)
      real(kind=8) :: zztxinb(ien,jen)  ! taux (file)
      real(kind=8) :: zztyinb(ien,jen)  ! tauy (file)
      real(kind=8) :: zzwinb(ien,jen)   ! pme  (file)
      real(kind=8) :: ztinb(ien,jen)    ! T    (boundary)
      real(kind=8) :: ziinb(ien,jen)    ! ice  (boundary)
      real(kind=8) :: ztxinb(ien,jen)   ! taux (boundary)
      real(kind=8) :: ztyinb(ien,jen)   ! tauy (boundary)
      real(kind=8) :: zwinb(ien,jen)    ! pme  (boundary)
      real         :: ztinp(ien,jen) =0.! T    (plasim)
      real         :: zhflp(ien,jen) =0.! hfl  (plasim)
      real         :: ziinp(ien,jen) =0.! ice  (plasim)
      real         :: ztxinp(ien,jen)=0.! taux (plasim)
      real         :: ztyinp(ien,jen)=0.! tauy (plasim)
      real         :: zwinp(ien,jen) =0.! pme  (plasim)
      real         :: zmld(ien,jen)  =0.! mixed layer depth (lsg)
!
!     factors controling the 'strength' of the coupling
!     ( X = z*X(plasim)+(1.-z)*X(boundary))
!     factors are provided by *sub* cpfac (in coupling interface cpl.f90)
!
      real :: zcpfact    ! coupling strength for T 
      real :: zcpfacice  ! coupling strength for ice
      real :: zcpfactau  ! coupling strength for tau
      real :: zcpfacpme  ! coupling strength for pme
!
      kmonth=mdatim(2)
!
!*    1. get boundary conditions
!
!     get seaice
!
!     a) plasim
!
      call cgeticeo(ziinp)
!
!     b) boundary condition from file
!
      call readbf(-13,no70,kmonth,zziinb)
      ziinb(:,:)=4.*zziinb(:,:)*wet(:,:,1)
!
!     get surface temperature 
!
!     a) PlaSim:
!
      if(ncoupling < 2) then
!
!     coupling via sst
!
!!       zmld(:,:)=dw(1) !! +zeta(:,:)
!
       call cgetssto(ztinp)
      else
!
!     coupling via hfl
!
       call cgethflo(zhflp)
       ztinp(:,:)=t(:,:,1)+wet(:,:,1)*zhflp(:,:)*dt/(dw(1)*rhonul*cp)
      endif
!
!     b) boundary condition from file
!
      call readbf(-2,no71,kmonth,zztinb)
      ztinb(:,:)=zztinb(:,:)-tkelvin
!
!     get wind stress components.
!
!     a) plasim
!
      call cgettauo(ztxinp,ztyinp)
!
!     b) boundary condition from file
!
      call readbf(52,no76,kmonth,ztxinb)
      call readbf(53,no75,kmonth,ztyinb)
!
!     get freshwater flux
!
!     a) plasim
!
      call cgetpmeo(zwinp)
!
!     b) boundary condition from file
!
      call readbf(-65,nowat,kmonth,zzwinb)
      zwinb(:,:)=zzwinb(:,:)
!
!*    compute final runoff and temperature according to the uncoupled lsg
!
      call runoff(zwinb,ztinb)
      call mktin(ztinb,ziinb)
!
!*    diagnose differences (may be used to compute flux corrections)
!
      fldsst(:,:)=ztinb(:,:)-ztinp(:,:)
      fldpme(:,:)=zwinb(:,:)-zwinp(:,:)
      fldice(:,:)=ziinb(:,:)-ziinp(:,:)
      fldtaux(:,:)=ztxinb(:,:)-ztxinp(:,:)
      fldtauy(:,:)=ztyinb(:,:)-ztyinp(:,:)
!
!*    2. flux correction for plasim input
!
      call cflukoa(ztinp,ztxinp,ztyinp,ziinp,zwinp)
!
!*    make final boundary conditions according to the coupling
!
      call cpfac(zcpfact,zcpfactau,zcpfacpme,zcpfacice)
      temin(:,:)=zcpfact*ztinp(:,:)+(1.-zcpfact)*ztinb(:,:)
      taxin(:,:)=zcpfactau*ztxinp(:,:)+(1.-zcpfactau)*ztxinb(:,:)
      tayin(:,:)=zcpfactau*ztyinp(:,:)+(1.-zcpfactau)*ztyinb(:,:)
      sicein(:,:)=zcpfacice*ziinp(:,:)+(1.-zcpfacice)*ziinb(:,:)
      wflin(:,:)=zcpfacpme*zwinp(:,:)+(1.-zcpfacpme)*zwinb(:,:)
!PHCONS:
      if (debug_conservation) then
        sicein(:,:) = 0.0
      endif
!
!*    3. Write boundary fields in internal quantities.
!     
      do j=1,jen
       do i=1,ien
        tbound(i,j)=temin(i,j) 
!
!     Divide wind stress by reference density
!
        taux(i,j)=taxin(i,j)/rhonul
        tauy(i,j)=tayin(i,j)/rhonul
!PHCONS: take out salicc*(sicein-sice) and rescale with new 
!        wet thickness; means: consider gain/loss of seaice thickness
        s(i,j,1) = wet(i,j,1)*(                                         &
     &    (   s(i,j,1)*(dw(1)+zeta(i,j)-sice(i,j))                      &
     &      - salicc*(sicein(i,j)-sice(i,j)) ) /                        &
     &                 (dw(1)+zeta(i,j)-sicein(i,j)) )
!
        sice(i,j)=wet(i,j,1)*sicein(i,j)
        brine(i,j)=0.
!
        if (nsflu==0) cycle
        fluwat(i,j)=wflin(i,j)
       end do
      end do
!
      end subroutine readbou

!     =================
!     subroutine RUNOFF
!     =================

      subroutine runoff(pwinb,ptinb)
      use lsgvar
      implicit none
!
      real(kind=8) :: pwinb(ien,jen),ptinb(ien,jen) 
!
!     Simple runoff model.
!
      real (kind=8) :: rivn(ien,jen),hi1(ien,jen),hi2(ien,jen)
      real (kind=8) :: flualt(ien,jen),flualt2(ien,jen)

      integer, save :: imal0 = 3
      integer :: i,j,k,imal,il,ir,iter,jant

      real (kind=8) :: eps,fmal,fakwvd,traspo,zsum,sum2,sumneu,sumalt
      real (kind=8) :: fluant,precsum1,smax,precsum2,trj,tr,zsm
      real (kind=8) :: precsum3,sumfl,precsum4,sumfla
      real (kind=8) :: sumsto,sumriv,arland,sumh1,sumh2

      eps=1.e-20
      imal=imal0
      fmal=1/float(imal)
      fakwvd=10.*360.*86400.
      traspo=fmal*1.*dt/(86400.*30.)
!
!
!     Smoothe to avoid grid separation.
!
!
      do i=1,ien
        do j=1,jen
          hi1(i,j)=0.0
          hi2(i,j)=0.0
        end do
      end do
!
      do j=3,jen-2
        do i=1,ien
          hi1(i,j)=pwinb(i,j)
        end do
      end do
!
      do j=3,jen-2
        do i=1,ien
          hi2(i,j)=hi1(i,j)
          il=i-1
          if (il<1) il=ien
          ir=i+1
          if (ir>ien) ir=1
          hi2(i,j)=((hi1(il,j-1)+hi1(i,j))+(hi1(i,j-1)+hi1(i,j))        &
     &             +(hi1(i,j+1)+hi1(i,j))+(hi1(ir,j+1)+hi1(i,j)))       &
     &             *0.5/(4.+eps)
        end do
      end do
!
      do j=3,jen-2
        do i=1,ien
          pwinb(i,j)=hi2(i,j)
        end do
      end do
!
!     test Gerrit, to make sure that global P-E ist zero
!
      zsum=0.
      sum2=0.
!
      do j=3,jen-2
        do i=1,ien
          flualt(i,j)=pwinb(i,j)
          zsum=zsum+pwinb(i,j)*dl(j)*0.5*dphi
          sum2=sum2+dl(j)*0.5*dphi
          pwinb(i,j)=pwinb(i,j)+brine(i,j)*(35.-5.)/35.
        end do
      end do
!
      sumneu=0.
      do j=3,jen-2
        do i=1,ien
          pwinb(i,j)=pwinb(i,j)-zsum/sum2  ! to avoid drift because of sea ice
!                                          since brine is included in P-E
          sumneu=sumneu+(dt*pwinb(i,j)+riv(i,j)+stor(i,j))*dl(j)
        end do
      end do
!
      sumalt=sumneu
!
      do j=3,13
        do i=1,ien
          fluant=pwinb(i,j)
          pwinb(i,j)=pwinb(i,j)-fluant
          jant=jen-j+1
          pwinb(i,jant)=pwinb(i,jant)+fluant
        end do
      end do
!
      do j=3,jen-2
        do i=1,ien
          if (iriv(i,j)==0) then
            stor(i,j)=0.
          else
            stor(i,j)=stor(i,j)+pwinb(i,j)*dt
            pwinb(i,j)=0.
          end if
        end do
      end do
!
!     Vanya---------------------------------
      sumneu=0.
      do j=3,jen-2
        do i=1,ien
          sumneu=sumneu+(dt*pwinb(i,j)+riv(i,j)+stor(i,j))*dl(j)
        end do
      end do
!
!       write (6,*) "Precision 1:",(sumalt-sumneu)/sumneu
      precsum1=(sumalt-sumneu)/sumneu
      sumalt=sumneu
!     Vanya---------------------------------
!
      do j=3,jen-2
        do i=1,ien
          if (iriv(i,j)==0) cycle
          if (stor(i,j)>stomax) then
            riv(i,j)=riv(i,j)+(stor(i,j)-stomax)
            stor(i,j)=stomax
          end if
!
!         Protections for wrong P-E field
!
          if (riv(i,j)<0.and.stor(i,j)>0.) then
            smax=dmin1(-riv(i,j),stor(i,j))
            riv(i,j)=riv(i,j)+smax
            stor(i,j)=stor(i,j)-smax
          end if
          if (stor(i,j)<0..and.riv(i,j)>0.) then
            smax=dmin1(riv(i,j),-stor(i,j))
            stor(i,j)=stor(i,j)+smax
            riv(i,j)=riv(i,j)-smax
          end if
          if (stor(i,j)<-stomax) then
            riv(i,j)=riv(i,j)+stor(i,j)+stomax
            stor(i,j)=-stomax
          end if
!
        end do
      end do

!
!     ------------------------------------------------
      sumneu=0.
      do j=3,jen-2
        do i=1,ien
          sumneu=sumneu+(dt*pwinb(i,j)+riv(i,j)+stor(i,j))*dl(j)
        end do
      end do
!      write (6,*) "Precision 2:",(sumalt-sumneu)/sumneu
      precsum2=(sumalt-sumneu)/sumneu
      sumalt=sumneu
!     ---------------------------------------------------
!
      imal=1
      do iter=1,imal
        do j=3,jen-2
          do i=1,ien
            rivn(i,j)=riv(i,j)
          end do
        end do
!
        do j=3,jen-2
          trj=traspo*sqrt(dphi/dl(j))
          if (trj>1.) then
            imal0=imal0+1
            trj=0.95
          end if
          do i=1,ien
            if (iriv(i,j)==0) cycle
            il=i-1
            if (il<1) il=ien
            ir=i+1
            if (ir>ien) ir=1
            tr=trj
            if (ptinb(i,j)<0.) tr=tr/(1.+1.*log(1.-ptinb(i,j)))
!
            if (iriv(i,j)==1) then                 ! to the north, j-1
              rivn(i,j-1)=rivn(i,j-1)+tr*riv(i,j)*dl(j)/dl(j-1)
            else if (iriv(i,j)==2) then            ! to the east, j+1
              rivn(ir,j+1)=rivn(ir,j+1)+tr*riv(i,j)*dl(j)/dl(j+1)
            else if (iriv(i,j)==3) then            ! to the south, j+1
              rivn(i,j+1)=rivn(i,j+1)+tr*riv(i,j)*dl(j)/dl(j+1)
            else if (iriv(i,j)==4) then            ! to the west, j-1
              rivn(il,j-1)=rivn(il,j-1)+tr*riv(i,j)*dl(j)/dl(j-1)
            end if
            rivn(i,j)=rivn(i,j)-riv(i,j)*tr
          end do
        end do
!
        do j=3,jen-2
          do i=1,ien
            riv(i,j)=rivn(i,j)
            if (stor(i,j)<stomax.and.iriv(i,j)>0) then
              zsm=dmin1(stomax-stor(i,j),riv(i,j))*dt*fmal/fakwvd
              stor(i,j)=stor(i,j)+zsm
              riv(i,j)=riv(i,j)-zsm
            end if
          end do
        end do
!
      end do
!
!
!---------------------------------------------------------------------
      sumneu=0.
      do j=3,jen-2
        do i=1,ien
          sumneu=sumneu+(dt*pwinb(i,j)+riv(i,j)+stor(i,j))*dl(j)
        end do
      end do
!      write (6,*) "Precision 3:",(sumalt-sumneu)/sumneu
      precsum3=(sumalt-sumneu)/sumneu
      sumalt=sumneu
!-----------------------------------------------------------------------
!
!
!     Conversion into P-E field
!
      do j=3,jen-2
        do i=1,ien
          if (iriv(i,j)==0) then
            pwinb(i,j)=pwinb(i,j)+0.5*riv(i,j)/dt
            if (iriv(i,j+1)==0) then
              pwinb(i,j+1)=pwinb(i,j+1)+0.5*riv(i,j)*dl(j)            &
     &                      /(dl(j+1)*dt)
              goto 3101
            end if
            if (iriv(i,j-1)==0) then
              pwinb(i,j-1)=pwinb(i,j-1)+0.5*riv(i,j)*dl(j)            &
     &                      /(dl(j-1)*dt)
              goto 3101
            end if
            il=i-1
            if (il<1) il=ien
            if (iriv(il,j-1)==0) then
              pwinb(il,j-1)=pwinb(il,j-1)+.5*riv(i,j)*dl(j)           &
     &                       /(dl(j-1)*dt)
              goto 3101
            end if
            ir=i+1
            if (ir>ien) ir=1
            if (iriv(ir,j+1)==0) then
              pwinb(ir,j+1)=pwinb(ir,j+1)+.5*riv(i,j)*dl(j)           &
     &                       /(dl(j+1)*dt)
              goto 3101
            end if
 3101       riv(i,j)=0.
          end if
        end do
      end do
!
!     -------------------------------------------------------
      sumneu=0.
      sumfl=0.
      do j=3,jen-2
        do i=1,ien
          sumneu=sumneu+(dt*pwinb(i,j)+riv(i,j)+stor(i,j))*dl(j)
          sumfl=sumfl+dt*pwinb(i,j)*dl(j)
        end do
      end do
!     write (6,*) "Precision 4:",(sumalt-sumneu)/sumneu
      precsum4=(sumalt-sumneu)/sumneu
      return
!
      sumalt=sumneu
      sumfla=sumfl
!
!     Smoothe to avoid grid separation in sea level.
!
      write(6,*) "smoothe to avoid grid separation in sea level"
!
      do j=1,jen
        do i=1,ien
          hi1(i,j)=0.
          flualt2(i,j)=pwinb(i,j)
          hi1(i,j)=pwinb(i,j)
        end do
      end do
!
!
      do j=3,jen-2
        do i=1,ien
          hi2(i,j)=hi1(i,j)
          if (wet(i,j,1)<0.5) cycle
          il=i-1
          if (il<1) il=ien
          ir=i+1
          if (ir>ien) ir=1
          hi2(i,j)=(wet(il,j-1,1)*(hi1(il,j-1)+hi1(i,j))/2.+wet(i,j-1,1)&
     &             *(hi1(i,j-1)+hi1(i,j))/2.+wet(i,j+1,1)               &
     &             *(hi1(i,j+1)+hi1(i,j))/2.+wet(ir,j+1,1)              &
     &             *(hi1(ir,j+1)+hi1(i,j))                              &
     &             /2.+(4.-wet(il,j-1,1)-wet(i,j-1,1)-wet(i,j+1,1)      &
     &             -wet(ir,j+1,1)+eps)*hi1(i,j))/(4.+eps)
        end do
      end do
!
      sumsto=0.
      sumriv=0.
      arland=0.
      sumneu=0.
      sumfl=0.
      sumh1=0.
      sumh2=0.
      do j=3,jen-2
        do i=1,ien
          pwinb(i,j)=wet(i,j,1)*hi2(i,j)
          sumneu=sumneu+(dt*pwinb(i,j)+riv(i,j)+stor(i,j))*dl(j)
          sumfl=sumfl+dt*pwinb(i,j)*dl(j)
          sumsto=sumsto+dl(j)*(1.-wet(i,j,1))*stor(i,j)
          sumriv=sumriv+dl(j)*(1.-wet(i,j,1))*riv(i,j)
          arland=arland+dl(j)*(1.-wet(i,j,1))
          sumh1=sumh1+hi1(i,j)
          sumh2=sumh2+hi2(i,j)
        end do
      end do
!
      end subroutine runoff

!     ================
!     subroutine MKTIN
!     ================

      subroutine mktin(ptinb,piinb)
      use lsgvar
      implicit none
!
      real(kind=8) :: ptinb(ien,jen),piinb(ien,jen)  
!
!     ------------------------------------------------------------------
!
!**** subroutine *mktin*
!
!     by larry 04.2008
!     (after mixg from G. Lohmann, MPI, nov. 1998)
!
!**   Purpose.
!     --------
!     *mktin*  computes a boundary temperature from input data
!
!     input.
!     ------
!     piinb   sea ice
!
!     Input/output.
!     -------------
!     ptinb   boundary temperature 
!
!     Interface.
!     ----------
!     *call* *mktin*
!
!     ------------------------------------------------------------------

      integer :: i,j,k,ir,il
      real (kind=8) :: cc
      real (kind=8) :: gg1,ccg,gg2
      real (kind=8) :: zvolt2,zrlap
!
      real (kind=8) :: hi1(ien,jen)
      real (kind=8) :: ztinb(ien,jen)
!
!*    1. set constants (after Prange et al.)
!
      gg1=15.       ! in [W/m**2/K]
      gg2=2.E12     ! in [W/K]
!
!     2. inintialize h1 and ztinb
!
      do j=1,jen
       do i=1,ien
        hi1(i,j)=0.
        ztinb(i,j)=t(i,j,1)
       end do
      end do
!
!     3. compute final sst by using
!        prescribed boundary temperatures (plus diffusion)
!        as in old *mixg*
!
      do j=3,jen-2
       do i=1,ien
        ptinb(i,j)=dmax1(tfreez,ptinb(i,j))
        hi1(i,j)=ztinb(i,j)-ptinb(i,j)
       end do
      end do
!
!     horizontal diffusion *difh*.
!     definition of values for the neighbouring cells.
!     r=right, l=left
!
       do j=3,jen-2
        do i=1,ien
         if (wet(i,j,1)<0.5) cycle
         il=i-1
         if (il<1) il=ien
         ir=i+1
         if (ir>ien) ir=1
         zvolt2=+wet(il,j-1,1)*(hi1(il,j-1)-hi1(i,j))+wet(i,j-1,1)      &
     &         *(hi1(i,j-1)-hi1(i,j))+wet(i,j+1,1)                      &
     &         *(hi1(i,j+1)-hi1(i,j))+wet(ir,j+1,1)                     &
     &         *(hi1(ir,j+1)-hi1(i,j))
         zrlap=-zvolt2/(dl(j)*dphi)
         ztinb(i,j)=ztinb(i,j)                                          &
     &             +(gg1*(ptinb(i,j)-ztinb(i,j))                        &
     &              -gg2*zrlap)*dt/(dw(1)*cp*rhonul)
         if (piinb(i,j)>0.1) then
          ztinb(i,j)=tfreez
         end if
         ztinb(i,j)=dmax1(tfreez,ztinb(i,j))
        end do
       end do
!
!      4. replace old boundary temperature
!
       ptinb(:,:)=ztinb(:,:)
!
      end subroutine mktin
!     ==================================================================
!     ------------------------------------------------------------------
!
      real (kind=8) function rho(s,tp,p)
      implicit none
      real (kind=8), intent(in) :: s,tp,p
!
!     ------------------------------------------------------------------
!
!**** *rho*
!
!     by E. Maier-Reimer.
!     Last modified by U. Mikolajewicz 10/87.
!
!     Purpose.
!     --------
!     The function *rho* computes the density
!     in SI-units according to the UNESCO 83 formula.
!
!**   Input (all via parameter)
!     -------------------------
!     t     pot temperature .
!     s     salinity        .
!     p     pressure        .
!
!     Output.
!     -------
!      Density in SI-units.
!
!
      real (kind=8), parameter :: a0  = 999.842594
      real (kind=8), parameter :: a1  =   6.793952e-2
      real (kind=8), parameter :: a2  =  -9.095290e-3
      real (kind=8), parameter :: a3  =   1.001685e-4
      real (kind=8), parameter :: a4  =  -1.120083e-6
      real (kind=8), parameter :: a5  =   6.536332e-9

      real (kind=8), parameter :: b0  =   8.24493e-1
      real (kind=8), parameter :: b1  =  -4.0899e-3
      real (kind=8), parameter :: b2  =   7.6438e-5
      real (kind=8), parameter :: b3  =  -8.2467e-7
      real (kind=8), parameter :: b4  =   5.3875e-9

      real (kind=8), parameter :: c0  =  -5.72466e-3
      real (kind=8), parameter :: c1  =   1.0227e-4
      real (kind=8), parameter :: c2  =  -1.6546e-6

      real (kind=8), parameter :: d0  =   4.8314e-4

      real (kind=8), parameter :: e0  = 19652.21
      real (kind=8), parameter :: e1  =   148.4206
      real (kind=8), parameter :: e2  =    -2.327105
      real (kind=8), parameter :: e3  =     1.360477e-2
      real (kind=8), parameter :: e4  =    -5.155288e-5

      real (kind=8), parameter :: f0  =  54.6746
      real (kind=8), parameter :: f1  =  -0.603459
      real (kind=8), parameter :: f2  =   1.09987e-2
      real (kind=8), parameter :: f3  =  -6.1670e-5

      real (kind=8), parameter :: g0  =   7.944e-2
      real (kind=8), parameter :: g1  =   1.6483e-2
      real (kind=8), parameter :: g2  =  -5.3009e-4

      real (kind=8), parameter :: h0  =   3.239908
      real (kind=8), parameter :: h1  =   1.43713e-3
      real (kind=8), parameter :: h2  =   1.16092e-4
      real (kind=8), parameter :: h3  =  -5.77905e-7

      real (kind=8), parameter :: ai0 =   2.2838e-3
      real (kind=8), parameter :: ai1 =  -1.0981e-5
      real (kind=8), parameter :: ai2 =  -1.6078e-6

      real (kind=8), parameter :: aj0 =   1.91075e-4

      real (kind=8), parameter :: ak0 =   8.50935e-5
      real (kind=8), parameter :: ak1 =  -6.12293e-6
      real (kind=8), parameter :: ak2 =   5.2787e-8

      real (kind=8), parameter :: am0 =  -9.9348e-7
      real (kind=8), parameter :: am1 =   2.0816e-8
      real (kind=8), parameter :: am2 =   9.1697e-10
!
!     Computation of in-situ temperature
!
      real (kind=8), parameter :: aa1 =   3.6504e-4
      real (kind=8), parameter :: aa2 =   8.3198e-5
      real (kind=8), parameter :: aa3 =  -5.4065e-7
      real (kind=8), parameter :: aa4 =   4.0274e-9
      real (kind=8), parameter :: aa5 =   1.7439e-5
      real (kind=8), parameter :: aa6 =  -2.9778e-7
      real (kind=8), parameter :: aa7 =   8.9309e-7
      real (kind=8), parameter :: aa8 =  -3.1628e-8
      real (kind=8), parameter :: aa9 =   2.1987e-10
      real (kind=8), parameter :: aa10=   4.1057e-9
      real (kind=8), parameter :: aa11=  -1.6056e-10
      real (kind=8), parameter :: aa12=   5.0484e-12
!
      integer :: iter

      real (kind=8) :: t,t2,s3h,rhow,akw,aw,bw,b,akst0,a,akstp,rhst0
      t=tp
!
!     Preliminary in-situ temperature.
!
      t=tp+p*(aa1+aa2*t+p*(aa7+aa8*t+p*aa11))
      do iter=1,1
!
!       Final temperature.
!
        t2=tp+p*(aa1+(aa2+(aa3+aa4*t)*t)*t)+p*(s-35.)*(aa5+aa6*t)       &
     &     +p*p*(aa7+(aa8+aa9*t)*t)-aa10*(s-35.)*p*p+p*p*p*(aa11+aa12*t)
        t=t2
      end do
!
!     Computation of density.
!
      s3h=sqrt(s**3)
      rhow=a0+t*(a1+t*(a2+t*(a3+t*(a4+t*a5))))
      akw=e0+t*(e1+t*(e2+t*(e3+t*e4)))
      aw=h0+t*(h1+t*(h2+t*h3))
      bw=ak0+t*(ak1+t*ak2)
      b=bw+s*(am0+t*(am1+t*am2))
      a=aw+s*(ai0+t*(ai1+ai2*t))+aj0*s3h
      akst0=akw+s*(f0+t*(f1+t*(f2+t*f3)))+s3h*(g0+t*(g1+g2*t))
      akstp=akst0+p*(a+b*p)
      rhst0=rhow+s*(b0+t*(b1+t*(b2+t*(b3+t*b4))))                       &
     &      +d0*s**2+s3h*(c0+t*(c1+c2*t))
      rho=rhst0/(1.-p/akstp)
      return
      end function rho
!     ==================================================================
!     ------------------------------------------------------------------
!
      subroutine rhof1(tp,s,p,rh,kdim)
      implicit none
!
!     ------------------------------------------------------------------
!
!**** *rhof1*
!
!     by E. Maier-Reimer.
!     Last modified by U. Mikolajewicz 10/87.
!
!     Purpose.
!     --------
!     *rhof1* computes density *rh*.
!     the densities are computed according to the
!     UNESCO-1983 formula.
!
!**   Input.
!     ------
!     tp     pot temperature.            )
!     s      salinity.                   ) parameter.
!     p      pressure (bar).             )
!     kdim   array dimension (ien*jen)
!
!     Output.
!     -------
!     rh     density in SI-units.  parameter.
!
!     Interface.
!     ----------
!     *call* *rhof1*(*tp*,*s *,*p *,*rh *).
!
!
!     Dimensions of parameters.
!     tp(ien,jen), s(ien,jen), rh(ien,jen)
!
!     ------------------------------------------------------------------
!
!     Declaration of parameters.
!     --------------------------
!
      integer       :: kdim
      real (kind=8) :: p,s(kdim),tp(kdim),rh(kdim)
!
!     Declaration of local variables.
!
      real (kind=8) :: t(kdim),t2(kdim),s3h(kdim)
!
      real (kind=8), parameter :: a0  = 999.842594
      real (kind=8), parameter :: a1  =   6.793952e-2
      real (kind=8), parameter :: a2  =  -9.095290e-3
      real (kind=8), parameter :: a3  =   1.001685e-4
      real (kind=8), parameter :: a4  =  -1.120083e-6
      real (kind=8), parameter :: a5  =   6.536332e-9

      real (kind=8), parameter :: b0  =   8.24493e-1
      real (kind=8), parameter :: b1  =  -4.0899e-3
      real (kind=8), parameter :: b2  =   7.6438e-5
      real (kind=8), parameter :: b3  =  -8.2467e-7
      real (kind=8), parameter :: b4  =   5.3875e-9

      real (kind=8), parameter :: c0  =  -5.72466e-3
      real (kind=8), parameter :: c1  =   1.0227e-4
      real (kind=8), parameter :: c2  =  -1.6546e-6

      real (kind=8), parameter :: d0  =   4.8314e-4

      real (kind=8), parameter :: e0  = 19652.21
      real (kind=8), parameter :: e1  =   148.4206
      real (kind=8), parameter :: e2  =    -2.327105
      real (kind=8), parameter :: e3  =     1.360477e-2
      real (kind=8), parameter :: e4  =    -5.155288e-5

      real (kind=8), parameter :: f0  =  54.6746
      real (kind=8), parameter :: f1  =  -0.603459
      real (kind=8), parameter :: f2  =   1.09987e-2
      real (kind=8), parameter :: f3  =  -6.1670e-5

      real (kind=8), parameter :: g0  =   7.944e-2
      real (kind=8), parameter :: g1  =   1.6483e-2
      real (kind=8), parameter :: g2  =  -5.3009e-4

      real (kind=8), parameter :: h0  =   3.239908
      real (kind=8), parameter :: h1  =   1.43713e-3
      real (kind=8), parameter :: h2  =   1.16092e-4
      real (kind=8), parameter :: h3  =  -5.77905e-7

      real (kind=8), parameter :: ai0 =   2.2838e-3
      real (kind=8), parameter :: ai1 =  -1.0981e-5
      real (kind=8), parameter :: ai2 =  -1.6078e-6

      real (kind=8), parameter :: aj0 =   1.91075e-4

      real (kind=8), parameter :: ak0 =   8.50935e-5
      real (kind=8), parameter :: ak1 =  -6.12293e-6
      real (kind=8), parameter :: ak2 =   5.2787e-8

      real (kind=8), parameter :: am0 =  -9.9348e-7
      real (kind=8), parameter :: am1 =   2.0816e-8
      real (kind=8), parameter :: am2 =   9.1697e-10
!
!     Computation of in-situ temperature
!
      real (kind=8), parameter :: aa1 =   3.6504e-4
      real (kind=8), parameter :: aa2 =   8.3198e-5
      real (kind=8), parameter :: aa3 =  -5.4065e-7
      real (kind=8), parameter :: aa4 =   4.0274e-9
      real (kind=8), parameter :: aa5 =   1.7439e-5
      real (kind=8), parameter :: aa6 =  -2.9778e-7
      real (kind=8), parameter :: aa7 =   8.9309e-7
      real (kind=8), parameter :: aa8 =  -3.1628e-8
      real (kind=8), parameter :: aa9 =   2.1987e-10
      real (kind=8), parameter :: aa10=   4.1057e-9
      real (kind=8), parameter :: aa11=  -1.6056e-10
      real (kind=8), parameter :: aa12=   5.0484e-12
!
!     integer :: iter
!
!     Preliminary in-situ temperature.
!
      t(:)=tp(:)+p*(aa1+aa2*tp(:)+p*(aa7+aa8*tp(:)+p*aa11))
!
!     do iter=1,1
!
!       Final temperature.
!
        t2(:)=tp(:)+p*(aa1+(aa2+(aa3+aa4*t(:))*t(:))*t(:))              &
     &       +p*(s(:)-35.)*(aa5+aa6*t(:))                               &
     &       +p*p*(aa7+(aa8+aa9*t(:))*t(:))-aa10*(s(:)-35.)             &
     &       *p*p+p*p*p*(aa11+aa12*t(:))
        t(:) = t2(:)
!     end do ! iter
!
!     Computation of the density.
!
      s3h(:)=sqrt(s(:)**3)
      rh(:)=0.
      where ((tp(:).ne.0.).and.(s(:).ne.0.))
        rh(:)=(a0+t(:)*(a1+t(:)*(a2+t(:)*(a3+t(:)*(a4+t(:)*a5))))       &
     &       +s(:)*(b0+t(:)*(b1+t(:)*(b2+t(:)*(b3+t(:)*b4))))           &
     &       +d0*s(:)**2+s3h(:)*(c0+t(:)*(c1+c2*t(:))))                 &
     &       /(1.-p/(p*(h0+t(:)*(h1+t(:)*(h2+t(:)*h3))+s(:)             &
     &       *(ai0+t(:)*(ai1+ai2*t(:)))+aj0*s3h(:)                      &
     &       +(ak0+t(:)*(ak1+t(:)*ak2)+s(:)                             &
     &       *(am0+t(:)*(am1+t(:)*am2)))*p)+e0+t(:)                     &
     &       *(e1+t(:)*(e2+t(:)*(e3+t(:)*e4)))+s(:)                     &
     &       *(f0+t(:)*(f1+t(:)*(f2+t(:)*f3)))+s3h(:)                   &
     &       *(g0+t(:)*(g1+g2*t(:)))))
      endwhere
      return
      end subroutine rhof1
!     ==================================================================
!     ------------------------------------------------------------------
!
      subroutine set_quick
      use lsgvar
      implicit none
!
!     ------------------------------------------------------------------
!
!**** *set_quick*.
!
!     Last modified by A. Paul 11/2000.
!
!     Purpose.
!     --------
!     *set_quick* initializes the common block /lsgadv/
!     containing the number of advection-diffusion sub-cycles *kiterm*,
!     the horizontal diffusion coefficient *a_th*, the vertical
!     diffusion coefficient *a_tv* and the coefficients for the 3rd
!     order advective scheme. The latter depend on the grid as well as
!     on the topography.
!
!     Called from *inp* (after call to *initop*).
!
!**   Input.
!     ------
!     Number of gridpoints in meridional (*jen*) and zonal (*ien*)
!     direction, number of layers (*ken*) from namelist *param*
!     and individual layer thickness in t-points (*ddz*) from
!     commonblock *lsgtrd*.
!
!     Output.
!     -------
!     Commonblock .
!
!     Interface.
!     ----------
!     *call* *set_quick*.
!
!     ------------------------------------------------------------------
!
!     Descriptor.
!     -----------
!
!     ------------------------------------------------------------------
!
!     Local variables.
!     ----------------
!
      integer :: i,j,k,kp2,kp1,km1

      real (kind=8) :: ash,abh,astar,zstar,arange,lambda
      real (kind=8) :: two,pi
      real (kind=8) :: dzt(ken)
!
      write(6,*) " Quick-like advective scheme invoked."
      write(6,*) " Advection of tracers after farrow and stevens (1995):"
      write(6,*) " Explicit in time, third order in space."
!
      first_call=.true.
      debug_advect=.false.
      no_melting=.true.
!
!     Set initial values and constants.
!     ---------------------------------
!
 8810 format (5x,a3,i2,a4,f8.2,a6,i2,a2,g15.6)
 8811 format (8x,a,e20.10)
!
      two=2.0
      pi=4.0*atan(1.0)
!
!     Set depth-dependent horizontal diffusion coefficient (m**2/s).
!     ------------------------------------------------------------------
!
!     Following Bryan and Lewis (1979), the depth dependence of
!     horizontal diffusivity is taken to reflect the ocean"s tendency
!     to diffuse more rapidly at the surface than at depth.
!
      ash=1000.0
      abh=500.0
      write(6,*) ' '
      write(6,8811)  'ash    = ',ash
      write(6,8811)  'abh    = ',abh
      do k=1,ken
        a_th(k)=abh+(ash-abh)*exp(-du(k)/500.0)
        write (6,8810) 'du(',k,') = ',du(k),' a_th(',k,') = ',a_th(k)
      end do
!
!     Set depth-dependent vertical diffusion coefficient (m**2/s).
!     ----------------------------------------------------------------
!
!     Following Bryan and Lewis (1979), the diffusion coefficient
!     increases with depth from 0.3e-04 m**2/s in the upper kilometer
!     to 1.3e-04 (m**2)/s in the deep ocean.
!
!     "astar" is the vertical diffusion coefficient at a reference depth
!     --> Bryan and Lewis (1979): astar  = 0.80e-04
!
!     "arange" defines the range from the surface to the bottom
!     --> Bryan and Lewis (1979): arange = 1.05e-04/pi
!
!     "zstar" is the reference depth
!     --> Bryan and Lewis (1979): zstar  = 2500.0
!
!     "lambda" defines the rate at which the vertical diffusion
!      coefficient varies with depth, particularly near the depth "zstar"
!     --> Bryan and Lewis (1979): lambda = 4.50e-03
!

!     recommended values according to 14C simulations
      if (ken==11) then
        astar=1.80e-04  !1.90e-04
       arange=1.01e-04 !1.00e-04
       zstar=3000.
       lambda=3.40e-03 !3.50e-03
      else if (ken==22) then
        astar=1.44e-04           !!FL 1.3E-4   !!ORI 1.44e-04
        arange=0.78e-04          !!FL 0.3E-4   !!ORI 0.78e-04
        zstar=2800.              !!FL 2800.    !!ORI 2800.
        lambda=3.60e-03          !!FL 4.5E-3   !!ORI 3.60e-03
      end if

      do k=1,ken
        a_tv(k)=astar+arange*atan(lambda*(dw(k)-zstar))
!
!       write (unit=*,fmt="(a3,i2,a4,f8.2,a6,i2,a4,g20.6)") "dw(",k,
!    &         ") = ",dw(k)," a_tv(",k,") = ",a_tv(k)
!
      end do
!
      write(6,*) ' '
      write(6,8811)  'pi     = ',pi
      write(6,8811)  'astar  = ',astar
      write(6,8811)  'zstar  = ',zstar
      write(6,8811)  'arange = ',arange
      write(6,8811)  'lambda = ',lambda
      do k=1,ken
        write(6,8810) 'dw(',k,') = ',dw(k),' a_tv(',k,') = ',a_tv(k)
      end do
!
!     Set number of advection-diffusion sub-cycles to be done.
!     --------------------------------------------------------
!
!     Recommendation: best results are obtained with
!     dt = 3.75 days and kiterm = 8 for ken = 11, and with
!     dt = 2 days and kiterm=15 for ken = 22, respectively.
!
!     ....taken from the namelist
!
      write(6,*) '  Number of advection-diffusion sub-cycles:'
      write(6,'(8x,a,i13)') 'kiterm = ',kiterm
!
!     Set coefficients for 3rd order advective scheme.
!     ------------------------------------------------
!
!     compute coefficients for zonal advection
!
      do i=1,ien
        quick_x(i,1)=1.0
        quick_x(i,2)=1.0
!
        curv_xp(i,1)=0.25
        curv_xp(i,2)=-0.5
        curv_xp(i,3)=0.25
!
        curv_xn(i,1)=0.25
        curv_xn(i,2)=-0.5
        curv_xn(i,3)=0.25
      end do
!
!     compute coefficients for meridional advection
!
      do j=1,jen
        quick_y(j,1)=1.0
        quick_y(j,2)=1.0
!
        curv_yp(j,1)=0.25
        curv_yp(j,2)=-0.5
        curv_yp(j,3)=0.25
!
        curv_yn(j,1)=0.25
        curv_yn(j,2)=-0.5
        curv_yn(j,3)=0.25
      end do
!
!     compute coefficients for vertical advection
!
      do k=1,ken
        kp2=min(k+2,ken)
        kp1=min(k+1,ken)
        km1=max(k-1,1)
        do j=1,jen
          do i=1,ien
!           ! if the grid point is wet, then use the individual layer
!           ! thickness, else use the vertical grid spacing
            if (wet(i,j,kp2)>0.5) then
              dzt(kp2)=ddz(i,j,kp2)
            else
              dzt(kp2)=ddw(kp2)
            end if
            if (wet(i,j,kp1)>0.5) then
              dzt(kp1)=ddz(i,j,kp1)
            else
              dzt(kp1)=ddw(kp1)
            end if
            if (wet(i,j,k)>0.5) then
              dzt(k)=ddz(i,j,k)
            else
              dzt(k)=ddw(k)
            end if
            if (wet(i,j,km1)>0.5) then
              dzt(km1)=ddz(i,j,km1)
            else
              dzt(km1)=ddw(km1)
            end if
            quick_z(i,j,k,1)=two*dzt(kp1)/(dzt(kp1)+dzt(k))
            quick_z(i,j,k,2)=two*dzt(k)/(dzt(kp1)+dzt(k))
!
            curv_zp(i,j,k,1)=two*dzt(k)*dzt(kp1)                        &
     &                       /((dzt(km1)+2.0*dzt(k)+dzt(kp1))           &
     &                       *(dzt(k)+dzt(kp1)))
            curv_zp(i,j,k,2)=-two*dzt(k)*dzt(kp1)                       &
     &                       /((dzt(k)+dzt(kp1))*(dzt(km1)+dzt(k)))
            curv_zp(i,j,k,3)=two*dzt(k)*dzt(kp1)                        &
     &                       /((dzt(km1)+2.0*dzt(k)+dzt(kp1))           &
     &                       *(dzt(km1)+dzt(k)))
!
            curv_zn(i,j,k,1)=two*dzt(k)*dzt(kp1)                        &
     &                       /((dzt(k)+2.0*dzt(kp1)+dzt(kp2))           &
     &                       *(dzt(kp1)+dzt(kp2)))
            curv_zn(i,j,k,2)=-two*dzt(k)*dzt(kp1)                       &
     &                       /((dzt(kp1)+dzt(kp2))*(dzt(k)+dzt(kp1)))
            curv_zn(i,j,k,3)=two*dzt(k)*dzt(kp1)                        &
     &                       /((dzt(k)+2.0*dzt(kp1)+dzt(kp2))           &
     &                       *(dzt(k)+dzt(kp1)))
          end do
        end do
      end do
!
!
      end subroutine set_quick
!     ==================================================================
!     ------------------------------------------------------------------
!
      subroutine start
      use lsgvar
      implicit none
!
!     ------------------------------------------------------------------
!
!**** *start*.
!
!     by E. Maier-Reimer.
!     Last modified by U. Mikolajewicz 9/87.
!     modified by Larry 08/08.
!
!     Purpose.
!     --------
!     *start* controls initialisation of the model.
!     Calls subroutines *inp*,*densin* and *matrix*,*inival2*
!        *readbou* and *outback*.
!
!**   Input.
!     ------
!     itop   number of tapeunit with topographic data.
!     idat                           initial data.
!     itapin                         inputparameters.
!        all via commonblock /lsgtap/.
!        all type integer.
!
!     Output.
!     -------
!     Initialized common blocks.
!
!     Interface.
!     ----------
!     *call* *start*.
!
      character(len=7) file
      integer :: nta,jn
!
!*    1. Input of start data.
!     
      call inp    ! controls initial input
!
!     Decides, from which file the start data must be read.
!
!     iswit=filswit "kleiswi" (dummy) contains the unit of last backup:
!
      open (iswit,file=filswit,form="formatted")
      rewind iswit
      read (iswit,*,err=1101,end=1101) idin,inorm
!
      write(no6,9145) filswit,iswit,idin,inorm
 9145 format(' START: reading file ',a,' from channel no',i3/           &
     &       '        where the last backup is written on:',            &
     &       ' unit idin=',i3,' inorm=',i3)
      goto 1110
!
 1101 continue
      write (no6,*) "Error in reading file ",filswit
      write(6,*) "idin,inorm",idin,inorm
      idin=ibac1                ! ibac1 filbac1 create new restart file.
      inorm=1
      write(6,*) "idin,inorm now ",idin,inorm
      rewind iswit
      write (iswit) idin,inorm
!
 1110 continue
      close (iswit)
!

      if (idin==idat) then      !   ibac1,ibac2,iswit/   84,85,86/
        file=filauf             !   take old restart file "kleiauf" idat=52
      else if (idin==ibac2) then
        file=filbac2            !   take backup restart file "kleiin2"
      else
        file=filbac1            !  "kleiin1", new restart file
      end if
!
!     Reads the start data from file and if no restart( *inorm*=0)
!     writes data on old restart file.
!
      if (inorm==0) then        !   copy input data on filauf
        write (no6,9141) idin,file,inorm
        call inival2(idin,file)
        file=filauf
        write (no6,9142) idat,file,iswit
        call outback(idat,iswit,1,file)
      else                      !   read input data ???
        file=filauf
        write (no6,9141) idat,file,inorm
        call inival2(idat,file)
        call densin
      end if
!
!     write (no6,*) " now in reading file ",idin," ",file," start at ",
!    &              inorm
!
!
 9141 format(' START (INIVAL2): reading restart from channel',i3,       &
     &       '  file=',a,'  inorm=',i3)
 9142 format(' START (OUTBACK): writing restart  on  channel',i3,       &
     &       '  file=',a,'  iswit=',i3)
!     ------------------------------------------------------------------
!
!
!
!*    2.        Computation of the density field.
!     ---------------------------------
      call densin
!
!*    3.        Computation of the matrix for barotropic velocities.
!     ----------------------------------------------------
!     Only if nsve=2.
!
!FL   compute dimensions for matrix and allocate space
!
      call mkmatrx
!
      if (nsve==2) call matrix
!
!
!*    4.        Input of data for restart.
!     --------------------------
      if (inorm==0) goto 5000
      if (idin==idat) then
        file=filauf
      else if (idin==ibac2) then
        file=filbac2
      else
        file=filbac1
      end if
!
!
      write (no6,9141) idin,file,inorm
      call inival2(idin,file)       ! read input data for inorm<>0
      call densin
!
!
!     5.        Prepare for run.
!     ----------------
 5000 continue
      nt=nt0
!
      end subroutine start
!     ==================================================================
!     ------------------------------------------------------------------
!
      subroutine step
      use lsgvar
      implicit none
!
!     ------------------------------------------------------------------
!
!**** *step*.
!
!     by E. Maier-Reimer.
!     last modified by U. Mikolajewicz 6/88.
!
!     purpose.
!     --------
!     *step* controls execution of a timestep.
!
!     Calls subroutines *dens* (computes density distribution).
!                       *diva* (computes velocities, surface elevation
!                               and barotropic stream function).
!                       *readbou* (reads boundary values).
!                       *advect*  (computes the  advection)
!                       and depending on *nsmix*:
!                       *nsmix*=1 : *mix*    (computes changes the
!                                            surface layer, without
!                                            ice, must be driven by
!                                            water temperatures.)
!                       *nsmix*=2 : *mixatt* (computes changes of *t*
!                                            and *s* in the surface
!                                            layer, with ice, must be
!                                            driven by air temperatures)
!                        nsmix*=3 : *mixg    (computes changes of *t* and *s*,
!                                            without ice, driven by
!                                            atmospheric model)
!
!
!      Some of these subroutines call other subroutines (see there).
!
!     ------------------------------------------------------------------
!
      integer :: i,j,k
      real (kind=8) :: zetsum,arrsum,zetm,depthick,dep,ssum,svol
      real (kind=8) :: volum,sal0,thick,tsum
!
!     diagnostics
!
      real (kind=8) :: ztold(ien,jen,ken)
!
      ztold(:,:,:)=t(:,:,:)
      call mktote(t,toteu,toted)
!
!*    1.        Read new boundary values.
!     -------------------------
!
!PHCONS:
      if (debug_conservation) then
        ssum=0.
        svol=0.
        tsum=0.
        volum=0.
        do k=1,ken
          do j=3,jen-2
            do i=1,ien
             if (wet(i,j,k)<0.5) cycle
              if (k==1) then
                volum=(zeta(i,j)-sice(i,j)+ddz(i,j,k))*dlh(i,j)*dphi
              else
                volum=ddz(i,j,k)*dlh(i,j)*dphi
              end if
             ssum=ssum+volum*s(i,j,k)
             tsum=tsum+volum*t(i,j,k)
             svol=svol+volum
            end do
          end do
        end do
        ssum = 0.5*ssum
        svol = 0.5*svol
        tsum = 0.5*tsum
        write(6,'('' BEFORE RBOU  total S, volume, mean S: '',2e18.10,  &
     &               2f18.13)') ssum,svol,ssum/svol,tsum/svol
      endif
      call readbou
!PHCONS:
      if (debug_conservation) then
        ssum=0.
        svol=0.
        tsum=0.
        volum=0.
        do k=1,ken
          do j=3,jen-2
            do i=1,ien
             if (wet(i,j,k)<0.5) cycle
              if (k==1) then
                volum=(zeta(i,j)-sice(i,j)+ddz(i,j,k))*dlh(i,j)*dphi
              else
                volum=ddz(i,j,k)*dlh(i,j)*dphi
              end if
             ssum=ssum+volum*s(i,j,k)
             tsum=tsum+volum*t(i,j,k)
             svol=svol+volum
            end do
          end do
        end do
        ssum = 0.5*ssum
        svol = 0.5*svol
        tsum = 0.5*tsum
        write(6,'(''  AFTER RBOU  total S, volume, mean S: '',2e18.10,  &
     &               2f18.13)') ssum,svol,ssum/svol,tsum/svol
      endif
!
!  run-time plausibility test of slm (and results) in debug mode only (2005/10)
!
!      open (88,file="testwet",form="formatted")
!      rewind 88
!      write(88,*) 'WET:'
!      do j=3,jen-2
!        write (88,'(72i1)') (int(wet(i,j,1)+0.1),i=1,ien)
!      end do
!      close (88)

      call dbgoprt(tbound,ien,jen,ndbox,ndboy,'i.step   / tbound')
!
      zetsum=0.
      arrsum=0.
      do j=3,jen-2
        do i=1,ien
          if (wet(i,j,1)<0.5) cycle
          zetsum=zetsum+zeta(i,j)*0.5*dl(j)*dphi
          arrsum=arrsum+0.5*dl(j)*dphi
        end do
      end do
!
      zetm=zetsum/arrsum
!
!*    2.        Compute mixing at the surface.
!     ------------------------------
!PHCONS:
      if (debug_conservation) then
        ssum=0.
        svol=0.
        tsum=0.
        volum=0.
        do k=1,ken
          do j=3,jen-2
            do i=1,ien
             if (wet(i,j,k)<0.5) cycle
              if (k==1) then
                volum=(zeta(i,j)-sice(i,j)+ddz(i,j,k))*dlh(i,j)*dphi
              else
                volum=ddz(i,j,k)*dlh(i,j)*dphi
              end if
             ssum=ssum+volum*s(i,j,k)
             tsum=tsum+volum*t(i,j,k)
             svol=svol+volum
            end do
          end do
        end do
        ssum = 0.5*ssum
        svol = 0.5*svol
        tsum = 0.5*tsum
        write(6,'('' BEFORE MIXG  total S, volume, mean S: '',2e18.10,  &
     &               2f18.13)') ssum,svol,ssum/svol,tsum/svol
      endif
      if (nsmix==3) call mixg
!PHCONS:
      if (debug_conservation) then
        ssum=0.
        svol=0.
        tsum=0.
        volum=0.
        do k=1,ken
          do j=3,jen-2
            do i=1,ien
             if (wet(i,j,k)<0.5) cycle
              if (k==1) then
                volum=(zeta(i,j)-sice(i,j)+ddz(i,j,k))*dlh(i,j)*dphi
              else
                volum=ddz(i,j,k)*dlh(i,j)*dphi
              end if
             ssum=ssum+volum*s(i,j,k)
             tsum=tsum+volum*t(i,j,k)
             svol=svol+volum
            end do
          end do
        end do
        ssum = 0.5*ssum
        svol = 0.5*svol
        tsum = 0.5*tsum
        write(6,'(''  AFTER MIXG  total S, volume, mean S: '',2e18.10,  &
     &              2f18.13)') ssum,svol,ssum/svol,tsum/svol
      endif
!
!
!*    3.        Compute distribution of density.
!     --------------------------------
!
      call densin
      do j=1,jen
        do i=1,ien
          convef(i,j)=0.
        end do
      end do
!
      do k=1,ken
        if (k==1) then
          dep=0.
        else
          dep=dw(k-1)
        end if
!
        do j=1,jen
          do i=1,ien
            if (k==1) then
              thick=ddz(i,j,k)+zeta(i,j)
            else
              thick=ddz(i,j,k)
            end if
            convef(i,j)=convef(i,j)-wet(i,j,k)*rhonul*g*rh(i,j,k)       &
     &                  *thick*(dep+thick*0.5)/dt
          end do
        end do
      end do
!
!PHCONS:
      if (debug_conservation) then
        ssum=0.
        svol=0.
        tsum=0.
        volum=0.
        do k=1,ken
          do j=3,jen-2
            do i=1,ien
             if (wet(i,j,k)<0.5) cycle
              if (k==1) then
                volum=(zeta(i,j)-sice(i,j)+ddz(i,j,k))*dlh(i,j)*dphi
              else
                volum=ddz(i,j,k)*dlh(i,j)*dphi
              end if
             ssum=ssum+volum*s(i,j,k)
             tsum=tsum+volum*t(i,j,k)
             svol=svol+volum
            end do
          end do
        end do
        ssum = 0.5*ssum
        svol = 0.5*svol
        tsum = 0.5*tsum
        write(6,'('' BEFORE DENS  total S, volume, mean S: '',2e18.10,  &
     &               2f18.13)') ssum,svol,ssum/svol,tsum/svol
      endif
      call dens
!PHCONS:
      if (debug_conservation) then
        ssum=0.
        svol=0.
        tsum=0.
        volum=0.
        do k=1,ken
          do j=3,jen-2
            do i=1,ien
             if (wet(i,j,k)<0.5) cycle
              if (k==1) then
                volum=(zeta(i,j)-sice(i,j)+ddz(i,j,k))*dlh(i,j)*dphi
              else
                volum=ddz(i,j,k)*dlh(i,j)*dphi
              end if
             ssum=ssum+volum*s(i,j,k)
             tsum=tsum+volum*t(i,j,k)
             svol=svol+volum
            end do
          end do
        end do
        ssum = 0.5*ssum
        svol = 0.5*svol
        tsum = 0.5*tsum
        write(6,'(''  AFTER DENS  total S, volume, mean S: '',2e18.10,  &
     &               2f18.13)') ssum,svol,ssum/svol,tsum/svol
      endif
!
      do k=1,ken
        if (k==1) then
          dep=0.
        else
          dep=dw(k-1)
        end if
!
        do j=1,jen
          do i=1,ien
            if (k==1) then
              thick=ddz(i,j,k)+zeta(i,j)
            else
              thick=ddz(i,j,k)
            end if
            convef(i,j)=convef(i,j)+wet(i,j,k)*rhonul*g*rh(i,j,k)       &
     &                  *thick*(dep+thick*0.5)/dt
            convad(i,j,1)=convef(i,j)
          end do
        end do
      end do
!
!
!*    4.        Compute velocities and surface elevation.
!     -----------------------------------------
!
!PHCONS:
      if (debug_conservation) then
        ssum=0.
        svol=0.
        tsum=0.
        volum=0.
        do k=1,ken
          do j=3,jen-2
            do i=1,ien
              if (wet(i,j,k)<0.5) cycle
              if (k==1) then
                volum=(zeta(i,j)-sice(i,j)+ddz(i,j,k))*dlh(i,j)*dphi
              else
                volum=ddz(i,j,k)*dlh(i,j)*dphi
              end if
              ssum=ssum+volum*s(i,j,k)
              tsum=tsum+volum*t(i,j,k)
              svol=svol+volum
            end do
          end do
        end do
        ssum = 0.5*ssum
        svol = 0.5*svol
        tsum = 0.5*tsum
        write(6,'('' BEFORE DIVA  total S, volume, mean S: '',2e18.10,  &
     &               2f18.13)') ssum,svol,ssum/svol,tsum/svol
      endif
      call diva2
!PHCONS:
      if (debug_conservation) then
        ssum=0.
        svol=0.
        tsum=0.
        volum=0.
        do k=1,ken
          do j=3,jen-2
            do i=1,ien
              if (wet(i,j,k)<0.5) cycle
              if (k==1) then
                volum=(zeta(i,j)-sice(i,j)+ddz(i,j,k))*dlh(i,j)*dphi
              else
                volum=ddz(i,j,k)*dlh(i,j)*dphi
              end if
              ssum=ssum+volum*s(i,j,k)
              tsum=tsum+volum*t(i,j,k)
              svol=svol+volum
            end do
          end do
        end do
        ssum = 0.5*ssum
        svol = 0.5*svol
        tsum = 0.5*tsum
        write(6,'(''  AFTER DIVA  total S, volume, mean S: '',2e18.10,  &
     &               2f18.13)') ssum,svol,ssum/svol,tsum/svol
        write(6,'('' BEFORE QUICK total S, volume, mean S: '',2e18.10,  &
     &               2f18.13)') ssum,svol,ssum/svol,tsum/svol
      endif
      call adv_quick
!
!     fix a potential salinity drift in long-term integrations
!     due to small imbalances in the hydrological cycle
!PHCONS: don't do that
      ssum=0.
      svol=0.
      tsum=0.
      volum=0.

      sal0=35.

      do k=1,ken
        do j=3,jen-2
          do i=1,ien
            if (wet(i,j,k)<0.5) cycle
            if (k==1) then
              volum=(zeta(i,j)-sice(i,j)+ddz(i,j,k))*dlh(i,j)*dphi
            else
              volum=ddz(i,j,k)*dlh(i,j)*dphi
            end if
           ssum=ssum+volum*s(i,j,k)
           tsum=tsum+volum*t(i,j,k)
           svol=svol+volum
          end do
        end do
      end do
      ssum = 0.5*ssum
      tsum = 0.5*tsum
      svol = 0.5*svol
!PHCONS: don't tweak salinity
!!      s=s-ssum/svol+sal0
      if (debug_conservation) then
        write(6,'(''  AFTER QUICK total S, volume, mean S: '',2e18.10,  &
     &               2f18.13)') ssum,svol,ssum/svol,tsum/svol
        write(6,*) 'von FL: '
      endif
!
! old coupled mode
!     Now do the atmospheric step:
!     call system("sh puma_asy.run")
!     call system("sh atm_forc.run")
!     ------------------------------------------------------------------
! # endif # comment ???
!
!
!     diagnostics
!
      dtdt(:,:,:,1)=(t(:,:,:)-ztold(:,:,:))/dt
!
      write(6,*) 'End of step: '
      call prfldiag
      call inventory
!
      return
      end subroutine step
!     ==================================================================
!     ------------------------------------------------------------------
!
      subroutine uvtrop2
      use lsgvar
      implicit none
!
!     ------------------------------------------------------------------
!
!**** *uvtrop2*.
!
!     by E. Maier-Reimer.
!     last modified by U. Mikolajewicz 9/87.
!     optimised by E. Kirk 21-Jan-2009
!
!     Purpose.
!     --------
!     *uvtrop2* computes barotropic velocities and stream function
!     according to the equation for the barotropic velocities.
!
!**   Input.
!     ------
!     skal     scaling factors.           )
!     elim     elimination factors.       )   via common /lsgmat/.
!     trisys   tridiagonalsystem.         )
!     zeta     old surface elevation.       (common /lsgsur/).
!
!     Output.
!     -------
!     ub       ) barotropic velocities  (common /lsgfie/).
!     vb       )
!
!     Interface.
!     ----------
!     *call* *uvtrop2*.
!
!
!     Declaration of local variables.
!     -------------------------------
      integer :: i,j,k,l,m

      real (kind=8) :: xr
      real (kind=8) :: zdli,zx,zy
      real (kind=8) :: b(matot)
      real (kind=8) :: x(matot)
!
!
!
!*    1.        Set initial values and constants.
!     ---------------------------------
!
!
      x(:) = 0.0

!*    1.        Computation of the terms.
!     -------------------------
!
!     Computation of the terms containing windstress and surfaceslope
!     and of the horizontal pressure differences in the entire water mass.
!
      m = 0
      do j = 3 , jen-2
        zdli = dli(j)
        do i = 1 , ien
          if (lsmv(i,j)) then
            l = mod(i,ien) + 1
            zx=taux(i,j)+g*zdli*depth(i,j)*(zeta(i,j  )-zeta(l,j  ))
            zy=tauy(i,j)+g*dpin*depth(i,j)*(zeta(l,j+1)-zeta(i,j-1))
            do k = 1 , ken
              zx = zx - zdli * (p(l,j  ,k)-p(i,j  ,k)) * delta(i,j,k)
              zy = zy - dpin * (p(i,j-1,k)-p(l,j+1,k)) * delta(i,j,k)
            enddo
            b(m+1) = zx / depth(i,j) * skal(m+1) * dt
            b(m+2) = zy / depth(i,j) * skal(m+2) * dt
            m = m + 2
          endif
        enddo
      enddo
      b(m+1:matot) = 0.0
!
!
!*   2.        Comput. of the right side of the triangular matrix.
!     ---------------------------------------------------
!
!
!     Elimination.
!
      do j=1,matrx-1
        b(j+1:j+kb) = b(j+1:j+kb) - elim(:,j) * b(j)
      enddo
!
!
!*    3.        Computation of the barotropic velocities.
!     -----------------------------------------
!
!     Solving the matrix equation.
!
      x(matrx) = b(matrx) / trisys(1,matrx)
      do l = matrx-1 , 1 , -1
        xr = dot_product(x(l+1:l+kb),trisys(2:kb+1,l))
        x(l) = (b(l)-xr) / trisys(1,l)
      enddo
!
!
!     Computation of *ub* and *vb*.
!
      ub(:,:) = 0.0
      vb(:,:) = 0.0
      m = 0
      do j = 3 , jen-2
        do i = 1 , ien
          if (lsmv(i,j)) then
             ub(i,j) = x(m+1)
             vb(i,j) = x(m+2)
             m = m + 2
          endif
        enddo
      enddo

      return
      end subroutine uvtrop2
!     ==================================================================
!     ------------------------------------------------------------------
!
      subroutine zetapsi
      use lsgvar
      implicit none
!
!     ------------------------------------------------------------------
!
!**** *zetapsi*.
!
!     by E. Maier-Reimer.
!     last modified by U. Mikolajewicz 9/87.
!     optimised by E. Kirk 21-Jan-2009
!
!     Purpose.
!     --------
!     *zetapsi* computes the surface elevation *zeta* and its
!     time derivative *zetado* from the depth-integrated equation
!     of continuity.
!
!**   Input.
!     ------
!     zeta     old surface elevation.       (common /lsgsur/).
!     ub
!     vb
!
!     Output.
!     -------
!     psi      barotropic stream function.
!     zeta     surface elevation.       (common /lsgsur/).
!     zetado   time derivative of *zeta*.(common /lsgsur/).
!
!     Interface.
!     ----------
!     *call* *zetapsi*.
!
!
!     Declaration of local variables.
!     -------------------------------
      integer :: i,j,k,l,m
      integer :: mindij(2)

      real (kind=8) :: dti,aree,areo,zetame,zetamo,dwe
      real (kind=8) :: zetas(ien,jen),wetdlh(ien,jen)
!!FL
      real (kind=8) :: dzetadt(ien,jen)
      integer :: ip
!!FL
!
!
!*    1.        Computation of the barotropic stream function.
!     ----------------------------------------------
!     with zero at Antarctica.
!
      psi(:,jen-2:jen) = 0.0

      do j = jen-2 , 1 , -1
        do i = 1 , ien
          l = mod(i,ien) + 1
          psi(i,j) = psi(l,j+2) + depth(i,j+1)*dphi*ub(i,j+1)*1.e-6
        enddo
      enddo
!
      psimax = maxval(abs(psi(:,:)))
      mindij = maxloc(abs(psi(:,:)))
      mindi  = mindij(1)
      mindj  = mindij(2)
!
!
!*    2.        New surface elevation.
!     ----------------------
!
!     Computes surface elevation.
!
!FL:  move computation of new zeta to subroutine adv_quick
!     to be more consistent (conservative)
!
      zetas(:,:) = zeta(:,:) ! Save zeta for later use
!
!      do j = 3 , jen-2
!        do i = 1 , ien
!          l = mod(i+ien-2,ien) + 1 ! i-1
!          zeta(i,j) = zeta(i,j) + dt * dlih(i,j)                        &
!     &               * (ub(l,j) * depth(l,j) - ub(i,j) * depth(i,j)     &
!     &               + dpin * (vb(i,j+1) * depth(i,j+1) * dlh(i,j+1)    &
!                     - vb(l,j-1) * depth(l,j-1) * dlh(l,j-1)))
!        enddo
!      enddo
!
!
!     Surface elevation at the poles.
!
!     zeta(:,    1:  2) = 0.0
!     zeta(:,jen-1:jen) = 0.0
!
!     Compensation for cutoff errors of zeta.
!
!     wetdlh(:,:) = wet(:,:,1) * dlh(:,:)
!
!     areo   = sum(wetdlh(:,3:jen-3:2))
!     aree   = sum(wetdlh(:,4:jen-2:2))
!     zetamo = sum(wetdlh(:,3:jen-3:2) * zeta(:,3:jen-3:2))/areo/120.0
!     zetame = sum(wetdlh(:,4:jen-2:2) * zeta(:,4:jen-2:2))/aree/120.0
!
!     zeta(:,3:jen-3:2) = wet(:,3:jen-3:2,1)*(zeta(:,3:jen-3:2)-zetamo)
!     zeta(:,4:jen-2:2) = wet(:,4:jen-2:2,1)*(zeta(:,4:jen-2:2)-zetame)
!
!     Computes time derivative of surface elevation.
!
!     dti = 1.0 / dt
!     dwe = dw(1) / 4.0
!     zetado(:,:) = min(max(dti*(zeta(:,:)-zetas(:,:)),-dwe),dwe)
!
      return
      end subroutine zetapsi

      subroutine inventory
      use lsgvar
      implicit none
      real (kind=8) :: zssum,zsvol,ztsum,zvolum,zsurf,zsurfa,ztdeep
      real (kind=8) :: zvolt,ztup,zvolu,zeu,zed
      integer :: i,j,k
!
      zssum=0.
      zsvol=0.
      ztsum=0.
      zvolum=0.
      zsurf=0.
      zsurfa=0.
      ztdeep=0.
      zvolt=0.
      ztup=0.
      zvolu=0.
      do k=1,ken
       do j=3,jen-2
        do i=1,ien
         if (wet(i,j,k)<0.5) cycle
          if (k==1) then
           zvolum=(zeta(i,j)-sice(i,j)+ddz(i,j,k))*dlh(i,j)*dphi
          else
           zvolum=ddz(i,j,k)*dlh(i,j)*dphi
          end if
          zssum=zssum+zvolum*s(i,j,k)
          ztsum=ztsum+zvolum*(t(i,j,k)+tkelvin)
          zsvol=zsvol+zvolum
          if(k==1) then
           zsurf=zsurf+zeta(i,j)*dlh(i,j)*dphi
           zsurfa=zsurfa+dlh(i,j)*dphi
           if(zeta(i,j) > 0.) then
            ztdeep=ztdeep+(t(i,j,1)+tkelvin)*zeta(i,j)*dlh(i,j)*dphi
            zvolt=zvolt+dlh(i,j)*dphi*zeta(i,j)
            ztup=ztup+(t(i,j,1)+tkelvin)*(ddz(i,j,1)-sice(i,j))         &
     &               *dlh(i,j)*dphi
            zvolu=zvolu+dlh(i,j)*dphi*(ddz(i,j,1)-sice(i,j))
           else
            ztup=ztup+(t(i,j,1)+tkelvin)*(ddz(i,j,1)+zeta(i,j)-sice(i,j))&
     &               *dlh(i,j)*dphi
            zvolu=zvolu+dlh(i,j)*dphi*(ddz(i,j,1)+zeta(i,j)-sice(i,j))
           endif
          elseif(k==2) then
           if(zeta(i,j) < 0.) then
            ztdeep=ztdeep+(t(i,j,2)+tkelvin)*(ddz(i,j,2)+zeta(i,j))     &
     &                   *dlh(i,j)*dphi
            zvolt=zvolt+dlh(i,j)*dphi*(ddz(i,j,2)+zeta(i,j))
            ztup=ztup+(t(i,j,2)+tkelvin)*zeta(i,j)                      &
     &               *dlh(i,j)*dphi
            zvolu=zvolu+dlh(i,j)*dphi*zeta(i,j)
           else
            ztdeep=ztdeep+(t(i,j,2)+tkelvin)*ddz(i,j,2)                 &
     &                   *dlh(i,j)*dphi
            zvolt=zvolt+dlh(i,j)*dphi*ddz(i,j,2)
           endif
          else
           ztdeep=ztdeep+(t(i,j,k)+tkelvin)*ddz(i,j,k)                  &
     &                   *dlh(i,j)*dphi
           zvolt=zvolt+dlh(i,j)*dphi*ddz(i,j,k)
          endif
        end do
       end do
      end do
      zssum = 0.5*zssum
      zsvol = 0.5*zsvol
      ztsum = 0.5*ztsum
      write(6,'(A16,E18.10,A4,E18.10,A6,E18.10)')                       &
     &  'Inventories: S= ',zssum,' T= ',ztsum,' Vol= ',zsvol
      write(6,'(A13,f18.13,A4,f18.13)')                                 &
     &  'Averages: S= ',zssum/zsvol,' T= ',ztsum/zsvol
      write(6,*) 'surface: ',zsurf/zsurfa
      write(6,*) 'upper ocean T: ',ztup/zvolu
      write(6,*) 'deep ocean T: ',ztdeep/zvolt
      call mktote(t,zeu,zed)
      write(6,*) 'upper ocean dE: ',(zeu-toteu)/dt
      write(6,*) 'deep ocean  dE: ',(zed-toted)/dt
      toted=zed
      toteu=zeu
!
      return
      end subroutine inventory

      subroutine mkglmean(pfield,zsum2)
      use lsgvar
      implicit none
      real (kind=8) :: pfield(ien,jen)
      real (kind=8) :: zsum,zarea,zsum2
      integer :: i,j
!
      zsum=0.
      zarea=0.
      do j=3,jen-2
       do i=1,ien
        if (wet(i,j,1)<0.5) cycle
        zsum=zsum+pfield(i,j)*dlh(i,j)*dphi
        zarea=zarea+dlh(i,j)*dphi
       end do
      end do
      zsum=0.5*zsum
      zarea=0.5*zarea
      zsum2=zsum/zarea
      write(6,*) 'sum= ',zsum,' avg= ',zsum2,' area= ',zarea
!
      return
      end subroutine mkglmean

      subroutine mktote(pt,peu,ped)
      use lsgvar
      implicit none
      real (kind=8) :: pt(ien,jen,ken),peu,ped
     
      integer :: i,j,k
!
      peu=0.
      ped=0.
      do k=3,ken
       ped=ped+SUM((pt(:,:,k)+tkelvin)*ddz(:,:,k)*dlh(:,:)*dphi*wet(:,:,k))
      enddo
      do j=3,jen-2
      do i=1,ien
       if(zeta(i,j) > 0) then
        peu=peu+(pt(i,j,1)+tkelvin)*(ddz(i,j,1)-sice(i,j))              &
     &         *dlh(i,j)*dphi*wet(i,j,1)
        ped=ped+(pt(i,j,1)+tkelvin)*zeta(i,j)*dlh(i,j)*dphi*wet(i,j,1)  &
     &         +(pt(i,j,2)+tkelvin)*ddz(i,j,2)*dlh(i,j)*dphi*wet(i,j,2)
       else
        peu=peu+(pt(i,j,1)+tkelvin)*(ddz(i,j,1)-sice(i,j)+zeta(i,j))    &
     &         *dlh(i,j)*dphi*wet(i,j,1)                                &
     &         -(pt(i,j,2)+tkelvin)*zeta(i,j)*dlh(i,j)*dphi*wet(i,j,2)   
        ped=ped+(pt(i,j,2)+tkelvin)*(ddz(i,j,2)+zeta(i,j))              &
     &         *dlh(i,j)*dphi*wet(i,j,2)
       endif
      enddo
      enddo 
      peu=cp*rhonul*peu/SUM(dlh(:,:)*dphi*wet(:,:,1))
      ped=cp*rhonul*ped/SUM(dlh(:,:)*dphi*wet(:,:,1))
!
      return
      end subroutine mktote
     
      subroutine mkmatrx
      use lsgvar
!
!     compute the dimension of the matrix 
!     and set/allocate parameters/arrays accordingly 
!
      integer :: l,ls
! 
!     l gives the number of lines of the matrix.
!
      l=0
      ls=0
      do jla=3,jen-2
       do ilo=1,ien
!
!       only sea-point.
!
        if (numx(ilo,jla) == 0) ls=ls+1
        if (numx(ilo,jla) == 0) cycle
        do kvar=1,2
         l = l + 1
        end do ! kvar
       end do ! ilo
      end do ! jla
!
!     set parameters
!
      matrx = l
      matot = matrx + kb
!
!     allocate
!
      allocate(elim(kb,matrx))
      allocate(trisys(km,matrx))
      allocate(skal(matrx))
!
      write(6,*) 'FLs matrx = ',l,ls
!
      return
      end subroutine mkmatrx

      subroutine prfldiag
      use lsgvar
!
!     debug diagnostic by FL
!
       write(6,*) 'Max s: ',maxval(s,MASK=(wet > 0.5))                  &
     &       ,' at ',maxloc(s,MASK=(wet > 0.5))
       write(6,*) 'Min s: ',minval(s,MASK=(wet > 0.5))                  &
     &       ,' at ',minloc(s,MASK=(wet > 0.5))
       write(6,*) 'Max t: ',maxval(t,MASK=(wet > 0.5))                  &
     &       ,' at ',maxloc(t,MASK=(wet > 0.5))
       write(6,*) 'Min t: ',minval(t,MASK=(wet > 0.5))                  &
     &       ,' at ',minloc(t,MASK=(wet > 0.5))
       write(6,*) 'Max u: ',maxval(utot,MASK=(wetvec > 0.5))            &
     &       ,' at ',maxloc(utot,MASK=(wetvec > 0.5))
       write(6,*) 'Min u: ',minval(utot,MASK=(wetvec > 0.5))            &
     &       ,' at ',minloc(utot,MASK=(wetvec > 0.5))
       write(6,*) 'Max v: ',maxval(vtot,MASK=(wetvec > 0.5))            &
     &       ,' at ',maxloc(vtot,MASK=(wetvec > 0.5))
       write(6,*) 'Min v: ',minval(vtot,MASK=(wetvec > 0.5))            &
     &       ,' at ',minloc(vtot,MASK=(wetvec > 0.5))
       write(6,*) 'Max w: ',maxval(w,MASK=(wet > 0.5))                  &
     &       ,' at ',maxloc(w,MASK=(wet > 0.5))
       write(6,*) 'Min w: ',minval(w,MASK=(wet > 0.5))                  &
     &       ,' at ',minloc(w,MASK=(wet > 0.5))
       write(6,*) 'Max z: ',maxval(zeta,MASK=(wet(:,:,1) > 0.5))        &
     &       ,' at ',maxloc(zeta,MASK=(wet(:,:,1) > 0.5))
       write(6,*) 'Min z: ',minval(zeta,MASK=(wet(:,:,1) > 0.5))        &
     &       ,' at ',minloc(zeta,MASK=(wet(:,:,1) > 0.5))
       write(6,*) 'Max ice: ',maxval(sice,MASK=(wet(:,:,1) > 0.5))      &
     &       ,' at ',maxloc(sice,MASK=(wet(:,:,1) > 0.5))
!
      return
      end subroutine prfldiag

      subroutine mkdzeta(pu,pv,pw,pzeta,pdt,pdzdt)
      use lsgvar
      implicit none
!
!     Declaration of local variables.
!     -------------------------------
!
      integer :: i,j,k,l,m
      integer :: kswitch = 1 
      real(kind=8) :: pu(ien,jen,ken),pv(ien,jen,ken),pw(ien,jen,ken)
      real(kind=8) :: pzeta(ien,jen)
      real(kind=8) :: pdzdt(ien,jen),pdt
      real(kind=8) :: zub(ien,jen),zvb(ien,jen),zzeta(ien,jen)
!
      real (kind=8) :: dti,aree,areo,zetame,zetamo,dwe
      real (kind=8) :: wetdlh(ien,jen)
!
!     compute new zeta
!
      if(kswitch == 0) then
!
!     a) from barotropic velocities
!

       zub(:,:)=0.
       zvb(:,:)=0.
       do j=3,jen-2
        do i=1,ien
         if(wetvec(i,j,1) < 0.5) cycle
         zub(i,j)=SUM(pu(i,j,:)*delta(i,j,:))/depth(i,j)
         zvb(i,j)=SUM(pv(i,j,:)*delta(i,j,:))/depth(i,j)
        enddo
       enddo
       do j=3,jen-2
        do i=1,ien
         l=mod(i+ien-2,ien)+1 ! i-1
         zzeta(i,j)=pzeta(i,j)+pdt*dlih(i,j)                            &
     &             *(zub(l,j)*depth(l,j)-zub(i,j)*depth(i,j)            &
     &              +dpin*(zvb(i,j+1)*depth(i,j+1)*dlh(i,j+1)           &
                          -zvb(l,j-1)*depth(l,j-1)*dlh(l,j-1)))               
        enddo
       enddo
      else
!
!     or b) from upper level divergenz only
!
       do j=3,jen-2
        do i=1,ien
         l=mod(i+ien-2,ien)+1 ! i-1
         zzeta(i,j)=pzeta(i,j)+pdt*(pw(i,j,1)                           &
     &       +dlih(i,j)*(pu(l,j,1)*delta(l,j,1)                         &
     &                  -pu(i,j,1)*delta(i,j,1)                         &
     &                  +dpin*(pv(i,j+1,1)*delta(i,j+1,1)*dlh(i,j+1)    &
     &                        -pv(l,j-1,1)*delta(l,j-1,1)*dlh(l,j-1))))
        enddo
       enddo
      endif
!
!     Surface elevation at the poles.
!
      zzeta(:,    1:  2) = 0.0
      zzeta(:,jen-1:jen) = 0.0
!
!     Compensation for cutoff errors of zeta.
!
      wetdlh(:,:) = wet(:,:,1) * dlh(:,:)
!
      areo   = sum(wetdlh(:,3:jen-3:2))
      aree   = sum(wetdlh(:,4:jen-2:2))
      zetamo = sum(wetdlh(:,3:jen-3:2)*zzeta(:,3:jen-3:2))/areo/120.0
      zetame = sum(wetdlh(:,4:jen-2:2)*zzeta(:,4:jen-2:2))/aree/120.0
!
      zzeta(:,3:jen-3:2)=wet(:,3:jen-3:2,1)*(zzeta(:,3:jen-3:2)-zetamo)
      zzeta(:,4:jen-2:2)=wet(:,4:jen-2:2,1)*(zzeta(:,4:jen-2:2)-zetame)
!
!     Computes time derivative of surface elevation.
!     (limit zeta)
!     
      dti = 1.0/pdt
      dwe = dw(1) / 4.0
      pdzdt(:,:) = min(max(dti*(zzeta(:,:)-pzeta(:,:)),-dwe),dwe)
!
      return
      end subroutine mkdzeta
!!FL
      subroutine mkns(pin,pn,ps)
      use lsgvar
!
      real (kind=8) :: pin(ien,jen),pn,ps
!
      pn=0.
      ps=0.
      do j=1,jen/2
       pn=pn+SUM(pin(:,j)**2)
       ps=ps+SUM(pin(:,jen+1-j)**2)
      enddo
!
      return
      end 
!!FL    

 
