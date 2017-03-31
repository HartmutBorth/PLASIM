! ****************
! * MODUL CATMOD *
! ****************
module catmod

! *********************************************
! *                                           *
! *                 CAT                       *
! *                                           *
! *        Computer Aided Turbulence          *
! *                                           *
! *        Version  1.1    July 2016          *
! *                                           *
! *********************************************
! *  Theoretical Meteorology - KlimaCampus    *
! *         University of Hamburg             *
! *********************************************
! *     Hartmut Borth - Edilbert Kirk         *
! *                                           *
! *           Valerio Lucarini                *
! *********************************************

! *********************************************
! * The latest version of the code and the    *
! * the User's Guide can be downloaded from   *
! * the github repository                     *
! *                                           *
! *  https://github.com/HartmutBorth/QUAD     *
! *********************************************


! *********************************************
! * For more details on parameters, variables *
! * and default values, see the User's Guide, *
! * which is included in the distribution of  *
! * CAT.                                      *
! *********************************************

! *****************
! * Model control *
! *****************
character (256) :: catversion = "July 2016, Version 0.1"

integer :: nshutdown = 0       ! Flag to stop program
logical :: lrsttest  = .false. ! flag set to .true. for restart test,
                               ! parameter nsteps is not written
                               ! to restart file. lrsttest is set to
                               ! true if a file RSTTEST is found in the
                               ! run directory. Call <make rsttest> to
                               ! do the restart test. To clean the run after
                               ! a restart test call <make cleanresttest>


! ***************************
! * Physics and mathematics *
! ***************************
real(8), parameter :: pi    = 4.0d0 * atan(1.0d0)
real(8), parameter :: twopi = pi+pi

complex(8), parameter :: ci = (0.0d0,1.0d0)  ! complex unit


! ****************
! * Input/output *
! ****************

!--- flags and switches
logical :: lrst    = .false.   ! true if <cat_rstini> exists
logical :: lcatnl  = .false.   ! true if <cat_namelist> exists

integer :: ios_catnl  = 1      ! 0 if <cat_namelist> is readable


!--- i/o units
integer, parameter :: nucatnl     = 10  ! cat namelist
integer, parameter :: nuin        = 15  ! initial conditions
integer, parameter :: nutseri     = 20  ! time series
integer, parameter :: nucfl       = 25  ! time series
integer, parameter :: nugp        = 30  ! gridpoint output fields
integer, parameter :: nusp        = 35  ! spectral output fields
integer, parameter :: nurstini    = 40  ! restart for reading initial state
integer, parameter :: nurstfin    = 45  ! restart for writing final state
integer, parameter :: nudiag      = 50  ! statistics of model run

integer, parameter :: nutmp       = 999 ! temporary unit


!--- i/o file names
character (256) :: cat_namelist  = "cat_namelist"
character (256) :: cat_tseri     = "cat_tseri"
character (256) :: cat_cfl       = "cat_cfl"
character (256) :: cat_gp        = "cat_gp"
character (256) :: cat_sp        = "cat_sp"
character (256) :: cat_rstini    = "cat_rstini"
character (256) :: cat_rstfin    = "cat_rstfin"
character (256) :: cat_diag      = "cat_diag"

!--- header for output of service format
integer :: ihead(8)


!--- codes of variables
integer, parameter :: qcde       = 138     ! potential vorticity code
integer, parameter :: qfrccde    = 1138    ! potential vorticity forcing code


!--- codes of variables to be read at initialization
integer, parameter :: ningp  = 5
integer, parameter :: ninsp  = 5
integer :: ingp(ningp) = &
   (/ qcde       , &
      qfrccde    , &
       0         , &
       0         , &
       0           &
   /)
integer :: insp(ninsp) = &
   (/  qcde      , &
       qfrccde   , &
       0         , &
       0         , &
       0           &
   /)

!--- codes of variables to be written to cat_gp or cat_sp
integer, parameter :: noutgp = 5
integer, parameter :: noutsp = 5
integer :: outgp(noutgp) = &
   (/ qcde   , &
       0     , &
       0     , &
       0     , &
       0       &
   /)
integer :: outsp(noutsp) = &
   (/  0  , &
       0  , &
       0  , &
       0  , &
       0    &
   /)


! *******************************
! * Basic diagnostic parameters *
! *******************************

integer :: tsps          ! time steps per second (cpu-time)

real    :: tmstart       ! time at start of cpu-time keeping
real    :: tmstop        ! time at stop of cpu-time keeping
real    :: tmrun         ! total cputime


! *************************************
! * Basic parameters of model physics *
! *************************************

!--- non-dimensional size of fluid domain
real(8) :: lx = twopi   ! in x-direction (scale L_x = X/2pi)
real(8) :: ly = twopi   ! in y-direction (scale L_x = X/2pi)

!--- parameters of evolution equations
real(8) :: alpha = 0.0  ! alpha = 1/R_hat**2 [1/m**2]
                        ! with the non-dimensional Rossby
                        ! radius Ro_hat = Ro/L_x

real(8) :: beta  = 0.0  ! ambient vorticity gradient [1/m*s]

!---------------------!
! dissipation methods !
!---------------------!
integer :: diss_mthd  = 1  ! dissipation method
                           ! 1: Laplacian viscosity and friction
                           !    characterized by the coefficients
                           !    sig and lam and the powers psig and
                           !    plam
                           ! 2: Laplacian viscosity and friction
                           !    characterized by the reciprocal
                           !    damping time-scales rtsig and rtlam,
                           !    the cut-off wave-numbers ksig and klam
                           !    and the powers psig and plam


!----------------------------------!
! Laplacian viscosity and friction !
!----------------------------------!

integer, parameter :: nsig = 2 ! maximum number of terms of Laplacian
                                   ! dissipation on small scales

integer, parameter :: nlam = 2 ! maximum number of terms of Laplacian
                                   ! dissipation on large scales

!--- definition of dissipation on small scales
real(8) :: psig (nsig) =      & ! powers of Laplacian parameterizing
       (/ 1.d0,               & ! small-scale dissipation psig > 0
          4.d0                &
       /)
real(8) :: sig(nsig)   =      & ! coefficients of different powers
       (/ 9.765625d-05,       & ! of the Laplacian [m^2*psig/s]
          9.094947d-11        &
       /)
real(8) :: rtsig(nsig) =      & ! 1/time-scale of different powers
       (/ 1.d-1,              & ! of the Laplacian [1/s]
          1.d+2               &
       /)
integer :: ksig (nsig) =      & ! lower "cut-off" wave number of different
       (/ 320,                & ! powers of Laplacian
          100                 &
       /)

!--- definition of dissipation on large scales
real(8) :: plam (nlam) =      & ! powers of Laplacian modelling viscosity
       (/ -1.d0,              & ! on large scales plam <= 0
          -4.d0               &
       /)
real(8) :: lam(nlam)   =      & ! coefficients of different powers
       (/ 0.d0,               & ! of the Laplacian [m^(-2*psig)/s]
          0.d0                &
       /)
real(8) :: rtlam(nlam) =      & ! 1/time-scale of different powers
       (/ 0.d0,               & ! of the Laplacian [1/s]
          0.d0                &
       /)
integer :: klam (nlam) =      & ! lower "cut-off" wave number of different
       (/ 1,                  & ! powers of Laplacian
          1                   &
       /)

!--- Random number generation
integer, parameter  :: mxseedlen         = 8 ! maximum seed length
integer             :: nseedlen          = 0 ! seed length
integer             :: myseed(mxseedlen) = 0 ! seed given by namelist
integer,allocatable :: seed(:)               ! seed defined by clock

!--- Perturbation of initial conditions
integer :: npert = 0 ! initial perturbation switch
                     ! all perturbations have zero mean-vorticity
                     !
                     ! 0 = no perturbation
                     ! 1 = white noise perturbation
                     ! 2 = white noise symmetric perturbation about the axis
                     !     y = pi (zero mean vorticity in the upper and
                     !     lower part of the fluid domain separately)
                     ! 3 = white noise anti-symmetricy perturbation about
                     !     the axis y = pi (zero mean vorticity in the upper
                     !     and lower power separately)
                     ! 4 = white noise symmetric perturbation about the axis
                     !     x = pi (zero mean vorticity in the upper and
                     !     lower part of the fluid domain separately)
                     ! 5 = white noise anti-symmetricy perturbation about
                     !     the axis x = pi (zero mean vorticity in the upper
                     !     and lower power separately)

real(8) :: apert = 1.0d-5 ! amplitude of perturbation


!--- Forcing
real(8) :: aforc = 0.015     ! amplitude of forcing
real(8) :: tforc = 0.001     ! memory time scale of forcing [s]

integer :: in(1000,2)
integer :: nk
integer :: itau

real(8)    :: ampcoeff

complex(8) :: phi

integer :: nforc = 0 ! forcing switch
                     ! 0 = no forcing
                     !
                     ! 1 = constant frocing defined by external field
                     ! 2 = spectral forcing with constant amplitudes
                     !     forcing ring and random phases
                     ! 3 = markov chain
                     ! 4 = forcing with constant wave numbers
                     !     but random uncorrelated phases (white noise)

integer :: kfmin  = 0   ! min radius = sqrt(kfmin) of spectr. forcing
integer :: kfmax  = 8   ! max radius = sqrt(kfmax) of spectr. forcing

!--- Scaling factors for different parts of the evolution equation
real(8) :: jac_scale = 1.0 ! scaling factor for jacobian


! ***************************************!
! * Basic parameters of numerical scheme !
! ***************************************!

!--- grid in physical space
integer :: ngx = 64    ! number of grid points (x-direction)
integer :: ngy = 64    ! number of grid points (y-direction)

integer :: nxy  = 4096                      ! ngx*ngy
real(8) :: rnx  = 1.56250000000000000E-002  ! 1/ngx
real(8) :: rny  = 1.56250000000000000E-002  ! 1/ngy
real(8) :: rnxy = 2.44140625000000000E-004  ! 1/nxy


!--- grid in spectral space
integer :: nkx  = 21  ! ngx/3 (max. wave number in x-direction)
integer :: nky  = 21  ! ngy/3 (max. wave number in y-direction)
integer :: nfx  = 43  ! 2*nkx+1 (x-dimension in fourier domain)
integer :: nfy  = 42  ! 2*nky   (y-dimension in fourier domain)


!--- Jacobian
integer :: jacmthd  = 1       ! approximation method of Jacobian
                              ! 0 : no Jacobian
                              ! 1 : divergence form (uq)_x + (vq)_y
                              ! 2 : advection form   uq_x  + vq_y
                              ! 3 : invariant form

!--- time integration
integer :: nsteps    = 10000  ! number of time steps to be integrated
integer :: tstep     = 0      ! current time step (since start of runs)
integer :: tstop     = 0      ! last time step of run

integer :: tstp_mthd = 1      ! time stepping method
                              ! 1 = third order Adams-Bashford

integer :: nstdout   = 2500   ! time steps between messages to
                              ! standard output (-1 no output)


integer :: ndiag     = 500    ! time steps between test outputs into
                              ! out_diag (-1 no test outputs)

integer :: ngp       = 100     ! time steps between output of grid point
                               ! fields (-1 = no output)
integer :: nsp       = 100     ! time steps between output of spectral
                               ! fields (-1 = no output)
integer :: ncfl      = 100     ! time steps between cfl-check
                               ! (-1 no cfl check)

integer :: ntseri    = 100     ! time steps between time-series output
                               ! (-1 no time-series)

real(8) :: dt        = 1.d-3   ! length of time step [s]

!--------------------------------------------!
! variables in physical/gridpoint space (GP) !
!--------------------------------------------!

!--- basic model variables
real(8), allocatable :: gq(:,:)    ! vorticity            [1/s]
real(8), allocatable :: gpsi(:,:)  ! stream function      [m^2/s]
real(8), allocatable :: gu(:,:)    ! velocity x-direction [m/s]
real(8), allocatable :: gv(:,:)    ! velocity y-direction [m/s]


!--- jacobian (utility variables)

!--- case 1
real(8), allocatable :: guq(:,:)    ! u*q  [m/s^2]
real(8), allocatable :: gvq(:,:)    ! v*q  [m/s^2]

!--- case 2 (!! variables contain different fields depending on context)
real(8), allocatable :: gqx (:,:) ! dq/dx [1/ms]
real(8), allocatable :: gqy (:,:) ! dq/dy [1/ms]
real(8), allocatable :: gjac(:,:) ! jacobian [s^-2]

!--- case 3
real(8), allocatable :: gv2mu2(:,:) ! v^2-u^2 [m^2/s^2]
real(8), allocatable :: guv(:,:)    ! u*v     [m^2/s^2]

!--- gui
real(4), allocatable :: ggui(:,:) ! single precision transfer


!----------------------------------------------------------!
! variables in spectral space, representation as imaginary !
! and real part (prefix f) used for fourier transformation !
! and complex representation (prefix c) used for time      !
! evolution                                                !
!----------------------------------------------------------!

!--- basic model variables (format used for fourier transform)
real(8), allocatable :: fpsi(:,:) ! stream function
real(8), allocatable :: fu(:,:)   ! velocity x-direction
real(8), allocatable :: fv(:,:)   ! velocity y-direction

!--- jacobian (utility variables)

!--- case 1
real(8), allocatable :: fuq    (:,:) ! u*q  [m/s^2]
real(8), allocatable :: fvq    (:,:) ! v*q  [m/s^2]

!--- case 2 (!! variables contain different fields depending on context)
real(8), allocatable :: fqx(:,:) ! u*dq/dx [s^-2] or dq/dx [1/ms]
real(8), allocatable :: fqy(:,:) ! v*dq/dy [s^-2] or dq/dy [1/ms]

!--- case 3
real(8), allocatable :: fv2mu2 (:,:) ! v^2-u^2 [m^2/s^2]
real(8), allocatable :: fuv    (:,:) ! u*v     [m^2/s^2]


!--- basic model variables (format used for time propagation)
complex(8), allocatable :: cq(:,:)     ! vorticity
complex(4), allocatable :: c4(:,:)     ! vorticity

complex(8), allocatable :: cjac0(:,:)  ! Jacobian at time  0
complex(8), allocatable :: cjac1(:,:)  ! Jacobian at time -1
complex(8), allocatable :: cjac2(:,:)  ! Jacobian at time -2

complex(8), allocatable :: cqfrc(:,:) ! external vorticity forcing

!--- operators in Fourier Space
integer   , allocatable :: kx(:),ky(:)
real(8)   , allocatable :: kx2(:), ky2(:)
real(8)   , allocatable :: k2(:,:), k2pa(:,:), rk2pa(:,:)
real(8)   , allocatable :: kxrk2pa(:,:),kyrk2pa(:,:)
real(8)   , allocatable :: kx2mky2(:,:), kxky(:,:)


!--- linear time propagation (dissipation and beta-term)
complex(8), allocatable :: cli(:,:)    ! linear time propagation


! *******************
! * External moduls *
! *******************

!--- gui communication (guimod)
integer :: ngui    = 1   ! global switch and parameter
                         ! ngui > 0 GUI = on
                         ! ngui specifies moreover the number of
                         ! time-steps between two calls of gui_transfer
                         ! which exchanges data between GUI and CAT

integer :: nguidbg = 0   ! GUI debug mode

integer :: ndatim(6) = 1 ! date/time display
real(4) :: parc(5) = 0.0 ! timeseries display


!--- code of predefined simulations (simmod)
integer :: nsim  = 0 ! 0 no predefined simulation is specified
                     ! nsim > 0 predefined simulation is run
                     !
                     ! available predefined simulations (experiments):
                     ! -----------------------------------------------
                     !    code       name
                     !
                     !     1    decaying_jet01
                     !    51    top-hat wind-forcing
                     !    52    top-hat wind-forcing
                     ! -----------------------------------------------
                     !
                     ! Predefined simulations are given in the dat
                     ! directory as sim_XXXX.nl with XXXX the code
                     ! To activate a given simulation one has to copy
                     ! a given sim_XXXX.nl to sim_namelist in the run
                     ! directory.
                     ! -----------------------------------------------

!--- usermod (usermod)
integer :: nuser = 0 ! 1/0 user mode is switched on/off.
                     ! In user mode it is possible to introduce
                     ! new code into CAT (see chapter <Modifying
                     ! CAT> in the User's Guide).


!--- postprocessing (included in cat)
integer :: npost = 0 ! 1/0 post processing is switched on/off


!--- grid point fields
real(8), allocatable :: gqxm(:)   ! mean vorticity allong x-dir [1/s]
real(8), allocatable :: gqym(:)   ! mean vorticity allong y-dir [1/s]

real(8), allocatable :: gpsixm(:) ! mean stream function allong x-dir [m^2/s]
real(8), allocatable :: gpsiym(:) ! mean stream function allong y-dir [m^2/s]

real(8), allocatable :: guxm(:)   ! mean x-velocity allong x-dir [m/s]
real(8), allocatable :: guym(:)   ! mean x-velocity allong y-dir [m/s]

real(8), allocatable :: gvxm(:)   ! mean y-velocity allong x-dir [m/s]
real(8), allocatable :: gvym(:)   ! mean y-velocity allong y-dir [m/s]

real(4), allocatable :: gguixm(:) ! single precision mean in x-dir for gui
real(4), allocatable :: gguiym(:) ! single precision mean in y-dir for gui

real(4), allocatable :: ggxm(:,:) ! single precision mean in y-dir for gui

end module catmod


! *******
! * CAT *
! *******
program cat
use catmod
implicit none

call cat_prolog

call cat_master

call cat_epilog

stop
end program cat


! #####################################
! #      SUBROUTINES & FUNCTIONIS     #
! #####################################


! *************************
! * SUBROUTINE CAT_PROLOG *
! *************************
subroutine cat_prolog
use catmod
implicit none

call cpu_time(tmstart)
inquire(file="RSTTEST",exist=lrsttest)
call init_rst
call init_diag
call read_resol
call init_pars
call cat_alloc
call init_ops
if (lrst) then
   call check_rst
   call read_rst
endif
call cat_readnl
call cat_open
if (nsim > 0) call simstart
if (nuser > 0) call userstart
call cat_input
call init_jacobian
call init_lintstep
call init_rand
call init_forc
call add_pert
call init_tstep
call init_post
if (ngui > 0) call guistart

return
end subroutine cat_prolog


! *************************
! * SUBROUTINE READ_RESOL *
! *************************
subroutine read_resol
use catmod

character (80) :: argtmp

call get_command_argument(1,argtmp)
read(argtmp,*) ngx

return
end subroutine read_resol


! ************************
! * SUBROUTINE INIT_DIAG *
! ************************
subroutine init_diag
use catmod
implicit none

open(nudiag,file=cat_diag)

write(nudiag, &
 '(" *************************************************")')
write(nudiag, &
 '(" * CAT ",a40," *")') trim(catversion)
write(nudiag, &
   '(" *************************************************",/)')


if (lrsttest) then
   write(nudiag, &
    '(" *************************************************")')
   write(nudiag, &
    '(" * !!! Run CAT in restart test mode !!!")')
   write(nudiag, &
    '(" *************************************************",/)')
endif

return
end subroutine init_diag


! ***********************
! * SUBROUTINE READ_RST *
! ***********************
subroutine read_rst
use catmod

write(nudiag, &
 '(/," ************************************************")')
write(nudiag, &
 '("  Reading parameters from restart file <cat_rstini>")')
write(nudiag, &
 '(" ************************************************",/)')


!--- namelist parameters

if (.not. lrsttest) call get_restart_iarray('nsteps',nsteps,1,1)
call get_restart_iarray('ngp',ngp,1,1)
call get_restart_iarray('nsp',nsp,1,1)
call get_restart_iarray('ingp',ingp,ningp,1)
call get_restart_iarray('insp',insp,ninsp,1)
call get_restart_iarray('outgp',outgp,noutgp,1)
call get_restart_iarray('outsp',outsp,noutsp,1)
call get_restart_iarray('tstp_mthd',tstp_mthd,1,1)
call get_restart_iarray('ncfl',ncfl,1,1)
call get_restart_array ('dt',dt,1,1)
call get_restart_array ('alpha',alpha,1,1)
call get_restart_array ('beta',beta,1,1)
call get_restart_array ('lx',lx,1,1)
call get_restart_array ('ly',ly,1,1)
call get_restart_array ('sig',sig,nsig,1)
call get_restart_array ('psig',psig,nsig,1)
call get_restart_array ('rtsig',rtsig,nsig,1)
call get_restart_iarray('ksig',ksig,nsig,1)
call get_restart_array ('lam',lam,nlam,1)
call get_restart_array ('plam',plam,nlam,1)
call get_restart_array ('rtlam',rtlam,nlam,1)
call get_restart_iarray('klam',klam,nlam,1)
call get_restart_iarray('diss_mthd',diss_mthd,1,1)
call get_restart_iarray('npert',npert,1,1)
call get_restart_array ('apert',apert,1,1)
call get_restart_iarray('nforc',nforc,1,1)
call get_restart_iarray('kfmin',kfmin,1,1)
call get_restart_iarray('kfmax',kfmax,1,1)
call get_restart_array ('aforc',aforc,1,1)
call get_restart_array ('tforc',tforc,1,1)
call get_restart_iarray('myseed',myseed,mxseedlen,1)
call get_restart_iarray('ntseri',ntseri,1,1)
call get_restart_iarray('nstdout',nstdout,1,1)
call get_restart_iarray('jacmthd',jacmthd,1,1)
call get_restart_iarray('ndiag',ndiag,1,1)
call get_restart_array ('jac_scale',jac_scale,1,1)

!--- additional parameters
call get_restart_iarray('ngx',ngx,1,1)
call get_restart_iarray('tstep',tstep,1,1)
call get_restart_iarray('seed',seed,nseedlen,1)

!--- fluid state
call get_restart_carray('cq',cq,nkx+1,nfy+1)
call get_restart_carray('cjac1',cjac1,nkx+1,nfy+1)
call get_restart_carray('cjac2',cjac2,nkx+1,nfy+1)

return
end subroutine read_rst


! *************************
! * SUBROUTINE CAT_READNL *
! *************************
subroutine cat_readnl
use catmod
implicit none


namelist /cat_nl/ nsteps    ,ngp       ,ngui    ,ingp    ,insp       , &
                  ncfl      ,dt        ,alpha   ,beta     ,lx        , &
                  ly        ,sig       ,psig    ,rtsig    ,ksig      , &
                  lam       ,plam      ,rtlam   ,klam     ,diss_mthd , &
                  nforc     ,kfmin     ,kfmax   ,aforc    ,tforc     , &
                  myseed    ,ntseri    ,nstdout ,jacmthd  ,ndiag     , &
                  jac_scale ,nsp       ,outgp   ,outsp    ,tstp_mthd , &
                  nsim      ,nuser     ,npost   ,npert    ,apert     , &
                  nguidbg

inquire(file=cat_namelist,exist=lcatnl)

if (lcatnl) then
  open(nucatnl,file=cat_namelist,iostat=ios_catnl)
  read(nucatnl,cat_nl)
endif

!--- check consistency of namelist parameters
if (jacmthd .lt. 0 .or. jacmthd .gt. 3) then
   write(nudiag, &
   '(/," ************************************************")')
   write(nudiag,*) "jacmthd = ", jacmthd ," does not exist"
   write(nudiag,*) "use default instead (jacmthd = 1)"
   write(nudiag,'("                         ")')
   write(nudiag, &
   '(" ************************************************",/)')
   jacmthd = 1
endif


return
end subroutine cat_readnl


! ***********************
! * SUBROUTINE CAT_OPEN *
! ***********************
subroutine cat_open
use catmod
implicit none


open(nudiag,file=cat_diag)

if (ntseri .ge. 0)  open(nutseri,file=cat_tseri)
if (ncfl .gt. 0)    open(nucfl,file=cat_cfl)
if (ngp .gt. 0)     open(nugp,file=cat_gp,form='unformatted')
if (nsp .gt. 0)     open(nusp,file=cat_sp,form='unformatted')

return
end subroutine cat_open


! ************************
! * SUBROUTINE INIT_PARS *
! ************************
subroutine init_pars
use catmod
implicit none

ngy  = ngx           ! restrict model to square domains

nxy  = ngx * ngy     ! # of gridpoints
nkx  = ngx / 3       ! max wavenumber in x
nky  = ngy / 3       ! max wavenumber in y
nfx  = nkx * 2 + 1   ! x dimension in fourier domain
nfy  = nky * 2       ! y dimension in fourier domain
rnx  = 1.0d0 / ngx   ! reverse of ngx
rny  = 1.0d0 / ngy   ! reverse of ngy
rnxy = 1.0d0 / nxy   ! reverse of gridpoint number


write(nudiag, &
 '(" *************************************************")')
write(nudiag, &
 '(" * Values of basic parameters                    *")')
write(nudiag, &
 '(" *************************************************")')
write (nudiag,*) "nkx  = ", ngx/3
write (nudiag,*) "nky  = ", ngy/3
write (nudiag,*) "nfx  = ", nfx
write (nudiag,*) "nfy  = ", nfy
write (nudiag,*) "rnx  = ", rnx
write (nudiag,*) "rny  = ", rny
write (nudiag,*) "rnxy = ", rnxy
write(nudiag, &
 '(" *************************************************",/)')

return
end subroutine init_pars


! ************************
! * SUBROUTINE CAT_ALLOC *
! ************************
subroutine cat_alloc
use catmod
implicit none

integer :: j

!-----------------------------
! 1D real and complex vectors
!-----------------------------
allocate(kx(0:nkx))   ; kx(:)   = [(j,j=0,nkx)]
allocate(ky(0:nfy))   ; ky(:)   = [(j,j=0,nky),(j,j=-nky,-1)]
allocate(kx2(0:nkx))  ; kx2(:)  = 0.0
allocate(ky2(0:nfy))  ; ky2(:)  = 0.0


!----------------------------
! 2D real and complex arrays
!----------------------------

!--- grid point space
allocate(gq(1:ngx,1:ngy))   ; gq(:,:)   = 0.0  ! vorticity
allocate(gu(1:ngx,1:ngy))   ; gu(:,:)   = 0.0  ! velocity in x-dir.
allocate(gv(1:ngx,1:ngy))   ; gv(:,:)   = 0.0  ! velocity in y-dir.

allocate(gpsi(1:ngx,1:ngy)) ; gpsi(:,:) = 0.0  ! stream function

allocate(ggui(1:ngx,1:ngy)) ; ggui(:,:)  = 0.0  ! GUI transfer

!--- spetral space
allocate(fu(0:nfx,0:nfy))   ; fu(:,:)   = 0.0 ! u
allocate(fv(0:nfx,0:nfy))   ; fv(:,:)   = 0.0 ! v

allocate(k2(0:nkx,0:nfy))      ; k2(:,:)      = 0.0 ! Laplacian
allocate(k2pa(0:nkx,0:nfy))    ; k2pa(:,:)    = 0.0 ! modified Laplacian
allocate(rk2pa(0:nkx,0:nfy))   ; rk2pa(:,:)   = 0.0 ! modified Laplacian^-1
allocate(kxrk2pa(0:nkx,0:nfy)) ; kxrk2pa(:,:) = 0.0 ! q --> v
allocate(kyrk2pa(0:nkx,0:nfy)) ; kyrk2pa(:,:) = 0.0 ! q --> u
allocate(kx2mky2(0:nkx,0:nfy)) ; kx2mky2(:,:)  = 0.0 ! utility for Jacobian 3
allocate(kxky  (0:nkx,0:nfy))  ; kxky  (:,:)  = 0.0 ! utility for Jacobian 3



allocate(cli(0:nkx,0:nfy))  ; cli(:,:)   = (0.0,0.0) ! linear time propagator

allocate(cq(0:nkx,0:nfy))   ; cq(:,:)    = (0.0,0.0) ! vorticity
allocate(c4(0:nkx,0:nfy))   ; c4(:,:)    = (0.0,0.0) ! vorticity

allocate(cqfrc(0:nkx,0:nfy)); cqfrc(:,:) = (0.0,0.0) ! ext. vorticity force
allocate(cjac0(0:nkx,0:nfy)); cjac0(:,:) = (0.0,0.0) ! Jacobian at time level  0
allocate(cjac1(0:nkx,0:nfy)); cjac1(:,:) = (0.0,0.0) ! Jacobian at time level -1
allocate(cjac2(0:nkx,0:nfy)); cjac2(:,:) = (0.0,0.0) ! Jacobian at time level -2


!allocate(aptmp(0:nkx,0:nfy)); aptmp(:,:) = 0.0 ! amplitude of complex field
!allocate(patmp(0:nkx,0:nfy)); patmp(:,:) = 0.0 ! phase of complex field

return
end subroutine cat_alloc


! ****************************
! * SUBROUTINE INIT_JACOBIAN *
! ****************************
subroutine init_jacobian
use catmod
implicit none
select case(jacmthd)
case(0)
   print *, "init_jacobian case 0"
case(1)
   !--- allocate fields used by jacobian

   !--- grid point fields
   allocate(guq(1:ngx,1:ngy)) ; guq(:,:) = 0.0  ! u*q
   allocate(gvq(1:ngx,1:ngy)) ; gvq(:,:) = 0.0  ! v*q

   !--- spectral fields
   allocate(fuq(0:nfx,0:nfy)) ; fuq(:,:) = 0.0  ! u*q
   allocate(fvq(0:nfx,0:nfy)) ; fvq(:,:) = 0.0  ! v*q
case(2)
   !--- allocate fields used by jacobian

   !--- grid point fields
   allocate(gqx(1:ngx,1:ngy))  ; gqx(:,:)  = 0.0  ! qx
   allocate(gqy(1:ngx,1:ngy))  ; gqy(:,:)  = 0.0  ! qy
   allocate(gjac(1:ngx,1:ngy)) ; gjac(:,:) = 0.0  ! jacobian in grid-space

   !--- spectral fields
   allocate(fqx(0:nfx,0:nfy)) ; fqx(:,:) = 0.0  ! u*qx or qx
   allocate(fqy(0:nfx,0:nfy)) ; fqy(:,:) = 0.0  ! v*qy or qy
case(3)
   !--- allocate fields used by jacobian
   !--- grid point fields
   allocate(gv2mu2(1:ngx,1:ngy)) ; gv2mu2(:,:) = 0.0  ! v^2-u^2
   allocate(guv   (1:ngx,1:ngy)) ; guv   (:,:) = 0.0  ! u*v
   !--- spetral fields
   allocate(fv2mu2(0:nfx,0:nfy)) ; fv2mu2(:,:) = 0.0  ! v^2-u^2
   allocate(fuv   (0:nfx,0:nfy)) ; fuv   (:,:) = 0.0  ! u*v
   print *, "init_jacobian case 3"
end select

return
end subroutine init_jacobian


! ****************************
! * SUBROUTINE STOP_JACOBIAN *
! ****************************
subroutine stop_jacobian
use catmod
implicit none

select case (jacmthd)

case (1)
   !--- deallocate fields used by jacobian

   !--- grid point fields
   deallocate(guq)
   deallocate(gvq)

   !--- spetral fields
   deallocate(fuq)
   deallocate(fvq)

case (2)
   !--- deallocate fields used by jacobian

   !--- grid point fields
   deallocate(gqx)
   deallocate(gqy)
   deallocate(gjac)

   !--- spetral fields
   deallocate(fqx)
   deallocate(fqy)

case (3)
   !--- deallocate fields used by jacobian

   !--- grid point fields
   deallocate(gv2mu2)
   deallocate(guv)

   !--- spetral fields
   deallocate(fv2mu2)
   deallocate(fuv)

end select

return
end subroutine stop_jacobian


! ************************
! * SUBROUTINE CAT_INPUT *
! ************************
subroutine cat_input
use catmod
implicit none

logical         :: lexist
integer         :: kcode,jj,kk,kkmax
real(8)         :: gtmp(1:ngx,1:ngy)
real(8)         :: ftmp(0:nfx,0:nfy) ! temporary field for i/o

character(2)    :: gtp
character(256)  :: fname

integer, parameter      :: ngtp         = 2
character(2), parameter :: gtp_ar(ngtp) = (/ "GP" , "SP" /)


if (.not. (tstep .eq. 0) ) return

do jj = 1,ngtp
   gtp = gtp_ar(jj)
   select case (gtp)
   case ("GP")
      kkmax = ningp
   case ("SP")
      kkmax = ninsp
   end select
   do kk = 1,kkmax
      select case (gtp)
      case ("GP")
         kcode = ingp(kk)
      case ("SP")
         kcode = insp(kk)
      end select
      if (kcode .gt. 0) then
         call checkvar(kcode,gtp,lexist)
         if (lexist) then
            select case (kcode)
            case (qcde)
               select case (gtp)
               case ("GP")
                  call read_gp(kcode,gtp,gq)
                  call grid_to_fourier(gq,cq,nfx,nfy,ngx,ngy)
               case ("SP")
                  call read_sp(kcode,gtp,ftmp)
                  call f2c(ftmp,cq)
               end select
               cq(0,0) = (0.0,0.0)
            case (qfrccde)
               select case (gtp)
               case ("GP")
                  call read_gp(kcode,gtp,gtmp)
                  call grid_to_fourier(gtmp,cqfrc,nfx,nfy,ngx,ngy)
               case ("SP")
                  call read_sp(kcode,gtp,ftmp)
                  call f2c(ftmp,cqfrc)
               end select
               cqfrc(0,0) = (0.0,0.0)
            end select
         else
            call cat_fname(gtp,kcode,fname)
            write(nudiag, &
            '(" *************************************************")')
            write(nudiag, &
            '(" *   File ", a ," not found, use default")') trim(fname)
            write(nudiag, &
            '(" *************************************************",/)')
         endif
      endif
   enddo
enddo

return
end subroutine cat_input


! ************************
! * SUBROUTINE CAT_FNAME *
! ************************
subroutine cat_fname(gtp,kcode,fname)
use catmod
implicit none

character(2)   :: gtp
character(256) :: fname
integer :: kcode

if (ngx < 100) then
  write(fname,'(a2,i2.2,"_var",i4.4,".srv")') gtp,ngx,kcode
elseif (ngx < 1000) then
  write(fname,'(a2,i3.3,"_var",i4.4,".srv")') gtp,ngx,kcode
else
  write(fname,'(a2,i4.4,"_var",i4.4,".srv")') gtp,ngx,kcode
endif

fname = trim(fname)

return
end subroutine cat_fname


! ***********************
! * SUBROUTINE CHECKVAR *
! ***********************
subroutine checkvar(kcode,gtp,lexist)
use catmod

logical :: lexist
character(2)   :: gtp
character(256) :: fname
integer :: kcode

call cat_fname(gtp,kcode,fname)
inquire(file=trim(fname),exist=lexist)

return
end subroutine checkvar


! **********************
! * SUBROUTINE READ_GP *
! **********************
subroutine read_gp(kcode,gtp,gpvar)
use catmod
implicit none

character(2)   :: gtp
character(256) :: fname
integer        :: kcode
real(8)        :: gpvar(1:ngx,1:ngy)

call cat_fname(gtp,kcode,fname)

open(nuin,file=fname,form='unformatted')

write(nudiag, &
'(" *************************************************")')
write(nudiag,'(" * Reading var",i4.4 " from file " a)') &
      kcode,trim(fname)
write(nudiag, &
'(" *************************************************",/)')

read (nuin) ihead
read (nuin) gpvar(:,:)
close(nuin)

return
end subroutine read_gp


! **********************
! * SUBROUTINE READ_SP *
! **********************
subroutine read_sp(kcode,gtp,spvar)
use catmod
implicit none

character(2)   :: gtp
character(256) :: fname
integer        :: kcode
real(8)        :: spvar(0:nfx,0:nfy)

call cat_fname(gtp,kcode,fname)

open(nuin,file=fname,form='unformatted')

write(nudiag, &
'(" *************************************************")')
write(nudiag,'(" * Reading var",i4.4 " from file " a)') &
      kcode,trim(fname)
write(nudiag, &
'(" *************************************************",/)')

read (nuin) ihead
read (nuin) spvar(:,:)
close(nuin)

return
end subroutine read_sp


! ******************
! * SUBROUTINE F2C *
! ******************
subroutine f2c(fvar,cvar)
use catmod
implicit none

integer     :: i,j
real(8)     :: fvar(0:nfx,0:nfy)
complex(8)  :: cvar(0:nkx,0:nfy)

do j = 0, nfy
   do i = 0, nkx
      cvar(i,j) = cmplx(fvar(i+i,j),fvar(i+i+1,j))
   enddo
enddo

return
end subroutine f2c


! ******************
! * SUBROUTINE C2F *
! ******************
subroutine c2f(cvar,fvar)
use catmod
implicit none

integer     :: i,j
real(8)     :: fvar(0:nfx,0:nfy)
complex(8)  :: cvar(0:nkx,0:nfy)



do j = 0, nfy
   do i = 0, nkx
      fvar(i+i  ,j) = real (cvar(i,j))
      fvar(i+i+1,j) = aimag(cvar(i,j))
   enddo
enddo

return
end subroutine c2f


! ****************************
! * SUBROUTINE INIT_LINTSTEP *
! ****************************
subroutine init_lintstep
use catmod
implicit none

real(8)    :: diss_sig,diss_lam,b_term
complex(8) :: arg
integer    :: jx,jy,n

!--------------------------------------------------------------
! Time propagator of linear part of differential equation
!
! cli = exp(arg)                      with
!
! arg = -dt*[diss_sig + diss_lam + beta_term]
!
!  with
!     diss_sig  = -sum_j=1,nsig sig(j)*k^(2*psig(j))
!     diss_lam  = -sum_j=1,nlam lam(j)*k^(2*plam(j))
!     b_term    = i*beta*kx/(k^2 + alpha)
!--------------------------------------------------------------

do jy = 0, nfy
   do jx = 0, nkx
      diss_sig = 0.0
      diss_lam = 0.0
      do n = 1,nsig
         diss_sig = diss_sig - sig(n)*k2(jx,jy)**psig(n)
      enddo
      do n = 1,nlam
         diss_lam = diss_lam - lam(n)*k2(jx,jy)**plam(n)
      enddo
      b_term    = beta * kxrk2pa(jx,jy)
      arg       = dt * cmplx(diss_sig+diss_lam,b_term)
      cli(jx,jy)  = exp(arg)
   enddo
enddo

cli(0,0) = (0.0,0.0)

return
end subroutine init_lintstep


! ************************
! * SUBROUTINE INIT_RAND *
! ************************
subroutine init_rand
use catmod
implicit none

integer :: k, clock

call random_seed(size=nseedlen)
allocate(seed(nseedlen))

if (myseed(1) /= 0) then
   seed(:) = 0
   k = nseedlen
   if (k .gt. mxseedlen) k = mxseedlen
   seed(1:k) = myseed(1:k)
else
   call system_clock(count=clock)
   seed(:) = clock + 37 * (/(k,k=1,nseedlen)/)
endif

call random_seed(put=seed)

return
end subroutine init_rand


! ************************
! * SUBROUTINE INIT_FORC *
! ************************
subroutine init_forc
use catmod
implicit none

real(8) :: k2tmp,ent
integer :: jx,jy

if (nforc .eq. 0) return

!---
ent=0.d0
do jy = 0, nfy
  do jx = 0, nkx
!   k2 = kx2(i) + ky2(j)
    k2tmp = k2(jx,jy)
    if (jy.ne.0.or.jx.ne.0) then
      if (sqrt(k2tmp).ge.kfmin.and.sqrt(k2tmp).le.kfmax) then
        nk = nk + 1
        in(nk,1) = jx
        in(nk,2) = jy
!       ent=ent+k2+alpha
        ent=ent+k2pa(jx,jy)
      endif
    endif
  enddo
enddo

ampcoeff=sqrt(aforc*dt/ent)*nxy

if (tforc.le.dt) then
  itau=1
  else
    itau=tforc/dt
endif

print *, "the forcing ring contains ",nk," wavevectors."
return


end subroutine init_forc


! *************************
! * SUBROUTINE INIT_FORC1 *
! *************************
subroutine init_forc1
use catmod
implicit none

if (nforc .eq. 0) return



!--- introduce default forcing in init_forc1
!--- can be overwritten by input file

!--- introduce markovian time-scale

!--- forcing is specified by input files

!--- examples are produced in simmod.f90

!--- produce masks

!--- normalize forcing


return
end subroutine init_forc1


! ***********************
! * SUBROUTINE ADD_PERT *
! ***********************
subroutine add_pert
use catmod
implicit none


integer :: jx,jy                ! loop index
real(8) :: anglex,delx          ! angle and grid distance in x-direction
real(8) :: angley,dely          ! angle and grid distance in y-direction
real(8) :: gqpert(1:ngx,1:ngy)  ! vorticity perturbation [1/s]
real(8) :: xtmp(1:ngx)          ! temporary vector (x-direction)
real(8) :: ytmp(1:ngy)          ! temporary vector (y-direction)

complex(8) :: cqpert(0:nkx,0:nfy)  ! spectral vorticiy perturbation [1/s]

if ( .not. (tstep .eq. 0) .and. npert .eq. 0) return

delx = 4.0*asin(1.0)/ngx
dely = 4.0*asin(1.0)/ngy

select case(npert)
   !--- zero-mean grid point white noise perturbations
   case(1)
      !--- without symmetry
      call random_number(gqpert)
      gqpert = 2.0*apert*(gqpert - sum(gqpert)*rnxy)
      call grid_to_fourier(gqpert,cqpert,nfx,nfy,ngx,ngy)
      cqpert(0,0) = (0.0d0,0.0d0)
      cq = cq + cqpert
   case(2)
      !--- symmetric about channel center
      do jy = 1,ngy/2
         call random_number(xtmp)
         xtmp = 2.0*apert*(xtmp - sum(xtmp)*rnx)
         gqpert(:,jy)       = xtmp(:)
         gqpert(:,ngy+1-jy) = xtmp(:)
      enddo
      call grid_to_fourier(gqpert,cqpert,nfx,nfy,ngx,ngy)
      cqpert(0,0) = (0.0,0.0)
      cq = cq + cqpert
   case(3)
      !--- anti-symmetric about channel center
      do jy = 1,ngy/2
         call random_number(xtmp)
         xtmp = 2.0*apert*(xtmp - sum(xtmp)*rnx)
         gqpert(:,jy)       =  xtmp(:)
         gqpert(:,ngy+1-jy) = -xtmp(:)
      enddo
      call grid_to_fourier(gqpert,cqpert,nfx,nfy,ngx,ngy)
      cqpert(0,0) = (0.0,0.0)
      cq = cq + cqpert
   case(4)
      !--- symmetric about x = pi
      do jx = 1,ngx/2
         call random_number(ytmp)
         ytmp = 2.0*apert*(ytmp - sum(ytmp)*rny)
         gqpert(jx,:)       = ytmp(:)
         gqpert(ngx+1-jx,:) = ytmp(:)
      enddo
      call grid_to_fourier(gqpert,cqpert,nfx,nfy,ngx,ngy)
      cqpert(0,0) = (0.0,0.0)
      cq = cq + cqpert
   case(5)
      !--- anti-symmetric about x = pi
!     do jx = 1,ngx/2
!        call random_number(ytmp)
!        ytmp = 2.0*apert*(ytmp - sum(ytmp)*rny)
!        gqpert(jx,:)       =  ytmp(:)
!        gqpert(ngx+1-jx,:) = -ytmp(:)
!     enddo
      do jy = 1,ngy/2
         call random_number(xtmp)
         xtmp = 2.0*apert*(xtmp - sum(xtmp)*rnx)
         gqpert(:,jy)       =  xtmp(:)
         gqpert(:,ngy+1-jy) = -xtmp(:)
      enddo
      gqpert = transpose(gqpert)
      call grid_to_fourier(gqpert,cqpert,nfx,nfy,ngx,ngy)
      cqpert(0,0) = (0.0,0.0)
      cq = cq + cqpert
   case(6)
      cqpert = (0.0d0,0.0d0)
      cqpert(1,nky/2) = 0.5d0*apert*cmplx(1.0d0,0.0d0)
      cqpert(nkx/2+1,1) = 0.5d0*apert*cmplx(1.0d0,0.0d0)
      cq = cq + cqpert
end select

return
end subroutine add_pert


! *********************
! * SUBROUTINE CQ2GUV *
! *********************
subroutine cq2guv
use catmod
implicit none

call cq2fuv
call fourier_to_grid(fu,gu,nfx,nfy,ngx,ngy)
call fourier_to_grid(fv,gv,nfx,nfy,ngx,ngy)

return
end subroutine cq2guv


! **********************
! * SUBROUTINE CQ2GQUV *
! **********************
subroutine cq2gquv
use catmod
implicit none

call cq2fuv
call fourier_to_grid(cq,gq,nfx,nfy,ngx,ngy)
call fourier_to_grid(fu,gu,nfx,nfy,ngx,ngy)
call fourier_to_grid(fv,gv,nfx,nfy,ngx,ngy)

return
end subroutine cq2gquv

! **********************
! * SUBROUTINE CQ2GQUV *
! **********************
subroutine cq2gqxqy
use catmod
implicit none

call cq2fqxqy
call fourier_to_grid(fqx,gqx,nfx,nfy,ngx,ngy)
call fourier_to_grid(fqy,gqy,nfx,nfy,ngx,ngy)

return
end subroutine cq2gqxqy

! *************************
! * SUBROUTINE CAT_WRTOUT *
! *************************
subroutine cat_wrtout
use catmod
implicit none

real(8)  :: ftmp(0:nfx,0:nfy) ! temporary field for i/o


integer :: kk,kcode

if((ngp .gt. 0) .and. mod(tstep,ngp).eq.0)  then
   do kk = 1, noutgp
      kcode = outgp(kk)
      select case (kcode)
      case (qcde)
         call cat_wrtgp(nugp,gq,qcde,0)
      end select
   enddo
endif

if((nsp .gt. 0) .and. mod(tstep,nsp).eq.0)  then
   do kk = 1, noutsp
      kcode = outsp(kk)
      select case (kcode)
      case (qcde)
         call c2f(cq,ftmp)
         call cat_wrtsp(nusp,ftmp,qcde,0)
      end select
   enddo
endif

if((ntseri .gt. 0) .and. mod(tstep,ntseri).eq.0) call write_tseri

if((ncfl .gt. 0) .and. mod(tstep,ncfl).eq.0)   call write_cfl

return
end subroutine cat_wrtout


! ***************************
! * SUBROUTINE GUI_TRANSFER *
! ***************************
subroutine gui_transfer
use catmod
implicit none

ggui(:,:) = gq(:,:) ! double precision -> single
call guiput("GQ" // char(0), ggui, ngx, ngy, 1)

! double -> single and center about ky = 0
c4(:,0:nky-1)     = cq(:,nky+1:nfy)  
c4(:,nky:nfy) = cq(:,0:nky)
call guiput("C4" // char(0), c4, nkx+1, nfy+1, 2)   ! fc

if (npost > 0) then
   !--- zonal means
   gguixm(:) = gqxm(:)   ! double -> single
   call guiput("GQXM" // char(0), gguixm, ngy, 1, 1)   ! q

   gguixm(:) = gpsixm(:) ! double -> single
   call guiput("GPSIXM" // char(0), gguixm, ngy, 1, 1) ! psi

!  gguixm(:) = guxm(:) ! double -> single
!  call guiput("GUXM" // char(0), gguixm, ngy, 1, 1)   ! u

!  gguixm(:) = gvxm(:) ! double -> single
!  call guiput("GVXM" // char(0), gguixm, ngy, 1, 1)   ! v

   ggxm(:,1) = gvym(:)
!  ggxm(:,2) = gvxm(:)
   ggxm(:,2) = gqym(:)
   call guiput("GXM" // char(0), ggxm, ngy, 2, 1)   ! u & v


   !--- meridional means
   gguiym(:) = gqym(:)   ! double -> single
   call guiput("GQYM" // char(0), gguiym, ngx, 1, 1) ! q

   gguiym(:) = gpsiym(:) ! double -> single
   call guiput("GPSIYM" // char(0), gguiym, ngx, 1, 1) ! psi

   gguiym(:) = guym(:) ! double -> single
   call guiput("GUYM" // char(0), gguiym, ngx, 1, 1)   ! u

   gguiym(:) = gvym(:) ! double -> single
   call guiput("GVYM" // char(0), gguiym, 1, ngx, 1)   ! v

endif

return
end subroutine gui_transfer


! ************************
! * SUBROUTINE CAT_WRTGP *
! ************************
subroutine cat_wrtgp(ku,gpvar,kcode,klev)
use catmod
implicit none

integer :: ku,kcode,klev
integer :: yy,mo,dd,hh,mm,ii
real(8) :: gpvar(ngx,ngy)

! Build a header for service format

mm = mod(tstep,60)
ii = tstep / 60
hh = mod(ii,24)
ii = ii / 24
dd = mod(ii,30) + 1
ii = ii / 30
mo = mod(ii,12) + 1
yy = ii / 12

ihead(1) = kcode
ihead(2) = klev
ihead(3) = dd + 100 * mo + 10000 * yy
ihead(4) = mm + 100 * hh
ihead(5) = ngx
ihead(6) = ngy
ihead(7) = 0
ihead(8) = 0

write (ku) ihead
write (ku) gpvar(:,:)

return
end subroutine cat_wrtgp


! ************************
! * SUBROUTINE CAT_WRTSP *
! ************************
subroutine cat_wrtsp(ku,spvar,kcode,klev)
use catmod
implicit none

integer :: ku,kcode,klev
integer :: yy,mo,dd,hh,mm,ii
real(8) :: spvar(0:nfx,0:nfy)

mm = mod(tstep,60)
ii = tstep / 60
hh = mod(ii,24)
ii = ii / 24
dd = mod(ii,30) + 1
ii = ii / 30
mo = mod(ii,12) + 1
yy = ii / 12

ihead(1) = kcode
ihead(2) = klev
ihead(3) = dd + 100 * mo + 10000 * yy
ihead(4) = mm + 100 * hh
ihead(5) = nfx+1
ihead(6) = nfy+1
ihead(7) = 0
ihead(8) = 0

write (ku) ihead
write (ku) spvar(:,:)

return
end subroutine cat_wrtsp


! *************************
! * SUBROUTINE INIT_TSTEP *
! *************************
subroutine init_tstep
use catmod
implicit none

real(8) :: dt2


dt2 = dt/2

!--- determine tstop and return if tstep > 0
if (tstep .eq. 0) then
   tstop = tstep + nsteps
else
   tstop = tstep - 1 + nsteps
   return
endif

select case (tstp_mthd)
case (1)
   call cq2gquv
   call jacobian
   cjac1(:,:) = cjac0(:,:)
   cjac2(:,:) = cjac1(:,:)

   call cat_wrtout
   !--- euler-step with dt/2
   cq(:,:) = cli(:,:)*(cq(:,:)+dt2*cjac0(:,:))

   call cq2gquv
   call jacobian
   !--- euler-step with dt
   cq(:,:) = cli(:,:)*(cq(:,:)+dt*cjac0(:,:))
   tstep=tstep+1

   call cq2gquv
   call jacobian
   call cat_wrtout

   !--- adams-bashford method 2nd order
   cq(:,:) = cli(:,:) * (cq(:,:) + dt2 * (3.0 * cjac0(:,:) -            &
             cli(:,:) * cjac1(:,:)))
   cq(0,0) = (0.0,0.0)

   if (nforc .ge. 1) call add_forc

   cjac1(:,:) = cjac0(:,:)

   tstep=tstep+1
end select

return
end subroutine init_tstep

! ************************
! * SUBROUTINE INIT_POST *
! ************************
subroutine init_post
use catmod
implicit none

allocate(gqxm(1:ngy))   ; gqxm(:)   = 0.0 ! q-mean in x-dir [1/s]
allocate(gqym(1:ngx))   ; gqym(:)   = 0.0 ! q-mean in y-dir [1/s]

allocate(gpsixm(1:ngy)) ; gpsixm(:) = 0.0 ! psi-mean in x-dir [m^2/s]
allocate(gpsiym(1:ngx)) ; gpsiym(:) = 0.0 ! psi-mean in y-dir [m^2/s]

allocate(guxm(1:ngy))   ; guxm(:)   = 0.0 ! u-mean in x-dir [m/s]
allocate(guym(1:ngx))   ; guym(:)   = 0.0 ! u-mean in y-dir [m/s]

allocate(gvxm(1:ngy))   ; gvxm(:)   = 0.0 ! v-mean in x-dir [m/s]
allocate(gvym(1:ngx))   ; gvym(:)   = 0.0 ! v-mean in y-dir [m/s]

allocate(gguixm(1:ngy)) ; gguixm(:) = 0.0 ! single precision mean in x-dir
allocate(gguiym(1:ngx)) ; gguiym(:) = 0.0 ! single precision mean in y-dir

allocate(ggxm(1:ngx,4)) ; ggxm(:,:) = 0.0 ! single precision mean in x-dir

return
end subroutine init_post


! ************************
! * SUBROUTINE STOP_POST *
! ************************
subroutine stop_post
use catmod
implicit none

if (npost > 0) then
   deallocate(gqxm)
   deallocate(gqym)

   deallocate(gpsixm)
   deallocate(gpsiym)

   deallocate(guxm)
   deallocate(guym)

   deallocate(gvxm)
   deallocate(gvym)

   deallocate(gguixm)
   deallocate(gguiym)
endif

return
end subroutine stop_post


! *************************
! * SUBROUTINE CAT_MASTER *
! *************************
subroutine cat_master
use catmod
implicit none

do while (tstep <= tstop)
   call cq2gquv
   call jacobian
   call post
   call cat_wrtout
   if (ngui > 0 .and. mod(tstep,ngui) == 0) then
      call gui_transfer
      call guistep_cat
   endif
   call step
   if (nforc .ge. 1) call add_forc
   if (nsim > 0) call simstep
   if (nuser > 0) call userstep
   tstep = tstep + 1
   if (nstdout.gt.0 .and. ngui == 0 .and. mod(tstep,nstdout) == 0) then
      write(*,*)' time step ',tstep
   endif
   if (nshutdown > 0) return
enddo

return
end subroutine cat_master


! *******************
! * SUBROUTINE POST *
! *******************
subroutine post
use catmod
implicit none

!--- calculate zonal means
gqxm   = sum(gq  ,1)*rnx
gpsixm = sum(gpsi,1)*rnx
guxm   = sum(gu  ,1)*rnx
gvxm   = sum(gv  ,1)*rnx

!--- calculate meridional means
gqym   = sum(gq  ,2)*rny
gpsiym = sum(gpsi,2)*rny
guym   = sum(gu  ,2)*rny
gvym   = sum(gv  ,2)*rny

return
end subroutine post


! *******************
! * SUBROUTINE STEP *
! *******************
subroutine step
use catmod
implicit none

real(8) :: c,c1


select case (tstp_mthd)
case (1)
   c  = dt/12.0
   c1 = 23.0 * c

   !--- adams-bashford 3rd order
   cq(:,:) = cli(:,:) * (cq(:,:) + c1 * cjac0(:,:) + c * cli(:,:) *  &
             (-16.0 * cjac1(:,:) + 5.0 * cli(:,:)*cjac2(:,:)))
   cq(0,0) = (0.0,0.0)

   !--- shift time-levels (pointers are faster)
   cjac2(:,:) = cjac1(:,:)
   cjac1(:,:) = cjac0(:,:)

end select

return
end subroutine step


! *************************
! * SUBROUTINE CAT_EPILOG *
! *************************
subroutine cat_epilog
use catmod
implicit none

if (ngui > 0) call guistop

call write_rst

call cpu_time(tmstop)
tmrun = tmstop - tmstart
tsps    = nint(nsteps / tmrun)

write(nudiag, &
 '(" *************************************************")')
write(nudiag,'("  Total time in seconds: ",f15.2)') tmrun
write(nudiag,'("  Time steps per second: ",i12)') tsps
write(nudiag, &
 '(" *************************************************")')

call stop_jacobian
call stop_post
if (nsim > 0) call simstop
if (nuser > 0) call userstop
call close_files

return
end subroutine cat_epilog


! ************************
! * SUBROUTINE WRITE_RST *
! ************************
subroutine write_rst
use catmod
implicit none


!--- namelist parameters
if (.not. lrsttest) call put_restart_iarray('nsteps',nsteps,1,1)
call put_restart_iarray('ngp',ngp,1,1)
call put_restart_iarray('nsp',nsp,1,1)
call put_restart_iarray('ingp',ingp,ningp,1)
call put_restart_iarray('insp',insp,ninsp,1)
call put_restart_iarray('outgp',outgp,noutgp,1)
call put_restart_iarray('outsp',outsp,noutsp,1)
call put_restart_iarray('tstp_mthd',tstp_mthd,1,1)
call put_restart_iarray('ncfl',ncfl,1,1)
call put_restart_array ('dt',dt,1,1)
call put_restart_array ('alpha',alpha,1,1)
call put_restart_array ('beta',beta,1,1)
call put_restart_array ('lx',lx,1,1)
call put_restart_array ('ly',ly,1,1)
call put_restart_array ('sig',sig,nsig,1)
call put_restart_array ('psig',psig,nsig,1)
call put_restart_array ('rtsig',rtsig,nsig,1)
call put_restart_iarray('ksig',ksig,nsig,1)
call put_restart_array ('lam',lam,nlam,1)
call put_restart_array ('plam',plam,nlam,1)
call put_restart_array ('rtlam',rtlam,nlam,1)
call put_restart_iarray('klam',klam,nlam,1)
call put_restart_iarray('diss_mthd',diss_mthd,1,1)
call put_restart_iarray('npert',npert,1,1)
call put_restart_array ('apert',apert,1,1)
call put_restart_iarray('nforc',nforc,1,1)
call put_restart_iarray('kfmin',kfmin,1,1)
call put_restart_iarray('kfmax',kfmax,1,1)
call put_restart_array ('aforc',aforc,1,1)
call put_restart_array ('tforc',tforc,1,1)
call put_restart_iarray('myseed',myseed,mxseedlen,1)
call put_restart_iarray('ntseri',ntseri,1,1)
call put_restart_iarray('nstdout',nstdout,1,1)
call put_restart_iarray('jacmthd',jacmthd,1,1)
call put_restart_iarray('ndiag',ndiag,1,1)
call put_restart_array ('jac_scale',jac_scale,1,1)


!--- additional parameters
call put_restart_iarray('ngx',ngx,1,1)
call put_restart_iarray('tstep',tstep,1,1)

call random_seed(get=seed)
call put_restart_iarray('seed',seed,nseedlen,1)

!--- fluid state
call put_restart_carray('cq',cq,nkx+1,nfy+1)
call put_restart_carray('cjac1',cjac1,nkx+1,nfy+1)
call put_restart_carray('cjac2',cjac2,nkx+1,nfy+1)

return
end subroutine write_rst


! **************************
! * SUBROUTINE CLOSE_FILES *
! **************************
subroutine close_files
use catmod

if (lrst)            close(nurstini)
if (lcatnl)          close(nucatnl)
close(nurstfin)
close(nudiag)
if (ntseri .gt. 0)   close(nutseri)
if (ncfl .gt. 0)     close(nucfl)
if (ngp .gt. 0)      close(nugp)
if (nsp .gt. 0)      close(nusp)

return
end subroutine close_files


! ***********************
! * SUBROUTINE JACOBIAN *
! ***********************
subroutine jacobian
use catmod
implicit none

integer    :: jx

select case (jacmthd)
case (0)
   !--- Jacobian switched off
   cjac0(:,:) = cmplx(0.0d0,0.0d0)

case (1)
   !--- Jacobian in divergence form (second representation in CAT-UG)
   guq = gu*gq
   gvq = gv*gq

   call grid_to_fourier(guq,fuq,nfx,nfy,ngx,ngy)
   call grid_to_fourier(gvq,fvq,nfx,nfy,ngx,ngy)

   !--- (!!! Jacobian has different sign than in CAT-UG !!!)
   do jx = 0, nkx
      cjac0(jx,:) = &
       cmplx( kx(jx)*fuq(jx+jx+1,:)+ky(:)*fvq(jx+jx+1,:), &
             -kx(jx)*fuq(jx+jx  ,:)-ky(:)*fvq(jx+jx  ,:))
   enddo

case (2)
   !--- Jacobian in advection form (first representation in CAT-UG)
   call cq2gqxqy
   gjac = gu*gqx + gv*gqy
   call grid_to_fourier(gjac,cjac0,nfx,nfy,ngx,ngy)

case (3)
   !--- Jacobian in invariant form (third representation in CAT-UG)
   gv2mu2 = gv**2 - gu**2
   guv    = gu*gv

   call grid_to_fourier(gv2mu2,fv2mu2,nfx,nfy,ngx,ngy)
   call grid_to_fourier(guv,fuv,nfx,nfy,ngx,ngy)

   do jx = 0, nkx
      cjac0(jx,:) = &
       cmplx(kxky(jx,:)*fv2mu2(2*jx  ,:) + kx2mky2(jx,:)*fuv(2*jx  ,:), &
             kxky(jx,:)*fv2mu2(2*jx+1,:) + kx2mky2(jx,:)*fuv(2*jx+1,:))
   enddo

end select

if (jac_scale .ne. 1.0 .or. jac_scale .ne. 0.0)  &
                           cjac0(:,:) = jac_scale*cjac0(:,:)

return
end subroutine jacobian

! *********************
! * SUBROUTINE CQ2FUV *
! *********************
subroutine cq2fuv
use catmod
implicit none

integer :: jx

do jx = 0, nkx
   fu(2*jx  ,:) = -aimag(cq(jx,:)) * kyrk2pa(jx,:)
   fu(2*jx+1,:) =  real (cq(jx,:)) * kyrk2pa(jx,:)
   fv(2*jx  ,:) =  aimag(cq(jx,:)) * kxrk2pa(jx,:)
   fv(2*jx+1,:) = -real (cq(jx,:)) * kxrk2pa(jx,:)
enddo

return
end subroutine cq2fuv


! ***********************
! * SUBROUTINE CQ2FQXQZ *
! ***********************
subroutine cq2fqxqy
use catmod
implicit none

integer :: jx, jy

do jx = 0, nkx
   do jy = 0, nfy
      fqx(2*jx  ,jy) =  aimag(cq(jx,jy)) * kx(jx)
      fqx(2*jx+1,jy) = -real (cq(jx,jy)) * kx(jx)

      fqy(2*jx  ,jy) =  aimag(cq(jx,jy)) * ky(jy)
      fqy(2*jx+1,jy) = -real (cq(jx,jy)) * ky(jy)
   enddo
enddo

return
end subroutine cq2fqxqy


! ***********************
! * SUBROUTINE ADD_FORC *
! ***********************
subroutine add_forc
use catmod
implicit none
integer :: i,j,ic,ifk
real(8) :: k2tmp
complex(8) :: psif
real(8) :: eni,enf,ran4

eni=0.d0
enf=0.d0

select case (nforc)
   case (1)
      cq = cq + cqfrc

   case (2)
      do ifk = 1,nk
         i = in(ifk,1)
         j = in(ifk,2)
!        psif = ampcoeff*exp(ci*twopi*0.1525125)
         psif = ampcoeff*exp(ci*twopi)
         k2tmp = k2pa(i,j)
         eni = eni+(cq(i,j)*cq(i,j))/k2tmp
         cq(i,j) = cq(i,j)-k2tmp*psif*rnxy
         enf = enf+(cq(i,j)*cq(i,j))/k2tmp
      enddo

   case (3)
      ic=tstep+1
      if (mod(ic,itau).eq.0.) then
         call random_number(ran4)
         phi=ran4*twopi*ci
      endif
      psif=ampcoeff*exp(phi)
      ! add the forcing to the spectral vorticity q -> q+ k2tmp *psif
      do ifk=1,nk
         i = in(ifk,1)
         j = in(ifk,2)
         k2tmp = k2pa(i,j)
         eni=eni+(cq(i,j)*cq(i,j))/k2tmp
         cq(i,j) = cq(i,j)-k2tmp*psif*rnxy
         enf=enf+(cq(i,j)*cq(i,j))/k2tmp
      enddo

   case (4)
      do ifk=1,nk
         i = in(ifk,1)
         j = in(ifk,2)
         call random_number(ran4)
         psif = ampcoeff*exp(ci*twopi*ran4)
         k2tmp = k2pa(i,j)
         eni=eni+(cq(i,j)*cq(i,j))/k2tmp
         cq(i,j) = cq(i,j)-k2tmp*psif*rnxy
         enf=enf+(cq(i,j)*cq(i,j))/k2tmp
      enddo

end select

return
end subroutine add_forc


! ************************
! * SUBROUTINE ADD_FORC1 *
! ************************
subroutine add_forc1
use catmod
implicit none


select case (nforc)
   case(1)

   case(2)

   case(3)

   case default
      return
end select

return
end subroutine add_forc1



! **************************
! * SUBROUTINE WRITE_TSERI *
! **************************
subroutine write_tseri
use catmod
implicit none

real(8) :: ener,enst,qint


qint = sum(gq(:,:))
ener = 0.5 * rnxy * (sum(gu(:,:) * gu(:,:)) + sum(gv(:,:) * gv(:,:)))
enst = 0.5 * rnxy *  sum(gq(:,:) * gq(:,:))

write(nutseri,*) tstep,qint,ener,enst

return
end subroutine write_tseri


! ************************
! * SUBROUTINE WRITE_CFL *
! ************************
subroutine write_cfl
use catmod
implicit none

real(8) :: cfl
real(8) :: maxu,maxv


!--- check courant number
maxu = maxval(abs(gu(:,:))) * dt * ngx / lx
maxv = maxval(abs(gv(:,:))) * dt * ngy / ly

cfl = max(maxu,maxv)

write(nucfl,*) tstep,cfl

return
end subroutine write_cfl


! ***********************
! * SUBROUTINE INIT_OPS *
! ***********************
subroutine init_ops
use catmod
implicit none

integer :: jx,jy

kx2    = kx*kx
kx2(0) = 1.0

ky2    = ky*ky

do jy = 0, nfy
   do jx = 0, nkx
      k2     (jx,jy) = (kx2(jx)+ky2(jy))        ! - Laplacian      psi --> -q
      k2pa   (jx,jy) = (kx2(jx)+ky2(jy))+alpha  ! - mod Laplacian  psi --> -q
      rk2pa  (jx,jy) = 1.0d0/(k2(jx,jy)+alpha)  ! 1/mod Laplacian    q --> -psi
      kxrk2pa(jx,jy) = kx(jx)*rk2pa(jx,jy)      ! q ---> v/i
      kyrk2pa(jx,jy) = ky(jy)*rk2pa(jx,jy)      ! q ---> u/i
      kx2mky2(jx,jy) = kx2(jx)-ky2(jy)          ! utility Jacobian 3
      kxky   (jx,jy) = kx(jx)*ky(jy)            ! utility Jacobian 3
   enddo
enddo

return
end subroutine init_ops


! ***********************
! * SUBROUTINE INIT_RST *
! ***********************
subroutine init_rst
use catmod
implicit none

inquire(file=cat_rstini,exist=lrst)
if (lrst) then
  open(nurstini,file=cat_rstini,form='unformatted')
endif

open(nurstfin,file=cat_rstfin,form='unformatted')

return
end subroutine init_rst



! *************************
! * SUBROUTINE CHECK_DIFF *
! *************************
subroutine check_diff(cold,cnew,ytext)
use catmod
complex(8) :: cold(0:nkx,0:nfy)
complex(8) :: cnew(0:nkx,0:nfy)
character(*) :: ytext

write (nudiag,*) ytext
write (nudiag,*) "max old  = ",abs(maxval(real(cold)))+abs(maxval(aimag(cold)))
write (nudiag,*) "max new  = ",abs(maxval(real(cnew)))+abs(maxval(aimag(cnew)))
write (nudiag,*) "max diff = ",abs(maxval(real (cnew-cold(0:nkx,0:nfy)))) + &
                               abs(maxval(aimag(cnew-cold(0:nkx,0:nfy))))
write (nudiag,*) "min old  = ",abs(minval(real(cold)))+abs(minval(aimag(cold)))
write (nudiag,*) "min new  = ",abs(minval(real(cnew)))+abs(minval(aimag(cnew)))
write (nudiag,*) "min diff = ",abs(minval(real (cnew-cold(0:nkx,0:nfy)))) + &
                               abs(minval(aimag(cnew-cold(0:nkx,0:nfy))))


return
end
