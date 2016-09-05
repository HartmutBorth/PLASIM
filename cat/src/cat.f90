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


!--- array with grid types
integer, parameter      :: ngtp         = 2
character(2), parameter :: gtp_ar(ngtp) = (/ "GP" , "SP" /)

!--- codes of variables to be read at initialization
integer, parameter :: ningp  = 5
integer, parameter :: ninsp  = 5
integer :: ingp(ningp) = &
   (/ 138, &
       0 , &     
       0 , &     
       0 , &     
       0   &    
   /) 
integer :: insp(ninsp) = &
   (/  0 , &
       0 , &     
       0 , &     
       0 , &     
       0   &    
   /) 

!--- codes of variables to be written to cat_gp or cat_sp 
integer, parameter :: noutgp = 5
integer, parameter :: noutsp = 5
integer :: outgp(noutgp) = &
   (/ 138, &
       0 , &
       0 , &
       0 , &
       0   &
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
       (/ 32,                 & ! powers of Laplacian
          32                  &
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

!--- Forcing
real(8) :: aforc = 0.015     ! amplitude of forcing
real(8) :: tforc = 0.001     ! memory time scale of forcing [s]

integer :: idum
integer :: in(1000,2)
integer :: nk
integer :: itau

real(8)    :: ampcoeff

complex(8) :: phi

integer :: nforc = 0 ! forcing switch
                     ! 0 = no forcing
                     ! 1 = constant forcing with
                     ! spectral forcing ring
                     ! 2 = markov chain
                     ! 3 = forcing with constant wave numbers
                     !     but random uncorrelated phases (white noise)

real(8) :: kfmin  = 4.0+1.0 ! min radius = sqrt(kfmin) of spectr. forcing
real(8) :: kfmax  = 4.0-1.0 ! max radius = sqrt(kfmax) of spectr. forcing

!--- Scaling
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


!--- grid in Fourier space
integer :: nkx  = 21  ! ngx/3 (max. wave number in x-direction)
integer :: nky  = 21  ! ngy/3 (max. wave number in y-direction)
integer :: nfx  = 43  ! 2*nkx+1 (x-dimension in fourier domain)
integer :: nfy  = 42  ! 2*nky   (y-dimension in fourier domain)


!--- Jacobian
integer :: jac_mthd  = 1      ! approximation method of Jacobian
                              ! 0 : no Jacobian
                              ! 1 : divergence form

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

!--- variables in physical/gridpoint (GP) space
real(8), allocatable :: gq(:,:)   ! vorticity            [1/s]
real(8), allocatable :: gpsi(:,:) ! stream function      [m^2/s]
real(8), allocatable :: gu(:,:)   ! velocity x-direction [m/s]
real(8), allocatable :: gv(:,:)   ! velocity y-direction [m/s]

real(8), allocatable :: guq(:,:)  ! u*vorticity
real(8), allocatable :: gvq(:,:)  ! v*vorticity

real(4), allocatable :: ggui(:,:) ! single precision transfer

!--- variables in Fourierspace, real representation (F)
real(8), allocatable :: fpsi(:,:) ! stream function
real(8), allocatable :: fu(:,:)   ! velocity x-direction
real(8), allocatable :: fv(:,:)   ! velocity y-direction

real(8), allocatable :: fuq(:,:)  ! u*vorticity
real(8), allocatable :: fvq(:,:)  ! v*vorticity

real(8), allocatable :: ftmp(:,:) ! temporary field for, e.g. i/o


!--- variables in Fourierspace, complex representation (
complex(8), allocatable :: cq(:,:)     ! vorticity
complex(8), allocatable :: cjac0(:,:)  ! Jacobian at time 0
complex(8), allocatable :: cjac1(:,:)  ! Jacobian at time -1
complex(8), allocatable :: cjac2(:,:)  ! Jaconbian at time -2

!--- operators in Fourier Space
integer   , allocatable :: ki(:),kj(:)  
real(8)   , allocatable :: ki2(:), kj2(:) 
real(8)   , allocatable :: k2n(:,:), rk2an(:,:)
real(8)   , allocatable :: kirk2an(:,:),kjrk2an(:,:)
complex(8), allocatable :: cli(:,:)                  ! linear time propagation


! *******************
! * External moduls *
! *******************

!--- gui communication (guimod)
integer :: ngui    = 1   ! global switch 1 = on
integer :: nguidbg = 0   ! GUI debug mode

integer :: ndatim(6) = 1 ! date/time display
real(4) :: parc(5) = 0.0 ! timeseries display


!--- predefined simulations (simmod)
integer :: nsim  = 0 ! 1/0 predefined simulations on/off.
                     ! Predefined simulations are specified in 
                     ! <sim_namelist> of simmod. 

!--- usermod (usermod)
integer :: nuser = 0 ! 1/0 user mode is switched on/off.
                     ! In user mode it is possible to introduce 
                     ! new code into CAT (see chapter <Modifying 
                     ! CAT> in the User's Guide).


!--- postprocessing (postmod)
integer :: npost = 0 ! 1/0 post processing is switched on/off  


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
call init_ltprop
call init_rand
call init_forc
call init_step

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
call get_restart_iarray('nforc',nforc,1,1)
call get_restart_iarray('kfmin',kfmin,1,1)
call get_restart_iarray('kfmax',kfmax,1,1)
call get_restart_array ('aforc',aforc,1,1)
call get_restart_array ('tforc',tforc,1,1)
call get_restart_iarray('myseed',myseed,mxseedlen,1)
call get_restart_iarray('ntseri',ntseri,1,1)
call get_restart_iarray('nstdout',nstdout,1,1)
call get_restart_iarray('jac_mthd',jac_mthd,1,1)
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
                  myseed    ,ntseri    ,nstdout ,jac_mthd ,ndiag     , &
                  jac_scale ,nsp       ,outgp   ,outsp    ,tstp_mthd , &
                  nsim      ,nuser     ,npost

inquire(file=cat_namelist,exist=lcatnl)

if (lcatnl) then
  open(nucatnl,file=cat_namelist,iostat=ios_catnl)
  read(nucatnl,cat_nl) 
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

! *******************************
! * 1D real and complex vectors *
! *******************************
allocate(ki(0:nkx))   ; ki(:)   = [(j,j=0,nkx)]
allocate(kj(0:nfy))   ; kj(:)   = [(j,j=0,nky),(j,j=-nky,-1)]
allocate(ki2(0:nkx))  ; ki2(:)  = 0.0
allocate(kj2(0:nfy))  ; kj2(:)  = 0.0



! ******************************
! * 2D real and complex arrays *
! ******************************

!--- grid point space
allocate(gq(1:ngx,1:ngy))   ; gq(:,:)   = 0.0  ! vorticity
allocate(gu(1:ngx,1:ngy))   ; gu(:,:)   = 0.0  ! velocity in x-dir.
allocate(gv(1:ngx,1:ngy))   ; gv(:,:)   = 0.0  ! velocity in y-dir.

allocate(gpsi(1:ngx,1:ngy)) ; gpsi(:,:) = 0.0  ! stream function
allocate(guq(1:ngx,1:ngy))  ; guq(:,:)  = 0.0  ! u*vorticity
allocate(gvq(1:ngx,1:ngy))  ; gvq(:,:)  = 0.0  ! v*vorticity

allocate(ggui(1:ngx,1:ngy)) ; ggui(:,:)  = 0.0  ! GUI transfer

!--- spetral space 
allocate(fu(0:nfx,0:nfy))   ; fu(:,:)   = 0.0 ! u 
allocate(fv(0:nfx,0:nfy))   ; fv(:,:)   = 0.0 ! v 
allocate(fuq(0:nfx,0:nfy))  ; fuq(:,:)  = 0.0 ! u*q
allocate(fvq(0:nfx,0:nfy))  ; fvq(:,:)  = 0.0 ! v*q
allocate(ftmp(0:nfx,0:nfy)) ; ftmp(:,:) = 0.0 ! temporary field for, e.g. i/o

allocate(k2n(0:nkx,0:nfy))     ;k2n(:,:)     = 0.0 ! Laplacian
allocate(rk2an(0:nkx,0:nfy))   ;rk2an(:,:)   = 0.0 ! modified Laplacian^-1
allocate(kirk2an(0:nkx,0:nfy)) ;kirk2an(:,:) = 0.0 ! q --> v
allocate(kjrk2an(0:nkx,0:nfy)) ;kjrk2an(:,:) = 0.0 ! q --> u

allocate(cli(0:nkx,0:nfy))  ; cli(:,:)   = (0.0,0.0) ! linear time propagator 
allocate(cq(0:nkx,0:nfy))   ; cq(:,:)    = (0.0,0.0) ! vorticity
allocate(cjac0(0:nkx,0:nfy)); cjac0(:,:) = (0.0,0.0) ! Jacobian at time level  0
allocate(cjac1(0:nkx,0:nfy)); cjac1(:,:) = (0.0,0.0) ! Jacobian at time level -1
allocate(cjac2(0:nkx,0:nfy)); cjac2(:,:) = (0.0,0.0) ! Jacobian at time level -2

return
end subroutine cat_alloc


! ************************
! * SUBROUTINE CAT_INPUT *
! ************************
subroutine cat_input
use catmod
implicit none

logical         :: lexist
integer         :: kcode,jj,kk,kkmax
character(2)    :: gtp
character(256)  :: fname

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
            case (138)
               select case (gtp)
               case ("GP")
                  call read_gp(kcode,gtp,gq)
                  call grid_to_fourier(gq,cq,nfx,nfy,ngx,ngy)
               case ("SP")
                  call read_sp(kcode,gtp,ftmp)
                  call f2c(ftmp,cq)
               end select
               cq(0,0) = (0.0,0.0)
            end select
         else
            call cat_fname(kcode,gtp,fname)
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
subroutine cat_fname(kcode,gtp,fname)
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

call cat_fname(kcode,gtp,fname)
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

call cat_fname(kcode,gtp,fname)

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

call cat_fname(kcode,gtp,fname)

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
! *******************
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


! **************************
! * SUBROUTINE INIT_LTPROP *
! **************************
subroutine init_ltprop
use catmod
implicit none

real(8)    :: k2,diss_sig,diss_lam,b_term
complex(8) :: arg
integer    :: i,j,n

!--------------------------------------------------------------
! Time propagator of linear part of differential equation
!
! cli = exp(arg)                      with
!
! arg = -dt*[diss_sig*k^2 + diss_lam*k^2 + beta_term]/(k^2+alpha)
!
!  with 
!     diss_sig  = sum_j=1,nsig sig(j)*k^(2*psig(j))
!     diss_lam  = sum_j=1,nlam lam(j)*k^(2*plam(j))
!     b_term    = i*beta*kx
!--------------------------------------------------------------
do j = 0, nfy
   do i = 0, nkx
      diss_sig = 0.0
      diss_lam = 0.0
      k2       = ki2(i)+kj2(j)
      do n = 1,nsig
         diss_sig = diss_sig - sig(n)*k2**psig(n)*k2/(k2+alpha)
      enddo
      do n = 1,nlam
         diss_lam = diss_lam - lam(n)*k2**plam(n)*k2/(k2+alpha)
      enddo
      b_term    = cmplx(0.0,beta * ki(i) / (k2+alpha))

      arg       = dt*(diss_sig + diss_lam + b_term)

      cli(i,j)  = exp(arg)
   enddo
enddo

cli(0,0) = (0.0,0.0)

return
end subroutine init_ltprop


! *************
! * INIT_RAND *
! *************
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


! *************
! * INIT_FORC *
! *************
subroutine init_forc
use catmod
implicit none

real(8) :: k2,ent
integer :: i,j

if (nforc .eq. 0) return

!---     
ent=0.d0
do j = 0, nfy
  do i = 0, nkx
    k2 = ki2(i) + kj2(j)
    if (j.ne.0.or.i.ne.0) then
      if (sqrt(k2).gt.kfmin.and.sqrt(k2).lt.kfmax) then
        nk = nk + 1
        in(nk,1) = i
        in(nk,2) = j
        ent=ent+k2+alpha
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

idum=-12895673

print *, "the forcing ring contains ",nk," wavevectors."
return
end subroutine init_forc


! ********
! * Q2GUV *
! *********
subroutine q2guv
use catmod
implicit none

call q2uv
call fourier_to_grid(fu,gu,nfx,nfy,ngx,ngy)
call fourier_to_grid(fv,gv,nfx,nfy,ngx,ngy)

return
end subroutine q2guv


! **********
! * Q2GQUV *
! **********
subroutine q2gquv
use catmod
implicit none

call q2uv
call fourier_to_grid(cq,gq,nfx,nfy,ngx,ngy)
call fourier_to_grid(fu,gu,nfx,nfy,ngx,ngy)
call fourier_to_grid(fv,gv,nfx,nfy,ngx,ngy)

return
end subroutine q2gquv


! ****************
! * WRITE_OUTPUT *
! ****************
subroutine write_output
use catmod
implicit none

integer :: kk,kcode

if((ngp .gt. 0) .and. mod(tstep,ngp).eq.0)  then
   do kk = 1, noutgp
      kcode = outgp(kk)
      select case (kcode)
      case (138)
         call write_gp(nugp,gq,138,0)
      end select
   enddo 
endif

if((nsp .gt. 0) .and. mod(tstep,nsp).eq.0)  then
   do kk = 1, noutsp
      kcode = outsp(kk)
      select case (kcode)
      case (138)
         call c2f(cq,ftmp)
         call write_sp(nusp,ftmp,138,0)
      end select
   enddo
endif

if((ntseri .gt. 0) .and. mod(tstep,ntseri).eq.0) call write_tseri

if((ncfl .gt. 0) .and. mod(tstep,ncfl).eq.0)   call write_cfl

return
end subroutine write_output


! ****************
! * GUI TRANSFER *
! ****************
subroutine gui_transfer
use catmod
implicit none

ggui(:,:) = gq(:,:) ! double precision -> single
call guiput("GQ" // char(0), ggui, ngx, ngy, 1)

return
end subroutine gui_transfer


! ************
! * WRITE_GP *
! ************
subroutine write_gp(ku,gpvar,kcode,klev)
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
end subroutine write_gp


! ************
! * WRITE_SP *
! ************
subroutine write_sp(ku,spvar,kcode,klev)
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
end subroutine write_sp


! ************************
! * SUBROUTINE INIT_STEP *
! ************************
subroutine init_step
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
   call q2gquv
   call jacobian
   cjac1(:,:) = cjac0(:,:)
   cjac2(:,:) = cjac1(:,:)

   call write_output
   !--- euler-step with dt/2
   cq(:,:) = cli(:,:)*(cq(:,:)+dt2*cjac0(:,:))

   call q2gquv
   call jacobian
   !--- euler-step with dt
   cq(:,:) = cli(:,:)*(cq(:,:)+dt*cjac0(:,:))
   tstep=tstep+1

   call q2gquv
   call jacobian
   call write_output

   !--- adams-bashford method 2nd order
   cq(:,:) = cli(:,:) * (cq(:,:) + dt2 * (3.0 * cjac0(:,:) -            &
             cli(:,:) * cjac1(:,:)))

   cq(0,0) = (0.0,0.0)

   if (nforc .ge. 1) call add_forc
   cjac1(:,:) = cjac0(:,:)
   tstep=tstep+1
end select

return
end subroutine init_step


! *************************
! * SUBROUTINE CAT_MASTER *
! *************************
subroutine cat_master
use catmod
implicit none

if (nshutdown > 0) return   ! if an error occured so far


do while (tstep <= tstop)
   call q2gquv
   call jacobian
   call write_output
   if (ngui > 0 .and. mod(tstep,ngui) == 0) then
      call gui_transfer
      call guistep_cat
   endif
   call step
   if (nforc .ge. 1) call add_forc
   if (nsim > 0) call simstep
   if (nuser > 0) call userstep
   tstep = tstep + 1
   if (nstdout.ge.0 .and. mod(tstep,nstdout) == 0) then
      write(*,*)' time step ',tstep
   endif
enddo

return
end subroutine cat_master


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
call put_restart_iarray('nforc',nforc,1,1)
call put_restart_iarray('kfmin',kfmin,1,1)
call put_restart_iarray('kfmax',kfmax,1,1) 
call put_restart_array ('aforc',aforc,1,1)
call put_restart_array ('tforc',tforc,1,1)
call put_restart_iarray('myseed',myseed,mxseedlen,1)
call put_restart_iarray('ntseri',ntseri,1,1)
call put_restart_iarray('nstdout',nstdout,1,1)
call put_restart_iarray('jac_mthd',jac_mthd,1,1)
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
if (ngp .gt. 0)      close (nugp)
if (nsp .gt. 0)      close (nusp)

return
end subroutine close_files


! ************
! * JACOBIAN *
! ************
subroutine jacobian
use catmod
implicit none

integer    :: i,j

select case (jac_mthd)

case (1)
   guq = gu*gq
   gvq = gv*gq

   call grid_to_fourier(guq,fuq,nfx,nfy,ngx,ngy)
   call grid_to_fourier(gvq,fvq,nfx,nfy,ngx,ngy)

   do j = 0, nfy
      do i = 0, nkx
         cjac0(i,j) = cmplx(  ki(i) * fuq(i+i+1,j) + kj(j)*fvq(i+i+1,j),  &
                            - ki(i) * fuq(i+i  ,j) - kj(j)*fvq(i+i  ,j) )
      enddo
   enddo
end select

if (jac_scale .ne. 1.0) cjac0(:,:) = jac_scale*cjac0(:,:)

return
end subroutine jacobian


! ********
! * Q2UV *
! ********
subroutine q2uv
use catmod
implicit none

integer :: i,j

do j = 0, nfy
   do i = 0, nkx
      fu(i+i  ,j) = -aimag(cq(i,j)) * kjrk2an(i,j)
      fu(i+i+1,j) =  real (cq(i,j)) * kjrk2an(i,j)
      fv(i+i  ,j) =  aimag(cq(i,j)) * kirk2an(i,j)
      fv(i+i+1,j) = -real (cq(i,j)) * kirk2an(i,j)
   enddo
enddo

return
end subroutine q2uv


! ***************
! * ADD_FORCING *
! ***************
subroutine add_forc
use catmod
implicit none
integer :: i,j,ic,ifk
real(8) :: k2
complex(8) :: psif
real(8) :: eni,enf,ran4

eni=0.d0
enf=0.d0

select case (nforc)
   case (1)
      do ifk = 1,nk
         i = in(ifk,1)
         j = in(ifk,2)
!        psif = ampcoeff*exp(ci*twopi*0.1525125)
         psif = ampcoeff*exp(ci*twopi)
         k2 = ki2(i)+kj2(j)+alpha
         eni = eni+(cq(i,j)*cq(i,j))/k2
         cq(i,j) = cq(i,j)-k2*psif*rnxy
         enf = enf+(cq(i,j)*cq(i,j))/k2
      enddo

   case (2)
      ic=tstep+1
      if (mod(ic,itau).eq.0.) then
         call random_number(ran4)
         phi=ran4*twopi*ci
      endif
      psif=ampcoeff*exp(phi)
      ! add the forcing to the spectral vorticity q -> q+ k2 *psif
      do ifk=1,nk
         i = in(ifk,1)
         j = in(ifk,2)
         k2 = ki2(i)+kj2(j)+alpha
         eni=eni+(cq(i,j)*cq(i,j))/k2
         cq(i,j) = cq(i,j)-k2*psif*rnxy
         enf=enf+(cq(i,j)*cq(i,j))/k2
      enddo

   case (3)
      do ifk=1,nk
         i = in(ifk,1)
         j = in(ifk,2)
         call random_number(ran4)
         psif = ampcoeff*exp(ci*twopi*ran4)
         k2 = ki2(i)+kj2(j)+alpha
         eni=eni+(cq(i,j)*cq(i,j))/k2
         cq(i,j) = cq(i,j)-k2*psif*rnxy
         enf=enf+(cq(i,j)*cq(i,j))/k2
      enddo

end select

return
end subroutine add_forc


! ***************
! * WRITE_TSERI *
! ***************
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


! *************
! * WRITE_CFL *
! *************
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


! ************
! * INIT_OPS *
! ************
subroutine init_ops
use catmod
implicit none

integer :: i,j    
 
ki2    = ki*ki
ki2(0) = 1.0    


kj2    = kj*kj

do j = 0, nfy
   do i = 0, nkx
      k2n(i,j)     = (ki2(i)+kj2(j))
      rk2an(i,j)   = 1.0d0/(k2n(i,j)+alpha)
      kirk2an(i,j) = ki(i)*rk2an(i,j)
      kjrk2an(i,j) = kj(j)*rk2an(i,j)
   enddo
enddo

return
end subroutine init_ops


! ************
! * INIT_RST *
! ************
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



! **************
! * CHECK_DIFF *
! **************
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
