! =======================
! VEGETATION MODULE SIMBA              ! by Axel Kleidon & Pablo Paiewonsky
! =======================

module vegmod
use landmod
implicit none

integer :: jhor

! parameter for soil component

real, parameter :: vws_crit  = 0.25             ! critical soil wetness

! parameter for vegetation component
! all values are in m-k-s units

real, parameter :: q10       = 2.0              ! q10 value for resp

!     * parameters for land surface parameters

real, parameter :: zlaimax   =  9.0               ! max LAI
real, parameter :: zlaimin   =  0.1               ! min LAI for bare soil
real, parameter :: valb_min  =  0.12              ! min albedo for veg
real, parameter :: valb_max  =  0.30              ! max albedo for bare soil
real, parameter :: vsalb_min =  0.20              ! snow albedo for veg
real, parameter :: vwmax_max =  0.5               ! max SWHC for veg
real, parameter :: vwmax_min =  0.05              ! SWHC for bare soil
real, parameter :: vz0_max   =  2.0               ! max roughness for veg 
! not used! Uses pz0_max from landmod instead, which is a namelist parameter (PP)
real, parameter :: vz0_min   =  0.05              ! min roughness for bare soil
real, parameter :: cveg_a    =  2.9               ! forest cover conversion
real, parameter :: cveg_b    =  1.0               ! forest cover conversion
real, parameter :: cveg_c    =  1.570796326794897 ! forest cover conversion
real, parameter :: cveg_d    =  -1.9              ! new forest cover parameter
real, parameter :: cveg_f    =  3.0              ! new soil cover parameter
real, parameter :: cveg_k    =  0.5               ! conversion LAI
real, parameter :: cveg_l    =  0.25              ! LAI-biomass relationship
real, parameter :: ct_crit   =  5.0               ! for temp limit function

!* beta factor defined as GPP  = (1 + beta ln(C/C0))
!  beta  = 0.3-0.6, C0 = 360ppm - Harvey 1989

real, parameter :: co2_ref   = 360.     ! co2 reference concentration in ppm
real, parameter :: co2_comp  = 0.       ! light compensation point in ppm
real, parameter :: co2_sens  = 0.3      ! beta factor

real            :: cveg_e               ! forest cover conversion
real            :: cveg_g               ! soil cover conversion

!     * dynamically calculated land surface parameters

real,dimension(NHOR) :: zforest=   0.0         ! forest cover = dforest(:)
real,dimension(NHOR) :: zwmax  =   0.0         ! bucket depth = dwmax(:)
real,dimension(NHOR) :: dvz0   =   vz0_min     ! roughness length
real,dimension(NHOR) :: dvalb  =   valb_max    ! albedo
real,dimension(NHOR) :: dvrhs  =   0.          ! surface conductance
real,dimension(NHOR) :: dvsoil =   0.          ! soil cover
real :: zalbsn =   0.          ! snow albedo

!     * output of cumulative fluxes

integer :: ibiomass
real    :: zlaim                 ! maximum lai
real    :: zvegm                 ! maximum veg cover
real    :: zft
real    :: zsvp                  ! saturation vapor pressure
real    :: zvpd                  ! vapor pressure deficit
real    :: zbeta
real    :: zglacfree             ! glacier free fraction

end module vegmod


!========!
! VEGINI !
!========!

subroutine vegini
use vegmod
implicit none

if (nrestart == 0) then
   where (dls(:) == 0) dforest(:) = 0.0
endif

end subroutine vegini


!=========!
! VEGSTOP !
!=========!

subroutine vegstop
use vegmod
implicit none

return
end subroutine vegstop


!=========!
! VEGSTEP !
!=========!

subroutine vegstep
use vegmod
implicit none

! Initialization and copy of puma fields needed

cveg_e = atan(cveg_d) ! forest cover
cveg_g = atan(-cveg_f)   ! soil cover
zbeta  = max(0.,1.+co2_sens*log((co2-co2_comp)/(co2_ref-co2_comp)))

! Following plasim arrays are used but not modified
! -------------------------------------------------
! dwatc(:) = soil wetness [m]
! dgppl(:) = GPP light limited
! dswfl(:) = short wave radiation [W/m2]
! devap(:) = surface evaporation (negative)

! Make local copies of some plasim arrays
! They are copied back if NBIOME is set to 1 (interactive)

zforest(:) = dforest(:) ! Forest cover (0.0 - 1.0)
zwmax(:)   = dwmax(:)   ! Bucket depth

! Initialize arrays (declared in plasimmod)

dgpp(:)   = 0.0  ! Gross primary production [kg C m-2 s-1]
dgppl(:)  = 0.0  ! GPP light limited        [kg C m-2 s-1]
dgppw(:)  = 0.0  ! GPP water limited        [kg C m-2 s-1]
dnpp(:)   = 0.0  ! Net primary production   [kg C m-2 s-1]
dveg(:)   = 0.0  ! Vegetation cover         [0..1]

do jhor = 1 , NHOR
  if (dls(jhor) > 0.0 .and. dglac(jhor) < 0.90) then ! land cell with < 90 % glacier

    zglacfree = 1.0 - dglac(jhor) ! glacier free fraction
    ! zsvp = saturation vapor pressure [Pa] (Magnus-Teten)
    ! zvpd = vapor pressure deficit    [Pa] (set to >= 0.01 Pa)

    zsvp = ra1*exp(ra2*(dt(jhor,NLEP)-TMELT)/(dt(jhor,NLEP)-ra4))
    zvpd = max(0.01, zsvp - dp(jhor)*dq(jhor,NLEV)/0.622)

    ! structurally limited leaf area index (eq. 6.6 PlaSim RM)
    ! zlaimax : set to 9.0 - local constant
    ! zlaimin : set to 0.1 - minimal LAI for bare soil

!    zlaim = zlaimax * zforest(jhor) + zlaimin * (1.0 - zforest(jhor))
!    zlaim = zlaimin + 0.63661977236758 * zlaimax * atan(cveg_l * dcveg(jhor))
     zlaim = zlaimin + 0.63661977236758 * zlaimax * atan(cveg_l * 2.0 * plai(jhor) * dcveg(jhor))
! where 0.6366... = 2/pi

    ! structurally limited vegetation cover (eq. 6.5 PlaSim RM)
    ! cveg_k : set to 0.5 (local constant)

    zvegm = 1.0 - exp(-cveg_k * zlaim)

    ! calculate water stress factor (eq. 6.3 PlaSim RM)
    ! water stress = soil wetness / (field capacity * critical soil wetness)

    if (zwmax(jhor) > 0.0) then
       dveg(jhor) = min(1.0,zvegm,max(0.0,dwatc(jhor) / (zwmax(jhor) * vws_crit)))
    endif

    ! calculate gross primary productivity

    zft = min(1.0,max(0.0,(dt(jhor,NLEP) - TMELT) / ct_crit))

    ! light limited gpp [kg C / m2 /s]

    dgppl(jhor) = rlue * zbeta * zft * dveg(jhor) * dfd(jhor,NLEP)

    ! water limited gpp [kg C / m2 /s]
    ! co2    : initialized with 360.0 [ppmv] - member of namelist $RADPAR
    ! co2conv:                               - member of namelist $LANDPAR

    dgppw(jhor) = co2conv * dp(jhor) * dveg(jhor) * (-devap(jhor)) / zvpd * (0.3 * co2)
    dgppw(jhor) = min(1.0,max(0.0,dgppw(jhor))) ! limit to range (0..1)

    ! combine light & water limitations (eq. 6.1 PlaSim RM)

    dgpp(jhor) = min(dgppl(jhor),dgppw(jhor))

    ! net primary production [kg C / m2 /s]
    ! pgrow : initialized to 1.0 from variable forgrow from $LANDPAR

    dnpp(jhor) = pgrow(jhor) * 0.5 * dgpp(jhor) * zglacfree

    ! litterfall [kg C / m2 /s]
    ! tau_veg : initialized to 10 [years] - member of namelist $LANDPAR

    dlitter(jhor) = dcveg(jhor) / tau_veg

    ! heterotrophic respiration [kg C / m2 /s]
    ! q10 : set to 2.0 (local constant)
    ! tau_soil : initialized to 42 [years] - member of namelist $LANDPAR

    dres(jhor) = q10**((dt(jhor,NLEP)-TMELT-10.)/10.) * dcsoil(jhor)/tau_soil

    ! no growth allocation [kg C / m2 /s]

    dnogrow(jhor)  = (1.0 - pgrow(jhor)) * 0.5 * dgpp(jhor)

    ! carbon stored in biomass (ncveg = time accelerator)

    dcveg(jhor) = dcveg(jhor) + (dnpp(jhor) - dlitter(jhor)) * deltsec * ncveg  !! FL

    ! carbon stored in soil

    dcsoil(jhor) = dcsoil(jhor) + (dlitter(jhor) - dres(jhor)) * deltsec * ncveg !! FL

    ! forest cover
    ! plai : above ground growth factor : initialized to 0.5 in LANDMOD

!   zforest(jhor) = atan((2.0*plai(jhor)*dcveg(jhor)-cveg_a) /cveg_b) /cveg_c +cveg_d
!   zforest(jhor) = min(1.0,max(0.0,cveg_e * zforest(jhor)))
    zforest(jhor) = (atan(dcveg(jhor) - cveg_a) - cveg_e) / (cveg_c - cveg_e)
    zforest(jhor) = min(1.0,max(0.0,zforest(jhor)))

    ! soil cover

!   dvsoil(jhor) = atan((2.0*(1.0-plai(jhor))*dcveg(jhor)-cveg_a)/cveg_b)/cveg_c+cveg_d
!   dvsoil(jhor) = min(1.0,max(0.0,cveg_e * dvsoil(jhor)))
    dvsoil(jhor) = (atan(dcveg(jhor) - cveg_f) - cveg_g) / (cveg_c - cveg_g)
    dvsoil(jhor) = min(1.0,max(0.0,dvsoil(jhor)))

    ! discretization of vegetation state

    if (rnbiocats >= 2.0) then
      ibiomass      = zforest(jhor) * rnbiocats
      zforest(jhor) = min(1.0, real(ibiomass)/(rnbiocats-1.0))
      ibiomass      = dvsoil(jhor) * rnbiocats
      dvsoil(jhor)  = min(1.0, real(ibiomass)/(rnbiocats-1.0))
    endif

    ! derivation of land surface parameters

    dlai(jhor)  = -log(1.0 - dveg(jhor))/cveg_k
    dvz0(jhor)  = pz0_max(jhor) * zforest(jhor)+vz0_min*(1.0-zforest(jhor))
    dvalb(jhor) = valb_min * dveg(jhor) + valb_max * (1.0 - dveg(jhor))
    zwmax(jhor) = vwmax_max * dvsoil(jhor) + vwmax_min * (1.0 - dvsoil(jhor))
    dvrhs(jhor) = pgs(jhor) * min(1.0,max(0.0,dwatc(jhor)/(zwmax(jhor)*vws_crit)))

    ! modify albedo and surface conductance due to snow

    if (dsnow(jhor) > 0.0) then

      ! copy PUMA-2 formulation of snow albedo and mix it with forested albedo

      zalbsn = (albsmax - albsmin) * (dt(jhor,NLEP)-263.16) / (TMELT-263.16)
      zalbsn = max(albsmin,min(albsmax,albsmax-zalbsn))
      dvalb(jhor)=dvalb(jhor)+(zalbsn-dvalb(jhor))*dsnow(jhor)/(dsnow(jhor)+0.01)
      dvalb(jhor)=vsalb_min * zforest(jhor) + dvalb(jhor) * (1.-zforest(jhor))
  
      dvrhs(jhor) = 1.0
    endif

    ! interactive coupling

    if (nbiome == 1) then
      dz0(jhor)     = sqrt(dvz0(jhor)*dvz0(jhor)+dz0climo(jhor)*dz0climo(jhor))
      dwmax(jhor)   = zwmax(jhor)
      drhs(jhor)    = dvrhs(jhor)
      dalb(jhor)    = dvalb(jhor)
      dforest(jhor) = zforest(jhor)
    endif

  endif ! (dls(jhor) > 0.0 .and. dglac(jhor) < 0.9)
enddo ! jhor

! Send some arrays to GUI

call guihor("DCVEG"   // char(0),dcveg  ,1,1000.0,0.0)
call guihor("DVSOIL"  // char(0),dvsoil ,1,1000.0,0.0)
call guihor("ZFOREST" // char(0),zforest,1,1000.0,0.0)

return
end subroutine vegstep
