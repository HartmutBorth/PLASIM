!**********************************!
! Mars module for Planet Simulator !
!**********************************!

subroutine planet_ini
use radmod           ! Includes pumamod

logical :: lex

! Preset parameters may be changed using the namelist

namelist /planet_nl/ nfixorb, eccen, mvelp, obliq &
                , rotspd, siderial_day, solar_day &
                , akap, alr, gascon, ra1, ra2, ra4 &
                , pnu, ga, plarad &
                , gsol0 &
                , yplanet

yplanet = "Mars"     ! Planet name
nplanet = 4          ! Planet index
mars    = 1          ! Switch

! *********
! Astronomy
! *********

eccen        =     0.09341233 ! Eccentricity (AMIP-II value)
mvelp        =   336.04084    ! Longitude of perihelion
nfixorb      =     1          ! Don't use Berger orbits
obliq        =    25.19       ! Obliquity [deg] (AMIP-II)
rotspd       =     1.0        ! Rotation speed (factor)
siderial_day = 88642.0        ! 24h 37m 22s
solar_day    = 88775.0        ! 24h 39m 35s

! **********
! Atmosphere
! **********

akap    =   0.2273    ! Kappa
alr     =   0.0025    ! Lapse rate
gascon  = 188.9       ! Gas constant
psurf   = 636.0       ! Mean surface pressure [Pa]
ra1     = 610.66      ! Parameter for Magnus-Teten-Formula
ra2     =  21.875     ! for saturation vapor pressure
ra4     =   7.65      ! over liquid water
tgr     = 210.0       ! mean ground temperature

! ********
! Calendar
! ********

n_days_per_month =  55 ! Martian days - not earth days
n_days_per_year  = 660 ! simplified calendar 12 x 55 days

! ********
! Numerics
! ********

ndivdamp = 100         ! Initial start divergence damping
pnu      = 0.15        ! Time filter constant
oroscale = 1.0         ! Scale orography

! *******
! Physics
! *******

ga            =       3.728 ! Gravity
plarad        = 3396200.0  ! Radius

! *********
! Radiation
! *********

gsol0         =  595.0     ! Solar constant
no3           =    0.0     ! No ozone

! ********
! Namelist
! ********

if (mypid == NROOT) then
   inquire(file=planet_namelist,exist=lex)
   if (lex) then
      open (nut,file=planet_namelist)
      read (nut,planet_nl)
      close(nut)
   endif
   write(nud,'(/,"****************************************")')
   write(nud,'("* planet_nl from    <",a16,"> *")') planet_namelist
   write(nud,'("****************************************")')
   write(nud,planet_nl)

endif ! (mypid == NROOT)

return
end


subroutine print_planet
use radmod

p_mass        =    0.6419  ! [10^24 kg]
p_volume      =   16.318   ! [10^10 km3]
p_radius_eq   = 3393.0     ! Equatorial radius
p_radius_po   = 3373.0     ! Polar radius
p_radius_me   = 3390.0     ! Mean radius
p_ellipticity =    0.0065  ! Ellipticity
p_density     = 3933.0     ! [kg/m3]
p_albedo      =    0.16    ! Bond albedo
p_blackt      =  216.6     ! Black body temperature
p_sidorbit    =  686.98    ! Siderial orbit period
p_sidrot      =   24.6229  ! Sidereal rotation period
p_inclination =   23.98    ! Equatorial inclination
p_perihelion  =  206.6     ! Perihelion [10^6 km]
p_aphelion    =  249.2     ! Aphelion [10^6 km]

write(nud,4000)
write(nud,1000)
write(nud,1100) 'Simulating:',yplanet
write(nud,1000)
write(nud,2000) 'Parameter','Units','Value'
write(nud,1000)
write(nud,3000) 'Mass'             ,'[10^24 kg]'  ,p_mass
write(nud,3000) 'Volume'           ,'[10^10 km3]' ,p_volume
write(nud,3000) 'Equatorial radius','[km]'        ,p_radius_eq
write(nud,3000) 'Polar radius'     ,'[km]'        ,p_radius_po
write(nud,3000) 'Mean radius'      ,'[km]'        ,p_radius_me
write(nud,3000) 'Ellipticity'      ,' '           ,p_ellipticity
write(nud,3000) 'Mean density'     ,'[kg/m3]'     ,p_density
write(nud,3000) 'Surface gravity'  ,'[m/s2]'      ,ga
write(nud,3000) 'Bond albedo'      ,' '           ,p_albedo
write(nud,3000) 'Solar irradiance' ,'[W/m2]'      ,gsol0
write(nud,3000) 'Black-body temperature','[K]'    ,p_blackt
write(nud,3000) 'Sidereal orbit period' ,'[days]' ,p_sidorbit
write(nud,3000) 'Sidereal rotation period','[hrs]',p_sidrot
write(nud,3000) 'Equatorial inclination'  ,'[deg]',p_inclination
write(nud,3000) 'Perihelion'       ,'[10^6 km]'   ,p_perihelion
write(nud,3000) 'Aphelion'         ,'[10^6 km]'   ,p_aphelion
write(nud,3000) 'Orbit eccentricity'    ,' '      ,eccen
write(nud,3000) 'Gas constant'          ,' '      ,gascon
write(nud,3000) 'Mean surface pressure' ,'[Pa]'   ,psurf 
write(nud,1000)
write(nud,4000)

return

 1000 format(50('*'))
 1100 format('* ',a24,1x,a21,' *')
 2000 format('* ',a24,1x,a11,a10,' *')
 3000 format('* ',a24,1x,a11,f10.4,' *')
 4000 format(/)
      end  


!          ===========
!          PLANET_STEP
!          ===========

subroutine planet_step
return
end

!          ===========
!          PLANET_STOP
!          ===========

subroutine planet_stop
return
end

