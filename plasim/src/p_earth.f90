!***********************************!
! Earth module for Planet Simulator !
!***********************************!

subroutine planet_ini
use radmod           ! Includes pumamod

logical :: lex

! Preset parameters may be changed using the namelist

namelist /planet_nl/ nfixorb, eccen, mvelp, obliq  &
                , rotspd, sidereal_day, solar_day  &
                , sidereal_year, tropical_year     &
                , akap, alr, gascon, ra1, ra2, ra4 &
                , pnu, ga, plarad &
                , gsol0 &
                , yplanet

yplanet = "Earth"    ! Planet name
nplanet = 3          ! 3rd. stone from the sun

! *********
! Astronomy
! *********

eccen         =     0.016715  ! Eccentricity (AMIP-II value)
mvelp         =   102.7       ! Longitude of perihelion
obliq         =    23.441     ! Obliquity [deg] (AMIP-II)
rotspd        =     1.0       ! Rotation speed (factor)
sidereal_day  =    86164.0916 !      23h 56m 04s
solar_day     =    86400.0    !      24h 00m 00s
sidereal_year = 31558149.0    ! 365d 06h 09m 09s
tropical_year = 31556956.0    ! 365d 05h 49m 16s

! **********
! Atmosphere
! **********

akap    =   0.286     ! Kappa (Poisson constant R/Cp)
alr     =   0.0065    ! Lapse rate
gascon  = 287.0       ! Gas constant
ra1     = 610.78      ! Parameter for Magnus-Teten-Formula
ra2     =  17.2693882 ! for saturation vapor pressure
ra4     =  35.86      ! over liquid water

! ********
! Numerics
! ********

pnu     = 0.1        ! Time filter constant

! *******
! Physics
! *******

ga            =       9.80665 ! Gravity (mean on NN)
plarad        = 6371220.0     ! Radius

! *********
! Radiation
! *********

gsol0         = 1365.0     ! Solar constant

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

p_mass        =    5.9736  ! [10^24 kg]
p_volume      =  108.321   ! [10^10 km3]
p_radius_eq   = 6378.0     ! Equatorial radius
p_radius_po   = 6356.0     ! Polar radius
p_ellipticity =    0.0034  ! Ellipticity
p_density     = 5520.0     ! [kg/m3]
p_albedo      =    0.385   ! Bond albedo
p_blackt      =  247.3     ! Black body temperature
p_perihelion  =  147.1     ! Perihelion [10^6 km]
p_aphelion    =  152.1     ! Aphelion [10^6 km]
p_sidorbit    =  sidereal_year / sidereal_day ! Sidereal orbit period

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
write(nud,3000) 'Mean radius'      ,'[km]'        ,plarad/1000.0
write(nud,3000) 'Ellipticity'      ,' '           ,p_ellipticity
write(nud,3000) 'Mean density'     ,'[kg/m3]'     ,p_density
write(nud,3000) 'Surface gravity'  ,'[m/s2]'      ,ga
write(nud,3000) 'Bond albedo'      ,' '           ,p_albedo
write(nud,3000) 'Solar irradiance' ,'[W/m2]'      ,gsol0
write(nud,3000) 'Black-body temperature','[K]'    ,p_blackt
write(nud,3000) 'Sidereal orbit period' ,'[days]' ,p_sidorbit
write(nud,3000) 'Sidereal rotation period','[h]'  ,sidereal_day/3600.0
write(nud,3000) 'Perihelion'       ,'[10^6 km]'   ,p_perihelion
write(nud,3000) 'Aphelion'         ,'[10^6 km]'   ,p_aphelion

if (nfixorb /= 0) then
   write(nud,3010) 'Using fixed orbit'       ,' '
   write(nud,3000) 'Longitude of perihelion' ,'[deg]' ,mvelp
   write(nud,3000) 'Equatorial inclination'  ,'[deg]' ,obliq
   write(nud,3000) 'Orbit eccentricity'      ,' '     ,eccen
else
   write(nud,3010) 'Using Berger orbit' ,'nfixorb=0'
endif

write(nud,3000) 'Rotation factor'       ,' '      ,rotspd
write(nud,3000) 'Gas constant'          ,' '      ,gascon
write(nud,1000)
write(nud,4000)

return

 1000 format(50('*'))
 1100 format('* ',a24,1x,a21,' *')
 2000 format('* ',a24,1x,a11,a10,' *')
 3000 format('* ',a24,1x,a11,f10.4,' *')
 3010 format('* ',a24,1x,a11,10x  ,' *')
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

