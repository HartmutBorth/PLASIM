*burn.sh v1.0*
--------------
Afterburner cmorizer
Extracts variables from PlaSim output using CMOR conventions.
Adjusts variable names, signs, units and computes derived variables.

Usage: burn.sh [options] infile outfile var
Options:
        -h                Prints this help
        -s                Extract 3D variables on sigma levels
        -a                Afterburner variable (do not cmorize)
        -p lev1,lev2,...  Specify pressure levels (in hPa)

Default pressure levels are  20,30,50,70,100,150,200,250,300,400,500,600,700,850,925,1000 hPa

Variables which can be extracted:
-------------------------------------------------------------------------------------
| 2D CMOR2 variables:                                                               |
-------------------------------------------------------------------------------------
|     Name |                                Long Name |      Units | From org. vars |
-------------------------------------------------------------------------------------
|       ps |                     Surface Air Pressure |         Pa |             pl |
|      psl |                       Sea Level Pressure |         Pa |            psl |
|      tas |             Near-Surface Air Temperature |          K |            tas |
|   tasmin |Daily Minimum Near-Surface Air Temperature|          K |         tasmin |
|   tasmax |Daily Maximum Near-Surface Air Temperature|          K |         tasmax |
|       ts |                      Surface Temperature |          K |             ts |
|      uas |               Eastward Near-Surface Wind |      m s-1 |            uas |
|      vas |              Northward Near-Surface Wind |      m s-1 |            vas |
|  evspsbl |                              Evaporation | kg m-2 s-1 |           evap |
|       pr |                            Precipitation | kg m-2 s-1 |             pr |
|      prc |                 Convective Precipitation | kg m-2 s-1 |            prc |
|      prl |                Large Scale Precipitation | kg m-2 s-1 |            prl |
|     prsn |                            Snowfall Flux | kg m-2 s-1 |           prsn |
|    mrros |                           Surface Runoff | kg m-2 s-1 |           mrro |
|     mrso |              Total Soil Moisture Content | kg m-2 s-1 |           mrso |
|      tsl |                      Temperature of Soil |          K |            tso |
|     rsdt |         TOA Incident Shortwave Radiation |      W m-2 |       rsut rst |
|     rsut |         TOA Outgoing Shortwave Radiation |      W m-2 |           rsut |
|     rlut |          TOA Outgoing Longwave Radiation |      W m-2 |           rlut |
|     rsds |  Surface Downwelling Shortwave Radiation |      W m-2 |       ssru rss |
|     rsus |    Surface Upwelling Shortwave Radiation |      W m-2 |           ssru |
|     rlds |   Surface Downwelling Longwave Radiation |      W m-2 |       stru rls |
|     rlus |     Surface Upwelling Longwave Radiation |      W m-2 |           stru |
|     hfls |          Surface Upward Latent Heat Flux |      W m-2 |           hfls |
|     hfss |        Surface Upward Sensible Heat Flux |      W m-2 |           hfss |
|      clt |                     Total Cloud Fraction |          % |            clt |
|      prw |                         Water Vapor Path |     kg m-2 |            prw |
|     huss |           Near-Surface Specific Humidity |          1 |            hus |
|     hurs |           Near-Surface Specific Humidity |          1 |       tas td2m |
|     orog |                         Surface Altitude |          m |             sg |
|      snd |                               Snow Depth |          m |            snd |
|      snm |                        Surface Snow Melt | kg m-2 s-1 |            snm |
|     tauu |    Surface Downward Eastward Wind Stress |         Pa |           tauu |
|     tauv |   Surface Downward Northward Wind Stress |         Pa |           tauv |
|      sit |                        Sea Ice Thickness |          m |            sit |
|      sic |                            Sea Ice Cover |          % |            sic |
-------------------------------------------------------------------------------------
| 3D CMOR2 variables:                                                               |
-------------------------------------------------------------------------------------
|     Name |                                Long Name |      Units | From org. vars |
-------------------------------------------------------------------------------------
|      hus |                        Specific Humidity |          1 |            hus |
|      hur |                        Relative Humidity |          % |            hur |
|       zg |                      Geopotential Height |          1 |             zg |
|       ta |                          Air Temperature |          K |             ta |
|       ua |                            Eastward Wind |      m s-1 |             ua |
|       va |                           Northward Wind |      m s-1 |             va |
|       wa |                              Upward Wind |      m s-1 |             wa |
|      wap |                           omega (=dp/dt) |     Pa s-1 |            wap |
|       cl |                      Cloud Area Fraction |          % |             cl |
|      clw |      Mass Fraction of Cloud Liquid Water |          1 |            clw |
-------------------------------------------------------------------------------------
| Not standard variables:                                                           |
-------------------------------------------------------------------------------------
|     Name |                                Long Name |      Units | From org. vars |
-------------------------------------------------------------------------------------
|     td2m |        Near-Surface Dewpoint Temperature |          K |           td2m |
|      rst |              TOA Net Shortwave Radiation |      W m-2 |            rst |
|      rls |           Surface Net Longwave Radiation |      W m-2 |            rls |
|      alb |                           Surface Albedo |          1 |            alb |
|      prl |                Large Scale Precipitation | kg m-2 s-1 |            prl |
|     sndc |                Surface Snow Depth Change |      m s-1 |           sndc |
|      stf |                          Stream Function |          m |            stf |
|      lsm |                            Land Sea Mask |          1 |            lsm |
|      mld |              Ocean Mixed Layer Thickness |          m |            mld |
-------------------------------------------------------------------------------------
Variables not in this list are extracted directly using afterburner with no change.
Examples: zeta,d,bld etc.
Note: uas and vas are not available in PlaSim output, derived from ua and va at sigma=1
Note: psl (sea level pressure) is not avaiable in PlaSim output,
      derived from huss (hus at sigma=1), tas and ps using the hypsometric formula
