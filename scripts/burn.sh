#!/bin/bash
## burn.sh
## Afterburner cmorizer
## (2018) by Jost von Hardenberg <j.vonhardenberg AT isac.cnr.it>
VERSION=v1.0

#set -ex
BURNER=/work/users/jost/plasim/postprocessor/burn7.x
BURNER=/Users/jost/plasim/plasimnew/postprocessor/burn7.x
PLEVS="20,30,50,70,100,150,200,250,300,400,500,600,700,850,925,1000"
VARS3D="wap wa zeta stf psi d zg hur ta ua va hus clw cl spd psl"
GACC=9.80665
HUM_A=6.116441
HUM_m=7.591386 
HUM_Tn=240.7263

function makevar {
	outvar=$1
	no3d=false;
	case $outvar in
		"ps")      invar="ps";   standard_name="surface_air_pressure" ;  long_name="Surface Air Pressure"; units="Pa";         expr="ps=ps*100";;
		"psl")      invar="psl";   standard_name="air_pressure_at_sea_level" ;  long_name="Sea Level Pressure"; units="Pa";     expr="psl=psl*100";; 
		"tas")     invar="tas";  standard_name="air_temperature" ;long_name="Near-Surface Air Temperature"; units="K"; expr="tas=tas";;
		"tasmin")     invar="tasmin";  standard_name="air_temperature" ;long_name="Daily Minimum Near-Surface Air Temperature"; units="K"; expr="tasmin=tasmin";;
		"tasmax")     invar="tasmax";  standard_name="air_temperature" ;long_name="Daily Maximum Near-Surface Air Temperature"; units="K"; expr="tasmax=tasmax";;
		"ts")     invar="ts";    standard_name="surface_temperature" ;long_name="Surface Temperature"; units="K"; expr="ts=ts";;
		"uas")     invar="ua";    standard_name="eastward_wind" ;long_name="Eastward Near-Surface Wind"; units="m s-1"; expr="uas=ua"; no3d=true;;
		"vas")     invar="va";    standard_name="northward_wind" ;long_name="Northward Near-Surface Wind"; units="m s-1"; expr="vas=va"; no3d=true;;
		"evspsbl") invar="evap"; standard_name="water_evaporation_flux" ;long_name="Evaporation";          units="kg m-2 s-1"; expr="evspsbl=-evap*1000";;
		"pr")     invar="pr";    standard_name="precipitation_flux" ;long_name="Precipitation"; units="kg m-2 s-1"; expr="pr=pr*1000";;
		"prc")     invar="prc";    standard_name="convective_precipitation_flux" ;long_name="Convective Precipitation"; units="kg m-2 s-1"; expr="prc=prc*1000";;
		"prsn")     invar="prsn";    standard_name="snowfall_flux" ;long_name="Snowfall Flux"; units="kg m-2 s-1"; expr="prsn=prsn*1000";;
		"mrros")     invar="mrro";    standard_name="surface_runoff_flux" ;long_name="Surface Runoff"; units="kg m-2 s-1"; expr="mrros=mrro*1000";;
		"mrso")     invar="mrso";    standard_name="soil_moisture_content" ;long_name="Total Soil Moisture Content"; units="kg m-2 s-1"; expr="mrso=mrso*1000";;
		"tsl")     invar="tso";    standard_name="soil_temperature" ;long_name="Temperature of Soil"; units="K"; expr="tsl=tso";;
		"rsdt")     invar="rsut rst";    standard_name="toa_incoming_shortwave_flux" ;long_name="TOA Incident Shortwave Radiation"; units="W m-2"; expr="rsdt=rst-rsut";;
		"rsut")     invar="rsut";    standard_name="toa_outgoing_shortwave_flux" ;long_name="TOA Outgoing Shortwave Radiation"; units="W m-2"; expr="rsut=-rsut";;
		"rlut")     invar="rlut";    standard_name="toa_outgoing_longwave_flux" ;long_name="TOA Outgoing Longwave Radiation"; units="W m-2"; expr="rlut=-rlut";;
		"rsds")     invar="ssru rss";    standard_name="surface_downwelling_shortwave_flux_in_air" ;long_name="Surface Downwelling Shortwave Radiation"; units="W m-2"; expr="rsds=rss-ssru";;
		"rsus")     invar="ssru";    standard_name="surface_upwelling_shortwave_flux_in_air" ;long_name="Surface Upwelling Shortwave Radiation"; units="W m-2"; expr="rsus=-ssru";;
		"rlds")     invar="stru rls";    standard_name="surface_downwelling_longwave_flux_in_air" ;long_name="Surface Downwelling Longwave Radiation"; units="W m-2"; expr="rlds=rls-stru";;
		"rlus")     invar="stru";    standard_name="surface_upwelling_longwave_flux_in_air" ;long_name="Surface Upwelling Longwave Radiation"; units="W m-2"; expr="rlus=-stru";;
		"hfls")     invar="hfls";    standard_name="surface_upward_latent_heat_flux" ;long_name="Surface Upward Latent Heat Flux"; units="W m-2"; expr="hfls=-hfls";;
		"hfss")     invar="hfss";    standard_name="surface_upward_sensible_heat_flux" ;long_name="Surface Upward Sensible Heat Flux"; units="W m-2"; expr="hfss=-hfss";;
		"clt")     invar="clt";  standard_name="cloud_area_fraction" ;   long_name="Total Cloud Fraction"; units="%"; expr="clt=clt*100";;
		"prw")     invar="prw";    standard_name="atmosphere_water_vapor_content" ;long_name="Water Vapor Path"; units="kg m-2"; expr="prw=prw";;
		"huss")     invar="hus";    standard_name="specific_humidity" ;long_name="Near-Surface Specific Humidity"; units="1"; expr="huss=hus"; no3d=true;;
		"hurs")     invar="tas td2m";    standard_name="specific_humidity" ;long_name="Near-Surface Specific Humidity"; units="1"; expr="hurs=100*10^($HUM_m*( td2m/(td2m + $HUM_Tn )-tas/(tas + $HUM_Tn )   ))	"; no3d=true;;
		"clw")     invar="clw";    standard_name="mass_fraction_of_cloud_liquid_water_in_air" ;long_name="Mass Fraction of Cloud Liquid Water"; units="1"; expr="clw=clw";;
		"hus")     invar="hus";    standard_name="specific_humidity" ;long_name="Specific Humidity"; units="1"; expr="hus=hus";;
		"hur")     invar="hur";    standard_name="relative_humidity" ;long_name="Relative Humidity"; units="%"; expr="hur=hur";;
		"zg")     invar="zg";    standard_name="geopotential_height" ;long_name="Geopotential Height"; units="1"; expr="zg=zg";;
		"ta")     invar="ta";    standard_name="air_temperature" ;long_name="Air Temperature"; units="K"; expr="ta=ta";;
		"ua")     invar="ua";    standard_name="eastward_wind" ;long_name="Eastward Wind"; units="m s-1"; expr="ua=ua";;
		"va")     invar="va";    standard_name="northward_wind" ;long_name="Northward Wind"; units="m s-1"; expr="va=va";;
		"wa")     invar="wa";    standard_name="upward_wind" ;long_name="Upward Wind"; units="m s-1"; expr="wa=wa";;
		"wap")     invar="wap";    standard_name="lagrangian_tendency_of_air_pressure" ;long_name="omega (=dp/dt)"; units="Pa s-1"; expr="wap=wap";;
		"cl")     invar="cl";    standard_name="cloud_area_fraction_in_atmosphere_layer" ;long_name="Cloud Area Fraction"; units="%"; expr="cl=cl*100";;
		"orog")     invar="sg";    standard_name="surface_altitude" ;long_name="Surface Altitude"; units="m"; expr="orog=sg/$GACC";;
		"snd")     invar="snd";    standard_name="surface_snow_thickness" ;long_name="Snow Depth"; units="m"; expr="snd=snd";;
		"snm")     invar="snm";    standard_name="surface_snow_melt_flux" ;long_name="Surface Snow Melt"; units="kg m-2 s-1"; expr="snm=snm*1000";;
		"tauu")     invar="tauu";    standard_name="surface_downward_eastward_stress" ;long_name="Surface Downward Eastward Wind Stress"; units="Pa"; expr="tauu=tauu";;
		"tauv")     invar="tauv";    standard_name="surface_downward_northward_stress" ;long_name="Surface Downward Northward Wind Stress"; units="Pa"; expr="tauv=tauv";;

		"td2m")     invar="td2m";    standard_name="dewpoint_temperature" ;long_name="Near-Surface Dewpoint Temperature"; units="K"; expr="td2m=td2m";;
		"rst")     invar="rst";    standard_name="toa_net_shortwave_flux" ;long_name="TOA Net Shortwave Radiation"; units="W m-2"; expr="rst=rst";;
		"rss")     invar="rss";    standard_name="surface_net_shortwave_flux_in_air" ;long_name="Surface Net Shortwave Radiation"; units="W m-2"; expr="rss=rss";;
		"rls")     invar="rls";    standard_name="surface_net_longwave_flux_in_air" ;long_name="Surface Net Longwave Radiation"; units="W m-2"; expr="rls=rls";;
		"alb")     invar="alb";    standard_name="surface_albedo" ;long_name="Surface Albedo"; units="1"; expr="alb=alb";;
		"prl")     invar="prl";    standard_name="large_scale_precipitation_flux" ;long_name="Large Scale Precipitation"; units="kg m-2 s-1"; expr="prl=prl*1000";;
		"sndc")     invar="sndc";    standard_name="surface_snow_depth_change" ;long_name="Surface Snow Depth Change"; units="m s-1"; expr="sndc=sndc";;
		"stf")     invar="stf";    standard_name="stream_function" ;long_name="Stream Function"; units="m"; expr="stf=stf";;
		"lsm")     invar="lsm";    standard_name="land_sea_mask" ;long_name="Land Sea Mask"; units="1"; expr="lsm=lsm";;

		"sit")     invar="sit";    standard_name="sea_ice_thickness" ;long_name="Sea Ice Thickness"; units="m"; expr="sit=sit";;
		"sic")     invar="sic";    standard_name="sea_ice_cover" ;long_name="Sea Ice Cover"; units="%"; expr="sic=sic*100";;
		"mld")     invar="mld";    standard_name="ocean_mixed_layer_thickness" ;long_name="Ocean Mixed Layer Thickness"; units="m"; expr="mld=mld";;
		*)         invar=$outvar ;    expr="none";;
	esac
}

function help {
HELPVARS2D="ps psl tas tasmin tasmax ts uas vas evspsbl pr prc prl prsn mrros mrso tsl rsdt rsut rlut rsds rsus rlds rlus hfls hfss clt prw huss hurs orog snd snm tauu tauv sit sic"
HELPVARS3D="hus hur zg ta ua va wa wap cl clw"
HELPEXTRA="td2m rst rls alb prl sndc stf lsm mld"

echo "burn.sh $VERSION"
echo "Afterburner cmorizer"
echo "Extracts variables from PlaSim output using CMOR conventions."
echo "Adjusts variable names, signs, units and computes derived variables."
echo
echo Usage: burn.sh [options] infile outfile var 
echo Options:
echo "        -h                Prints this help"
echo "        -s                Extract 3D variables on sigma levels"
echo "        -a                Afterburner variable (do not cmorize)"
echo "        -p lev1,lev2,...  Specify pressure levels (in hPa)"
echo
echo "Default pressure levels are " $PLEVS "hPa"
echo
echo Variables which can be extracted:
echo -------------------------------------------------------------------------------------
echo "| 2D CMOR2 variables:                                                               |"
echo -------------------------------------------------------------------------------------
printf "| % 8s | % 40s | % 10s | % 14s |\n" Name "Long Name" Units "From org. vars"
echo -------------------------------------------------------------------------------------
for v in $HELPVARS2D; do
	makevar $v
	printf "| % 8s | % 40s | % 10s | % 14s |\n" $v "$long_name" "$units" "$invar" 
done
echo -------------------------------------------------------------------------------------
echo "| 3D CMOR2 variables:                                                               |"
echo -------------------------------------------------------------------------------------
printf "| % 8s | % 40s | % 10s | % 14s |\n" Name "Long Name" Units "From org. vars"
echo -------------------------------------------------------------------------------------
for v in $HELPVARS3D; do
	makevar $v
	printf "| % 8s | % 40s | % 10s | % 14s |\n" $v "$long_name" "$units" "$invar" 
done
echo -------------------------------------------------------------------------------------
echo "| Not standard variables:                                                           |"
echo -------------------------------------------------------------------------------------
printf "| % 8s | % 40s | % 10s | % 14s |\n" Name "Long Name" Units "From org. vars"
echo -------------------------------------------------------------------------------------
for v in $HELPEXTRA; do
	makevar $v
	printf "| % 8s | % 40s | % 10s | % 14s |\n" $v "$long_name" "$units" "$invar" 
done
echo -------------------------------------------------------------------------------------

echo Variables not in this list are extracted directly using afterburner with no change. 
echo Examples: zeta,d,bld etc.
echo Note: uas and vas are not available in PlaSim output, derived from ua and va at sigma=1
echo "Note: psl (sea level pressure) is not avaiable in PlaSim output," 
echo "      derived from huss (hus at sigma=1), tas and ps using the hypsometric formula"
}

function fixvar {
   local infile=$1

   if [ "$expr" != "none" ] && [ $after != 1 ]; then 
	  cdo setattribute,$outvar@standard_name="$standard_name" -setattribute,$outvar@long_name="$long_name" -setattribute,$outvar@units="$units" -expr,"$expr" $infile $2
    else
	  cp $infile $2
   fi
}

function getvar {
  local infile=$1
  local outfile=$2
  local code=$3
  local lev1=1
  rm -f burn.log
  if [ "$no3d" == "true" ]; then
	 lev1=10
  fi 
  if [[ " $VARS3D " =~ .*\ $code\ .* ]] && [ "$no3d" != "true" ]; then
	if [ $sigma == 0 ]; then
	   cat > burn$$.nl <<EOF
	   code=$code,vtype=p,htype=g,mean=0,netcdf=1
	   hpa=$PLEVS
EOF
        else
	   cat > burn$$.nl <<EOF
	   code=$code,vtype=s,htype=g,mean=0,netcdf=1
EOF
        fi
  else
	cat > burn$$.nl <<EOF2
	code=$code,vtype=s,htype=g,mean=0,netcdf=1
	modlev=$lev1
EOF2
  fi
  $BURNER $infile $outfile < ./burn$$.nl >> burn.log
  rm ./burn$$.nl
}

##################
## Main 
##################

sigma=0
after=0
while getopts "h?sap:" opt; do
   case "$opt" in
        h|\?) help ;;
        s)  sigma=1 ;;
        a)  after=1 ;;
        p)  PLEVS=$OPTARG;;
   esac
done
shift $((OPTIND-1))

if [ $# -lt 3 ] ; then
help
exit
fi

infile=$1 
outfile=$2
outvar=$3

makevar $outvar

rm -f temp$$.*.nc
for var in $invar
do
	getvar $infile temp$$.$var.nc $var
done
cdo merge temp$$.*.nc temp$$.all.nc 2>> burn.log 

fixvar temp$$.all.nc $outfile $outvar 2>> burn.log
rm -f temp$$.*.nc
