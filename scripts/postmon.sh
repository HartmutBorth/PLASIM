#!/bin/bash
# Extracts main variables from PlaSim output and computes monthly means
# Edit $POSTDIR and $INDIR below
#set -ex

VARS2D="ps psl tas tasmin tasmax ts uas vas evspsbl pr prc prl prsn mrros mrso tsl rsdt rsut rlut rsds rsus rlds rlus hfls hfss clt prw huss snd snm tauu tauv sit sic"
VARS3D="hus hur zg ta ua va wa wap cl clw"
VARSEXTRA="rst rls prl sndc stf mld"
VARSFX="orog lsm"
PLEVS="20,30,50,70,100,150,200,250,300,400,500,600,700,850,925,1000"

BURNSH="/home/jost/work/plasim/scripts/burn.sh"
#BURNSH="/home/jost/work/plasim/scripts/burn.sh -p $PLEVS" # uncomment if you want to specify different levels

if [ $# -lt 3 ] ; then
        echo Usage: post.sh exp year1 year2	
	exit
fi

exp=$1
y1=$2
y2=$3
fn=MOST

POSTDIR="$SCRATCH/plasim/exps/$exp/post"
INDIR="$WORK/plasim/plasim/exps/$exp"

mkdir -p $POSTDIR
y1f=$(printf "%03d" $y1)

mkdir -p $POSTDIR/fx
for var in $VARSFX
do
	echo "Processing $var (fixed)"
     $BURNSH $INDIR/$fn.$y1f $POSTDIR/fx/temp$$.nc $var
     cdo -s seltimestep,1 $POSTDIR/fx/temp$$.nc $POSTDIR/fx/${exp}_${var}.nc 
     rm -f $POSTDIR/fx/temp$$.nc 
done

for var in $VARS2D $VARS3D $VARSEXTRA
do
  for (( y=$y1 ; y<=$y2; y++ )); do
     yf=$(printf "%03d" $y)
     echo "Processing $var year $yf"
     mkdir -p $POSTDIR/$yf
     $BURNSH $INDIR/$fn.$yf $POSTDIR/$yf/temp$$.nc $var
     cdo -s monmean $POSTDIR/$yf/temp$$.nc $POSTDIR/$yf/${exp}_${yf}_${var}.nc 
     rm -f $POSTDIR/$yf/temp$$.nc 
  done
done
