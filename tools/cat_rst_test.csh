#!/bin/csh
unset noclobber

#--- initialize test mode
touch RSTTEST

#--- remove files from possible previous test
rm -f cat_rstini
rm -f one*
rm -f two*
rm -f full*

#--- set namelist parameter nsteps in cat_namelist to 10000 
sed -i -e 's/nsteps.*/nsteps = 10000/' cat_namelist

#--- deactivate the gui in cat_namelist
sed -i -e 's/ngui.*/ngui = 0/' cat_namelist

#--- compile and run first part of run
make
make run
cp cat_diag   one_cat_diag
cp cat_gp     one_cat_gp
cp cat_rstfin one_cat_rstfin
cp cat_tseri  one_cat_tseri
mv cat_rstfin cat_rstini
#--- run second part of run
make run
cp cat_diag   two_cat_diag
cp cat_gp     two_cat_gp
cp cat_rstfin two_cat_rstfin
cp cat_tseri  two_cat_tseri

#--- prepare namelist file for total run
sed -i -e 's/nsteps.*/nsteps = 20000/' cat_namelist

#--- prepare and do total run
rm cat_rstini
make run
cp cat_diag   full_cat_diag
cp cat_gp     full_cat_gp
cp cat_rstfin full_cat_rstfin
cp cat_tseri  full_cat_tseri

#--- compare if final restart files differ
diff two_cat_rstfin full_cat_rstfin

if ($? == 0) echo "Restart test passed"
