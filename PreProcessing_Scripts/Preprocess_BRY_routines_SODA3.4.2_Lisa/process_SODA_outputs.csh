#!/bin/csh

# This script preprocesses SODA Reanalyses outputs
# to make them in the DODS classical Format (as ECCO/ORCA format when downloaded by openDAP procedure)
# format : YnMn
# This is the first step to create BRY file in Roms input format
#
# Step 2: Add bartropic currents (Ubar,Vbar) with add_UVbar.m
# 
# Step 3: Create bry files with make_OGCM_soda.m
#
########################################################
#  Define files to process
########################################################
clear
#
# SODA File names to process (outputs SODA)
# =========================================
set path_in='/media/maillard/HelloWorld/CROCO/SODA/'
set path_out='/media/maillard/HelloWorld/CROCO/SODA/COUPLED_SIMU/'
set MODEL=soda3.4.2_5dy_ocean_reg_
#
# Time period
# ============
set NY_START=2014
set NY_END=2019
set NM_START=1
set NM_END=12
#
########################################################
########################################################

@ NY_END++
@ NM_END++

set NY=$NY_START
while ($NY != $NY_END)
  if ($NY == $NY_START) then
    set NM=$NM_START
  else
    set NM=1
  endif
  set MY_YEAR=$NY
  @ MY_YEAR++
  if ($MY_YEAR == $NY_END) then
    set MONTH_END=$NM_END
  else
    set MONTH_END=13
  endif
   while ($NM != $MONTH_END)

	if ($NM < 10) then
      set TIME=Y${NY}M${NM}
	  set month=0${NM}
	else
      set TIME=Y${NY}M${NM}
	  set month=${NM}
	  
	endif

      echo "Concatene : $TIME"
ncrcat -v temp,salt,u,v,ssh,taux,tauy $path_in/${MODEL}${NY}_${month}*.nc $path_in/${MODEL}${TIME}.nc
      echo "Extract : $TIME"
#ncks -d xt_ocean,349,584 -d yt_ocean,109,230 -d xu_ocean,349,584 -d yu_ocean,109,230 $path_in/${MODEL}${TIME}.nc $path_out/${MODEL}EPAC_${TIME}.nc
#ncks -d xt_ocean,347,586 -d yt_ocean,107,232 -d xu_ocean,347,587 -d yu_ocean,107,232 $path_in/${MODEL}${TIME}.nc $path_out/${MODEL}EPAC_${TIME}.nc
ncks -d xt_ocean,347,586 -d yt_ocean,115,193 -d xu_ocean,347,587 -d yu_ocean,115,193 $path_in/${MODEL}${TIME}.nc $path_out/${MODEL}npac12_${TIME}.nc


rm -f $path_in/${MODEL}${TIME}.nc

 @ NM++
  end
  @ NY++
end

























            
