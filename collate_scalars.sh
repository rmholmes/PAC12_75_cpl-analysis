#!/bin/bash

nam="PAC12_75_cpl"
exp="exp07"

mkdir tmp_sed
rm tmp_sed/out_all.txt

for f in /g/data/e14/rmh561/croco/archive/${nam}/${nam}_${exp}/jobs_${nam}_exp*/20*/croco.log
do
    echo $f
    sed -n '/^ STEP/,${p;/^ CREATE RST/q}' $f > tmp_sed/out_tmp.txt
    sed -i '1,1d' tmp_sed/out_tmp.txt
    sed -i '/ CREATE RST/d' tmp_sed/out_tmp.txt
    sed -i '/      ONLINE_BULK/d' tmp_sed/out_tmp.txt
    sed -i '/ Open Meteo/d' tmp_sed/out_tmp.txt
    sed -i '/      GET_BRY/d' tmp_sed/out_tmp.txt

    cat tmp_sed/out_tmp.txt >> tmp_sed/out_all.txt
    rm tmp_sed/out_tmp.txt
done

mv tmp_sed/out_all.txt ${nam}_${exp}_scalarvars.txt
rm -r tmp_sed/

         
# Extract all lines from a given file:
# sed -n '/^ STEP/,${p;/^ CREATE RST/q}' croco.log > out.txt

# Remove top and bottom line:
# sed '1,1d' out.txt | head -n -1 > out2.txt

# Remove bad lines:
# sed -i '/      ONLINE_BULK/d' out2.txt > out3.txt
# sed '/ Open Meteo/d' out3.txt > out4.txt
# sed '/      GET_BRY/d' out4.txt > out5.txt


