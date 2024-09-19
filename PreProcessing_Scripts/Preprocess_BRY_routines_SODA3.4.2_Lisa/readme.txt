L. Maillard
19/02/2021
-----------

Preprocess BRY files from SODA data.

1 - download SODA data: ex for regrided soda 3.4.2 product 
        yr=1990
        fn=soda3.4.2_5dy_ocean_reg_${yr}*.nc
        wget -r -l1 --no-parent --progress=bar -nd -A ${fn} https://dsrs.atmos.umd.edu/DATA/soda3.4.2/REGRIDED/ocean/

2 - preprocesses SODA Reanalyses outputs to make them in the DODS classical Format (as ECCO/ORCA format when downloaded by openDAP procedure)
       chose time span and wanted spatial boundaries (to put in the ncks function, should be a bit wider that what you really need) in process_SODA_outputs.csh and execute the file

3 - Add bartropic currents (Ubar,Vbar)
       chose time span in add_UVbar.m
       in matlab, execute start.m then add_UVbar.m

4 - Create BRY files 
       chose time span, select bry and ini creation in crocotools_param.m 
       in matlab, execute start.m then make_OGCM_soda.m
 
 
