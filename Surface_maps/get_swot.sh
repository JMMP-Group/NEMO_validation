#!/bin/bash

start_date=20230328
final_date=20250112

user=XXXX # filled by user
pass=XXXX # filled by user

 #1. Loop over the dates
date=$(date +%Y%m%d -d $start_date)
while test  ${date} -lt ${final_date};  do
    YYYY=$(date +%Y -d $date)
    MM=$(date +%m -d $date)
    DD=$(date +%d -d $date)
    echo "$YYYY-$MM-$DD"
	
    out_path=/gws/ssde/j25b/jmmp/CO9_AMM15_validation/SWOT/
    fn_out=${out_path}dt_global_allsat_phy_l4_$YYYY$MM$DD.nc 
    dwnld_str="https://tds-odatis.aviso.altimetry.fr/thredds/ncss/grid/dataset-duacs-experimental-dt-phy-grids-nadirs-and-wide-swath/v3_0/miost/dt_global_allsat_phy_l4_${YYYY}${MM}${DD}_20250112.nc?var=ugos&var=vgos&north=65&west=-30&east=15&south=40&horizStride=1&time_start=${YYYY}-${MM}-${DD}T00:00:00Z&time_end=${YYYY}-${MM}-${DD}T00:00:00Z&&&accept=netcdf4-classic&addLatLon=true"
    wget -O $fn_out $dwnld_str --user=user --password=pass
    date=$(date +%Y%m%d -d $date' +1 day')
done

