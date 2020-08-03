#!/bin/bash
dev=png
if [ psi.grd -nt psi.nc ]; then
    cdo -f nc import_binary psi.ctl psi.nc
fi
ncl -nQ dev=\"${dev}\" plot_inv_qgpv.ncl
ncl -nQ dev=\"${dev}\" plot_slp.ncl
ncl -nQ dev=\"${dev}\" plot_zuv.ncl