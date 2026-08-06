#!/bin/tcsh

# Usage:
# ./make_extrapol_map.sh <addition> <timepoint>

set add = $1
set tim = $2

# ===========================================
# Input files
# ===========================================

set light_file = ${tim}_scaled.hkl
set diff_file  = ${tim}_diff.phs
set phase_file = dark_phase.hkl

# ===========================================
# Output
# ===========================================

set outdir = add${add}

if (! -d $outdir) then
    echo "Creating directory $outdir"
    mkdir -p $outdir
endif

# ===========================================
# Crystal information
# ===========================================

set cell = "164.690 164.690 433.020 90.00 90.00 90.00"
set symm = 91
set grid = "240 240 576"

# ===========================================
# Resolution
# ===========================================

set resmin = 30.0
set resmax = 2.60

# ===========================================
# Extrapolate dark structure factors
# ===========================================

python extrapolate.py --light ${light_file} --diff ${diff_file} --phase ${phase_file} --scale ${add} --output ${outdir}/${tim}_darkFC_extrap.hkl


f2mtz HKLIN ${outdir}/${tim}_darkFC_extrap.hkl HKLOUT ${outdir}/temp_dummy.mtz << END-f2mtz
CELL $cell
SYMM $symm
LABOUT H   K  L   F_${tim}  SIGF_${tim}  PHI_${tim} 
CTYPE  H   H  H   F      Q                P
END
END-f2mtz

cad HKLIN1 ${outdir}/temp_dummy.mtz HKLOUT ${outdir}/temp_${tim}_extra_dum.mtz << END-cad
LABIN FILE 1 E1=F_${tim} E2=SIGF_${tim} E3=PHI_${tim} 
CTYPE  FILE 1 E1=F E2=Q E3=P
RESO FILE 1 $resmin $resmax
END-cad


freerflag HKLIN ${outdir}/temp_${tim}_extra_dum.mtz HKLOUT ${outdir}/${tim}_extra${add}.mtz << END-free
freerfrac 0.05
END-free

rm $outdir/temp*

