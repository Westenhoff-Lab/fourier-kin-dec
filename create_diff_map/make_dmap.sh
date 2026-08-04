#!/bin/tcsh

# Usage:
# ./make_dmap.sh <light_label>
# The workflow for calculating weighted difference electron density maps follows the workflow proposed by Schmidt, 2023

set bin_nam = $1

# ===========================================
# Input files (EDIT THESE)
# ===========================================

set dark_model = dark.pdb
set model_F    = dark.mtz
set dark_obs   = FOBS_dark.mtz
set light_obs  = FOBS_${bin_nam}.mtz

# ===========================================
# Crystal information (EDIT THESE)
# ===========================================

set cell = "164.690 164.690 433.020 90.00 90.00 90.00"
set symm = 91
set grid = "240 240 576"

# =============================
# MTZ column labels (EDIT THESE)
# =============================

# Model MTZ (dark.mtz)
set model_F_col   = F
set model_PHI_col = PHIC

# Dark observed MTZ (FOBS_dark.mtz)
set dark_F_col    = F_DARK
set dark_SIGF_col = SIGF_DARK

# Light observed MTZ (FOBS_<time>.mtz)
set light_F_col   = F_${bin_nam}
set light_SIGF_col = SIGF_${bin_nam}

# =============================

set resmax = 2.60
set scalmin = 30.0
set mapmin = 30.0
set phasmin = 30.0


# =============================

# phase file must be calculated from refined model in dark
# folder
# convert the .mtz file to a readable .hkl and add column with 1.0 as sigmas

mtz2various HKLIN $model_F  HKLOUT model_phs.hkl << mtz_phs
     LABIN FP=$model_F_col PHIC=$model_PHI_col
     OUTPUT USER '(3I5,F12.3,'  1.00  ',F12.3)'
mtz_phs


# get the phase file with added sig column into the system

f2mtz HKLIN model_phs.hkl HKLOUT FC_dark.mtz << f2m_phs 
# SKIP 5  Miss out 5 lines of header
CELL $cell
SYMM $symm
LABOUT H   K  L   FC_D SIG_FC_D PHI_D
CTYPE  H   H  H   F     Q        P
f2m_phs


# cad things together

cad \
HKLIN1 FC_dark.mtz    \
HKLIN2 $dark_obs \
HKLIN3 $light_obs     \
HKLOUT all.mtz \
<< END-cad

LABIN FILE 1 E1=FC_D E2=SIG_FC_D E3=PHI_D
CTYP  FILE 1 E1=F E2=Q E3=P 
LABIN FILE 2 E1=$dark_F_col E2=$dark_SIGF_col
CTYP  FILE 2 E1=F E2=Q
LABIN FILE 3 E1=$light_F_col E2=$light_SIGF_col
CTYP  FILE 3 E1=F E2=Q

END
END-cad


# scale the datasets
# 1 scale dark to FC dark
# 2 scale light to dark


echo " SCALEIT NUMBER 1 "


scaleit \
HKLIN all.mtz    \
HKLOUT all_sc1.mtz    \
<< END-scaleit1
TITLE FPHs scaled to FP
reso $scalmin $resmax      # Usually better to exclude lowest resolution data
WEIGHT
#Exclude FP data if: FP < 5*SIGFP & if FMAX > 1000000
EXCLUDE FP SIG 3 FMAX 10000000
REFINE ANISOTROPIC 
LABIN FP=FC_D SIGFP=SIG_FC_D  -
  FPH1=$dark_F_col SIGFPH1=$dark_SIGF_col -
  FPH2=$light_F_col SIGFPH2=$light_SIGF_col
CONV ABS 0.0001 TOLR  0.000000001 NCYC 150
END
END-scaleit1

echo " SCALEIT NUMBER 2 "

scaleit \
HKLIN all_sc1.mtz    \
HKLOUT all_sc2.mtz    \
<< END-scaleit2
TITLE FPHs scaled to FP
reso $scalmin $resmax      # Usually better to exclude lowest resolution data
WEIGHT
#Exclude FP data if: FP < 5*SIGFP & if FMAX > 1000000
REFINE ANISOTROPIC 
EXCLUDE FP SIG 3 FMAX 10000000
LABIN FP=$dark_F_col SIGFP=$dark_SIGF_col -
  FPH1=$light_F_col SIGFPH1=$light_SIGF_col
CONV ABS 0.0001 TOLR  0.000000001 NCYC 150
END
END-scaleit2

goto free

# here fhscal
fhscal:

fhscal \
hklin all_sc1.mtz \
hklout all_sc2.mtz <<END-fhscal
TITLE scale  by Kraut method
BIAS 1 ! if we trust the standard deviations
LABIN FP=$dark_F_col SIGFP=$dark_SIGF_col FPH=$light_F_col SIGFPH=$light_SIGF_col
AUTO
END
END-fhscal

free:
freerflag HKLIN all_sc2.mtz HKLOUT all_sc2_free.mtz <<+
freerfrac 0.05
+


echo "Calculate unweighted maps"


maps:

fft HKLIN all_sc2.mtz MAPOUT ${bin_nam}_nonw.map << endfft
  RESO $mapmin  $resmax
  GRID $grid
  BINMAPOUT
  LABIN F1=$light_F_col SIG1=$light_SIGF_col F2=$dark_F_col SIG2=$dark_SIGF_col PHI=PHI_D
endfft

# split the scaled file into:
# light_scaled file with amplitudes and sigmas
# dark_scaled file with amplitudes and sigmas
# dark_phase file with amplitudes, sigmas and phases
# calculate the weighted map from these files

mtz2various HKLIN all_sc2.mtz  HKLOUT ${bin_nam}_scaled.hkl << end_mtzv1
     LABIN FP=$light_F_col SIGFP=$light_SIGF_col
     OUTPUT USER '(3I5,2F12.3)'
     RESOLUTION 60.0 $resmax 
end_mtzv1


mtz2various HKLIN all_sc2.mtz  HKLOUT dark_scaled.hkl << end_mtzv2
     LABIN FP=$dark_F_col SIGFP=$dark_SIGF_col
     OUTPUT USER '(3I5,2F12.3)'
     RESOLUTION 33.0 $resmax 
end_mtzv2


mtz2various HKLIN all_sc2.mtz  HKLOUT dark_phase.hkl << end_mtzv3
     LABIN FP=FC_D SIGFP=SIG_FC_D PHIC=PHI_D
     OUTPUT USER '(3I5,3F12.3)'
     RESOLUTION 33.0 $resmax 
end_mtzv3

# the python script will produce a difference structure factor file
# h k l DF weight Phase
# according to Ursby and Bourgeois weighting

echo "Ursby and Bourgeois weighting"

python weight.py --light ${bin_nam}_scaled.hkl --dark dark_scaled.hkl --phase dark_phase.hkl --output ${bin_nam}_diff.phs

# convert .phs into .mtz

f2mtz HKLIN ${bin_nam}_diff.phs HKLOUT ${bin_nam}_dwt.mtz << end_weight
CELL $cell
SYMM $symm
LABOUT H   K  L   DOBS_${bin_nam}  FOM_${bin_nam}  PHI
CTYPE  H   H  H   F      W   P
END
end_weight

#calculate weighted difference map

fft HKLIN ${bin_nam}_dwt.mtz MAPOUT ${bin_nam}_wd.map << END-wfft
  RESO $mapmin  $resmax 
  GRID $grid
  BINMAPOUT
  LABIN F1=DOBS_${bin_nam} W=FOM_${bin_nam} PHI=PHI
END-wfft

mapmask mapin ${bin_nam}_wd.map mapout ${bin_nam}_wdex.map xyzin $dark_model << ee
extend xtal
border 0.0
ee


#cp {$bin_nam}_wdex.map 30ps/dmaps/${bin_nam}.wdex.map


