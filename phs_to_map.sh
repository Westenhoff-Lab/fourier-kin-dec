#!/bin/tcsh

# While executing, ex:./make_dmap.sh 731fs

set bin_nam = $1
set resmax = 1.80
set scalmin = 30.0
set mapmin = 30.0
set phasmin = 30.0

set cell = "66.830 66.830 40.920 90.00 90.00 120.00"
set symm = 173
set grid = "120 120 80"

f2mtz HKLIN ${bin_nam}.phs HKLOUT ${bin_nam}_dwt.mtz << end_weight
CELL $cell
SYMM $symm
LABOUT H   K  L   DOBS_${bin_nam}  FOM_${bin_nam}  PHI
CTYPE  H   H  H   F      W   P
END
end_weight

fft HKLIN ${bin_nam}_dwt.mtz MAPOUT ${bin_nam}_wd.map << END-wfft
  RESO $mapmin $resmax
  GRID $grid
  BINMAPOUT
  LABIN F1=DOBS_${bin_nam} W=FOM_${bin_nam} PHI=PHI
END-wfft
