# fourier-kin-dec
fourier-kin-dec performs kinetic deconvolution of time-resolved difference structure factors directly in Fourier (reciprocal) space. Given a kinetic model (species occupancies over time) and a series of experimental difference structure factors + uncertainties, the program reconstructs the difference structure factors corresponding to the individual kinetic intermediates while propagating experimental uncertainties.

# Requirements
The script requires Python 3 together with the numpy and pandas packages. The remaining imported modules (argparse, pathlib and collections) are part of the Python standard library. CCP4 is only required to generate the input `.hkl` and `.phs` files from MTZ files (e.g. using `mtz2various`). Once the input files have been prepared, the kinetic deconvolution can be run without any CCP4 programs. If you want to visualize your results, you will need the `f2mtz` and `fft` program in CCP4. The CCP4 installation instructions can be found here: https://www.ccp4.ac.uk/download/doc/installation.html

# Input data 
The program requires:
- A **dark-state phase file** (`.hkl`) containing reflection phases.
- A **dark-state amplitude/sigma file** (`.hkl`).
- For each timepoint:
  - A **difference structure factor file** (`.phs`) containing amplitudes and phases.
  - A corresponding **sigma file** (`.hkl`) containing the experimental uncertainties.
- A **concentration matrix** (`.csv`) containing the occupancy of each kinetic species at every timepoint. The names of the timepoints need to match the `.phs` and `.hkl` files. 

The expected column layout for each input file is defined at the beginning of the Python script (`PHS_COLUMNS`, `HKL_COLUMNS`, and `DARK_PHASE_COLUMNS`). 
The required `.hkl` and `.phs` files can be written from MTZ files using `mtz2various`. Example command for exporting the **dark-state phase file** are shown below. 

`mtz2various HKLIN input.mtz HKLOUT dark_phase.hkl << EOF \\
LABIN FP=<Amplitude_column> SIGFP=<Amplitude_sigma_column> PHIC=<Phase_column> \\
OUTPUT USER '(3I5,3F12.3)' \\
RESOLUTION <min_resolution> <max_resolution> \\
EOF`

For the **dark-state amplitude/sigma file**, as well as the **difference structure factor file** with the corresponding **sigma file**, they can be exported in the same way with adjusting the LABIN lables.

# Simulated data - Example 
As described in the paper, difference structure factors were simulated for 17 timepoints containing different concentrations of four structural intermediates of the photoactive yellow protein (PYP). 

The provided example dataset contains:
- A dark-state phase file (`.hkl`)
- A dark-state amplitude/sigma file (`.hkl`)
- Difference structure factor files (`.phs`) for all timepoints with corresponding sigma files (`.hkl`). These files are in the folders `phs` and `hkl`. 
- The concentration matrix describing the population of each structural intermediate at each timepoint

These files can be used to reproduce the example analysis described in the paper. They also serve as a reference for the required input format, including the expected column organization and data labels.

# Usage
The script can be run on the simulated data with the following commands:

`python kinetic_deconvolution.py --concentrations concentrations.csv --phs-dir phs --hkl-dir hkl --dark-sigma dark_scaled.hkl --dark-phase dark_phase.hkl --out-prefix state`

The script first loads and displays the concentration matrix, including the number and names of the structural intermediates. It then performs phase correction and amplitude sign adjustment of the input difference structure factors. For the provided simulated dataset, no reflections require phase correction because the phases and amplitudes were already generated consistently.

The script filters reflections based on their presence across the dataset. By default, reflections present in at least 12 of the 17 timepoints are retained. For the ideal simulated dataset, all 9786 reflections pass this criterion.

The resulting kinetic modes are saved as `.phs` files containing the reconstructed difference structure factor amplitudes and phases, together with corresponding `.hkl` files containing the propagated uncertainties. These files can be used for further structure factor extrapolation.

The generated `.phs` files can be converted to `.mtz` format using `f2mtz` from CCP4. Electron density maps can then be calculated using `fft` in CCP4 for visualization. We provide a `phs_to_map.sh` CCP4 script, which can be used to convert the final `.phs` files to DED maps. It can be run by simply doing:
`bash phs_to_map.sh state_1_pr0` 
