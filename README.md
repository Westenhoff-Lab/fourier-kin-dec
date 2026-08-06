# Kinetic Deconvolution in Fourier Space 
This python code performs kinetic deconvolution of time-resolved difference structure factors directly in Fourier (reciprocal) space. Given a kinetic model (species occupancies over time) and a series of experimental difference structure factors $\Delta F$ + uncertainties $\sigma$, the program reconstructs the difference structure factors corresponding to the individual kinetic intermediates while propagating experimental uncertainties.

# Requirements
The script requires Python 3 together with the numpy and pandas packages. The remaining imported modules (argparse, pathlib and collections) are part of the Python standard library when you create a conda environment. We recommend to install CCP4, which is required to generate the required input `.hkl` and `.phs` files through the provided `create_diff_map` pipeline. Once the input files have been prepared, the kinetic deconvolution can be run without any other crystallographic programs. For visualization or conversion of your results, you can use the `f2mtz` and `fft` program in CCP4 or the equivalent tools in phenix. Furthermore, the provided `make_extrapol_map.sh` script also requires CCP4. The CCP4 installation instructions can be found here: https://www.ccp4.ac.uk/download/doc/installation.html


# Input data 
The program requires:
- A **dark-state phase file** (`.hkl`) containing reflection phases.
- A **dark-state amplitude/sigma file** (`.hkl`).
- For each timepoint:
  - A **difference structure factor file** (`.phs`) containing amplitudes and phases.
  - A corresponding **sigma file** (`.hkl`) containing the experimental uncertainties.
- A **concentration matrix** (`.csv`) containing the occupancy of each kinetic species at every timepoint. The names of the timepoints need to match the `.phs` and `.hkl` files. 

The expected column layout for each input file is defined at the beginning of the Python script (`PHS_COLUMNS`, `HKL_COLUMNS`, and `DARK_PHASE_COLUMNS`). 

## Generating the input data via our diff_map scripts 
The recommended way to generate the required `.phs/.hkl` input files is using the pipeline provided in the `create_diff_map` folder. The script `make_dmap.sh` generates weighted difference structure factors and difference maps from raw crystallographic data and exports all input files required for the kinetic deconvolution. Given raw structure factor amplitudes for the light data (`FOBS_<timepoint>.mtz`), raw amplitudes for the reference data (`FOBS_<reference>.mtz`), a refined reference model (`<reference>.pdb`), and the corresponding MTZ file containing reference phases (`<reference>.mtz`), the script runs the full CCP4/Python pipeline.

Before running the script, adjust the input file names, column labels, and crystal information at the beginning of the script using a text editor. The script can then be executed as:
`tcsh make_dmap.sh <timepoint>` 

## Alternative way to generate the input data 
Alternatively, the required `.hkl` and `.phs` files can be written from MTZ files using `mtz2various`. Example command for exporting the **dark-state phase file** are shown below. 

`mtz2various HKLIN input.mtz HKLOUT dark_phase.hkl << EOF \\
LABIN FP=<Amplitude_column> SIGFP=<Amplitude_sigma_column> PHIC=<Phase_column> \\
OUTPUT USER '(3I5,3F12.3)' \\
RESOLUTION <min_resolution> <max_resolution> \\
EOF`

The same approach can be used for the dark-state amplitude/sigma file, as well as the difference structure factor file with the corresponding sigma file by adjusting the `LABIN` labels accordingly.

# Usage
Before running the script, make sure that you have the corresponding `.phs` and `.hkl` file in two directories (in this case it is phs and hkl) and make sure that the labels in the **concentration matrix** `.csv` file, has the same naming as the files. 

The script can be run with the following commands:

`python kinetic_deconvolution.py --concentrations concentrations.csv --phs-dir phs --hkl-dir hkl --dark-sigma dark_scaled.hkl --dark-phase dark_phase.hkl --out-prefix state`

The script first loads and displays the concentration matrix, including the number and names of the structural intermediates. It then performs phase correction and amplitude sign adjustment of the input difference structure factors.

The script filters reflections based on their presence across the dataset. By default, reflections present in at least 12 of the 17 timepoints are retained.

The resulting kinetic modes are saved as `.phs` files containing the deconvoluted difference structure factor amplitudes and phases, together with corresponding `.hkl` files containing the propagated uncertainties. These files can be used for further structure factor extrapolation (see below `How to continue from here`).

The generated `.phs` files can be converted to `.mtz` format using `f2mtz` from CCP4. Electron density maps can then be calculated using `fft` in CCP4 for visualization. We provide a `phs_to_map.sh` CCP4 script, which can be used to convert the final `.phs` files to DED maps. It can be run in the terminal by simply doing:
`tcsh phs_to_map.sh <state_name>` 
You can then take a look at the deconvoluted maps in coot using the reference `.pdb`. 

# Simulated data - Example 
As described in the paper, difference structure factors were simulated for 17 timepoints containing different concentrations of four structural intermediates of the photoactive yellow protein (PYP). For the provided simulated dataset, no reflections require phase correction because the phases and amplitudes were already generated consistently. Similarly, all 9786 reflections pass the filtering criterion.

The provided example dataset contains:
- A dark-state phase file (`.hkl`)
- A dark-state amplitude/sigma file (`.hkl`)
- Difference structure factor files (`.phs`) for all timepoints with corresponding sigma files (`.hkl`). These files are in the folders `phs` and `hkl`. 
- The concentration matrix describing the population of each structural intermediate at each timepoint

These files can be used to reproduce the example analysis described in the paper. They also serve as a reference for the required input format, including the expected column organization and data labels.

# How to continue from here 
The deconvoluted weighted difference structure factors (`.phs` file), together with the propagated sigmas (`.hkl` file) can now be used to calculate extrapolated structure factors. For that, we provide a CCP4/python script, which takes the generated files together with the calculated structure factors fo the reference state (`dark_phase.hkl` file) and calculate **extrapolated structure factors** by
F<sub>ext</sub> = F<sub>calc</sub> + N × (ω × ΔF<sub>obs</sub>), where N is an extrapolation factor which corresponds to the photoactivation yield 2/occupancy. Several methods have been proposed to estimate N, and users should determine it using their method of choice.

Before running the script, adjust the input file names and crystal information at the beginning of the script using a text editor. The script can then be executed as:
`tcsh make_extrapol_map.sh <extrapolation_factor> <timepoint/state_name>` 

From the resulting extrapolated `.mtz` file, structures can be refined using either CCP4 or phenix. Important to note here is, that the refined structures will most likely result in high R<sub>work</sub> and R<sub>free</sub> values due to the phase error introduced when adding difference amplitudes to the dark state structure factors. The true light phase can be estimated, and a script for generating **phased extrapolated structure factors** will be made available soon. 

# Conventional workflow 
The `create_diff_map` and `create_extrapol_map` workflow can also be used for timepoints without using the deconvolution approach. This workflow is based on Schmidt (2023, Structural Dynamics). 

# Citation 
When using the provided scripts, please cite as follows:

- Grunewald, L. et al., Visualizing Reaction Pathways via Reciprocal Space Kinetic Decomposition, submitted (2026).

Please also cite the software and workflow on which these scripts are based:
- Agirre, J. et al., The CCP4 suite: integrative software for macromolecular crystallography, Acta Crystallogr. D. Biol. Crystallogr. 67, 235–242 (2011).
- Schmidt, M., Practical considerations for the analysis of time-resolved x-ray data, Struct. Dyn. 10, 044303 (2023).
