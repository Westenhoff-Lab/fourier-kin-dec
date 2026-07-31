# fourier-kin-dec
fourier-kin-dec performs kinetic deconvolution of time-resolved difference structure factors directly in Fourier (reciprocal) space. Given a kinetic model (species occupancies over time) and a series of experimental difference structure factors, the program reconstructs the difference structure factors corresponding to the individual kinetic intermediates while propagating experimental uncertainties.

# Requirements
The script requires Python 3 together with the numpy and pandas packages. The remaining imported modules (argparse, pathlib and collections) are part of the Python standard library. CCP4 is only required to generate the input .hkl and .phs files from MTZ files (e.g. using mtz2various). Once the input files have been prepared, the kinetic deconvolution can be run without any CCP4 programs. 
