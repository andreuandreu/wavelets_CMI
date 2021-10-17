# wavelets_CMI
Wavelet decomposition and Transfer Entropy/ Conditional Mutual Information package for univalued time series


## Table of contents: 
## Description:

Set of tools in python and julia
it consists of 3 main parts
- generate or read artificial signals (python)
- decompose signals using Continous Complex Wavelet analysis (python)
- Measure correlation and causality using julia (python)

Associated with these parts there are the corresponding plotting and data generating tools


## Getting Started

find the corresponding scripts and notebooks to do the analysis for each section in:

- Artificial sygnals:
    scripts/create_synthetic_signals.py -> in development
    scripts/read_surrogate_data.py
    scripts/rossler_tools_wavelets.py
    scripts/generate_colored_noises.py
- Wavelet decomposition
   scripts/wavelet_analysis.py 
   scripts/wavelets_experiment.py -> in development
- CMI/TE
   scripts/TE_functions.jl
   notebooks/rossler_TE-checkpoint.ipynb -> in development




### Prerequisites

Python 3, Julia 1.5 or above

### Installing

For the julia part of this project

0. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box.




