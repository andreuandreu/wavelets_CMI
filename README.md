# wavelets_CMI
Wavelet decomposition and Transfer Entropy/ Conditional Mutual Information package for univalued time series


## Table of contents: 

find the corresponding scripts and notebooks to do the analysis for each section in:

- Artificial signals:
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
## Description:

Set of tools in python and julia
it consists of 3 main parts
- generate or read artificial signals (python)
- decompose signals using Continous Complex Wavelet analysis (python)
- Measure correlation and causality using julia (python)

Associated with these parts there are the corresponding plotting and data generating tools


## Getting Started

# the code has four three main parts as of now:

0. create configuration file to define the options to run
1. read and generate data to be transformed into wavelets.
2. use the given data/wavelets to compute pairs of CMI TE using julia package
3. plot the results of the CMI TE computations

# the steps to run the code

1. place the data to be analized in ./data/input
2. got to the script ./scripts/wavelet_warping_calls.py and read the instructions there
   IMPORTANT the most significant things in this script are 
      1 select the 'frequencies'
      2 make sure the 'units' and signal are 

   - 1st provide a univaluated signal, chose one function from the provide_signals script
      1 chose a tag in the dictionary name_source
      2 if no tag exists, create it
      3 make sure the amplitude of the signal is renormalized to values vetween -2 to 2

   - 2nd settle in the units of the signal by choosing a value for sampling_dt
      1 edit the 'unit' variable with a string bearing the physical unit
      2 make sure the sampling_dt corresponds to the unit selected!!!!
         correct the code to include if statements to make sure, for each dataset!

   - 3rd do the fast fourier transfrom of the signal using the FFT(t, sig) function

   - 4th select a scaling or frequency regime to compute the wavelets in these frequencies, these can be made
      1 by hand creating an arrray
      or 2 by selectring a range of lin frequencies by seting :
        step_period
        max_period 
        min_period  
      or 3 by making it automatic using the last mode of the fourier transfer wc.scales_fourier(...)

   - 4th compute the wavelets using one of the two possible methods by selecting the tag in wav_method, or both
      wav_method = 'pywt' -> wc.pywt_compute_wavelets(...)
      wav_method = 'niko' -> wc.niko_compute_wavelets(...)


3. Set the config file at ./confs/config_embeding_char.ini
 there set 
   [folders]
   data_folder = ./data/output/
   input_folder #TO BE CHANGED WITH APPROPRIATE DATA NAME
   export_folder #TO BE CHANGED  WITH APPROPRIATE DATA NAME

   [names]
   in_data_tag #TO BE CHANGED WITH APPROPRIATE DATA NAME
   sufixes = _pha.npy,_amp.npy,_ska.npy
   pha_amp_com = _pha,_pha
   pha_amp_more = _pha,_amp

   [emb_par]
   bins = 150
   max_tau = 10
   jump_tau = 5
   embeding_dimension = 2,1,1,1,1
   period_range = 6,84

   [prob_est]
   name_tag = Prob-est_VisFreq_b
   prob_kind = VisFreq

   [surrogates]
   surr_kind = circular
   surr_num = 11

4.  at ./confs/config_embeding_char.ini  select or unselect the lines 
      bashCommand = "julia --project=. ./scripts/compute_TE.jl ./confs/config_embeding_char.ini"
      process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
      output, error = process.communicate()

   in order to run or not the julia code.


5. at ./scripts/plot_TE-corelation.py plot the apropiate plots nedded.






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

This will install all necessary packages for you to be able to run the scripts and everything should work out of the box.




