# NMR analysis utility functions and scripts

This is a collection of utility functions for analyzing NMR data. I wrote these during my PhD and MSc work in [Carl Michal's lab](https://www.phas.ubc.ca/~michal/) in the Physics department of the University of British Columbia in Vancouver, Canada.

These haven't been put into proper python modules, and I used absolute imports via the imp package:
```python
import imp
read_xnmr_bin = imp.load_source("read_xnmr_bin",<path to read_xnmr_bin.py>)
nmrtools = imp.load_source("nmrtools",<path to nmrtools.py>)
read_xnmr_params = imp.load_source("read_xnmr_params",<path to read_xnmr_params.py>)
```
This worked well for me since I was the only one using these functions. They should be modularized if they are used by multiple people.

## External libraries
All of these require SciPy, NumPy, and Matplotlib. Many of them require LMFIT. These packages can be installed using PIP. Other requirements will be indicated below in the descriptions of each package.

## Coding style
Except where necessary (eg. numerically integrating the Super-Lorentzian lineshape, performing Fourier Transforms), I haven't put very much emphasis on speed or efficiency. My later code is substantially better: more concise, more self-commenting, more vectorizing and use of built-in NumPy function, more use of exceptions, and faster.

Exception handling is ugly and not terribly informative. Sorry.

I performed all of my analysis using IPython, so I frequently import the IPython debugger in the scripts.

## Overview of the code
See the docstrings in the functions for more information.

#### ft_analysis.py
A collection of functions for calculating Fourier transforms. This uses the PyFFTW library, which allows for crazy fast FTs compared to NumPy's implementations. Includes:
  * A function to calculate the entire FID (real and imaginary parts) from a real spectrum. Useful for digital filtering.
  * Functions to calculate forward/reverse FTs, including the frequencies.
  * Functions to do auto-phasing of FIDs prior to taking the FT.

#### nmrtools.py
Despite its name, there are only two functions in here which would be of any use. The rest are used to read/write/analyze old Vnmr data files. However, now that the 400 MHz spectrometer runs VnmrJ, these aren't needed. Besides, this older code was some of the first Python I wrote, so it's likely horrendous anyway.

The useful functions here are:
  * A function which can convert a spectrum in Hz to ppm using an external reference (get_ppm())
  * A function for converting between tensor representations (convert_tensor())

#### plot_tools.py
Don't bother with this. It was some clunky code for plotting grids in Matplotlib. I think the newest Matplotlib has improve functions to do this.

#### read_xnmr_bin.py
Has one function only. It reads the binary Xnmr datafiles into a NumPy array. Can handle one and two dimensional data. Read the docstring before using, there are some quirks with how Xnmr saves raw and processed data.

#### read_xnmr_params.py
Has one function only. Reads Xnmr experiment parameters and returns them as a dictionary.

#### save_xnmr_bin.py
Saves an xnmr file. To be used after read_xnmr_bin() since it requires an existing Xnmr file directory structure. Includes a helper function to change the param file of an Xnmr file directory.

#### reg_nnls.py
Performs regularized NNLS fitting to CPMG decay curves, extracting the T2 distributions. See the docstring for more details.

#### lineshapes.py and lineshapes_c.c

##### lineshapes_c.c
This is a collection of c functions which calculate the Super-Lorentzian lineshapes in the time and frequency domain. Since these lineshapes require integration over many sub-spectra, doing the for-loops in c explcitly is faster than in Python (by my experience, about 2x faster). But how to call these from Python? ...There's a number of [ways to do this](https://www.scipy-lectures.org/advanced/interfacing_with_c/interfacing_with_c.html). I used the [ctypes module](https://www.scipy-lectures.org/advanced/interfacing_with_c/interfacing_with_c.html#id3)

In order to create a .so library file which can be loaded using ctypes, compile this file with
```bash
gcc -fPIC -shared -o lineshapes_cfuncs.so -O3 -ffast-math lineshapes_c.c
```

##### lineshapes.py
lineshapes.py contains a number of functions to calculate both FID and spectral lineshapes (ie lineshapes in the frequency and the time domains). This includes Gaussians, Lorentzians, and Super-Lorentzians. Functions which start with py_ indicate the lineshape is computed using python. Functions that start with c_ indicate the lineshape is calculated using c, relying upon the lineshapes_cfuncs.so library compiled from linshapes_c.c

See the docstrings for more information. The functions which call the c functions are set up to allow multi-threading: the arguments are a list.