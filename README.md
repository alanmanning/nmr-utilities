# NMR analysis utility functions and scripts

This is a collection of utility functions that I have written over the years. The coding style for some of the earlier scripts is pretty rough, but everything works.

These haven't been put into proper python modules, and I used absolute imports via the imp package:
```python
import imp
read_xnmr_bin = imp.load_source("read_xnmr_bin",drop_dir+"nmr-analysis/"+"read_xnmr_bin.py")
nmrtools = imp.load_source("nmrtools",drop_dir+"nmr-analysis/"+"nmrtools.py")
read_xnmr_params = imp.load_source("read_xnmr_params",drop_dir+"nmr-analysis/read_xnmr_params.py")
```
This worked well for me since I was the only one using these functions. They should be modularized if they are used by multiple people.

## Description of the utilities

### 