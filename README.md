# Forest Scaling workflow

[![DOI](https://zenodo.org/badge/234630889.svg)](https://zenodo.org/badge/latestdoi/234630889)

This is a repository containing all data and code necessary to reproduce the analyses in:

 John M. Grady, Quentin D. Read, Sydne Record, Nadja Rüger, Phoebe L. Zarnetske, Anthony I. Dell, Stephen P. Hubbell, Sean T. Michaletz, Brian J. Enquist. 2023. Life history scaling in a tropical forest. *Journal of Ecology*.

(This citation will be updated upon final acceptance of the manuscript.)

## Contents

- The root folder contains 19 numbered scripts that should be run in order to replicate all analysis presented in the main text and supplements of the manuscript.
- The `data` folder contains all raw data needed.
- The `model_scripts` folder contains the Stan model files.
- The `R_functions` folder contains the source code for the `forestscaling` package, needed to run the analysis, as well as an additional file containing needed functions.
- The `shell_scripts` folder contains files with shell commands needed to fit the Stan models.

## Instructions

- First, install the `forestscaling` R package which includes some functions and plotting themes needed to run the analysis R code.
The source file for the package is included in this repository. Install from source using this R command:

```
install.packages('R_functions/forestscaling_0.1.tar.gz', repos = NULL, type = 'source')
```

- Next, run scripts `01` through `06` to set up all needed files to run the Stan models. *Note that these scripts also process and create some data objects associated with light interception, leaf area, and crown volume. These do not appear in the analyses in the final published version of the manuscript but may be of interest.*
- Next, fit the Stan models using the Bash shell commands in the file `cmdstan_calls.sh` (these were run on a compute cluster with Slurm software).
- Finally, run scripts `07` through `19` to process, analyze, and visualize model output.

*This README last updated by QDR on 27 October 2023*
