# forestscalingworkflow

This is a repository with all the raw data and code needed to reproduce the analysis in the forest light scaling manuscript. 

The pipeline is complex and is split into 14 R scripts. Scripts 1 through 6 include the initial data processing pipeline in R. After Script 6, the models are fit using CmdStan (outside of R). Stan scripts are found in the `model_scripts` directory in this repository, and shell scripts used for Stan model fitting are in the `shell_scripts` directory. Scripts 7-14 are post-processing R scripts that turn the model fit output into tables of summary statistics and produce CSVs of observed and fitted data that are used for creating plots.

Many of the functions called in these scripts are found in an R package called [forestscaling](https://github.com/qdread/forestscaling), created exclusively to package some of the code needed for this manuscript.

Please contact the authors of the manuscript or open an issue on this repository if you have trouble reproducing any aspect of the analysis.

*Last updated by Quentin Read on 17 January 2020*
