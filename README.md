# forestscalingworkflow

This is a repository with all the raw data and code needed to reproduce the analysis in the manuscript "Life history scaling and the division of energy in forests," which will shortly be submitted to *Ecology Letters*. 

John M. Grady, Quentin D. Read, Sydne Record, Nadja Rüger, Phoebe L. Zarnetske, Anthony I. Dell, Stephen P. Hubbell, Sean T. Michaletz, Alexander Shenkin, and Brian J. Enquist. 2020. Life history scaling and the division of energy in forests. 

The workflow is split into 20 numbered R scripts. 

- Scripts 1 through 6 include the initial data processing pipeline in R. 
- After Script 6, the models are fit using CmdStan (outside of R). Stan scripts are found in the `model_scripts` directory in this repository, and shell scripts used for Stan model fitting are in the `shell_scripts` directory. 
- Scripts 7-15 are post-processing R scripts that turn the model fit output into tables of summary statistics and produce CSVs of observed and fitted data that are used for creating plots. 
- Scripts 16-19 are additional analyses that appear in the manuscript. (Note: as of 22 March 2021, the analysis in Script 19 no longer appears anywhere in the manuscript but is retained for completeness.)
- Script 20 generates all the main-text and supplemental figures found in the manuscript. (Note: as of 22 March 2021, there may be a discrepancy between the figures created in this script and those that appear in the manuscript, due to revisions.)

*Note*: If you would like to reproduce only the plots in the manuscript, without having to reproduce the data processing and model fitting steps, run script 20 only. The final outputs of the data processing and model fitting steps are included in this repository to enable easy reproduction of the figures.

Many of the functions called in these scripts are found in an R package called [forestscaling](https://github.com/qdread/forestscaling), created to package some of the code needed for this manuscript.

Please contact the authors of the manuscript or open an issue on this repository if you have trouble reproducing any aspect of the analysis.

*Last updated by Quentin Read on 22 March 2021*
