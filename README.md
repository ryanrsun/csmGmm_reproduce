# csmGmm_reproduce

Code to reproduce simulation results, figures, and real data analysis results from the paper "Testing a Large Number of Composite Null Hypotheses Using Conditionally Symmetric Multidimensional Gaussian Mixtures in Genome-Wide Studies" by Ryan Sun, Zachary McCaw, and Xihong Lin.

## Overview

- Please first install the csmGmm package, for example through the command devtools::install_github("ryanrsun/csmGmm").

- Please navigate to the Data/ folder and either run the download_data.sh script or download the other necessary files manually using the other_data_locations.txt file. These other files are too large to be bundled with the rest of the code, so we keep them separately in Dropbox.

- All code provided is designed to be run on a high-performance computing cluster. We provide batch submission scripts for clusters using LSF.

- In general, one folder is provided for each figure or each table. For example, all the code needed to reproduce Figure 1 is placed inside the Fig1 folder.

- When a folder does not exist for a figure, its code can be found in the preceding figure folder. For example, Supplemental Figures 21-22 rely on much of the same code, so there is only a SuppFig21 folder and no SuppFig22 folder.

- We use the 'here' package to refer to file locations within this project. When running this code, please make sure the working directory is set to somewhere within the csmGmm_reproduce folder to ensure that the here package works correctly.

- Each folder also contains example LSF batch submission scripts that end with the .lsf extension. These submission scripts are used to run the R scripts in parallel on a computing cluster, e.g. to run 200 simulations at once for Figure 1A. A few lines in these scripts will generally need to be slightly modified. For example, the job submission queues will have different names at different institutions. The number of jobs to run will be placed in a comment, e.g. "#Run 1-800" means that an array of jobs with IDs from 1 to 800 should be run.

- After making any necessary adjustments, all .lsf files in a folder should be run (e.g. with the command "bsub <run_Fig1A.lsf"), and then the user should wait until all jobs from the folder are finished running.

- Each folder also contains an R script that starts with the word "plot" such as "plot_Fig1.R." After all jobs in a folder have finished running, this script should be run to produce the final figure or table from the manuscript.

- An upper bound on the amount of time needed to execute the code in each folder is roughly equivalent to the largest -W argument given in the lsf files of each folder. For example, -W 3:00 indicates the code should take less than 3 hours. The actual time needed will depend on many factors such as the hardware of the computing cluster. If there is a waiting time for jobs to be run, that time will also need to be added on top of the -W time. We also provide an upper bound on time estimates here for convenience:

Fig1 - 3 hours

Fig2 - 3 hours

Fig3 - 7 hours

Fig4 - 4 hours

Tab1/Tab2 - 5 hours

SuppFig1 - 3 hours

SuppFig2 - 3 hours
SuppFig3 - 7 hours
SuppFig4 - 7 hours
SuppFig5 - 4 hours
SuppFig6-SuppFig35 - 3 hours
SuppTab1 - Instant
SuppTab2 - 23 hours
SuppTab3 - Instant
SuppTab4 - Instant
