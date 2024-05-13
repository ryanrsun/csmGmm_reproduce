# csmGmm_reproduce

Code to reproduce simulation results, figures, and real data analysis results from the paper "Testing a Large Number of Composite Null Hypotheses Using Conditionally Symmetric Multidimensional Gaussian Mixtures in Genome-Wide Studies" by Ryan Sun, Zachary McCaw, and Xihong Lin.

## Overview
- All code provided is designed to be run on a high-performance computing cluster. We provide batch submission scripts for clusters using LSF.

- In general, one folder is provided for each figure or each table. For example, all the code needed to reproduce Figure 1 is placed inside the Fig1 folder.

- When a folder does not exist for a figure, it's code can be found in the preceding figure folder. For example, Supplemental Figures 21-22 rely on much of the same code, so there is only a SuppFig21 folder and no SuppFig22 folder.

- Each folder contains _sim.R scripts needed to perform simulations or analysis. At the top of each script, there is a clearly marked box that identifies any code that may need to be changed. For example, the output directory will generally need to be changed to accomodate the varied directory structures of different users. Other than changing these few lines of code, it is not necessary to interact with these scripts.

- Each folder also contains LSF batch submission scripts that end with the .lsf extension. A few lines in these scripts will also generally need to be slightly modified. For example, the job submission queues will generally have different names at different institutions.

- After making any necessary adjustments, all .lsf files in a folder should be run (e.g. with the command "bsub <run_Fig1A.lsf"), and then the user should wait until all jobs from the folder are finished running.

- Each folder also contains an R script that starts with the word "plot" such as "plot_Fig1.R." After all jobs in a folder have finished running, this script should be run to produce the final figure or table from the manuscript.

