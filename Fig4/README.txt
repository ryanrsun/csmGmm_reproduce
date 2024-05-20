Here we detail the steps for reproducing Figure 4 and Tables 1 and 2:

1. Open each of the .R scripts and change the top few lines where necessary, such as to point to the correct output folder.

2. (optional) Run the script calculate_eqtl_stats_synthetic.R. We are not able to provide the original raw GTEx data for this step due to privacy concerns. There is a synthetic substitute synthDat.txt in the Data/ folder, but it will not reproduce the original findings. We have also provided the original output of this step in a folder name eqtl_data/. The location of eqtl_data/ can be found in the Data/ folder. Using that output in the next steps will reproduce our work.

3. Make necessary changes to .lsf files, such as using the correct queue names and runtimes.

4. Run the run_analyze_med.lsf and run_pleiotropy_analysis.lsf jobs.

5. After step 4 is completed, run the run_summarize_S1.lsf and run_summarize_S2.lsf jobs.

6. Run the plot_data_analysis.R script to plot Figure 4 and create Tables 1 and 2.

 






