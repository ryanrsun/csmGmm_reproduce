For Figure 4 and Tables 1 and 2, the order of operations is:


1. Make the necessary changes at the top of each .R script, for example to change the directory where the summary statistics are saved.

2. (Optional) Run the calculate_eqtl_stats_synthetic.R script. This step is optional because it requires the individual level GTEx data which we cannot share. We have provided synthetic data, but using the synthetic data will not reproduce the original results. We have also provided the original output of this step, which can be shared, in the Data folder. This original output can be used in the following steps to reproduce our findings.

3. Make the necessary changes to each batch submission script (e.g. queue names or runtimes). 

4. Run analyze_med.lsf and run_pleiotropy_analysis.lsf first and wait for the jobs to finish.

5. Run run_summarize_S1.lsf and run_summarize_S2.lsf after Step 4 has completed.

6. Run the plot_data_analysis.R script to reproduce Figure 4 and Tables 1 and 2.  




