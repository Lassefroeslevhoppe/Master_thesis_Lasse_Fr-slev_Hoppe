Pipeline for my analysis of data from core M0058A

Author: Lasse Frøslev Hoppe

Affiliation: Master students, Section for GeoGenetics, Globe Institute, University of Copenhagen
___________________________________________________________________________________________________
Overview

This study analysed 20 samples from sediment core m0058A collected at the IODP exp. 325.
The scripts provided here show how I edited, extracted and analysed the data from LCA-files.

The processing of fastq-files was conducted through a pipeline created by the Center of Geogenetics that is not yet published. I edited the pipeline to fit my data, but it was run by Dr Sara Behnamian. I cannot share the script here - therefore, see the master thesis itself for a description of the workflow. Further, the fastq and lca files are not shared as they have not been published either.

Filtering and editing of the LCA files are provided in scripts 1-4 in the depository here on GitHub. 
Data analysis is provided in scripts 5-11 and all scripts are numerized to show the order of which they are supposed to run.
____________________________________________________________________________________________________
Requirements
Python v 3.7.10

AdapterRemoval v 2.3.3

fastp v 0.24.0 

SGA v 0.10.15 

BBMAP v 39.10 

Bowtie2 v 2.4.2

metaDMG v 0.4.143 

Pandas v 2.2.1 

Matplotlib v 3.5.3 

NumPy v 2.2.0 

statsmodels v 0.12.2 

scikit-learn v 1.6.1 
___________________________________________________________________________________________________
Environment for running STARS (time series analysis)

Because I was not allowed to install any packages on the server Dandycomp, I created an environment to allow me to install the R_shift package. I chose to make a conda-environment as I was showed how by one of the people I was sharing a workspace with when I first started 
learning about coding: "conda create STARS_env". To this environment, I loaded Pandas, Nympy, subprocess, os and tempfile. This was done to handle the r-script that I included in the scripts for detecting Regime Shifts. 
Most inportantly for this script to run, the R_shift package must be installed in the envrionenment.  
___________________________________________________________________________________________________
Pipeline (the scripts must be run in order)

1_Appending_depth_and_age-py --> Adding columns with information of the sample age and the position where the sample was collected from the core. 

2_Highest_tax_level_to_excel.py --> The highest taxonomic level for each read was extracted and pooled (counting the reads) and put in an excel file to get an overview of the data.

3_counting_n_reads_per_fastq_file.py --> Because the fastq-files were so big, I counted the number of reads for each and put them in an excel file.

4_Plotting_percentage_of_fast_to_lca_superkingdoms_commented.py --> Here I counted the number of reads per LCA-file and used the number of reads from the excel created in script 3 to plot the difference in reads before and after filtering. 

5_Passing_QC_finalized.py --> This script plots the number of genera passing the multistep quality control I set up

6_Panel_plots.py --> This script is used by running it in python and typing in the genus to create panel plots for: eg. "6_panel_plots.py Porites". These plots were used to authenticate the genera reads as truely ancient. 

7_Heatmap_all_genera.py --> This script gathered all the reads for genera in the whitelist (excel file created manually), applied a CLR data transformation and plotted the genera present for each sample.

8_Heatmap_marine_genera.py --> Same as above, but for a restricted set of genera.

9_PCA_all_genera.py --> Created PCA biplots for both SS and DS data. 

10_PCA_finalized_marine_genera.py --> Same as above, but for the resticted set of marine genera. 

11_Time_series_analysis.py --> Created a time series analysis with STARS using the R_shift packgage and the envrionment described above. 

12_Time_sereies_analysis_marine_genera.py --> Same as above but only for the set of marine genera. 
____________________________________________________________________________________________________________________________
Additional files

In the repository, I have also included the output excel files (such as the n_reads of fastq files). Also, the whitelists for all genera passing damage QC and the marine genera are provided as excel files.  
____________________________________________________________________________________________________________________________
Final remarks

The scripts here are provided for my supervisor, co-supervisor and examiner/censor to show how I processed and analysed the data in my master thesis. This is meant as a supplement and 
and not as a stand-alone part of the project. For a more detailed description see the Master Thesis' methodology section. 
