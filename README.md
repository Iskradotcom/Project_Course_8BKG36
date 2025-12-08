# Project_Course_8BKG36
This repository contains the script &amp; data used during the Project Course in Bioinformatics

The pre-processing of data was done in Galaxy and it included the following steps:
- Firstly, FastQ which is used for quality control and checks whether the data has any initial issues <br/>
- Then, trimming with fastp was performed and this sees what has changed after the filtering/after FastQ, which results in a fastp HTML report <br/>
- Trimming with Cutadapt: removes adapter sequences, primers, polyA + other types of unwanted sequcnes that could cause errors in the downstream analysis (Galaxy, n.d.) <br/>
- MultiQC, which can combine several QC results into a single report <br/>
- Mapping with HISAT2 and STAR = they're both splice-aware aligners, which means that they detect + align the reads (specifically, HISAT2 sees how well the sequenced reads align with the reference genome) <br/>
- And finally, featureCounts: used to quantify the reads, ENTREZIDs are converted for easier biological interpretation <br/>


Steps performed in R: 
- Reading in of the data from Galaxy using featureCounts and annotateMyIDs <br/>
- Creating a metadata file that contains the information about the 2 treatment groups, i.e., Sensitive and Resistant <br/>
- Data normalisation <br/>
- PCA (Principle Component Analysis) <br/>
- Statisitcal tests, e.g., T-test for all genes <br/>
- Adjusting the generated p-values using Bonferroni correction <br/>
- Calculating the log fold change / log2FC <br/>
- Visualisation using different plots, those being the volcano, scatter and heatmap plots <br/>
- Lastly, pathway analysis using kegg 

