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

Explanation of graphs:<br/>
Graph 1 is a PCA plot and it shows the clear distance/difference between the sensitive and resistant groups, with PC1 (91.36%) being able to explain a high amount of variance in the data. <br/>
Graph 2 is a volcano plot and this plot has thresholds in place for the data - a y-intercept for the -log10(adjusted p-values) and a threshold for the fold change that between -1 and +1. Furthermore, the plot shows the up and down-regulated significant genes and also genes that aren't signficant. The type of plot is useful to have as it helps with visualising very large datasets easily so they're easy to view and the data points take a certain shape that aids in interpretation.<br/>
Graph 3a and 3b are heatmaps. Specifically, graph 3a is filtered by the top 20 lowest p-values and graph 3b by the top 20 genes with the highest FC values. From looking at the graphs, I could spot, for example, the gene CALB1 having one of the lowest p-values & also importantly a high logFC value. It has a different shade intensity between the sensitive and resistant samples, meaning that since calcium-binding protein is involved in signal pathways, this could mean that signalling pathways are affected differently in the sensitive and resistant grps (possibly the cancer cells have less signalling to their surrounding environment, to ensure their survival/so they're not detected by cancer-killing cells). <br/>
Graph 4 is the scatter plot and within it, we are able to see a clear divide between the upregulated genes (to the left, above the line) and the down regulated genes (more to the right, under the line). A scatterplot is a useful plot to have because one is able to see the distribution of the data points easily & the legend gives additional information about the data, using colours (similar to the heatmap) and shades. 
