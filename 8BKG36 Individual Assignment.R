
### DAY 1: 12.11.2025
###Individual assignment - Name: Iskra Gutesa -  Course code: 8BKG36

rm(list = ls())
setwd("C:/Users/Iskra/Documents/Bioinformatics Project Course")


library(magrittr)
library(stringr)
library(purrr)
library(dplyr)

df_files <- list.files(pattern = ".tabular")

names_input <- df_files %>%
  str_remove(".tabular")

for(i in 1:length(names_input)) {
  df12 <- read.delim(df_files[i])
  assign(names_input[i], df12)
}

# Alternative 2 (more efficient)
count_table <- list(SRR10416711, SRR10416712, SRR10416713, SRR10416714, SRR10416715, SRR10416716) %>%
  reduce(full_join, by = "Geneid")

------### DAY 2: 13.11.2025 ---------------------

# Requires: 1) a count table with 6 of the samples as well as a geneid column. 2) a dataframe containing the entrez ids and the symbol names (name geneid)

GeneNames=c(GeneID$SYMBOL)
count_table$GeneNames <- GeneNames

counts_updated <- merge(count_table, GeneID, by.x = "Geneid", by.y = "ENTREZID")
counts_updated$Geneid <- NULL 
 
# check if any nas are present in the table
table(is.na(counts_updated))

counts_updated <- na.omit(counts_updated) # remove rows with NAs

rownames(counts_updated) <- counts_updated$SYMBOL # put the symbol as rownames
counts_updated$SYMBOL <- NULL # remove the symbol column from the dataframe
counts_updated$GeneNames <- NULL


---# Metadata file:---------

metadata <- data.frame(sample_id = colnames(counts_updated))
#metadata <- head(metadata, -1)
metadata$treatment <- c("Sensitive", "Sensitive", "Sensitive", "Resistant", "Resistant", "Resistant")
rownames(metadata)<-metadata$sample_id
metadata$sample_id <- NULL


---## Normalisation of data ------

sequencing_depth <- colSums(counts_updated) 
counts_updated_CPM <- counts_updated 
for (n in 1:ncol(counts_updated)){ 
  counts_updated_CPM [,n] <-  
    1e6*counts_updated[,n]/sequencing_depth[n] 
}

counts_updated_logCPM <- log2(counts_updated_CPM + 1)



# needed to find a way to remove rows that are all zero. 
counts_updated_logCPM<- counts_updated_logCPM[rowSums(counts_updated_logCPM != 0) > 0, ]



------### Carry out a PCA -------------

library(ggfortify)
library(stats)
library(ggplot2)


#***PCA Graph (Graph 1)

count_pca <- prcomp(t(counts_updated_logCPM)) # t() = transpose
autoplot(count_pca, data=metadata, colour="treatment", size=4) +
  theme_classic()






---------### Statistical tests : t-tests or linear models------------------------------------------------

pval_counts_updated <- t.test(counts_updated_logCPM["AARS1", 1:3], 
                    counts_updated_logCPM["AARS1", 4:6], alternative= "two.sided") 
pval_counts_updated$p.value

# For example, here I am comparing the expression of the gene "AARS1" between the 6 SRRs from the 2 treatment groups
# and it presents a p-value of 0.0001867305, which could convey statistical significance/the observed difference
# between 2 treatment grps would be unlikely due to chance. 


p_value <- data.frame(p_value=rep(0,nrow(counts_updated_logCPM))) 
rownames(p_value) <- rownames(counts_updated_logCPM)  
for (i in 1:nrow(counts_updated_logCPM)){ 
  p_value$p_value[i] <- t.test(counts_updated_logCPM[i,1:3], 
                               counts_updated_logCPM[i,4:6], alternative = "two.sided")$p.value 
} 

# Now, a data frame called "p_value" has been created and it's a table that contains the varying p values of each of the 6 samples. 


count_lm <- counts_updated_logCPM %>% as.matrix() 
lm_results <- data.frame(GeneID=rownames(count_lm), p_value=NA) 
  sampleID.resistant <- rownames(metadata)[metadata$treatment=="Resistant"]
sampleID.sensitive <- rownames(metadata)[metadata$treatment=="Sensitive"]

for (i in 1:nrow(counts_updated_logCPM)){ 
  df <- data.frame(expr=count_lm[i,], 
                         sampleID=colnames(count_lm)) %>% na.omit() 
  if (nrow(metadata)>0){ 
    df$treatment <- NA 
    df$treatment[df$sampleID %in% sampleID.sensitive] <- "Sensitive" 
    df$treatment[df$sampleID %in% sampleID.resistant] <- "Resistant" 
    lmFit <- lm(expr ~ treatment, df) %>% 
      summary() %>% 
      coefficients() %>% 
      as.data.frame()
    lm_results[i,"p_value"] <- lmFit["treatmentSensitive","Pr(>|t|)"] 
  } 
}    
p_value <- lm_results %>% arrange(p_value) 
rownames(p_value) <- p_value$GeneID 
p_value$GeneID <- NULL 

## 

--------------### Adjusting the p-values: ---------------------------------------------------------

#sig_pval$ajusted_pval <- p.adjust(sig_pval$p_value, method = "BH")

# I used BH instead of Bonferroni, as BH is more commonly used for these type of analyses and Bonferroni can be very 
# conservative and problematic with high-dimensional data analysis. 


sig_pval <- p_value[p_value$p_value<0.05,,drop=FALSE]
# here, I want to just have the 'significant' genes, i.e., those that have a p-value below 0.05.

names_sig_genes <- rownames(sig_pval)
sig_genes <- 
  subset(counts_updated_logCPM,rownames(counts_updated_logCPM) %in% 
           names_sig_genes)

#* **changing data frames**

# I got hints from Microsoft CoPilot on how to do this. 
p_value$adjusted_pval <- p.adjust(p_value$p_value, method = "BH")




# with the sig_genes name, now I have a data frame only with the significant genes. 

-------## Calculating logFC values ---------------

FC_resistant <- rowMeans(counts_updated_logCPM[,4:6])
FC_sensitive <- rowMeans(counts_updated_logCPM[,1:3])
FC <- as.data.frame(FC_sensitive-FC_resistant)
colnames(FC) <- "FoldChange" 




save(df, file= "counts_updated_logCPM_IskraG.RData")



#•••••••••••• Visualisation ••••••••••••••••••••••••••••••••••••••••••••••••••

#** VOLCANO PLOT (Graph 2) : 


p_value$log2FC <- FC$FoldChange

# Standard volcano plot
p <- p_value %>%
  ggplot(mapping = aes(x=log2FC, y=-log10(adjusted_pval)))+
  geom_point()+
  theme_minimal()

p2 <- p +
  geom_hline(yintercept = -log10(0.05), col="red")

# Color the dots based on their differential expression
# Let's put a threshold on the Fold Change so that we only consider genes
# that have a significant pvalue and are above a certain threshold when
# it comes to the log2 fold change

p_value$diffexpressed[p_value$adjusted_pval > 0.05] <- "Not_significant"
p_value$diffexpressed[p_value$log2FC > -1 & p_value$log2FC < 1  & p_value$adjusted_pval < 0.05] <- "NO" 
p_value$diffexpressed[p_value$log2FC > 1 & p_value$adjusted_pval < 0.05] <- "UP"
p_value$diffexpressed[p_value$log2FC < -1 & p_value$adjusted_pval < 0.05] <- "DOWN"

p <- p_value %>%
  ggplot(mapping=aes(x=log2FC, y=-log10(adjusted_pval), col=diffexpressed))+
  geom_point()+
  theme_minimal()

p2 <- p +
  geom_hline(yintercept = -log10(0.05), col="red")+
  geom_vline(xintercept= c(-1,1), col="red") # add lines for log2FC -1 and 1

mycolors <- c("blue", "green3", "orchid4", "tomato3")
names(mycolors) <- c("DOWN", "UP", "NO")
p3 <- p2+
  scale_color_manual(values=mycolors)



###################################################################################################
#** HEATMAP (Graph 3a):

#_# log2FC

# Because there are so many genes, let's look at the top 20 with the highest log2FC

significant_DEGs_log2FC <- p_value[
  order(p_value$log2FC, decreasing = TRUE),
  ,
  drop = FALSE
]

significant_DEGs_subset_log2FC <- significant_DEGs_log2FC[1:20, ]
# ordering the top 20 w/ biggest log FC, I used the help of Chat-GPT here as it advised me 
# that my code before could be error prone and is messy, due to the manual selection of 
# the rows that I was doing before (I was subtracting from the total no. of entries). 

names_sig_subset_log2FC <- rownames(significant_DEGs_subset_log2FC)

sig_DEGs1 <- counts_updated_logCPM %>%
  subset(rownames(counts_updated_logCPM) %in% names_sig_subset_log2FC) %>%
  as.data.frame() # only includes genes of interest

sapply(sig_DEGs1, class) # check class of the variables in the data frame
matrix_sig_DEGs1 <- as.matrix(sig_DEGs1)  
heatmap(matrix_sig_DEGs1, scale="row", Colv=T, Rowv = F)# uses the heatmap function in base R, scale is used to normalize between the columns
#pheatmap(matrix_sig_DEGs1,scale = "row")

coul <- colorRampPalette(brewer.pal(11, "RdBu"))(256)
heatmap(matrix_sig_DEGs1, scale="row", Colv = T, Rowv = T, col= coul)




#pheatmap

pheatmap(matrix_sig_DEGs1) # this codes for the basic heatmap

pheatmap(matrix_sig_DEGs1, scale= "column", color = colorRampPalette(c("navy", "white", "red"))(50),
         cluster_rows = T, cluster_cols = T)

pheatmap(matrix_sig_DEGs1, scale= "column", color = colorRampPalette(c("navy", "white", "red"))(50),
         cluster_rows = T, cluster_cols = T)













#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

## p-value -> with the top 20 lowest p-values (Graph 3b)

significant_DEGs <- p_value[order(p_value$adjusted_pval),,drop=FALSE]

significant_DEGs_subset <- significant_DEGs[c(1:20),] # extracting those with lowest p-values


names_sig_subset <- rownames(significant_DEGs_subset) # extracting vector of the names of the significant genes



# Subset the gene expression matrix so that it only includes the genes of interest
sig_DEGs <- counts_updated_logCPM %>%
  subset(rownames(counts_updated_logCPM) %in% names_sig_subset) %>%
  as.data.frame()

sapply(sig_DEGs, class) # checking the class of the variables in the data frame, which is numeric
matrix_sig_DEGs <- as.matrix(sig_DEGs)
heatmap(matrix_sig_DEGs, scale="row", Colv=T, Rowv = F)# uses the heatmap function in base R, scale is used to normalize between the columns
pheatmap(matrix_sig_DEGs,scale = "row")
install.packages("RColorBrewer")
library(RColorBrewer)

coul <- colorRampPalette(brewer.pal(11, "RdBu"))(256)
heatmap(matrix_sig_DEGs, scale="column", Colv = T, Rowv = T, col= coul)

pheatmap(matrix_sig_DEGs) # basic heatmap

pheatmap(matrix_sig_DEGs, scale= "column", color = colorRampPalette(c("navy", "white", "red"))(50),
         cluster_rows = T, cluster_cols = T)

pheatmap(matrix_sig_DEGs, scale= "column", color = colorRampPalette(c("navy", "white", "red"))(50),
         cluster_rows = T, cluster_cols = T)

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# heatmap with logFC expression & the top 10 most upregulated

# upregulated_logFC <-p_value[order(p_value)]









#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#** SCATTERPLOT (Graph 4): 

library(ggplot2)
library(ggrepel)
library(pheatmap)

p_value # is my data frame containing all adjusted p-values and the logFC for all genes

sig_pvalue <- p_value[p_value$adjusted_pval<0.05,,drop=FALSE]

sig_pvalue <- sig_pvalue[order(sig_pvalue$log2FC),,drop=FALSE] # the order is based on logFC

top_genes_sp <- sig_pvalue[c(1:10, 7499:7508),] # here the top down-regulated and top-upregulated are extracted into the data frame called top_genes_sp

top_genes_sp_names <- rownames(top_genes_sp)

top_genes_expression <- counts_updated_logCPM[(rownames(counts_updated_logCPM) %in% top_genes_sp_names),]

top_genes_expression$SensitiveMean<-rowMeans(top_genes_expression[,1:3]) # adding in sensitive mean 
top_genes_expression$ResistantMean<-rowMeans(top_genes_expression[,4:6]) # adding in resistant mean

head(top_genes_expression)


top_genes_expression$genename<-rownames(top_genes_expression) # add a column with genename to use as labels for the plot


ggplot(top_genes_expression, aes(x=SensitiveMean, y=ResistantMean, label=genename))+
  geom_point(aes(color=ResistantMean))+
  geom_abline(intercept = 0)+
  theme_bw()+
  ggtitle(label="Top 10 up and downregulated genes")+
  scale_color_gradient2(low = "darkblue", mid = "white",
                        high = "red", space = "Lab" )+
  geom_text_repel()+
  theme(plot.title = element_text(hjust = 0.5, size=20))




##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Pathway analysis:

significant_DEGs$gene_symbol <- rownames(significant_DEGs)
background <- rownames(counts_updated_logCPM) # a background of all genes measured needs to be used when running the pathway analysis
significant_DEGs <- significant_DEGs[significant_DEGs$adjusted_pval < 0.05 & abs(significant_DEGs$log2FC)>=1, ]   #since I had some adjusted p-values that were over 0.05, here I am filtering them

install.packages("clusterProfiler")
library(clusterProfiler)

DEG_entrez<- bitr(significant_DEGs$gene_symbol, 
                  fromType = "SYMBOL",
                  toType = "ENTREZID", 
                  OrgDb = "org.Hs.eg.db") # adds a column of entrezID to the data frame

background <-  bitr(background, 
                    fromType = "SYMBOL",
                    toType = "ENTREZID", 
                    OrgDb = "org.Hs.eg.db")

ekegg <- enrichKEGG(
  as.character(DEG_entrez[,2]),
  organism = "hsa",
  keyType = "kegg",
  universe = background$ENTREZID,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH")

dotplot(ekegg, showCategory=20, x="count") + ggtitle("Pathway analysis")

# here only a singular data point is present, corresponings to cytokine-cytokine receptor interaction, and 
# according to a scientific article on Frontiers, gemcitabine can "influence the anti-tumor immune 
#responses via several mechanisms, such as... cytokine production", so this could provide some
# explanation. 