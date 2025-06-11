#Packages for the installation:

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2", dependencies = TRUE, force = TRUE)
BiocManager::install("ComplexHeatmap", dependencies = TRUE, force = TRUE)
BiocManager::install("clusterProfiler", dependencies = TRUE, force = TRUE)
library("BiocManager")
library("DESeq2")

if (!require("ggplot2", quietly = TRUE))
  install.packages("ggplot2")

if (!require("ggfortify", quietly = TRUE))
  install.packages("ggfortify")

if (!require("gridExtra", quietly = TRUE))
  install.packages("gridExtra")

if (!require("carData", quietly = TRUE))
  install.packages("carData")

if (!require("car", quietly = TRUE))
  install.packages("car")

if (!require("factoextra", quietly = TRUE))
  install.packages("factoextra")

if (!require("pcaMethods", quietly = TRUE))
  install.packages("pcaMethods")

if (!require("FactoMineR", quietly = TRUE))
  install.packages("FactoMineR")

if (!require("dplyr", quietly = TRUE))
  install.packages("dplyr")

if (!require("RColorBrewer", quietly = TRUE))
  install.packages("RColorBrewer")

if (!require("wordcloud", quietly = TRUE))
  install.packages("wordcloud")

library(ComplexHeatmap)
library(ggbiplot)
library("ggplot2")
library("ggfortify")
library("gridExtra")
library("carData")
library("car")
library("factoextra")
library('pcaMethods')
library(ggplot2)
library(dplyr)
library(clusterProfiler)
library(RColorBrewer)
library(wordcloud)
library(pathview)

# Read the .txt file after FeatureCounts
simple_counts <- read.table("simple_counts_STAR.txt", 
                            header = TRUE, sep="\t", na.strings = "NA", row.names = 1,
                            dec=".")

# Rename the columns 
colnames(simple_counts)[1] <- "Female_Oregon_1"
colnames(simple_counts)[2] <- "Female_Oregon_2"
colnames(simple_counts)[3] <- "Female_CHD1n_1"
colnames(simple_counts)[4] <- "Female_CHD1n_2"
colnames(simple_counts)[5] <- "Male_Oregon_1"
colnames(simple_counts)[6] <- "Male_Oregon_2"
colnames(simple_counts)[7] <- "Male_CHD1n_1"
colnames(simple_counts)[8] <- "Male_CHD1n_2"

# Evaluate the data quality
# Correlation between samples
correlation <- cor(simple_counts)

# Evaluate the square of Pearson correlation coefficient
cor.sq <- round(correlation**2, 3)
spearman <- cor(simple_counts, method = "spearman")

# Create table with describtion of experiment design

expr.matrix <- as.matrix(simple_counts)
is.matrix(simple_counts)
is.matrix(expr.matrix)
print(dim(expr.matrix))

# Create the design of experiment 

expr.design <- data.frame(row.names = colnames(simple_counts), condition = c("control", "control", "mutant", "mutant", "control", "control", "mutant", "mutant"))

expr.design <- data.frame(row.names = colnames(expr.matrix), # name of the raws are the name of samples
                          condition = c("control", "control",
                                        "mutant", "mutant",
                                        "control", "control",
                                        "mutant", "mutant"))

# Combine all prepeared data in common  variable

dds <- DESeqDataSetFromMatrix(countData = expr.matrix, colData = expr.design, design = ~ condition) 
dds <- DESeq(dds)
res <- results(dds)
head(res)

# Save obtained results with DEGs as a table in local directory

write.table(as.data.frame(res), "DESeq_Dm.txt", sep = "\t", col.name = T, row.names = T, quote = F, na = "NA")


# Normalize data on library size
dds <- estimateSizeFactors(dds)
sizeFactors(dds)

# Obtain normalized readcounts 
normalized.counts <- counts(dds, normalized = T)

# Save in table
write.table(normalized.counts, "Norm_Sac.txt", sep = "\t", col.name = T, row.names = T, quote = F, na = "NA")

# Estimate the dispersion 
plotDispEsts(dds)

# Visualize the results of data analysis
plotMA(dds, ylim = c(-2,2), main = "DESeq2")

# Evaluate statistically significant DEGs (FDR < 0.05)
sum(res$padj < 0.05, na.rm = T)
667

# Select in variable genes with FDR < 0.05
resSig <- subset(res, padj < 0.05)
head(res)
head(resSig)

### Select DEGs
# Divide into different variables genes with FC > 2 and FC < 0.05

resSigUp <- subset(resSig, log2FoldChange > 1)
resSigDown <- subset(resSig, log2FoldChange < -1)

sum(resSig$log2FoldChange > 1, na.rm = T)
249

sum(resSig$log2FoldChange < -1, na.rm = T)
209

# Save results in tables
write.table(as.data.frame(resSigUp), "DESeq2_Up.txt",
            sep = "\t", col.names = TRUE, row.names = TRUE,
            quote = FALSE, na = "NA")

write.table(as.data.frame(resSigDown), "DESeq2_Down.txt",
            sep = "\t", col.names = TRUE, row.names = TRUE,
            quote = FALSE, na = "NA")

#Heatmap

# Unit in common variable name of rows with DEGs
names <- c(row.names(resSigDown), row.names(resSigUp))

# Select the rows from table normalized.counts, where rows match to rows name in variable names
de.genes <- subset(normalized.counts, row.names(normalized.counts) %in% names)
Heatmap(de.genes, show_row_names = F)

# logarithmize the data to normalize the expression level
log.de.genes <- log10(de.genes)
Heatmap(log.de.genes, show_row_names = F)

# error (since our data contains genes that were not expressed in some samples) - add some number to all readcounts
de.genes <- de.genes +0.01
log.de.genes <- log10(de.genes)
Heatmap(log.de.genes, show_row_names = F)

# Scaled to improve visualization
scaled.row.de.genes <- t(scale(t(log.de.genes)))
Heatmap(scaled.row.de.genes, show_row_names = F)

#Volcano plot 
library(ggplot2)

# Save variable res in table 
de.results <- as.data.frame(res)

# Since the statistically significant values ​​of the adjusted p-value - padj that interest us are too small, we take their logarithm
ggplot(de.results, aes(log2FoldChange, -log(padj, 10))) +
  geom_point(size = 0.4)

library(dplyr)

# Let's reduce the y-axis and show which genes are growing, falling, etc.
de.results <- de.results %>%
  mutate(Expression = case_when(log2FoldChange >= log(2)
                                & padj < 0.05 ~ "Upregulated",
                                log2FoldChange <= -log(2)
                                & padj < 0.05 ~ "Downregulated",
                                T ~ "Unchanged"))


ggplot(de.results, aes(log2FoldChange, -log(padj, 10), color = Expression)) +
  geom_point(size = 0.4) +
  scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) + 
  theme_bw() +
  ylim(0, 20)

####Prepare Input

# reading in input from deseq2
df = read.csv("drosphila_example_de.csv", header=TRUE)

deseq_dm <- read.table("DESeq_Dm.txt", 
                       header = TRUE, sep="\t", na.strings = "NA", row.names = 1,
                       dec=".")

# we want the log2 fold change 
original_gene_list <- deseq_dm$log2FoldChange

# name the vector
names(original_gene_list) <- deseq_dm$X

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

# Exctract significant results (padj < 0.05)
sig_genes_df = subset(deseq_dm, padj < 0.05)

# From significant results, we want to filter on log2fold change
genes <- sig_genes_df$log2FoldChange

# Name the vector
names(genes) <- sig_genes_df$X

# omit NA values
genes <- na.omit(genes)

# filter on min log2fold change (log2FoldChange > 2)
genes <- names(genes)[abs(genes) > 2]

####Create enrichGO object

#Check which options are available with the keytypes command, for example 

keytypes(org.Dm.eg.db)

####Create the object

go_enrich <- enrichGO(gene = genes,
                      universe = names(gene_list),
                      OrgDb = organism, 
                      keyType = 'ALIAS',
                      readable = T,
                      ont = "BP",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)

###Barplot

barplot(go_enrich, 
        drop = TRUE, 
        showCategory = 10, 
        title = "GO Biological Pathways",
        font.size = 8)


###Dotplot

dotplot(go_enrich)






