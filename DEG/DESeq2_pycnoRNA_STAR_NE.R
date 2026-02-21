## Import the libraries that we're likely to need in this session
library(DESeq2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(ggpubr)
library(wesanderson)
library(vsn)  ### First: BiocManager::install("vsn") AND BiocManager::install("hexbin")

# Navigate to star_countsMatrix.txt directory

countsTable <- read.table("star_countsMatrix.txt", header=TRUE, row.names=1)
countsTable_NE <- countsTable[, 1:10]
head(countsTable_NE)
dim(countsTable_NE)
# 24184    10


#DESeq2 doesn't like decimals
countsTableRound <- round(countsTable_NE) 
head(countsTableRound)


#import metadata table (Naive and Exposed only "NE")
conds <- read.delim("pycno_samples_NE.txt", header=TRUE, stringsAsFactors = TRUE, row.names=1)
head(conds)

# How many reads we have from each sample?
colSums(countsTableRound) 
mean(colSums(countsTableRound))
barplot(colSums(countsTableRound), names.arg=colnames(countsTableRound),cex.names=1, las=3,ylim=c(0,25000000), main = "Number of mapped Reads")
abline(h=mean(colSums(countsTableRound)), col="red", lwd=2)

# the average number of counts per gene
rowSums(countsTableRound)
mean(rowSums(countsTableRound))
median(rowSums(countsTableRound))

#apply(countsTableRound,2,mean) # 2 in the apply function does the action across columns
apply(countsTableRound,1,mean) # 1 in the apply function does the action across rows

hist(apply(countsTableRound,1,mean), xlim=c(1,1000), ylim=c(0,10000),breaks=10000, xlab="avg count", main="mean count frequency")

#### Number of genes under with total < 10 counts per samples
filt_count_round <- countsTableRound %>%
  filter(rowSums(countsTableRound) < 100)

dim(filt_count_round) # 8425   10


hist(apply(filt_count_round,1,mean), xlim=c(0,20), ylim=c(0,10000),breaks=10, xlab="avg count", main="mean count frequency")



#==============================================================================# 
#### Create a DESeq object and define the experimental design ~
#==============================================================================# 

# the ~ '0 + ' designates no intercept 
dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData=conds, 
                                 design= ~ 0 + health_site_status)

#filter out reads with less than 10 reads per sample
dim(dds) # [1] 24184    10
dds <- dds[rowSums(counts(dds)) > 100] #100 because 10 samples ~avg more than 10 counts per sample
dim(dds) # [1] 15743    10


# Run the DESeq model to test for differential gene expression
dds <- DESeq(dds)
resultsNames(dds)
# [1] "health_site_statusexposed" "health_site_statusnaive"


#==============================================================================# 
#### PCA from dds object
#==============================================================================#

vsd <- vst(dds, blind=FALSE)

data <- plotPCA(vsd, intgroup = "health_site_status", returnData = TRUE)
percentVar <- round(100 * attr(data, "percentVar"))

color_vector <- c("skyblue", "lightcoral", "magenta")

ggplot(data, aes(PC1, PC2, color = health_site_status)) +
  geom_point(size = 6, alpha = 0.85) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  scale_color_manual(
    values = setNames(color_vector, unique(data$health_site_status)),
    breaks = unique(data$health_site_status),  # Ensure correct mapping
    labels = toupper(unique(data$health_site_status))  # Convert labels to uppercase
  ) +
  labs(color = "Health Status", title="A. Gene Expression") + 
  theme_light()
  #guides(color = T)  # This removes the color legend




#####################
# Order and summarize the results from specific contrasts
#####################

res<- results(dds, alpha=0.05)
res <- res[order(res$padj),]
head(re)  

# baseMean log2FoldChange     lfcSE      stat      pvalue        padj
# <numeric>      <numeric> <numeric> <numeric>   <numeric>   <numeric>
#   g14574 2100.7335        3.15223  0.313024  10.07025 7.47879e-24 1.13827e-19
# g9913   328.3688      -27.25972  3.026350  -9.00746 2.10887e-19 1.60485e-15
# g4570    92.6155      -25.26112  3.026709  -8.34607 7.05729e-17 3.58040e-13
# g11468   73.6005      -25.09949  3.026838  -8.29231 1.11066e-16 4.22606e-13
# g6421    66.9932      -24.96465  3.026900  -8.24760 1.61611e-16 4.91945e-13
# g1002    42.5237       -3.32000  0.420940  -7.88712 3.09242e-15 7.84444e-12


summary(res)


# out of 15743 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 743, 4.7%
# LFC < 0 (down)     : 1263, 8%
# outliers [1]       : 523, 3.3%
# low counts [2]     : 0, 0%
# (mean count < 8)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

## Remove NA's
res <- res[!is.na(res$padj),]

## grab deferentially expressed transcripts
degs <- row.names(res[res$padj < 0.05,])
length(degs)



# test for differences between naive and exposed groups. 
res_NE <- results(dds, contrast=c("health_site_status","exposed","naive"), alpha=0.05)
res_NE <- res_NE[order(res_NE$padj),]
head(res_NE)
#write.csv(res_NE, "res_NE.csv", row.names = TRUE)

summary(res_NE)

# out of 15743 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 743, 4.7%
# LFC < 0 (down)     : 1263, 8%
# outliers [1]       : 523, 3.3%
# low counts [2]     : 0, 0%
# (mean count < 8)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

res_NE <- res_NE[!is.na(res_NE$padj),]
degs_NE <- row.names(res_NE[res_NE$padj < 0.05,])


## ---------------- heat map:
library(pheatmap)

color_vector <- c(naive = "blue", exposed = "orange")

res_NE <- res_NE[!is.na(res_NE$padj),]
topgenes <- head(rownames(res_NE),50)
mat <- assay(vsd)[topgenes,]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(dds)['health_site_status'])

df$health_site_status <- factor(df$health_site_status)

heatmap <- pheatmap(mat, annotation_col = df, annotation_colors = list(health_site_status = color_vector), cluster_rows = TRUE,cluster_cols = F)




### plot individual genes
# Counts of specific top interaction gene! (important validatition that the normalization, model is working)
d <-plotCounts(dds, gene="g12090", intgroup = ("health_site_status"), returnData=TRUE)
d

ggplot(d, aes(x = health_site_status, y = count, fill = health_site_status)) +
  geom_boxplot() +
  labs(title = "Average Count by Health Site Status",
       x = "Health Site Status",
       y = "Average Count") +
  theme_minimal()




