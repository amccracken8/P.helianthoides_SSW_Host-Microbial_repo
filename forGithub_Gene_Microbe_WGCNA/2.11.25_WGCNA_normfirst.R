#WGCN code
# adapted from https://bioinformaticsworkbook.org/tutorials/wgcna.html#gsc.tab=0

library("tidyverse")
library("WGCNA")
library("DESeq2")
library("data.table")

#load in counts tables
meta <- read_table("pycno_samples_NE.txt")


# take table.lev7.csv and subset it to samples to compare:
microbes <- read.csv("level-7.csv")
micro_sub <- microbes[microbes$sampleID %in% meta$sample,]

#rename sampleID to match Genes
micro_sub <- micro_sub %>% mutate(sampleID = meta$sampleID)
#write.csv(micro_sub, "micro_sub.csv", row.names = F)


##############################################################
# Start Here
micro_all <- read.csv("micro_sub.csv")
genes_all <- read.csv("genes_all.csv")



#remove unwanted metadata from tables
metadata <- micro_all %>% dplyr::select(colnames(micro_all)[!grepl("^(d_)", colnames(micro_all))])
micro_all <- micro_all %>% dplyr::select(sampleID, colnames(micro_all)[grepl("^(d_)", colnames(micro_all))])
genes_all <- genes_all %>% dplyr::select(sampleID,colnames(genes_all)[grepl("^(g)", colnames(genes_all))])

# note: for WGCNA rows must be samples and columns must be genes
### make sampleID the rowname
micro_all <- micro_all %>% 
  column_to_rownames(var = "sampleID")

genes_all <- genes_all %>% 
  column_to_rownames(var = "sampleID")


# or remove counts where >70% of the samples have zero coutns. min of 3/10
micro_filt_0 <- micro_all %>%
  dplyr::select(where(~ is.numeric(.) && mean(. == 0) <= 0.70))
genes_filt_0 <- genes_all %>%
  dplyr::select(where(~ is.numeric(.) && mean(. == 0) <= 0.70  && sum(.) > 100))


###############################
### NORMALIZATION 
###############################

### DESEQ normalization BEFORE concatonization. 

## combined normalizing with DESeq2
# DESeq2 expects rows to be genes and cols to be conditions (opposite of WGCNA) so must transpose first and then back again after

t_micro_filt <- t(micro_filt_0)
t_genes_filt <- t(genes_filt_0)


meta <- read.table("pycno_samples_NE.txt", sep="\t", header=T, row.names=1)

dds_micro <- DESeqDataSetFromMatrix(countData = t_micro_filt , colData = meta, design = ~ 1)
dds_genes <- DESeqDataSetFromMatrix(countData = t_genes_filt , colData = meta, design = ~ 1)


# dds_micro <- DESeqDataSetFromMatrix(countData = t_micro_filt , colData = meta, design = ~ 0 + health_site_status)
# dds_genes <- DESeqDataSetFromMatrix(countData = t_genes_filt , colData = meta, design = ~ 0 + health_site_status)


dds_micro <- DESeq(dds_micro)
dds_genes <- DESeq(dds_genes)


normalized_micro <- counts(dds_micro, normalized = TRUE) # extract normalized counts
normalized_genes <- counts(dds_genes, normalized = TRUE) # extract normalized counts


df_combined <- rbind(normalized_micro, normalized_genes)


normalized_counts_log <- log2(df_combined + 1) #optional log normalization of deseq2 normalized counts


#######################################################
###### WGCNA ############
datExpr <- data.frame(t(normalized_counts_log)) # transpose back to rows = samples and cols=genes/taxa
#write.csv(datExpr, "dataExpr_normfirst.txt")

datExpr <- read.csv("dataExpr_normfirst.txt",row.names = 1)

# check for outliers and "bad" genes to exclude
gsg = goodSamplesGenes(datExpr, verbose = 3) # all are good

# Choose a set of soft-thresholding powers

powers = c(c(1:10), seq(from = 12, to = 20, by = 2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# Choose a power based on the scale-free topology criterion
softPower = sft$powerEstimate

### visualize threshhold values
par(mfrow = c(1,2));
cex1 = 0.9;

plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = paste("Scale independence")
)
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red"
)
abline(h = 0.90, col = "red")


plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity")
)
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers,
     cex = cex1, col = "red")

### Picked threshhold of 7
picked_power = 7

### create the network using the blockwiseModules command

temp_cor <- cor       
cor <- WGCNA::cor         # Force it to use WGCNA cor function (fix a namespace conflict issue)



#################################################################

### run WGCNA - make network and detect modules
adjacency = adjacency(datExpr, power = softPower, type = "signed")

#write.csv(adjacency,"adjacency.nf.csv")

TOM = TOMsimilarity(adjacency)
dissTOM = 1 - TOM
geneTree = hclust(as.dist(dissTOM), method = "average")

dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = 30)

#write_lines(dynamicMods,"dynmicMods_nf.txt")
moduleColors = labels2colors(dynamicMods)

#write_lines(moduleColors,"moduleColors_nf.txt")


# Plot the dendrogram and the module colors
sizeGrWindow(8, 6)
plotDendroAndColors(geneTree, moduleColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

#################################################################

#======================re-start here for plots ====================================================#
#==== files to load in from WGCNA run to make plots and analysis below =================#
datExpr<-read.csv("dataExpr_normfirst.txt", row.names=1)
moduleColors <- read_lines("moduleColors_nf.txt")
meta <- read.table("pycno_samples_NE.txt", sep="\t", header=T, row.names=1)



##################################################################
# modules and their associated genes
# Assuming 'datExpr' has samples as rows and genes as columns
geneNames <- colnames(datExpr)

# Create a data frame with gene names and their corresponding module colors
geneModuleData <- data.frame(Gene = geneNames, Module = moduleColors)


# Split the genes by their module colors
modules <- split(geneModuleData$Gene, geneModuleData$Module)


# Prepare the data for writing to a file
moduleList = do.call(rbind, lapply(names(modules), function(module) {
  data.frame(Module = module, Gene = modules[[module]])
}))

moduleGeneList <- aggregate(Gene ~ Module, data = moduleList, FUN = function(x) paste(x, collapse = ", "))

#write.table(moduleList,"module_list_long_nf.txt",sep="\t", row.names=F)
#write.table(moduleGeneList,"module_list_short_nf.txt", row.names=F)






########### Module correlations between traits ###################

# Get Module Eigengenes per cluster
MEs0 <- moduleEigengenes(datExpr, colors = moduleColors)$eigengenes


# Reorder modules so similar modules are next to each other
MEs0 <- orderMEs(MEs0)
module_order = names(MEs0) %>% gsub("ME","", .)

#write.csv(MEs0, "Module_Engine_MEs0_nf.csv")
#MEs0<- read.csv("Module_Engine_MEs0_nf.csv", row.names=1)


# Add treatment names
MEs0$treatment = row.names(MEs0)

# tidy & plot data
mME = MEs0 %>%
  pivot_longer(-treatment) %>%
  mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order)
  )

#plot

mME %>% ggplot(., aes(x=treatment, y=name, fill=value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1,1)) +
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Module-trait Relationships", y = "Modules", fill="corr")


########################################################
## KME (module membership (kME) values,) for each gene #
#remove $treatment from MES0

kME <- signedKME(datExpr, MEs0)
kME <- as.data.frame(kME)
colnames(kME) <- sub("kME", "", colnames(kME))
write.csv(kME, "gene_micro_moduleMembership_nf.csv")
#test <- read.csv("gene_micro_moduleMembership_nf.csv")

### create df with gene, module, and kME.

# Assuming 'datExpr' has samples as rows and genes as columns
geneNames <- colnames(datExpr)

# Create a data frame with gene names and their corresponding module colors
geneModule <- data.frame(Gene = geneNames, Module = moduleColors)

geneModule$EigengeneValue <- NA

# Loop through each row of geneModule and assign eigengene values based on module and gene
for (i in 1:nrow(geneModule)) {
  gene <- geneModule$Gene[i]
  module <- geneModule$Module[i]
  
  # Get the eigengene vector for the current module
  eigengene_vector <- kME[[module]]
  
  # Assign the eigengene value to the current gene
  gene_index <- which(rownames(kME) == gene)
  geneModule$EigengeneValue[i] <- eigengene_vector[gene_index]
}

#write.csv(geneModule,"gene_micro_eigns_nf.csv", row.names=F)



######### module corelation relationships #####################
# Convert the trait data to a numeric format if not already done

meta$status <- ifelse(meta$health_site_status == "naive", 0, 
                           ifelse(meta$health_site_status == "exposed", 1, NA))


# Compute the correlation between module between naive and exposed "traits" 
moduleTraitCor <- cor(MEs0, meta$status, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples = nrow(datExpr))

# Prepare the text for the heatmap
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)


# Plot the heatmap
sizeGrWindow(10, 6)
par(mar = c(10, 10, 4, 2) + 0.1)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(meta)[4], # Assuming the first column is sample names
               yLabels = names(MEs0),
               ySymbols = names(MEs0),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               cex.lab = 0.7, 
               zlim = c(-1, 1),
               main = paste("Module-health relationships"))




## -----------------------------------------------------------------------------

## pull out only significant cors
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples = nrow(datExpr))

# Filter the modules with p-values less than 0.05
significantModules <- moduleTraitPvalue[, 1] < 0.05

# Subset the module eigengenes and trait correlation matrices to include only significant modules
significantModuleTraitCor <- moduleTraitCor[significantModules, , drop = FALSE]
significantModuleTraitPvalue <- moduleTraitPvalue[significantModules, , drop = FALSE]

# Prepare the text for the heatmap
significantTextMatrix <- paste(signif(significantModuleTraitCor, 2), "\n(",
                               signif(significantModuleTraitPvalue, 1), ")", sep = "")
dim(significantTextMatrix) <- dim(significantModuleTraitCor)

# Set the margins to be larger to accommodate the full text

par(mar = c(5, 10, 4, 2) + 0.1)
# Plot the heatmap for the significant modules
labeledHeatmap(Matrix = significantModuleTraitCor,
               xLabels = "Mod Cor and Pval", # Trait name
               yLabels = rownames(significantModuleTraitCor), # Significant module eigengenes
               ySymbols = rownames(significantModuleTraitCor),
               colorLabels = FALSE,
               colors = colorRampPalette(c("lightgreen", "white", "purple"))(50),
               textMatrix = significantTextMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               cex.lab = 0.7, 
               zlim = c(-1, 1),
               main = paste("Significant Module-health relationships"))



sig_mod_list <- row.names(significantModuleTraitCor)
sig_mod_list <- sig_mod_list[!sapply(sig_mod_list, is.na)]
sig_mod_list <- substring(sig_mod_list, 3)

#write.csv(sig_mod_list, "sig_mod_list_nf.csv")

# plot corelations 
mME_sig <-  mME[mME$name %in% sig_mod_list, ]

mME_sig <- mME_sig %>%
  mutate(group = ifelse(startsWith(treatment, "A"), "Naive", 
                        ifelse(startsWith(treatment, "E"), "Exposed", NA)))


sigMods<-read.table("module_list_sig_nf.txt",header=T)
sigMods_gm <- sigMods %>% filter(grepl("^d", Gene))
mME_sig_micro <- mME_sig[mME_sig$name %in% sigMods_gm$Module,]

wgcna_plot <- mME_sig_micro %>% ggplot(., aes(x=treatment, y=name, fill=value)) +
  geom_tile() +
  theme_bw() +
  geom_vline(xintercept = 5.5, linetype = "dashed", color = "black", size = 1) +
  scale_fill_gradient2(
    low = "springgreen3",
    high = "purple3",
    mid = "white",
    midpoint = 0,
    limit = c(-1,1)) +
  theme(axis.text.x = element_text(angle=90)) +

  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12,face = "bold"),  # Horizontal labels, centered, size increased
    axis.text.y = element_text(size = 12,face = "bold"),  # Increase size of y-axis labels too if needed
    axis.title.x = element_text(size = 12,face = "bold"),
    axis.title.y = element_text(size = 12,face = "bold")
  )+
  
  labs(title = "Significant Module-trait Relationships", y = "Modules", x="Health-Group", fill="corr") +
  scale_x_discrete(labels = function(x) {
    # Manually create x-axis labels with only "Naive" for the first set and "Exposed" for the second
    labels <- rep("", length(x))  # Initialize empty labels
    labels[3] <- "Naive"          # Set "Naive" at the first treatment
    labels[length(x)-2] <- "Exposed" # Set "Exposed" at the last treatment
    return(labels)
  })
  






####### get list of genes in significant modules

moduleGeneList_sig <- moduleGeneList[moduleGeneList$Module %in% sig_mod_list, ]
#write.table(moduleGeneList_sig,"module_list_sig_nf.txt", row.names=F)

#test<-read.table("module_list_sig_nf.txt",header=T)

### Get list of genes/microbes in module

# Function to get genes by module using dplyr
get_genes_by_module <- function(data, module_name) {
  genes <- data %>%
    filter(Module == module_name) %>%
    pull(Gene)
  return(genes)
}

get_genes_by_module(moduleGeneList_sig, "bisque4")
