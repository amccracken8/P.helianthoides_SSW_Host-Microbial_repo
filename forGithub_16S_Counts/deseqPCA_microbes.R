#WGCN code
# adapted from https://bioinformaticsworkbook.org/tutorials/wgcna.html#gsc.tab=0

library("tidyverse")
library("WGCNA")
library("DESeq2")

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



#remove unwanted metadata from tables
metadata <- micro_all %>% select(colnames(micro_all)[!grepl("^(d_)", colnames(micro_all))])
micro_all <- micro_all %>% select(sampleID, colnames(micro_all)[grepl("^(d_)", colnames(micro_all))])

# note: for WGCNA rows must be treatments and columns must be genes
### make sampleID the rowname
micro_all <- micro_all %>% 
  column_to_rownames(var = "sampleID")

micro_filt_0 <- micro_all %>%
  select(where(~ is.numeric(.) && mean(. == 0) <= 0.70))


######## normalizaation

meta <- read.table("pycno_samples_NE.txt", sep="\t", header=T, row.names=1)

tdf.micro <- t(micro_filt_0)

dds_micro <- DESeqDataSetFromMatrix(countData = tdf.micro , colData = meta, design = ~ 0 + health_site_status)
dds_micro <- DESeq(dds_micro)
resultsNames((dds_micro))

#####################
# PCA Visualizations#
#####################

vsd <- varianceStabilizingTransformation(dds_micro, blind = FALSE)

data <- plotPCA(vsd, intgroup = "health_site_status", returnData = TRUE)
percentVar <- round(100 * attr(data, "percentVar"))

color_vector <- c("skyblue", "lightcoral", "magenta")

m <- ggplot(data, aes(PC1, PC2, color = health_site_status)) +
  geom_point(size = 6, alpha = 0.85) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  scale_color_manual(
    values = setNames(color_vector, unique(data$health_site_status)),
    breaks = unique(data$health_site_status), 
    labels = toupper(unique(data$health_site_status)) 
  ) +
  labs(color = "Health Status", title="B. Microbial Abundance") + 
  theme_light() + 
  theme(axis.title = element_text(face = "bold"))

m # + guides(color = T)


### generate PCA from Deseq2 on expression data 
library(patchwork)

e+m + plot_layout(guides = "collect", axis="collect") 


