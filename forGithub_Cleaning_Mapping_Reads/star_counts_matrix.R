# script for making a counts matrix from "ReadsPerGene.out.tab" output from STAR output
# 01/9/24


# navigate to where star-counts files are located

library(tidyverse)
library(dplyr)

AH1R10 <- read.table("AH1R10_ReadsPerGene.out.tab", header = FALSE, sep = "\t",skip = 4)
AH1R17 <- read.table("AH1R17_ReadsPerGene.out.tab", header = FALSE, sep = "\t",skip = 4)
AH1R3 <- read.table("AH1R3_ReadsPerGene.out.tab", header = FALSE, sep = "\t",skip = 4)
AH2R12 <- read.table("AH2R12_ReadsPerGene.out.tab", header = FALSE, sep = "\t",skip = 4)
AH2R5 <- read.table("AH2R5_ReadsPerGene.out.tab", header = FALSE, sep = "\t",skip = 4)

E12H3R <- read.table("E12H3R_ReadsPerGene.out.tab", header = FALSE, sep = "\t",skip = 4)
E12H4R <- read.table("E12H4R_ReadsPerGene.out.tab", header = FALSE, sep = "\t",skip = 4)
E13H8R <- read.table("E13H8R_ReadsPerGene.out.tab", header = FALSE, sep = "\t",skip = 4)
E13H9R <- read.table("E13H9R_ReadsPerGene.out.tab", header = FALSE, sep = "\t",skip = 4)
E15H3R <- read.table("E15H3R_ReadsPerGene.out.tab", header = FALSE, sep = "\t",skip = 4)


samples <- list('AH1R10','AH1R17','AH1R3','AH2R12','AH2R5','E12H3R','E12H4R','E13H8R','E13H9R','E15H3R')

df_list <- list(AH1R10,AH1R17,AH1R3,AH2R12,AH2R5,E12H3R,E12H4R,E13H8R,E13H9R,E15H3R)


# Rename columns to preserve sampleID incountrs matrix
colnames(AH1R10) <- c("Gene", "AH1R10", "F", "R")
colnames(AH1R17) <- c("Gene", "AH1R17", "F", "R")
colnames(AH1R3) <- c("Gene", "AH1R3", "F", "R")
colnames(AH2R12) <- c("Gene", "AH2R12", "F", "R")
colnames(AH2R5) <- c("Gene", "AH2R5", "F", "R")

colnames(E12H3R) <- c("Gene", "E12H3R", "F", "R")
colnames(E12H4R) <- c("Gene", "E12H4R", "F", "R")
colnames(E13H8R) <- c("Gene", "E13H8R", "F", "R")
colnames(E13H9R) <- c("Gene", "E13H9R", "F", "R")
colnames(E15H3R) <- c("Gene", "E15H3R", "F", "R")



# take 2nd row of each data frame and add it to new counts_df

counts_df <- data.frame(matrix(nrow = 24184, ncol = 0),stringsAsFactors = FALSE)

for (i in 1:10) {
    counts_df <- cbind(counts_df, df_list[[i]][, 2, drop = FALSE])
}

#colnames ORDER IS IMPORTANT must be the same order as df_list
colnames(counts_df) <- c('AH1R10','AH1R17','AH1R3','AH2R12','AH2R5','E12H3R','E12H4R','E13H8R','E13H9R','E15H3R')

test <- counts_df

rownames(test) <- AH1R10[, 1]

write.csv(test, "star_counts_df.csv", row.names = TRUE)

