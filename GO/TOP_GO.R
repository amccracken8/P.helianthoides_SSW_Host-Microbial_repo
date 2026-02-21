# R - TopGO gene ontology enrichment test for enrichment in Exposed and Naive groups in comparison to others. 
# Fisher exact test on if genes were significantly differential expressed or not 

# BiocManager::install("topGO")
library(topGO)
library(GO.db)
library(tidyverse)
library(GOfuncR)

pycno_go <- read.table("pycnoGO_NE.tsv", header=T)
NE_stat <- read.csv("res_NE.csv")
colnames(NE_stat)[1] <- "geneID"

dim(pycno_go[pycno_go$geneID %in% NE_stat$geneID,])[1] # confirm data frames match


# significant = 1 nonsignificant = 0 for Fisher Exact Test of TopGO Enrichment
NE_stat <- NE_stat %>%
  mutate(Group = case_when(
    log2FoldChange > 0 & padj < 0.05 ~ "up_E",
    log2FoldChange < 0 & padj < 0.05 ~ "up_N",
    TRUE ~ "nonsig"
  ))

sum(NE_stat$Group=="up_E") #1263
sum(NE_stat$Group=="up_N") #743


# Convert GO_terms column to a list
gene2GO <- strsplit(pycno_go$GO, ";")
names(gene2GO) <- pycno_go$geneID

geneList_N <- setNames(NE_stat$padj, NE_stat$geneID)
geneList_E <- setNames(NE_stat$padj, NE_stat$geneID)

### Define Target genes list for group enrichment. For example, if it is up in exposed, and significant (the up_E group) then we assign a value on 1 and pull out all genes with that value for the GOdata function below.  
geneList_N <- setNames(as.integer(NE_stat$Group == "up_N"), NE_stat$geneID)
geneList_E <- setNames(as.integer(NE_stat$Group == "up_E"), NE_stat$geneID)


# Create topGO data object
GOdata_N <- new("topGOdata",
              description = "GO enrichment analysis",
              ontology = "BP",  # Use "BP" for Biological Process
              allGenes = geneList_N,
              geneSel = function(x) x == 1,# Define how to select target group expressed genes
              nodeSize = 10, # define min number of terms for category to be considered
              annot = annFUN.gene2GO,
              gene2GO = gene2GO)

# Create topGO data object
GOdata_E <- new("topGOdata",
              description = "GO enrichment analysis",
              ontology = "BP",  
              allGenes = geneList_E,
              geneSel = function(x) x == 1,
              nodeSize = 10,
              annot = annFUN.gene2GO,
              gene2GO = gene2GO)




#run test
result_E <- runTest(GOdata_E, algorithm = "parentChild", statistic = "fisher", numChar=1000)
result_N <- runTest(GOdata_N, algorithm = "parentChild", statistic = "fisher", numChar=1000)


all_nodes_E <- usedGO(GOdata_E)
all_nodes_N <- usedGO(GOdata_N)


# Display top GO terms
E_pc_Res <- GenTable(GOdata_E, classicFisher = result_E, topNodes = length(all_nodes_E),numChar=1000)
N_pc_Res <- GenTable(GOdata_N, classicFisher = result_N, topNodes = length(all_nodes_N),numChar=1000)


# set values as numeric. 
E_pc_Res$classicFisher <- as.numeric(E_pc_Res$classicFisher)
E_pc_Res$Annotated <- as.numeric(E_pc_Res$Annotated)
E_pc_Res$Significant <- as.numeric(E_pc_Res$Significant)

N_pc_Res$classicFisher <- as.numeric(N_pc_Res$classicFisher)
N_pc_Res$Annotated <- as.numeric(N_pc_Res$Annotated)
N_pc_Res$Significant <- as.numeric(N_pc_Res$Significant)


E_pc_filt <- E_pc_Res[E_pc_Res$classicFisher < 0.05, ]
N_pc_filt <- N_pc_Res[N_pc_Res$classicFisher < 0.05, ]

#write.csv(N_pc_filt, "Naive_TOPGO_pc.csv", rownames=F)
#write.csv(E_pc_filt, "Exposed_TOPGO_pc.csv", rownames=F)


### Add group identifier before combining datasets. 
E_pc_filt$group <- "Exposed"
N_pc_filt$group <- "Naive"


# save dfs  
N_pc_filt <- read.csv("Naive_TOPGO_pc.csv")
E_pc_filt <- read.csv("Exposed_TOPGO_pc.csv")

N_pc_filt$gene_ratio <- N_pc_filt$Significant / N_pc_filt$Annotated
E_pc_filt$gene_ratio <- E_pc_filt$Significant / E_pc_filt$Annotated

#======================================================# 
# Choose N_pc_filt or E_pc_filt as df
#======================================================# 

df = N_pc_filt

###  removecategories that are too broad > 500 genes per category such as "metabolic process" 
df <- df[!df$Annotated >= 500, ]

# Initialize the vector for redundant terms
redundant_terms <- c()

# Only GO:0140694 caused errors and was removed (membraneless organelle assembly) 
GO_list <- df$GO.ID
items_to_remove <- c("GO:0140694")
GO_list <- GO_list[!GO_list %in% items_to_remove]
terms_list <- df$Term


for (i in GO_list) {
  term.i = df$Term[df$GO.ID==i] # term associated with term ID. 
  parents <- get_parent_nodes(i)$parent_name  # Extract parent terms
  parents <- parents[parents != term.i]
  
  # Check if any parent term is already in the list of terms
  if (any(parents %in% terms_list)) {
    redundant_terms <- c(redundant_terms, i)
  }
}


# add removed items to redundent terms
N_pc_final <- df[!df$GO.ID %in% redundant_terms, ]

df.p <- N_pc_final[order(E_pc_final$classicFisher), ]

#subset to top 30 or can list select significant GO terms.
df.top <- df.p %>% arrange(classicFisher)%>%
  slice_head(n = 30)


## plot Naive and Exposed Enrichment
N <- ggplot(df.top, aes(x=-log10(classicFisher), y=reorder(Term, -log10(classicFisher)), size=Significant, color=gene_ratio)) + 
  scale_color_continuous(low = "skyblue", high = "darkblue") +
  geom_point(alpha = 0.8) +
  xlab("-log10(P-value)") + 
  ylab("") + 
  labs(color = "DEG Ratio", size="# of DEGs") +  # Change the color legend label
  theme(legend.position = "right", axis.text.y = element_text(face = "bold"))



E <- ggplot(df.top, aes(x=-log10(classicFisher), y=reorder(Term, -log10(classicFisher)), size=Significant, color=gene_ratio)) + 
  scale_color_continuous(low = "pink", high = "firebrick4") +
  geom_point(alpha = 0.8) +
  xlab("-log10(P-value)") + 
  ylab("") + 
  labs(color = "DEG Ratio", size="# of DEGs") +  # Change the color legend label
  theme(legend.position = "right", axis.text.y = element_text(face = "bold"))



# combine plots
library(patchwork)
N|E

