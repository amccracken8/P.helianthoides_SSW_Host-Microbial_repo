# R-TopGO 

# BiocManager::install("topGO")
library(topGO)
library(GO.db)
library(tidyverse)
library(GOfuncR)

pycno_go <- read.table("pycnoGO_NE.tsv", header=T)
red_stat <- read.csv("red_PA_stat.csv")
purple_stat <- read.csv("purple_PA_stat.csv")

colnames(purple_stat)[1] <- "geneID"

dim(pycno_go[pycno_go$geneID %in% purple_stat$geneID,])[1]

# define groups to explore in topGO
sum(red_stat$sig==1) #770
sum(purple_stat$sig==1) #246


# Convert GO_terms column to a list
gene2GO <- strsplit(pycno_go$GO, ";")
names(gene2GO) <- pycno_go$geneID

# Create a named vector of log fold changes
geneList <- setNames(purple_stat$sig, purple_stat$geneID)

# Create topGO data object
GOdata_purple <- new("topGOdata",
              description = "GO enrichment analysis",
              ontology = "BP",  
              allGenes = geneList,
              geneSel = function(x) x == 1,# Define how to select target group expressed genes
              nodeSize = 10, # define min number of terms for category to be considered
              annot = annFUN.gene2GO,
              gene2GO = gene2GO)



#run test - elim method takes into account parental relationships
result_purple <- runTest(GOdata_purple, algorithm = "parentChild", statistic = "fisher", numChar=1000)

all_nodes_purple <- usedGO(GOdata_purple)


# Display top GO terms
pc_Res <- GenTable(GOdata_purple, classicFisher = result_purple, topNodes = length(all_nodes_purple),numChar=1000)


# set values as numeric. 
pc_Res$classicFisher <- as.numeric(pc_Res$classicFisher)
pc_Res$Annotated <- as.numeric(pc_Res$Annotated)
pc_Res$Significant <- as.numeric(pc_Res$Significant)
pc_Res_filt <- pc_Res[pc_Res$classicFisher < 0.05, ]


write.csv(pc_Res_filt, "purple_TOPGO_pc.csv", row.names=F)



# save dfs  
#======================================================# 

purple_pc_filt <- read.csv("purple_TOPGO_pc.csv")

red_pc_filt$gene_ratio <- red_pc_filt$Significant / red_pc_filt$Annotated
purple_pc_filt$gene_ratio <- purple_pc_filt$Significant / purple_pc_filt$Annotated


df = purple_pc_filt

# manually remove "error" causing go terms and catagories that are too broad > 1000 genes per catagory such as "metabolic process" 
df <- df[!df$Annotated >= 300, ]

# Initialize the vector for redundant terms
redundant_terms <- c()

# Get the list of terms from E_elim_filt
# remove super large catagories > 1000 terms
GO_list <- df$GO.ID
items_to_remove <- c("GO:0140694")
GO_list <- GO_list[!GO_list %in% items_to_remove]
terms_list <- df$Term

# Iterate over each term in the list
for (i in GO_list) {
  # Get parent terms for the current term
  
  term.i = df$Term[df$GO.ID==i] # term associated with term ID. 
  parents <- get_parent_nodes(i)$parent_name  # Extract parent terms
  parents <- parents[parents != term.i]
  
  # Check if any parent term is already in the list of terms
  if (any(parents %in% terms_list)) {
    redundant_terms <- c(redundant_terms, i)
  }
}



# add removed items to redundent terms
purple_pc_final <- df[!df$GO.ID %in% redundant_terms, ]

df.p <- purple_pc_final[order(purple_pc_final$classicFisher), ]


#subset to top 30 or can list select significant GO terms.
df.top <- df.p %>% arrange(classicFisher)%>%
  slice_head(n = 30)

#can add back in any terms of interest for visualization
## df.top <- rbind(df.top, E_pc_filt[E_pc_filt$GO.ID == "GO:0001666",]) # hypoxia


purple <- ggplot(df.top, aes(x=-log10(classicFisher), y=reorder(Term, -log10(classicFisher)), size=Significant, color=gene_ratio)) + 
  scale_color_continuous(low = "lightblue", high = "purple3") +
  geom_point(alpha = 0.8) +
  xlab("-log10(classicFisher)") + 
  ylab("") + 
  labs(color = "DEG Ratio", size="# of DEGs") +  # Change the color legend label
  theme(legend.position = "right", axis.text = element_text(face = "bold")) 


# combine plots
library(patchwork)
red|purple
