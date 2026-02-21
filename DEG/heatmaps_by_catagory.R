# heatmaps of Differentially expressed genes organized by catagory

library(tidyverse)
library(data.table)
library(patchwork)
library(dplyr)
library(stringr)
library(pheatmap)

#rename ALL vectors? 

###differentually abundant genes
degs <- read.csv("NE_DEG_annotated.csv")

###differentially abundant microbes

microbes <- read.csv("8.1.24_res.NvE_Sylva.csv") # 244 microbial families found passing the >10 reads/sample average 
dams <- microbes[microbes$diff_abn.site.animal.healthSH==TRUE,] # 39 differential abundance microbes


### Gene ontology mapping to DEGs
pGO <- fread("Pycno_GO_formatted.tsv")
colnames(pGO) <- c("gene_id","GO")
degs_GO <- left_join(degs,pGO, by="gene_id")


dataExpr <- fread("dataExpr_normfirst.txt")
dataExpr <- data.table(dataExpr)


# load in normalized DESeq2 counts data
#dataExpr <- read.csv("dataExpr.txt", header=T, row.names=1)
rownames(dataExpr)

dataExpr$group <- case_when(
  dataExpr$V1 %in% dataExpr$V1[1:5] ~ "Naive",
  dataExpr$V1 %in% dataExpr$V1[6:10] ~ "Exposed"
)


df <- dataExpr %>%
  pivot_longer(cols = -c(V1, group),  # Exclude 'V1' and 'group' columns
               names_to = "gene",
               values_to = "expression")

colnames(df)[1] <- "sample"
colnames(df)[3] <- "Gene"


# adding GO terms from gp.pi table which were previously mapped each go term to Immune, tissue, neuro, and stress catagories
go.pi <- read.csv("pixy_pi_go.csv")
  
cats <- go.pi %>% dplyr::select(Gene,immune_go,nerv_go,tissue_go,stress_go)


df.sets <- left_join(df,cats, by="Gene")

###subset to only those that are differentially expressed
df.sets <- df.sets[df.sets$Gene %in% degs$gene_id,]


# add annotations
ann.df <- read.csv("NE_ALL_annotated.csv")

ann.df <- ann.df %>% dplyr::select(gene_id,annote)
colnames(ann.df)[1] <- "Gene"

df.sets.ann <- left_join(df.sets,ann.df, by="Gene")

### include if you want all genes seperatly plotted, or average accross paraloges
df.sets.ann <- df.sets.ann %>%
  mutate(annote_clean_gene = paste(Gene, str_extract(annote, "(?<=\\s)[^\\[]+(?=\\s*\\[)"), sep = " "))

df.sets.ann <- df.sets.ann %>%
  mutate(annote_clean = str_extract(annote, "(?<=\\s)[^\\[]+(?=\\s*\\[)"))


############# Immune Subset heatmap ##################
immune_df <- df.sets.ann %>%
  filter(
    str_detect(annote, "ficolin") |
      str_detect(annote, "echinoidin") |
      str_detect(annote, "macrophage") |
      str_detect(annote, "TNF|tumor necrosis") |
      str_detect(annote, "lymphocyte") |
      str_detect(annote, "toll") | 
      str_detect(annote, "t-cells")|
      str_detect(annote, "interleukin")|
      str_detect(annote, "complement") |
      str_detect(annote, "scavenger") |
      str_detect(annote, "defense") |
      str_detect(annote, "cell death")
    # "deleted in malignant brain"
    # "HDD11"
    # "E3 ubiquitin-protein"
  )

complement_df <- df.sets.ann %>%
  filter(
    str_detect(annote, "ficolin") |
      str_detect(annote, "echinoidin") |
      str_detect(annote, "complement") |
      str_detect(annote, "scavenger") 
  )

inflam_df <- df.sets.ann %>%
  filter(
    str_detect(annote, "tumor necrosis factor receptor") |
      str_detect(annote, "TNF") |
      str_detect(annote, "interleukin") |
      str_detect(annote, "death") |
      str_detect(annote, "leukocyte") 
  )

recog_df <- df.sets.ann %>%
  filter(
    str_detect(annote, "macrophage") |
      str_detect(annote, "t-cell") |
      str_detect(annote, "toll") |
      str_detect(annote, "lymphocyte") |
      str_detect(annote, "CD9") |
      str_detect(annote, "HDD11") |
      str_detect(annote, "deleted in malignant brain")|
      str_detect(annote, "lectin lectoxin") |
      str_detect(annote, "lectin L6") 
  )
  

############# Tissue Subset heatmap ##################
tissue_df <- df.sets.ann %>%
  filter(
    str_detect(annote, "collagen") |
      str_detect(annote, "short-chain") |
      str_detect(annote, "fibrinogen") |
      str_detect(annote, "laminin") |
      str_detect(annote, "disintegrin") | 
      str_detect(annote, "fibroblast") |
      str_detect(annote, "integrin") | 
      str_detect(annote, "latrophilin") | 
      str_detect(annote, "ADAMTS") |
      str_detect(annote, "cadherin") |
      str_detect(annote, "footprint") |
      str_detect(annote, "extracellular matrix") |
      str_detect(annote, "fibronectin") |
      str_detect(annote, "proteoglycan") |
      str_detect(annote, "metalloproteinase") |
      str_detect(annote, "adhesion")
  )

ecm_df <- df.sets.ann %>%
  filter(
    str_detect(annote, "collagen") |
    str_detect(annote, "fibronectin") |
    str_detect(annote, "laminin") |
    str_detect(annote, "proteoglycan 4") |
    str_detect(annote, "extracellular matrix") 
    #str_detect(annote, "footprint") 
    
  )

rem.rep_df <- df.sets.ann %>%
  filter(
    str_detect(annote, "leukocyte elastase inhibitor") |
    str_detect(annote, "fibrinogen") |
    str_detect(annote, "vascular endothelial growth factor receptor") |
    str_detect(annote, "serine-rich adhesin for platelets-like") |
    str_detect(annote, "metalloproteinase") |
    str_detect(annote, "ADAMTS") 
  )

adh.sig_df <- df.sets.ann %>%
  filter(
    str_detect(annote, "integrin alpha-9") |
    str_detect(annote, "adhesion") |
    str_detect(annote, "cadherin") |
    str_detect(annote, "latrophilin") |
    str_detect(annote, "fibroblast growth factor receptor") 
    # str_detect(annote, "rho GTPase-activating protein") 
  )

############# neural Subset heatmap ##################
neural_df <- df.sets.ann %>%
  filter(
    str_detect(annote, "acetylcholine") |
      str_detect(annote, "synaptotagmin") |
      str_detect(annote, "synaptic") |
      str_detect(annote, "G-protein") |
      str_detect(annote, "Acetylcholinesterase") | 
      str_detect(annote, "Dopamine") |
      str_detect(annote, "GABA") | 
      str_detect(annote, "histamine") | 
      str_detect(annote, "neuroendocrine") | 
      str_detect(annote, "voltage-dependent") |
      str_detect(annote, "voltage-gated") |
      str_detect(annote, "calcium channel") |
      str_detect(annote, "NMDA") |
      str_detect(annote, "glutamate receptor") |
      str_detect(annote, "neuro") |
      str_detect(annote, "neurotrypsin") |
      str_detect(annote, "synapse") |
      str_detect(annote, "SNARE") |
      str_detect(annote, "neuronal ")
  )

neural_df_select <- df.sets.ann %>%
  filter(
      str_detect(annote, "synaptotagmin") |
      str_detect(annote, "synaptic") |
      str_detect(annote, "G-protein") |
      str_detect(annote, "GABA") | 
      str_detect(annote, "FRMF") | 
      str_detect(annote, "histamine") | 
      str_detect(annote, "calcium channel") |
      str_detect(annote, "NMDA") |
      str_detect(annote, "glutamate receptor") |
      str_detect(annote, "neurotrypsin") |
      str_detect(annote, "SNARE") 
  )

############# Stress and Detox Subset heatmap ##################
stress_df <- df.sets.ann %>%
  filter(
    str_detect(annote, "superoxide") | # Oxidative Stress and Antioxidant Defense
    str_detect(annote, "thioredoxin") |
    str_detect(annote, "NADH") |
    str_detect(annote, "NADPH oxidase") |
    str_detect(annote, "xanthine") |
    str_detect(annote, "peroxiredoxin") |
    str_detect(annote, "Glutathione S-transferase") | 
    str_detect(annote, "Glutaredoxin") |
    str_detect(annote, "Heat shock") | # Heat shock and chaperones
    str_detect(annote, "hsp") |
    str_detect(annote, "Proteasome") |
    str_detect(annote, "p450") | # Detoxification and Xenobiotic Metabolism
    str_detect(annote, "glucuronosyltransferase") |
    str_detect(annote, "sulfotransferase") |
    str_detect(annote, "carboxylesterase") |
    str_detect(annote, "ATP-binding cassette") |
    str_detect(annote, "repair") | # DNA Repair and Damage
    str_detect(annote, "E3 ubiquitin") 
    )

stress_df_select <- df.sets.ann %>%
  filter(
    str_detect(annote, "superoxide") | # Oxidative Stress and Antioxidant Defense
      str_detect(annote, "thioredoxin") |
      str_detect(annote, "NADPH oxidase") |
      str_detect(annote, "xanthine") |
      str_detect(annote, "peroxiredoxin") |
      str_detect(annote, "Glutathione S-transferase") | 
      str_detect(annote, "Glutaredoxin") |
      str_detect(annote, "Heat shock") | # Heat shock and chaperones
      str_detect(annote, "hsp") |
      str_detect(annote, "Proteasome") |
      str_detect(annote, "p450") | # Detoxification and Xenobiotic Metabolism
      str_detect(annote, "glucuronosyltransferase") |
      str_detect(annote, "sulfotransferase") |
      str_detect(annote, "carboxylesterase") |
      str_detect(annote, "repair") # DNA Repair and Damage
  )

########### Plot ##############
library(gridExtra)

heat.plot <- function(df,title,lfc_cut=0.0) {
  
  
  degs_cut <- degs[abs(degs$log2FoldChange) >= lfc_cut,]
  
  df <- df[df$Gene %in% degs_cut$gene_id,]
  
  ##########.
  
  df <- df %>%
    group_by(Gene) %>%
    #mutate(norm_expression = (expression - min(expression)) / (max(expression) - min(expression)))
    mutate(norm_expression = 2 * ((expression - min(expression)) / (max(expression) - min(expression))) - 1)
  
    
  # Pivot to wide format (genes as rows, samples as columns)
  df_wide <- df %>%
    dplyr::select(sample, Gene, norm_expression) %>%
    spread(key = sample, value = norm_expression)
  
  # Create a vector of custom labels for the rows (annote_clean)
  row_labels <- df %>%
    distinct(Gene, annote_clean) %>%
    arrange(match(Gene, df_wide$Gene)) %>%
    pull(annote_clean)
  
  # Create a vector for column annotation
  sample_labels <- ifelse(colnames(df_wide)[-1] %in% df$sample[df$group == "Naive"], "Naive", "Exposed")
  
  
  # Convert the expression data (wide format) to a matrix (excluding 'unique_label' column)
  df_matrix <- as.matrix(df_wide[, -1])  # Exclude 'unique_label' column
  
  col_anno <- data.frame(group = sample_labels)
  rownames(col_anno) <- colnames(df_matrix) 
  
  annotation_colors <- list(group = c("Naive" = "skyblue", "Exposed" = "salmon"))
  
  p <- pheatmap(df_matrix, 
                cellwidth = 15, 
                cellheight = 15, 
                border_color = NA,
                cluster_rows = TRUE,
                cluster_cols = FALSE,
                scale = "none",
                color = colorRampPalette(c("mediumseagreen","white","purple3"))(100),
                show_rownames = TRUE,
                show_colnames = FALSE,
                labels_row = row_labels,
                treeheight_row = 0,
                treeheight_col = 0,
                fontsize = 12,
                fontfamily = "Arial"
  )
  
  #library(gridExtra)
  
  #grid.arrange(i$gtable, t$gtable, ncol=2)
  
  return(p)
  
}


#immune plots
comp <- heat.plot(complement_df,"complement cascade")
reg <- heat.plot(inflam_df,"immune regulation and cell death")
rec <- heat.plot(recog_df,"pathogen recognition and defense")

grid.arrange(comp$gtable, rec$gtable,reg$gtable, ncol=1)



#tissue plots
ecm <- heat.plot(ecm_df,"ECM")
adh <- heat.plot(adh.sig_df,"singal and adhesion")
rem <- heat.plot(rem.rep_df,"remodeling genes")

grid.arrange(comp$gtable, reg$gtable, rec$gtable, ncol=3)






#########.

df <- complement_df
#### pick to separate or combine like-genes

### optional filter for lfc > set paramiter #####

lfc_cut = 0.0

degs_cut <- degs[abs(degs$log2FoldChange) >= lfc_cut,]

df <- df[df$Gene %in% degs_cut$gene_id,]

##########.

df <- df %>%
  group_by(Gene) %>%
  mutate(norm_expression = (expression - min(expression)) / (max(expression) - min(expression)))

# Pivot to wide format (genes as rows, samples as columns)
df_wide <- df %>%
  dplyr::select(sample, Gene, norm_expression) %>%
  spread(key = sample, value = norm_expression)

# Create a vector of custom labels for the rows (annote_clean)
row_labels <- df %>%
  distinct(Gene, annote_clean) %>%
  arrange(match(Gene, df_wide$Gene)) %>%
  pull(annote_clean)

# Create a vector for column annotation
sample_labels <- ifelse(colnames(df_wide)[-1] %in% df$sample[df$group == "Naive"], "Naive", "Exposed")


# Convert the expression data (wide format) to a matrix (excluding 'unique_label' column)
df_matrix <- as.matrix(df_wide[, -1])  # Exclude 'unique_label' column

col_anno <- data.frame(group = sample_labels)
rownames(col_anno) <- colnames(df_matrix) 

annotation_colors <- list(group = c("Naive" = "skyblue", "Exposed" = "salmon"))

pheatmap(df_matrix, 
         border_color = "grey20",
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         scale = "none",
         color = colorRampPalette(c("springgreen3", "white", "purple"))(100),
         show_rownames = TRUE,
         show_colnames = FALSE,
         main = "Immune Response",
         labels_row = row_labels,
         treeheight_row = 0,
         treeheight_col = 0,
         annotation_col = col_anno,
         annotation_colors = annotation_colors,
         fontsize=9
         )

library(gridExtra)

grid.arrange(i$gtable, t$gtable, ncol=2)



