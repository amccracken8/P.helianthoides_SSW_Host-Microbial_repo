# maing plots from GO MWU output of Salmon Module from WGNCA
#Top correlated genes in module

library(tidyverse)
library(data.table)
library(patchwork)

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
rownames(dataExpr)
dataExpr$group <- case_when(
  dataExpr$V1 %in% dataExpr$V1[1:5] ~ "Naive",
  dataExpr$V1 %in% dataExpr$V1[6:10] ~ "Exposed"
)

### select all genes in red module that are also DEGs for plotting
red_deg_genes <- read.csv("red_mod_degs_nf.csv")

redgenes <- read_lines("red_mod_genes_nf.txt")

dataExpr_reddeg <- dataExpr%>% dplyr::select(c("V1", "group", redgenes))



rename_vector <- setNames(degs_GO$annote, degs_GO$gene_id)
colnames(dataExpr_reddeg) <- rename_vector[colnames(dataExpr_reddeg)]
colnames(dataExpr_reddeg)<-make.unique(colnames(dataExpr_reddeg))



###### Plotting Module Genes #######

rgenes <- dataExpr_reddeg %>%
  dplyr::select(V1, contains("ficolin-2"), contains("echinoidin"),contains("collagenase"),contains("C-type lectin lectoxin"),contains("collagen alpha"),contains("fibrinogen"), group) %>%
  pivot_longer(cols = -c(V1, group),  # Exclude 'V1' and 'group' columns
               names_to = "gene",
               values_to = "expression") %>%
  mutate(gene_type = case_when(
    str_detect(gene, "ficolin") ~ "ficolin-2",
    str_detect(gene, "echinoidin") ~ "echinoidin",
    str_detect(gene,"collagen alpha") ~ "collagen alpha-1(XVII)",
    str_detect(gene,"collagenase") ~ "type IV collagenase",
    #str_detect(gene,"lymphocyte") ~ "lymphocyte antigen 6",
    str_detect(gene,"C-type lectin lectoxin-Phi1-like") ~ "C-type lectin lectoxin-Phi1",
    #str_detect(gene,"adhesion") ~ "adhesion G protein-coupled receptor",
    str_detect(gene,"fibrinogen") ~ "fibrinogen",
    #str_detect(gene,"ADAMTS") ~ "ADAMTS",
    #str_detect(gene,"disintegrin") ~ "ADAMTS",
    TRUE ~ "other"  # Optional: handle other cases if needed
  ))



rgenes$group <- factor(rgenes$group, levels = c("Naive", "Exposed"))
rgenes$gene_type <- stringr::str_wrap(rgenes$gene_type, width = 15) 


unique(rgenes$gene_type) # catch special invisible characters
rgenes$gene_type <- factor(rgenes$gene_type, levels = c("echinoidin", "ficolin-2", "C-type lectin\nlectoxin-Phi1", "collagen\nalpha-1(XVII)","type IV\ncollagenase","fibrinogen"))


rg <- ggplot(rgenes, aes(x = group, y = expression, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = c("Exposed" = "lightcoral", "Naive" = "skyblue2")) + 
  # Points for each sub-transcript mean
  stat_summary(aes(group = gene, alpha=0.5),
               fun = mean, geom = "point", size = 1.5, position = position_dodge(width = 0.5)) +
  stat_summary(aes(group = gene, alpha=0.5), 
               fun = mean, 
               geom = "line", 
               position = position_dodge(width = 0.5)) +
  facet_wrap(~ gene_type, scales = "free_x", ncol = 6, strip.position = "bottom") +
  labs(title = "Red Module Genes",
       x = "Group",
       y = "log2 normalized Expression") +
  theme_minimal() +
  theme(axis.text.x = element_blank(), 
        strip.placement = "outside",  
        strip.text.x = element_text(size = 10, angle = 45, hjust = 0.5, face = "bold"), 
        axis.title.x = element_blank(),
        panel.spacing = unit(0.5, "lines"), 
        panel.border = element_rect(color = "grey", fill = NA, size = 1)) + 
  guides(color = "none",alpha="none")









########### Plotting Microbes ##############


# pull out genes of interest to plot
red_mico <- dataExpr %>% dplyr::select("V1", "d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Enterobacterales.f__Vibrionaceae.g__Vibrio.s__","d__Bacteria.p__Bacteroidota.c__Bacteroidia.o__Flavobacteriales.f__Flavobacteriaceae.g__Flavobacterium.s__","d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Enterobacterales.f__Moritellaceae.g__Moritella.s__","d__Bacteria.p__Firmicutes.c__Clostridia.o__Peptostreptococcales.Tissierellales.f__Peptostreptococcales.Tissierellales.g__JTB215.s__","d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rhizobiales.f__Rhizobiaceae.g__Pseudahrensia.s__","d__Bacteria.p__Fusobacteriota.c__Fusobacteriia.o__Fusobacteriales.f__Fusobacteriaceae.g__Psychrilyobacter.s__","group")

# rename and shorten microbe names
colnames(red_mico) = c("V1", "Vibrio.spp","Flavobacterium.spp","Moritella.spp","JTB215.spp","Pseudahrensia.spp","Psychrilyobacter.spp","group")


# elongate df to plot
red_long <- pivot_longer(red_mico, cols = 2:7, names_to = "microbe", values_to = "Value")

red_long$group <- factor(red_long$group, levels = c("Naive", "Exposed"))


# wrap text of long names
red_long$Gene <- stringr::str_wrap(red_long$microbe, width = 15) 

red_long$microbe <- factor(red_long$microbe, levels = c("Vibrio.spp","Flavobacterium.spp","Moritella.spp","JTB215.spp","Pseudahrensia.spp","Psychrilyobacter.spp"))


# Create the boxplot
rm <- ggplot(red_long, aes(x = group, y = Value, fill = group)) +
  geom_boxplot() +
  facet_grid(. ~ microbe, switch = "x") +  # Move gene names to the bottom
  scale_fill_manual(values = c("Exposed" = "lightcoral", "Naive" = "skyblue2")) +  # Match exact values from group column
  labs(title = "Red Module Microbes",
       y = "log2 normalized counts") +  # Remove the x-axis title
  theme_minimal() +
  theme(axis.text.x = element_blank(),  # remove x axis lables
        strip.placement = "outside",  # Place facet labels outside
        strip.text.x = element_text(size = 10,angle =50,face = "bold"),  # Adjust facet label size
        axis.title.x = element_blank(),  # Remove x-axis title
        panel.spacing = unit(0.5, "lines"), 
        panel.border = element_rect(color = "grey", fill = NA, size = 1)) # Add black borders between facets
 

red.plot <- (rg/rm)+
  plot_layout(guides = "collect",axes = "collect") 





############## new plots #####################


library(ggplot2)
library(ggdist)   # for half-violin geom
library(dplyr)

g <- ggplot(rgenes, aes(x = group, y = expression, fill = group, color = group)) +
  # # Half violin (distribution)
  # ggdist::stat_halfeye(
  #   adjust = 0.6,
  #   width = 0.8,
  #   justification = -0.3,  # shift left/right
  #   .width = 0,
  #   point_colour = NA,
  #   alpha = 0.8
  # ) +
  # "Rain" points
  geom_jitter(
    aes(color = group),
    width = 0.15,
    alpha = 0.8,
    size = 1.5
  ) +
  # Optional boxplot overlay for summary
  geom_boxplot(
    width = 0.75,
    outlier.shape = NA,
    alpha = 0.70,
    color = "grey1",
    linewidth = 0.4
  ) +
  scale_fill_manual(values = c("Exposed" = "lightcoral", "Naive" = "skyblue2")) +
  scale_color_manual(values = c("Exposed" = "firebrick3", "Naive" = "steelblue4")) +
  facet_wrap(~ gene_type, scales = "free_x", ncol = 6, strip.position = "bottom") +
  labs(title = "Red Module Genes",
       x = "Group",
       y = "log2 normalized Expression") +
  theme_minimal() +
  theme(
    #panel.grid = element_blank(),      
    axis.text.x = element_blank(),
    strip.placement = "outside",
    strip.text.x = element_text(size = 10, angle = 0, hjust = 0.5, face = "bold"),
    axis.title.x = element_blank(),
    panel.spacing = unit(0.5, "lines"),
    panel.border = element_rect(color = "grey", fill = NA, size = 0.5),
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5)
  ) +
  guides(fill = guide_legend(position = "bottom"))
g 




r <- ggplot(red_long, aes(x = group, y = Value, fill = group, color = group)) +
  # Half violin for distribution
  # ggdist::stat_halfeye(
  #   adjust = 1.0,
  #   width = 1.0,
  #   justification = -0.3, 
  #   .width = 0,
  #   point_colour = NA,
  #   alpha = 0.8
  # ) +
  # "Rain" points
  geom_jitter(
    aes(color = group),
    width = 0.15,
    alpha = 0.8,
    size = 1.5
  ) +
  # Optional boxplot overlay
  geom_boxplot(
    width = 0.75,
    outlier.shape = NA,
    alpha = 0.70,
    color = "grey1",
    linewidth = 0.4
  ) +
  scale_fill_manual(values = c("Exposed" = "lightcoral", "Naive" = "skyblue2")) +
  scale_color_manual(values = c("Exposed" = "firebrick3", "Naive" = "steelblue4")) +
  facet_grid(. ~ microbe, switch = "x") +  # keep microbes along x-axis
  labs(title = "Red Module Microbes",
       y = "log2 normalized counts") +
  theme_minimal() +
  theme(
    #panel.grid = element_blank(),    
    axis.text.x = element_blank(),
    strip.placement = "outside",
    strip.text.x = element_text(size = 10, angle = 0, face = "bold.italic"),
    axis.title.x = element_blank(),
    panel.spacing = unit(0.5, "lines"),
    panel.border = element_rect(color = "grey", fill = NA, size = 0.5),
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5)
  ) +
  guides(fill = guide_legend(position = "bottom"))
r

red.plot <- (g / r) +
  plot_layout(guides = "collect", axes = "collect") &
  theme(legend.position = "bottom")
