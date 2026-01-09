### exploring close relationships in the RED module

setwd("C:/Users/andre/Desktop/GitHub/AK_Pycno_RNAxMetagen/16s_Sylva/qiime2_fullset/WGCNA/Networkplot")

library('igraph')
library('network')
library('networkD3')
library('intergraph')
library("tidyverse")

# adjacency matrix must be pre-loaded from WGCNA module before running code. Saving the file as a matrix corrupted it so using the raw object from the WGCNA code up to building the adjacency matrix is necessary before running code below. 
moduleColors <- read_lines("moduleColors_nf.txt")

RedGenes = names(datExpr)[moduleColors == "red"] 

adjacencyModule = adjacency[RedGenes, RedGenes]

# picking out target connections to vibrio directly within the "red" module
targetConnections = adjacencyModule["d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Enterobacterales.f__Vibrionaceae.g__Vibrio.s__", ]

topConnections = sort(targetConnections, decreasing = TRUE)

topcons = names(topConnections)[1:30]

ann.df <- read.csv("C:/Users/andre/Desktop/GitHub/AK_Pycno_RNAxMetagen/Pycno_NE_STAR/annotations/NE_ALL_annotated.csv")




###### Igrah top connectons########

adjacencySubset = adjacencyModule[topcons, topcons]

g <- graph.adjacency(as.matrix(adjacencySubset), mode = "undirected", weighted = TRUE, diag = FALSE)

# Node degree (number of connections, considering edge weights)
nodeDegree <- strength(g, mode = "all", weights = E(g)$weight)


# Color target gene red, others light blue
nodeColors <- ifelse(V(g)$name == "d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Enterobacterales.f__Vibrionaceae.g__Vibrio.s__", "red", "lightblue")

# Plot the network
plot(
  g,
  vertex.size = nodeDegree / max(nodeDegree) * 15 + 3,  # Scale size by degree
  vertex.label.cex = 0.7,
  vertex.color = nodeColors,
  edge.width = E(g)$weight * 5,  # Scale edge thickness by weight
  layout = layout_with_fr  # Fruchterman-Reingold layout (good for small networks)
)



################ genes of interest

#####
ann.df <- read.csv("C:/Users/andre/Desktop/GitHub/AK_Pycno_RNAxMetagen/Pycno_NE_STAR/annotations/NE_ALL_annotated.csv")

ann.df <- ann.df %>% select(gene_id,annote)

ann.df <- ann.df %>%
  mutate(
    ann_short = str_remove(annote, "^\\S+\\s"),             
    ann_short = str_remove(ann_short, "\\s\\[.*?\\]$")
  )

# there are far too many genes to plot and it is a mess. here is a subset of key terms pulled from the annotation column searching for specific genes of interest within the red module
gen.of.int <- ann.df %>%
  filter(
    str_detect(annote, "ficolin") |
      str_detect(annote, "echinoidin") |
      str_detect(annote, "complement") |
      str_detect(annote, "scavenger") |
      str_detect(annote, "tumor necrosis factor receptor") |
      str_detect(annote, "TNF") |
      str_detect(annote, "interleukin") |
      str_detect(annote, "death") |
      str_detect(annote, "leukocyte") |
      str_detect(annote, "macrophage") |
      str_detect(annote, "t-cell") |
      str_detect(annote, "toll") |
      str_detect(annote, "lymphocyte") |
      str_detect(annote, "CD9") |
      str_detect(annote, "HDD11") |
      str_detect(annote, "deleted in malignant brain")|
      str_detect(annote, "lectin lectoxin") |
      str_detect(annote, "lectin L6") |
      str_detect(annote, "techylectin") |

      str_detect(annote, "collagen") |
      str_detect(annote, "fibronectin") |
      str_detect(annote, "laminin") |
      str_detect(annote, "proteoglycan 4") |
      str_detect(annote, "extracellular matrix") |
      str_detect(annote, "leukocyte elastase inhibitor") |
      str_detect(annote, "fibrinogen") |
      str_detect(annote, "vascular endothelial growth factor receptor") |
      str_detect(annote, "serine-rich adhesin for platelets-like") |
      str_detect(annote, "metalloproteinase") |
      str_detect(annote, "ADAMTS") |
      str_detect(annote, "integrin alpha-9") |
      str_detect(annote, "adhesion") |
      str_detect(annote, "cadherin") |
      str_detect(annote, "latrophilin") |
      str_detect(annote, "footprint") |
      str_detect(annote, "fibroblast growth factor receptor") |

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
      str_detect(annote, "neuronal ") |

      str_detect(annote, "superoxide") | # Oxidative Stress and Antioxidant Defense
      str_detect(annote, "thioredoxin") |
      str_detect(annote, "NADH") |
      str_detect(annote, "NADPH oxidase") |
#      str_detect(annote, "xanthine") |
      str_detect(annote, "peroxiredoxin") |
      str_detect(annote, "Glutathione S-transferase") |
      str_detect(annote, "Glutaredoxin") |
      str_detect(annote, "Heat shock") | # Heat shock and chaperones
      str_detect(annote, "hsp") |
      str_detect(annote, "Proteasome") |
      str_detect(annote, "p450") | # Detoxification and Xenobiotic Metabolism
      str_detect(annote, "P450") |
      str_detect(annote, "glucuronosyltransferase") |
      str_detect(annote, "sulfotransferase") |
#      str_detect(annote, "carboxylesterase") |
      str_detect(annote, "ATP-binding cassette") |
      str_detect(annote, "repair") # DNA Repair and Damage
  )

#####

target_gene <- c("d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Enterobacterales.f__Vibrionaceae.g__Vibrio.s__")

genes_of_interest <- c(gen.of.int$gene_id)

filteredConnections <- targetConnections[names(targetConnections) %in% genes_of_interest]

filteredConnections <- sort(filteredConnections, decreasing = TRUE)

topcons <- unique(c(target_gene, names(filteredConnections)[1:15]))

id_to_annot <- setNames(ann.df$ann_short, ann.df$gene_id)


adjacencySubset <- adjacencyModule[topcons, topcons]
adjacencySubset[adjacencySubset <= 0.5] <- 0


rownames(adjacencySubset) <- ifelse(
  rownames(adjacencySubset) %in% names(id_to_annot),
  id_to_annot[rownames(adjacencySubset)],
  rownames(adjacencySubset)
)

colnames(adjacencySubset) <- ifelse(
  colnames(adjacencySubset) %in% names(id_to_annot),
  id_to_annot[colnames(adjacencySubset)],
  colnames(adjacencySubset)
)



g_subset <- graph.adjacency(as.matrix(adjacencySubset), mode = "undirected", weighted = TRUE, diag = FALSE)


nodeColors <- ifelse(
  V(g_subset)$name == target_gene, # This is the original Vibrio node name
  "red",
  "lightblue"
)

V(g_subset)$name[V(g_subset)$name == "d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Enterobacterales.f__Vibrionaceae.g__Vibrio.s__"] <- "Vibrio spp."


edgeColors <- ifelse(
  E(g_subset)$weight > 0.9, "red",
  ifelse(E(g_subset)$weight >= 0.7, "orange", "grey")
)

node_degree <- degree(g_subset)

#write_graph(g_subset, file = "network.graphml", format = "graphml")


plot(
  g_subset,
  vertex.size = node_degree*2,
  vertex.label.cex = 1.0,
  vertex.color = nodeColors,
  edge.width = E(g_subset)$weight * 5,
  edge.color = edgeColors,
  layout = layout_with_fr(g_subset, niter = 1000, area = 1.5 * vcount(g_subset)) 
)


tkid <- tkplot(
  g_subset,
  vertex.size = node_degree * 2,
  vertex.color = nodeColors,
  vertex.label.color = "black",
  vertex.label.cex = 1.2,
  vertex.label.font = 2,
  edge.width = E(g_subset)$weight * 5,
  edge.color = alpha(edgeColors, 0.5)
)

layout_coords <- tkplot.getcoords(tkid)

plot(
  g_subset,
  layout = layout_coords,
  vertex.size = node_degree * 2,
  vertex.color = nodeColors,
  vertex.label.color = "black",
  vertex.label.cex = 1.0,
  vertex.label.font = 2,
  edge.width = E(g_subset)$weight * 5,
  edge.color = alpha(edgeColors, 0.5)
)





#########  Switch names for  numbers to clean up plot

V(g_subset)$label <- seq_along(V(g_subset))

# Plot the graph with numeric labels
tkid<- tkplot(
  g_subset,
  vertex.size = node_degree * 2,
  vertex.color = nodeColors,
  vertex.label.color = "black",
  vertex.label.cex = 1.0,
  vertex.label.font = 2,
  edge.width = E(g_subset)$weight * 5,
  edge.color = alpha(edgeColors, 0.5),
  layout = layout_with_fr
  
  )

plot(
  g_subset,
  layout = layout_coords,
  vertex.size = node_degree * 2,
  vertex.color = nodeColors,
  vertex.label.color = "black",
  vertex.label.cex = 1.5,
  vertex.label.font = 2,
  edge.width = E(g_subset)$weight * 5,
  edge.color = alpha(edgeColors, 0.5)
)

plot(1, 1, type = "n", xlab = "", ylab = "", axes = FALSE, main = "")


legend_labels <- paste(seq_along(V(g_subset)), V(g_subset)$name, sep = ": ")  # Combine number and name
legend(
  "topright",                          # Position legend (adjust as necessary)
  legend = legend_labels,
  text.col = "black",
  bg = "white",
  cex = 1.0
) 




############ saving and outputing for python ##########

# Assign attributes to vertices
V(g_subset)$color <- nodeColors
V(g_subset)$label <- seq_along(V(g_subset))          # numeric labels
V(g_subset)$degree <- degree(g_subset)              # degree

# Now build a data frame
nodes_df <- data.frame(
  id = seq_along(V(g_subset)),        # numeric ID
  name = V(g_subset)$name,            # node names
  color = V(g_subset)$color,          # assigned color
  degree = V(g_subset)$degree         # node degree
)

E(g_subset)$color <- edgeColors  # assign colors to edges if not already

edges_df <- data.frame(
  from = as.numeric(head_of(g_subset, E(g_subset))),
  to = as.numeric(tail_of(g_subset, E(g_subset))),
  weight = E(g_subset)$weight,
  color = E(g_subset)$color
)

head(edges_df)


# Export to CSV for Python
write.csv(nodes_df, "nodes.csv", row.names = FALSE)
write.csv(edges_df, "edges.csv", row.names = FALSE)
