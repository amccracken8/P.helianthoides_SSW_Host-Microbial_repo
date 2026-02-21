library(microbiome)
library(phyloseq)
library(qiime2R)
library(ANCOMBC)
library(knitr)
library(kableExtra)


#Use table-lev7.qza for species and table.qza for ASV
physeq_fam<-qza_to_phyloseq(
  features="table-lev7.qza",
  tree=,
  taxonomy=,
  metadata = "pyc_sylva_manifest.txt")

# Read in manifest
meta <- read.table("pyc_sylva_manifest.txt", header=T)
Manifest <- read.table("pycno_samples_NE.txt", header = TRUE)

# suset metadata by desired manifest file
meta_sub <- meta[meta$sampleID %in% Manifest$manifestID,]
  

#Subset samples, just healthy animals (inclusive of Naive and Exposed samples)
pseq_healthy <- subset_samples(physeq_fam, animal.health == "healthy")

#Subset samples, just impacted sites (inclusive of Exposed and Wasting samples)
pseq_impacted <- subset_samples(physeq_fam, site.status == "impacted")

#subset to samples desired
pseq_sub <- subset_samples(physeq_fam, animalID %in% meta_sub$animalID)

#ANCOMBC for site-animal-health Naive vs Exposed (Taxa-lev7 : Species)
#Compared by Site-animal-health status as (HH=Naive, SH=Exposed)
out.NvE_Sylva <- ancombc(
  phyloseq = pseq_sub, 
  group = "site.animal.health",
  formula = "site.animal.health", #dots not hyphens for file name?
  p_adj_method = "fdr", 
  #zero_cut = 0.9, # by default prevalence filter of 10% is applied, this causes us to lose too many taxa that aren't in many of the samples
  lib_cut = 0,
  struc_zero = TRUE, 
  neg_lb = TRUE, 
  tol = 1e-5, 
  max_iter = 100, 
  conserve = TRUE, 
  alpha = 0.05, 
  global = TRUE
)
res.NvE_Sylva <- as.data.frame(out.NvE_Sylva$res)
#colnames(resSHtoSS) <- c("Exp_beta", "Exp_se", "Exp_W", "Exp_p_val", "Exp_q_val", "Exp_DA")

#write.csv(as.data.frame(res.NvE_Sylva),file="8.1.24_res.NvE_Sylva.csv")

