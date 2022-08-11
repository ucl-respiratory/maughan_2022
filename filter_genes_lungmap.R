# Load LungMAP genes
library(gdata)
library(pheatmap)
lungmap <- read.table('data/LungMAP/LMEX0000003691_normalized_counts.txt')
lungmap.pheno <- read.csv('data/LungMAP/LMEX0000003691_sample_metadata.csv')

# Tidy up data a bit
lungmap.pheno$Sample.Inventory.ID <- make.names(lungmap.pheno$Sample.Inventory.ID)
colnames(lungmap) <- gsub('_', '.', colnames(lungmap))
samps <- intersect(colnames(lungmap), lungmap.pheno$Sample.Inventory.ID)
lungmap.pheno <- lungmap.pheno[match(samps, lungmap.pheno$Sample.Inventory.ID),]
lungmap <- lungmap[,samps]
all(colnames(lungmap) == lungmap.pheno$Sample.Inventory.ID)

# Divide into age cohorts
lungmap.1 <- lungmap[,which(lungmap.pheno$Age.Cohort == '1. Neonate (up to 30 days)')]
lungmap.2 <- lungmap[,which(lungmap.pheno$Age.Cohort == '2. Infant (> 30 days and < 1 year)')]
lungmap.3 <- lungmap[,which(lungmap.pheno$Age.Cohort == '3. Child (>= 1 year and < 11 years)')]
lungmap.4 <- lungmap[,which(lungmap.pheno$Age.Cohort == '5. Adult (20+ years)')]

# Danaher-like functions:
sim <- function(x,y) {
  s = sum((x - mean(x)) * (y - mean(y))) / (((length(x) - 1)/2) * (var(x) + var(y)))
  return(s)
}

check.similarity <- function(candidates, dataset) {
  candidates <- candidates[which(candidates %in% rownames(dataset))]
  n <- length(candidates)
  m <- matrix(NA, nrow=n, ncol=n)
  for(i in 1:n) {
    for(j in 1:n) {
      m[i,j] <- sim(as.numeric(dataset[candidates[i],]), as.numeric(dataset[candidates[j],]))
    }
  }
  rownames(m) <- candidates
  colnames(m) <- candidates
  
  # Plot the heatmap and extract the obviously co-correlated cluster:
  # hm <- pheatmap(m, display_numbers = T, fontsize_number = 5, fontsize = 6)
  hm <- pheatmap(m, display_numbers = F,  fontsize_number = 5, fontsize = 6)
  
  return(hm)
}

# Gene lists to test:
load('data/travaglini_airway_marker_candidates.RData')

# Test them:

######################################################################
# Basal
######################################################################
hm <- check.similarity(genes.basal, lungmap.1)
c <- cutree(hm$tree_col, k=2)
genes.basal.filt.1 <- names(c[which(c==2)])
# Neonatal - proliferation markers dominate here

hm <- check.similarity(genes.basal, lungmap.2)
c <- cutree(hm$tree_col, k=3)
genes.basal.filt.2 <- names(c[which(c==1)])
a1 <- as.grob(hm)
# Infant - two clusters, proliferation and cell markers

hm <- check.similarity(genes.basal, lungmap.3)
c <- cutree(hm$tree_col, k=2)
genes.basal.filt.3 <- names(c[which(c==1)])
a2 <- as.grob(hm)
# Child - merged to give a larger cluster of proliferation and cell markers

hm <- check.similarity(genes.basal, lungmap.4)
c <- cutree(hm$tree_col, k=2)
genes.basal.filt.4 <- names(c[which(c==1)])
a3 <- as.grob(hm)
# Adult - two clusters again seen

# Overall - ignore the neonates as they seem quite different
# Merge others
genes.basal.filt <- intersect(
  genes.basal.filt.2,
  intersect(genes.basal.filt.3, genes.basal.filt.4)
)


######################################################################
# Secretory
######################################################################
hm <- check.similarity(genes.secretory, lungmap.1)
c <- cutree(hm$tree_col, k=2)
genes.secretory.filt.1 <- names(c[which(c==2)])
# Neonatal - proliferation markers dominate here

hm <- check.similarity(genes.secretory, lungmap.2)
c <- cutree(hm$tree_col, k=3)
genes.secretory.filt.2 <- names(c[which(c==1)])
# Infant - two clusters, proliferation and cell markers
b1 <- as.grob(hm)

hm <- check.similarity(genes.secretory, lungmap.3)
c <- cutree(hm$tree_col, k=2)
genes.secretory.filt.3 <- names(c[which(c==1)])
b2 <- as.grob(hm)
# Child - merged to give a larger cluster of proliferation and cell markers

hm <- check.similarity(genes.secretory, lungmap.4)
c <- cutree(hm$tree_col, k=2)
genes.secretory.filt.4 <- names(c[which(c==1)])
b3 <- as.grob(hm)
# Adult - two clusters again seen

# Overall - ignore the neonates as they seem quite different
# Merge others
genes.secretory.filt <- intersect(
  genes.secretory.filt.2,
  intersect(genes.secretory.filt.3, genes.secretory.filt.4)
)


######################################################################
# Ciliated
######################################################################
hm <- check.similarity(genes.ciliated, lungmap.1)
c <- cutree(hm$tree_col, k=3)
genes.ciliated.filt.1 <- names(c[which(c==3)])
# Neonatal - proliferation markers dominate here

hm <- check.similarity(genes.ciliated, lungmap.2)
c <- cutree(hm$tree_col, k=2)
genes.ciliated.filt.2 <- names(c[which(c==1)])
c1 <- as.grob(hm)
# Infant - two clusters, proliferation and cell markers

hm <- check.similarity(genes.ciliated, lungmap.3)
c <- cutree(hm$tree_col, k=2)
genes.ciliated.filt.3 <- names(c[which(c==1)])
c2 <- as.grob(hm)
# Child - merged to give a larger cluster of proliferation and cell markers

hm <- check.similarity(genes.ciliated, lungmap.4)
c <- cutree(hm$tree_col, k=2)
genes.ciliated.filt.4 <- names(c[which(c==1)])
c3 <- as.grob(hm)
# Adult - two clusters again seen

# Overall - ignore the neonates as they seem quite different
# Merge others
genes.ciliated.filt <- intersect(
  genes.ciliated.filt.2,
  intersect(genes.ciliated.filt.3, genes.ciliated.filt.4)
)


save(
  genes.basal.filt,
  genes.secretory.filt,
  genes.ciliated.filt,
  file = 'data/filtered_gene_lists.RData'
)
