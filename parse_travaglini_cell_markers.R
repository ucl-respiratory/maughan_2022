# Make gene lists for use in the Danaher co-correlation method
# Take from Travaglini et al 2020 - uses supp. table 4 from https://www.nature.com/articles/s41586-020-2922-4#Sec32
# Groups (from figure 1):
# - Basal: 4,5,6,7
# - Secretory: 1,8,9
# - Ciliatory: 2,3
#
# The paper uses both 10x and smart-seq 2 - use only those identified with BOTH methods (more stringent)
# Require significant p value AND positive fold change (we only care about upregulated genes as candidate markers)
library(gdata)
f <- 'data/Travaglini_tableS4.xls'
sheets <- sheetNames(f)
getCluster <- function(c, cutoff = 0.01) {
  print(paste('Reading cluster', c))
  c <- paste0('Cluster ', c)
  x1 <- read.xls(f, sheet = c, skip = 1, stringsAsFactors = F)
  genes1 <- x1$Gene[which(x1$p_val_adj < cutoff & x1$avg_logFC > 1)]
  
  s <- paste0(c, ' (SS2)')
  if(s %in% sheets) {
    x2 <- read.xls(f, sheet = s, skip = 1, stringsAsFactors = F)
    genes2 <- x2$Gene[which(x2$p_val_adj < cutoff & x2$avg_logFC > 1)]
  } else {
    genes2 <- genes1
  }
  
  
  return(unique(intersect(genes1, genes2)))
}

genes.basal <- unique(c(
  getCluster(4), getCluster(5),getCluster(6), getCluster(7)
))
genes.secretory <- unique(c(
  getCluster(1), getCluster(8),getCluster(9)
))
genes.ciliated <- unique(c(
  getCluster(2), getCluster(3)
))

save(genes.basal, genes.secretory, genes.ciliated, file = "data/travaglini_airway_marker_candidates.RData")
