#############################################################################
# Cell-intrinsic differences between human airway epithelial cells from adults and children
#############################################################################
#
# This analysis pipeline produces figures and supplementary tables for the above publication by Maughan et al.
# It is published to improve reproducibility and transparency.
#
# Differential analysis of two datasets:
# - Unsorted epithelial cells from adults and children
# - Basal cells sorted from adults and children
#
# Our key questions are:
# - Are there marked differences between adults and children? (unbiased analysis)
# - What pathways are involved?
# - Are there apparent differences in cellular composition?
# - Are there differences in a hypothesis-driven specific analysis of COVID-related genes?
#
# In each case, ask these questions first of unsorted epithelial cells then of sorted basal cells
#
# Figures:
#
# Supplementary figures:
#
# Supplementary data:
# - DE genes in unsorted and basal
# - DE pathways in unsorted and basal

#############################################################################
# Load required packages
#############################################################################
libs <- c(
  'ggplot2', 'cowplot', 'ggpubr', 'ggsignif', 'ggplotify', 'ggrepel',
  'pheatmap', 'gdata', 'DESeq2', 'tidyverse', 'fgsea',
  'WriteXLS', 'umap', 'org.Hs.eg.db'
)
for(lib in libs) {
  library(lib, character.only = T)
}
theme_set(theme_cowplot())

#############################################################################
# Analysis settings
#############################################################################
# Set output directory to save files:
version <- '1.2'
# opdir <- paste0("results/v", version, '/')

# For live creating in Dropbox:
opdir <- paste0("~/Dropbox/Lizzie_Paeds_Paper/Cell Reports Revision/adam_figures/v", version, "/")

dir.create(opdir, recursive = T, showWarnings = F)

# Set a seed for reproducibility
set.seed(20)

# Plotting settings
adults.col <- '#D53D20'
paeds.col <- '#14A1D5'
gg.col <- scale_color_manual(values = c('Adult' = adults.col, 'Paediatric' = paeds.col))
gg.fill <- scale_fill_manual(values = c('Adult' = adults.col, 'Paediatric' = paeds.col))
annot.colors <- list(group = c('Adult' = adults.col, 'Paediatric' = paeds.col))

# Figure width - assume double column = 183mm (see https://www.nature.com/nature/for-authors/formatting-guide)
# -> 7.204" (needed for save_plot base_width method)
fig.width <- 7.204
fig.maxheight <- 9.72


# When running in production, we save the versions of these packages for reproducibility:
vfile <- paste0(opdir, "package.versions.csv")
version.data <- installed.packages()[libs,]
write.csv(version.data, file = vfile)

fdr_cutoff <- 0.05
fc_cutoff <- 1.2

#############################################################################
# Load data
#############################################################################

# Load counts data
# This object should contain four objects: counts.unsorted, pheno.unsorted, counts.basal, pheno.basal
# These represent count matrices and metadata for the two datasets: unsorted epithelium and sorted basal cells
# Count matrices should be formatted with gene symbols as rownames and sample names as column names
load('data/counts.cache.newsamples.RData')


# Load some relevant gene sets
# Use Hallmark gene sets from msigdb for broad-brush analysis
# Download from https://www.gsea-msigdb.org/gsea/msigdb/download_file.jsp?filePath=/msigdb/release/7.1/h.all.v7.1.symbols.gmt
gsets.raw <- readLines("data/h.all.v7.1.symbols.gmt")
gsets.hm <- lapply(gsets.raw, function(x) {
  unlist(strsplit(x, '\t'))[-c(1:2)]
})
names(gsets.hm) <- sapply(gsets.raw, function(x) {
  unlist(strsplit(x, '\t'))[1]
})
# Viral response gene lists (expert curated)
genes.viral.raw <- read.xls('data/viral_gene_list.xlsx', sheet = 'Epithelial Viral Infection', stringsAsFactors = F)
genes.viral.raw$Gene.ID <- str_trim(genes.viral.raw$Gene.ID)
genes.viral <- genes.viral.raw[,1]
genes.covid.raw <- read.xls('data/viral_gene_list.xlsx', sheet = 'SARS-CoV-2 Interest', stringsAsFactors = F)
genes.covid.raw$Gene.ID <- str_trim(genes.covid.raw$Gene.ID)
genes.covid <- genes.covid.raw[,1]
gsets.viral <- list(
  'viral' = genes.viral,
  'covid' = genes.covid
)

# Basal, secretory and ciliated marker genes
# genes.basal = read.xls('data/RH_Airway Gene Lists.xlsx', sheet = 'Basal', header=F, stringsAsFactors=F)[,1]
# genes.ciliated = read.xls('data/RH_Airway Gene Lists.xlsx', sheet = 'Ciliated', header=F, stringsAsFactors=F)[,1]
# genes.secretory = read.xls('data/RH_Airway Gene Lists.xlsx', sheet = 'Secretory', header=F, stringsAsFactors=F)[,1]
load('data/travaglini_airway_marker_candidates.RData')

# Filter these genes based on the principle that markers of the same cell type should correlate.
# This method is sufficiently complex that we separate it into a new file.
# source('filter.celltype.genes.R')
# This contains the output of running filter_genes_lungmap.R
load('data/filtered_gene_lists.RData')

# Danaher breakdown:
source('do.danaher.R')

# Exclude genes on sex chromosomes due to gender imbalances:
allgenes <- data.frame(org.Hs.egSYMBOL)
x <- data.frame(org.Hs.egCHR)
allgenes <- merge(allgenes, x, by='gene_id')
allgenes <- allgenes[which(!duplicated(allgenes$symbol)),]
chr.lookup <- allgenes$chromosome
names(chr.lookup) <- allgenes$symbol

# counts.unsorted <- counts.unsorted[which(!(chr.lookup[rownames(counts.unsorted)] %in% c('X', 'Y'))),]
# counts.basal <- counts.basal[which(!(chr.lookup[rownames(counts.basal)] %in% c('X', 'Y'))),]
# counts.cultured <- counts.cultured[which(!(chr.lookup[rownames(counts.cultured)] %in% c('X', 'Y'))),]

# Sex analysis of GTEX data:
# Exclude sex chromosomes
load('data/gtex.lung.cache.RData')
gtex.lung.pheno$SEX <- factor(gtex.lung.pheno$SEX)
sel <- which(!(chr.lookup[rownames(gtex.lung.counts)] %in% c('X','Y')))
dds.sex <- DESeqDataSetFromMatrix(gtex.lung.counts[sel,], gtex.lung.pheno, design = ~SEX)
dds.sex <- DESeq(dds.sex)
resultsNames(dds.sex)
res.sex <- results(dds.sex, name = 'SEX_2_vs_1')
res.sex <- res.sex[order(abs(res.sex$stat), decreasing = T),]
res.sex$chr <- chr.lookup[rownames(res.sex)]
sex.sig.genes <- rownames(res.sex)[which(res.sex$padj < fdr_cutoff & abs(res.sex$log2FoldChange) > fc_cutoff)]


########################################################################
# Dimension reduction
########################################################################
# Consider all samples (unsorted, sorted, cultured) and try PCA and UMAP
# Do this on vst-transformed data
counts.vst <- vst(as.matrix(counts.all))
rownames(counts.vst) <- rownames(counts.all)

pca <- prcomp(t(counts.vst))
x <- pheno.all
x$PC1 <- pca$x[,1]
x$PC2 <- pca$x[,2]
fig <- ggplot(x, aes(x=PC1, y=PC2)) + 
  geom_point(aes(color=group, shape=type)) +
  gg.col
fig
save_plot(paste0(opdir, 'fig.X_PCA.pdf'), fig, base_width = fig.width, base_height = fig.width)
# Debugging plot to identify outlier samples:
ggplot(x, aes(x=PC1, y=PC2, label = sample)) + 
  geom_point(aes(color=group, shape=type)) +
  gg.col +
  geom_label_repel()

# What is PC1?
pc1 <- pca$rotation[,1]
names(pc1) <- rownames(counts.vst)
pc1 <- pc1[order(abs(pc1), decreasing = T)]
head(pc1)
pc2 <- pca$rotation[,2]
names(pc2) <- rownames(counts.vst)
pc2 <- pc1[order(abs(pc2), decreasing = T)]
head(pc2)


# conf <- umap.defaults
# conf$n_neighbors <- 5
u <- umap(t(counts.vst)) 
x$UMAP1 <- u$layout[,1]
x$UMAP2 <- u$layout[,2]
fig <- ggplot(x, aes(x=UMAP1, y=UMAP2)) + 
  geom_point(aes(color=group, shape=type), size=2.5) +
  gg.col
save_plot(paste0(opdir, 'fig.X_UMAP.pdf'), fig, base_width = fig.width, base_height = fig.width)

########################################################################
# Differential analysis of adults vs paeds (unbiased)
########################################################################
# Write a function for running DESeq on each of our count/pheno arrays:
do.deseq <- function(counts, pheno) {
  # Filter unexpressed genes
  # sel <- which(rowSums(counts) > 0)
  
  # Get a normalised counts matrix first
  # Require genes to be expressed in >1 sample:
  sel <- which(apply(counts, 1, function(x) {
    length(which(x > 0)) > 1
  }))
  dds <- DESeqDataSetFromMatrix(counts[sel,], pheno, ~group)
  dds <- estimateSizeFactors(dds)
  counts.norm <- counts(dds, normalized = T)
  
  
  # Exclude genes on X and Y chromosomes from further analysis
  counts <- counts[which(!(chr.lookup[rownames(counts)] %in% c('X', 'Y'))),]
  
  # Require genes to be expressed in >1 sample:
  sel <- which(apply(counts, 1, function(x) {
    length(which(x > 0)) > 1
  }))
  
  dds <- DESeqDataSetFromMatrix(counts[sel,], pheno, ~group)
  dds <- DESeq(dds)
  resultsNames(dds)
  res <- results(dds, name = 'group_Paediatric_vs_Adult')
  res <- res[order(abs(res$stat), decreasing = T),]
  # Note the NAs here represent genes expressed in a single sample
  # Remove these
  res <- res[which(!is.na(res$padj)),]
  res$sex.related <- rownames(res) %in% sex.sig.genes
  
  # Use these stats as 'ranks' to enter into GSEA
  ranks <- res$stat
  names(ranks) <- rownames(res)
  set.seed(20)
  gsea <- fgsea(gsets.hm, ranks, nperm = 1000, maxSize = 500)
  gsea <- gsea[which(gsea$padj < fdr_cutoff),]
  gsea <- data.frame(gsea)
  gsea$leadingEdge <- sapply(gsea$leadingEdge, function(x) {paste(x, collapse=',')})
  
  return(list(
    'dds' = dds,
    'res' = res,
    'gsea' = gsea,
    'norm' = counts.norm
  ))
}

# Run DEseq on each dataset:
x <- do.deseq(counts.unsorted, pheno.unsorted)
dds.unsorted <- x$dds
res.unsorted <- x$res
gsea.unsorted <- x$gsea
counts.unsorted.norm <- x$norm

x <- do.deseq(counts.basal, pheno.basal)
dds.basal <- x$dds
res.basal <- x$res
gsea.basal <- x$gsea
counts.basal.norm <- x$norm

x <- do.deseq(counts.cultured, pheno.cultured)
dds.cultured <- x$dds
res.cultured <- x$res
gsea.cultured <- x$gsea
counts.cultured.norm <- x$norm

# Normalised count matrices:
# counts.unsorted.norm <- counts(dds.unsorted, normalized=TRUE)
# counts.basal.norm <- counts(dds.basal, normalized=TRUE)
# counts.cultured.norm <- counts(dds.cultured, normalized=TRUE)

# Compare first unsorted cells, then sorted basal cells
# Store outputs in a supplementary table


# Annotations for heatmaps:
annot.unsorted <- data.frame(row.names = pheno.unsorted$sample, group = pheno.unsorted$group)
annot.basal <- data.frame(row.names = pheno.basal$sample, group = pheno.basal$group)
annot.cultured <- data.frame(row.names = pheno.cultured$sample, group = pheno.cultured$group)


# Plot a heatmap
x <- counts.unsorted.norm[rownames(res.unsorted)[which(res.unsorted$padj < fdr_cutoff & abs(res.unsorted$log2FoldChange) > fc_cutoff)],]
fig <- pheatmap(x, scale='row', annotation_col = annot.unsorted, fontsize = 6, annotation_colors = annot.colors)
# save_plot(paste0(opdir, "Fig.hm.p0.05.pdf"), as.ggplot(fig), base_width = 10, base_height = 10)

x <- counts.unsorted.norm[rownames(res.unsorted)[which(res.unsorted$padj < 0.01 & abs(res.unsorted$log2FoldChange) > 1.2)],]
fig.2a <- pheatmap(x, scale='row', annotation_col = annot.unsorted, annotation_colors = annot.colors)
# save_plot(paste0(opdir, "Fig.hm.p0.01.pdf"), as.ggplot(fig.2a), base_width = 10, base_height = 10)

########################################################################
# Repeat with sorted basal cells
########################################################################
# sel <- which(rowSums(counts.basal) > 0)
# dds.basal <- DESeqDataSetFromMatrix(counts.basal[sel,], pheno.basal, ~group)
# dds.basal <- estimateSizeFactors(dds.basal)
# dds.basal <- estimateDispersions(dds.basal)
# dds.basal <- nbinomWaldTest(dds.basal)
# 
# # dds.basal <- DESeq(dds.basal)
# resultsNames(dds.basal)
# res.basal <- results(dds.basal, name = 'group_Paediatric_vs_Adult')
# res.basal <- res.basal[order(abs(res.basal$stat), decreasing = T),]
# # Note the NAs here represent genes expressed in a single sample
# # Remove these
# res.basal <- res.basal[which(!is.na(res.basal$padj)),]
# res.basal$sex.related <- rownames(res.basal) %in% sex.sig.genes
# 
# # Use these stats as 'ranks' to enter into GSEA
# # Use KEGG gene sets from msigdb (C2)
# ranks <- res.basal$stat
# names(ranks) <- rownames(res.basal)
# gsea.basal <- fgsea(gsets.hm, ranks, nperm = 1000, maxSize = 500)
# gsea.basal <- gsea.basal[which(gsea.basal$padj < 0.05),]
# gsea.basal <- data.frame(gsea.basal)
# gsea.basal$leadingEdge <- sapply(gsea.basal$leadingEdge, function(x) {paste(x, collapse=',')})
# 
# counts.basal.norm <- counts(dds.basal, normalized=TRUE)



########################################################################
# Repeat with cultured basal cells
########################################################################
# sel <- which(rowSums(counts.cultured) > 0)
# dds.cultured <- DESeqDataSetFromMatrix(counts.cultured[sel,], pheno.cultured, ~group)
# dds.cultured <- DESeq(dds.cultured)
# resultsNames(dds.cultured)
# res.cultured <- results(dds.cultured, name = 'group_Paediatric_vs_Adult')
# res.cultured <- res.cultured[order(abs(res.cultured$stat), decreasing = T),]
# # Note the NAs here represent genes expressed in a single sample
# # Remove these
# res.cultured <- res.cultured[which(!is.na(res.cultured$padj)),]
# res.cultured$sex.related <- rownames(res.cultured) %in% sex.sig.genes
# 
# # Use these stats as 'ranks' to enter into GSEA
# # Use KEGG gene sets from msigdb (C2)
# ranks <- res.cultured$stat
# names(ranks) <- rownames(res.cultured)
# gsea.cultured <- fgsea(gsets.hm, ranks, nperm = 1000, maxSize = 500)
# gsea.cultured <- gsea.cultured[which(gsea.cultured$padj < 0.05),]
# gsea.cultured <- data.frame(gsea.cultured)
# gsea.cultured$leadingEdge <- sapply(gsea.cultured$leadingEdge, function(x) {paste(x, collapse=',')})
# 
# counts.cultured.norm <- counts(dds.cultured, normalized=TRUE)


########################################################################
# All sample analysis
########################################################################
# DEseq using all samples and correcting for type
sel <- which(rowSums(counts.all) > 0)
dds.all <- DESeqDataSetFromMatrix(counts.all[sel,], pheno.all, ~group + type)
dds.all <- DESeq(dds.all)
resultsNames(dds.all)
res.all <- results(dds.all, name = 'group_Paediatric_vs_Adult')
res.all <- res.all[order(abs(res.all$stat), decreasing = T),]
# Note the NAs here represent genes expressed in a single sample
# Remove these
res.all <- res.all[which(!is.na(res.all$padj)),]

counts.all.norm <- counts(dds.all, normalized=TRUE)
colnames(counts.all.norm) <- colnames(counts.all)

########################################################################
# Rank order plots
########################################################################
source('rankOrderPlot.R')
fig.rankorder <- plot_grid(
  rankOrderPlot(genes.basal.filt, counts.unsorted, "Unsorted: Basal genes"),
  rankOrderPlot(genes.secretory.filt, counts.unsorted, "Unsorted: Secretory genes"),
  rankOrderPlot(genes.ciliated.filt, counts.unsorted, "Unsorted: Ciliated genes"),
  
  rankOrderPlot(genes.basal.filt, counts.basal, "Sorted: Basal genes"),
  rankOrderPlot(genes.secretory.filt, counts.basal, "Sorted: Secretory genes"),
  rankOrderPlot(genes.ciliated.filt, counts.basal, "Sorted: Ciliated genes"),
  
  rankOrderPlot(genes.basal.filt, counts.cultured, "Cultured: Basal genes"),
  rankOrderPlot(genes.secretory.filt, counts.cultured, "Cultured: Secretory genes"),
  rankOrderPlot(genes.ciliated.filt, counts.cultured, "Cultured: Ciliated genes")
)
save_plot(paste0(opdir, "fig.X.rankorder_celltypes.pdf"), fig.rankorder, base_width = 1.5*fig.width, base_height = 1.5*fig.width)

# Repeat adults and kids:
# sel1 <- which(pheno.unsorted$group == 'Adult')
# sel2 <- which(pheno.basal$group == 'Adult')
# sel3 <- which(pheno.cultured$group == 'Adult')
# plot_grid(
#   rankOrderPlot(genes.basal.filt, counts.unsorted[,sel1], "Unsorted: Basal genes"),
#   rankOrderPlot(genes.secretory.filt, counts.unsorted[,sel1], "Unsorted: Secretory genes"),
#   rankOrderPlot(genes.ciliated.filt, counts.unsorted[,sel1], "Unsorted: Ciliated genes"),
#   
#   rankOrderPlot(genes.basal.filt, counts.basal[,sel2], "Sorted: Basal genes"),
#   rankOrderPlot(genes.secretory.filt, counts.basal[,sel2], "Sorted: Secretory genes"),
#   rankOrderPlot(genes.ciliated.filt, counts.basal[,sel2], "Sorted: Ciliated genes"),
#   
#   rankOrderPlot(genes.basal.filt, counts.cultured[,sel3], "Cultured: Basal genes"),
#   rankOrderPlot(genes.secretory.filt, counts.cultured[,sel3], "Cultured: Secretory genes"),
#   rankOrderPlot(genes.ciliated.filt, counts.cultured[,sel3], "Cultured: Ciliated genes")
# )
# sel1 <- which(pheno.unsorted$group == 'Paediatric')
# sel2 <- which(pheno.basal$group == 'Paediatric')
# sel3 <- which(pheno.cultured$group == 'Paediatric')
# plot_grid(
#   rankOrderPlot(genes.basal.filt, counts.unsorted[,sel1], "Unsorted: Basal genes"),
#   rankOrderPlot(genes.secretory.filt, counts.unsorted[,sel1], "Unsorted: Secretory genes"),
#   rankOrderPlot(genes.ciliated.filt, counts.unsorted[,sel1], "Unsorted: Ciliated genes"),
#   
#   rankOrderPlot(genes.basal.filt, counts.basal[,sel2], "Sorted: Basal genes"),
#   rankOrderPlot(genes.secretory.filt, counts.basal[,sel2], "Sorted: Secretory genes"),
#   rankOrderPlot(genes.ciliated.filt, counts.basal[,sel2], "Sorted: Ciliated genes"),
#   
#   rankOrderPlot(genes.basal.filt, counts.cultured[,sel3], "Cultured: Basal genes"),
#   rankOrderPlot(genes.secretory.filt, counts.cultured[,sel3], "Cultured: Secretory genes"),
#   rankOrderPlot(genes.ciliated.filt, counts.cultured[,sel3], "Cultured: Ciliated genes")
# )


# Viral and COVID:
fig.rankorder.covid <- plot_grid(
  rankOrderPlot(genes.viral, counts.unsorted, "Unsorted: Viral genes"),
  rankOrderPlot(genes.covid, counts.unsorted, "Unsorted: COVID genes"),
  
  rankOrderPlot(genes.viral, counts.basal, "Sorted: Viral genes"),
  rankOrderPlot(genes.covid, counts.basal, "Sorted: COVID genes"),
  
  rankOrderPlot(genes.viral, counts.cultured, "Cultured: Viral genes"),
  rankOrderPlot(genes.covid, counts.cultured, "Cultured: COVID genes"),
  
  ncol = 2
)
save_plot(paste0(opdir, "fig.X.rankorder_covid.pdf"), fig.rankorder.covid, base_width = 1.5*fig.width, base_height = 1.5*fig.width)

########################################################################
# Cell types
########################################################################
# Make a heatmap showing basal, secretory and ciliated gene expression in the unsorted dataset

annot.genes <- data.frame(
  row.names = c(genes.basal.filt, genes.secretory.filt, genes.ciliated.filt),
  celltype = c(rep('Basal', length(genes.basal.filt)), rep('Secretory', length(genes.secretory.filt)), rep('Ciliated', length(genes.ciliated.filt)))
)
genes <- intersect(rownames(annot.genes), rownames(counts.unsorted.norm))
pheatmap(counts.unsorted.norm[genes, ], annotation_row = annot.genes, annotation_col = annot.unsorted, scale='row', cluster_rows = F, annotation_colors = annot.colors)

genes <- intersect(rownames(annot.genes), rownames(counts.basal.norm))
pheatmap(counts.basal.norm[genes, ], annotation_row = annot.genes, annotation_col = annot.basal, scale='row', cluster_rows = F, annotation_colors = annot.colors)

genes <- intersect(rownames(annot.genes), rownames(counts.cultured.norm))
pheatmap(counts.cultured.norm[genes, ], annotation_row = annot.genes, annotation_col = annot.cultured, scale='row', cluster_rows = F, annotation_colors = annot.colors)

# Different way to display this:
df <- lapply(genes, function(gene) {
  x <- pheno.unsorted
  x$gene = gene
  x$val = as.numeric(counts.unsorted.norm[gene,])
  if(gene %in% genes.basal) {x$set = 'Basal'}
  if(gene %in% genes.secretory) {x$set = 'Goblet'}
  if(gene %in% genes.ciliated) {x$set = 'Ciliated'}
  x
})
df <- do.call('rbind', df)
df$gene <- factor(df$gene, levels = genes)
ggplot(df, aes(x = gene, y=log(val))) +
  geom_boxplot(aes(fill=set)) + geom_point() +
  facet_grid(group ~ .) +
  theme(axis.text.x = element_text(angle = 90))
# Collapse by sample using geomean:
source('geomean.R')
df2 <- aggregate(df[,'val'], by=list(df$sample, df$group, df$set), FUN=geomean)
colnames(df2)[which(colnames(df2) == 'Group.2')] <- 'Group'
# p.bas <- wilcox.test(log(x) ~ Group, data = df2[which(df2$Group.3 == 'Basal'),])$p.value
# p.sec <- wilcox.test(log(x) ~ Group, data = df2[which(df2$Group.3 == 'Secretory'),])$p.value
# p.cil <- wilcox.test(log(x) ~ Group, data = df2[which(df2$Group.3 == 'Ciliated'),])$p.value
# stat.test <- compare_means(x ~ Group, data = df2[which(df2$Group.3 == 'Basal'),], method='wilcox.test')
# 
# ggplot(df2, aes(x=Group.3, y=log(x), fill=Group)) +
#   geom_boxplot(outlier.shape = NA) +
#   geom_point(aes(group=Group), position = position_dodge(width=0.75)) +
#   theme(axis.title.x = element_blank()) +
#   ylab("Log normalised expression") +
#   stat_compare_means(aes(group = Group)) +
#   gg.fill
  # stat_pvalue_manual(stat.test, y.position = 10, xmin = 0.25, xmax = 0.75)
  # stat_pvalue_manual(data.frame(group1='Basal', group2='Secretory', p=p.bas, y.position=10, xmin=0.25, xmax=1.75))

fig.1c <- ggplot(df2, aes(x=Group, y=log(x), fill=Group)) +
  geom_boxplot() +
  geom_boxplot(outlier.shape = NA) +
  geom_point() +
  facet_wrap(~Group.3) +
  geom_signif(comparisons = list(c('Adult', 'Paediatric')), map_signif_level = F) +
  guides(fill=F) +
  ylab('Log normalised expression') +
  theme(axis.title.x = element_blank()) +
  gg.fill
  # stat_compare_means(aes(group=Group)) +
  # ylim(c(1.5, 11))

# save_plot(paste0(opdir, "Fig.1c.celltype_comp.pdf"), fig.1c, base_width = fig.width, base_height = fig.width)

# ggplot(pdata2, aes(x=name, y=val)) + 
#   geom_boxplot(aes(fill=Outcome), outlier.size = 0.5) +
#   geom_point(position = position_dodge(width=0.75), aes(group=Outcome), size=0.5) +


# Danaher breakdown:
dan <- do.danaher(counts.unsorted.norm)
dan <- dan[,which(apply(dan, 2, function(x) {!all(is.na(x))}))]
x <- pheno.unsorted
x$til.score <- dan[x$sample,]$til.score
ggplot(x, aes(x=group, y=til.score)) + 
  geom_boxplot(aes(fill=group)) + geom_point() +
  geom_signif(comparisons = list(c('Adult', 'Paediatric'))) +
  gg.fill

pheatmap(t(dan[,-14]), annotation_col = annot.unsorted, annotation_colors = annot.colors)

# Compare Danaher results adults vs paeds:
dan.comp <- sapply(colnames(dan), function(x) {
  p <- pheno.unsorted
  p$val <- as.numeric(dan[,x])
  wilcox.test(val ~ group, data = p)$p.value
})
dan.comp <- data.frame(
  celltype = names(dan.comp),
  p = as.numeric(dan.comp),
  p.adj = p.adjust(as.numeric(dan.comp))
)



# Repeat cell type analysis with all to show selection worked!
genes <- intersect(rownames(annot.genes), rownames(counts.all.norm))
annot.all <- data.frame(row.names = pheno.all$sample, group = pheno.all$group, type = pheno.all$type)
fig <- pheatmap(counts.all.norm[genes, ], annotation_row = annot.genes, annotation_col = annot.all, scale='row', cluster_rows = F, annotation_colors = annot.colors, show_colnames = T)
fig

save_plot(paste0(opdir, "fig.S6A.celltype.heatmap.pdf"), fig, base_width = fig.width, base_height = fig.maxheight)

########################################################################
# Viral infections
########################################################################
# Plot viral genes as above:
genes <- intersect(rownames(counts.unsorted.norm), genes.viral)
df <- lapply(genes, function(gene) {
  x <- pheno.unsorted
  x$gene = gene
  x$val = as.numeric(counts.unsorted.norm[gene,])
  x
})
df <- do.call('rbind', df)
df$gene <- factor(df$gene, levels = genes)
ggplot(df, aes(x = gene, y=log(val))) +
  geom_boxplot() + geom_point() +
  facet_grid(group ~ .) +
  theme(axis.text.x = element_text(angle = 90))
# Collapse by sample using geomean:
source('geomean.R')
df2 <- aggregate(df[,'val'], by=list(df$sample, df$group), FUN=geomean)
ggplot(df2, aes(x=Group.2, y=log(x))) +
  geom_boxplot(aes(fill=Group.2)) +
  geom_point() +
  geom_signif(comparisons = list(c('Adult', 'Paediatric'))) +
  ylab("log(mean of viral response genes)") +
  theme(axis.title.x = element_blank())



# Consider differential analysis and heatmaps specific to viral/COVID genes
genes <- unique(genes.viral)
missing <- genes[which(!(genes %in% rownames(counts.basal)))]
if(length(missing) > 0) {
  message(paste("Check gene names for: ", paste(missing, collapse=", ")))
}
# annot.row <- data.frame(row.names = genes, set = ifelse(genes %in% genes.viral, "Viral response", "COVID related"))

# counts.basal.norm <- counts(dds.basal, normalized=TRUE)

annot.unsorted <- data.frame(row.names = pheno.unsorted$sample, group = pheno.unsorted$group)
hm.viral.unsorted <- pheatmap(counts.unsorted.norm[intersect(genes, rownames(counts.unsorted.norm)),], annotation_col = annot.unsorted, scale='row', annotation_colors = annot.colors, fontsize = 8, annotation_legend = F)

annot.basal <- data.frame(row.names = pheno.basal$sample, group = pheno.basal$group)
hm.viral.basal <- pheatmap(counts.basal.norm[intersect(genes, rownames(counts.basal.norm)),], annotation_col = annot.basal,  scale='row', annotation_colors = annot.colors, fontsize = 8, annotation_legend = F)

annot.cultured <- data.frame(row.names = pheno.cultured$sample, group = pheno.cultured$group)
hm.viral.cultured <- pheatmap(counts.cultured.norm[intersect(genes, rownames(counts.cultured.norm)),], annotation_col = annot.cultured,  scale='row', annotation_colors = annot.colors, fontsize = 8, annotation_legend = F)

# Perform DE analysis on these gene sets alone
sel <- which(rownames(counts.unsorted) %in% genes & rowSums(counts.unsorted) != 0)
dds.unsorted.viral <- DESeqDataSetFromMatrix(counts.unsorted[sel,], pheno.unsorted, ~group)
dds.unsorted.viral <- DESeq(dds.unsorted.viral)
resultsNames(dds.unsorted.viral)
res.unsorted.viral <- results(dds.unsorted.viral, name = 'group_Paediatric_vs_Adult')
res.unsorted.viral <- res.unsorted.viral[order(abs(res.unsorted.viral$stat), decreasing = T),]
res.unsorted.viral <- res.unsorted.viral[which(!is.na(res.unsorted.viral$padj)),]

# Repeat for basal
sel <- which(rownames(counts.basal) %in% genes & rowSums(counts.basal) != 0)
dds.basal.viral <- DESeqDataSetFromMatrix(counts.basal[sel,], pheno.basal, ~group)
dds.basal.viral <- DESeq(dds.basal.viral)
resultsNames(dds.basal.viral)
res.basal.viral <- results(dds.basal.viral, name = 'group_Paediatric_vs_Adult')
res.basal.viral <- res.basal.viral[order(abs(res.basal.viral$stat), decreasing = T),]
res.basal.viral <- res.basal.viral[which(!is.na(res.basal.viral$padj)),]
# Result - none of these genes come up as significant.

# Repeat for cultured
sel <- which(rownames(counts.cultured) %in% genes & rowSums(counts.cultured) != 0)
dds.cultured.viral <- DESeqDataSetFromMatrix(counts.cultured[sel,], pheno.cultured, ~group)
dds.cultured.viral <- DESeq(dds.cultured.viral)
resultsNames(dds.cultured.viral)
res.cultured.viral <- results(dds.cultured.viral, name = 'group_Paediatric_vs_Adult')
res.cultured.viral <- res.cultured.viral[order(abs(res.cultured.viral$stat), decreasing = T),]
res.cultured.viral <- res.cultured.viral[which(!is.na(res.cultured.viral$padj)),]



# Add a boxplot of the COVID-specific genes:
gene.plot <- function(counts.norm, pheno, genes = genes.covid) {
  x <- pheno
  genes <- intersect(genes, rownames(counts.norm))
  for(gene in genes) {
    x[,gene] <- as.numeric(counts.norm[gene,])
  }
  df <- lapply(genes, function(gene) {
    y <- x
    y$gene <- gene
    y$val <- y[,gene]
    y
  })
  df <- do.call('rbind', df)
  df$gene <- factor(as.character(df$gene), levels = genes)
  
  # Calculate corrected p-values:
  ps <- sapply(genes, function(x) {
    wilcox.test(val ~ group, data = df[which(df$gene == x),])$p.value
  })
  names(ps) <- genes
  ps.adj <- p.adjust(ps, method = 'fdr')
  
  
  scaleFUN <- function(x) sprintf("%.1f", x)
  
  plots <- lapply(genes, function(gene) {
    sel <- which(df$gene == gene)
    m <- log(max(df$val[sel])) * 1.1
    p <- data.frame(group1='Adult', group2='Paediatric', p = signif(ps.adj[[gene]],2), y.position = m)
    ggplot(df[sel,], aes(x=group, y=log(val))) + 
      geom_boxplot(aes(fill=group)) + geom_point() + 
      stat_pvalue_manual(p, vjust = 1.5) +
      # geom_signif(comparisons = list(c('Adult', 'Paediatric')), vjust = 1.3, map_signif_level = F) +
      # facet_grid(gene ~ ., scales = 'free_y') +
      ylab('Log normalised expression') +
      theme(
        axis.title = element_blank(), axis.text.x = element_blank(), legend.direction = 'horizontal', legend.justification = 'center',
        axis.text.y = element_text()
      ) +
      gg.fill +
      guides(fill=F) +
      annotate("text", x=2.5, y=m * 0.95, label = gene, hjust="right", vjust = 0, size = 3) +
      scale_y_continuous(labels = scaleFUN)
      # ggtitle(gene)
  })
  plot_grid(plotlist = plots, ncol = 1)
}
# fig.covid2 <- gene.plot(counts.unsorted.norm, pheno.unsorted)
# fig.covid2


# x <- pheno.unsorted
# genes <- intersect(genes.covid, rownames(counts.unsorted.norm))
# for(gene in genes) {
#   x[,gene] <- as.numeric(counts.unsorted.norm[gene,])
# }
# df <- lapply(genes, function(gene) {
#   y <- x
#   y$gene <- gene
#   y$val <- y[,gene]
#   y
# })
# df <- do.call('rbind', df)
# df$gene <- factor(as.character(df$gene), levels = genes)
# fig.covid <- ggplot(df, aes(x=group, y=log(val))) + 
#   geom_boxplot(aes(fill=group)) + geom_point() + 
#   geom_signif(comparisons = list(c('Adult', 'Paediatric')), vjust = 1.3, map_signif_level = F) +
#   facet_grid(gene ~ ., scales = 'free_y') +
#   ylab('Log normalised expression') +
#   theme(axis.title.x = element_blank(), legend.direction = 'horizontal', legend.justification = 'center') +
#   gg.fill
# 
# leg <- get_legend(fig.covid)
# fig.covid <- fig.covid + guides(fill=F)
fig.covid <- gene.plot(counts.unsorted.norm, pheno.unsorted)

fig <- plot_grid(
  # plot_grid(
    as.ggplot(hm.viral.unsorted), 
    fig.covid, 
    ncol=2, rel_widths = c(2,1), labels = c('A', 'B')
  # ), 
  # leg,
  # ncol=1, 
  # rel_heights = c(20,1)
)
fig
save_plot(paste0(opdir, "fig.6.COVID.genes.unsorted.pdf"), fig, base_width = fig.width * 1.5, base_height = fig.maxheight * 1.5)

# Repeat for basal
# x <- pheno.basal
# genes <- intersect(genes.covid, rownames(counts.basal.norm))
# for(gene in genes) {
#   x[,gene] <- as.numeric(counts.basal.norm[gene,])
# }
# df <- lapply(genes, function(gene) {
#   y <- x
#   y$gene <- gene
#   y$val <- y[,gene]
#   y
# })
# df <- do.call('rbind', df)
# fig.covid.basal <- ggplot(df, aes(x=group, y=log(val))) + 
#   geom_boxplot(aes(fill=group)) + geom_point() + 
#   stat_compare_means(method = 'wilcox.test', method.args = list(alternative = 'two.sided')) +
#   # geom_signif(comparisons = list(c('Adult', 'Paediatric')), vjust = 1.3, map_signif_level = F) +
#   facet_grid(gene ~ ., scales = 'free_y') +
#   ylab('Log normalised expression') +
#   theme(axis.title.x = element_blank(), legend.direction = 'horizontal', legend.justification = 'center') +
#   gg.fill
# # Check the wilcox test p-values:
# ps <- sapply(unique(df$gene), function(x) {
#   wilcox.test(log(val) ~ group, data = df[which(df$gene == x),])$p.value
# })
# names(ps) <- unique(df$gene)
# 
# leg <- get_legend(fig.covid.basal)
# fig.covid.basal <- fig.covid.basal + guides(fill=F)
# fig.covid.basal

fig.covid.basal <- gene.plot(counts.basal.norm, pheno.basal)

fig <- plot_grid(
  # plot_grid(
    as.ggplot(hm.viral.basal), 
    fig.covid.basal, 
    ncol=2, rel_widths = c(2,1), labels = c('A', 'B')
  # ), 
  # leg,
  # ncol=1, 
  # rel_heights = c(20,1)
)
fig
save_plot(paste0(opdir, "fig.S6.COVID.genes.basal.pdf"), fig, base_width = fig.width * 1.5, base_height = fig.maxheight * 1.5)


# Repeat for cultured
# x <- pheno.cultured
# genes <- intersect(genes.covid, rownames(counts.cultured.norm))
# for(gene in genes) {
#   x[,gene] <- as.numeric(counts.cultured.norm[gene,])
# }
# df <- lapply(genes, function(gene) {
#   y <- x
#   y$gene <- gene
#   y$val <- y[,gene]
#   y
# })
# df <- do.call('rbind', df)
# fig.covid.cultured <- ggplot(df, aes(x=group, y=log(val))) + 
#   geom_boxplot(aes(fill=group)) + geom_point() + 
#   geom_signif(comparisons = list(c('Adult', 'Paediatric')), vjust = 1.3, map_signif_level = F) +
#   facet_grid(gene ~ ., scales = 'free_y') +
#   ylab('Log normalised expression') +
#   theme(axis.title.x = element_blank(), legend.direction = 'horizontal', legend.justification = 'center') +
#   gg.fill
# 
# leg <- get_legend(fig.covid.cultured)
# fig.covid.cultured <- fig.covid.cultured + guides(fill=F)

fig.covid.cultured <- gene.plot(counts.cultured.norm, pheno.cultured)

fig <- plot_grid(
  # plot_grid(
    as.ggplot(hm.viral.cultured), 
    fig.covid.cultured, 
    ncol=2, rel_widths = c(2,1), labels = c('A', 'B')
  # ), 
  # leg, 
  # ncol=1, 
  # rel_heights = c(20,1)
)

save_plot(paste0(opdir, "fig.S6.COVID.genes.cultured.pdf"), fig, base_width = fig.width * 1.5, base_height = fig.maxheight * 1.5)



# Similar for viral response genes:

######################################################
# Figures (built from above)
######################################################
x <- counts.unsorted.norm[rownames(res.unsorted)[which(res.unsorted$padj < 0.01 & abs(res.unsorted$log2FoldChange) > 1.2)],]
fig.2a <- pheatmap(x, scale='row', annotation_col = annot.unsorted, annotation_colors = annot.colors, annotation_legend = F, fontsize = 6)

pathway.plot <- function(gsea) {
  x <- gsea[which(gsea$padj < 0.05),]
  # x$pathway <- as.character(str_trunc(x$pathway, 30))
  x <- x[order(x$NES),]
  x$pathway <- gsub("HALLMARK_", "", x$pathway)
  x$pathway <- gsub("_", " ", x$pathway)
  x$pathway <- factor(x$pathway, levels = x$pathway)
  x$group <- ifelse(x$NES < 0, 'Adult', 'Paediatric')
  x$group <- factor(x$group, levels = c('Adult', 'Paediatric'))
  ggplot(x, aes(x=pathway, y=-NES)) +
    geom_col(aes(fill=group)) +
    coord_flip() +
    scale_x_discrete(position = "top") +
    theme(axis.text.y = element_text(size=6), axis.title.y = element_blank(), 
          axis.title.x = element_text(size = 8), legend.direction = 'horizontal', legend.justification = "center") + 
    ylab('GSEA Normalised Enrichment Score') +
    gg.fill
}

fig.2b <- pathway.plot(gsea.unsorted)

leg <- get_legend(fig.2b)
fig.2b <- fig.2b + guides(fill = F)

fig.2 <- plot_grid(as.ggplot(fig.2a), fig.2b, leg, ncol=1, rel_heights = c(20,10,1.5), labels = c('A','B'), axis = 't') #, align = 'h')
fig.2

save_plot(paste0(opdir, "fig.2.pdf"), fig.2, base_width = fig.width, base_height = fig.maxheight)


# Supplementary heatmap of basal genes
x <- counts.basal.norm[rownames(res.basal)[which(res.basal$padj < 0.01 & abs(res.basal$log2FoldChange) > 1.2)],]
fig.basal.hm <- pheatmap(x, scale='row', annotation_col = annot.basal, annotation_colors = annot.colors, annotation_legend = F, fontsize = 6, show_rownames = F, show_colnames = F)

# save_plot(paste0(opdir, "fig.5.basal.hm.pdf"), fig.basal.hm, base_width = fig.width, base_height = fig.maxheight)

# Supplementary heatmap of cultured genes
x <- counts.cultured.norm[rownames(res.cultured)[which(res.cultured$padj < 0.01 & abs(res.cultured$log2FoldChange) > 1.2)],]
fig.cultured.hm <- pheatmap(x, scale='row', annotation_col = annot.cultured, annotation_colors = annot.colors, annotation_legend = F, fontsize = 6, show_rownames = F, show_colnames = F)

# save_plot(paste0(opdir, "fig.X.cultured.hm.pdf"), fig.cultured.hm, base_width = fig.width, base_height = fig.maxheight)

# Figure 5: Basal + cultured + cultured pathways

fig.5b <- pathway.plot(gsea.basal) + guides(fill = F)

fig.5c <- pathway.plot(gsea.cultured)
leg <- get_legend(fig.5c)
fig.5c <- fig.5c + guides(fill = F)


fig.5 <- plot_grid(
  plot_grid(
    as.ggplot(fig.basal.hm), 
    plot_grid(
      NULL, fig.5b, NULL,
      ncol=1,
      rel_heights = c(1,2,1)
    ),
    as.ggplot(fig.cultured.hm),
    fig.5c,
    ncol=2
    # labels = c("A","B","C","D")
  ),
  leg,
  ncol = 1,
  rel_heights = c(20,1)
)
fig.5

save_plot(paste0(opdir, "fig.5.pdf"), fig.5, base_width = fig.width, base_height = fig.maxheight)

######################################################
# Supplementary table 1: Gene lists
######################################################
WriteXLS(list(
  "Epithelial Viral Infection" = genes.viral.raw,
  "SARS-CoV-2 Related" = genes.covid.raw,
  "Basal markers" = data.frame(
    Gene = genes.basal,
    Passed.Filter = ifelse(genes.basal %in% genes.basal.filt, "YES", "NO")
  ),
  "Secretory markers" = data.frame(
    Gene = genes.secretory,
    Passed.Filter = ifelse(genes.secretory %in% genes.secretory.filt, "YES", "NO")
  ),
  "Ciliated markers" = data.frame(
    Gene = genes.ciliated,
    Passed.Filter = ifelse(genes.ciliated %in% genes.ciliated.filt, "YES", "NO")
  )
), ExcelFileName = paste0(opdir, "SupTable1.xls"), AdjWidth = T)


######################################################
# Supplementary table 2: DE genes
######################################################
DEgenes.unsorted <- data.frame(res.unsorted[which(res.unsorted$padj < 0.01 & abs(res.unsorted$log2FoldChange) > 1.2),])
DEgenes.basal <- data.frame(res.basal[which(res.basal$padj < 0.01 & abs(res.basal$log2FoldChange) > 1.2),])
DEgenes.cultured <- data.frame(res.cultured[which(res.cultured$padj < 0.01 & abs(res.cultured$log2FoldChange) > 1.2),])

DEpaths.unsorted <- data.frame(gsea.unsorted)
DEpaths.basal <- data.frame(gsea.basal)
DEpaths.cultured <- data.frame(gsea.cultured)

# TODO: flip ES/NES

x <- c('DEgenes.unsorted', 'DEpaths.unsorted', 'DEgenes.basal', 'DEpaths.basal', 'DEgenes.cultured', 'DEpaths.cultured')
WriteXLS(x, ExcelFileName = paste0(opdir, 'SupTable2.xls'), row.names = T, AdjWidth = T)


# Replicate experimental plots:
source('replicate_plots.R')


######################################################
# Data for GEO
######################################################
# Submit as 3 separate datasets, to be wrapped in a superset
# Per-sample metadata is dumped to Excel here to be copied into a manually created template
geo.shared <- data.frame(
  V1= c(
    "title", 
    "summary",
    "overall design",
    "contributor",
    "contributor",
    "supplementary file", 
    "SRA_center_name_code"
  ),
  V2 = c(
    "Cell-intrinsic differences between human tracheal epithelial cells from children and adults",
    "The airway epithelium is a key protective barrier whose integrity is preserved by the self-renewal and differentiation of basal progenitor cells. Epithelial cells are central to the pathogenesis of multiple chronic lung diseases for which age is a principle risk factor. Children are also less susceptible to SARS-CoV-2 infection, suffer less severe symptoms than adults and have a lower rate of mortality. Few studies have addressed differences between airway epithelial cells in children and adults. Here, we perform bulk RNA sequencing studies in laser-captured whole epithelium, FACS-sorted basal cells and cultured basal cells, as well as in vitro cell proliferation experiments, to address the intrinsic molecular differences between paediatric and adult airway basal cells. We find that, while the cellular composition of the paediatric and adult tracheal epithelium is broadly similar, in cell culture, paediatric airway epithelial cells displayed higher colony forming ability, better in vitro growth and outcompeted adult cells in competitive proliferation assays. Although recurring differences between airway epithelial gene expression were seen between samples from children and adults, RNA sequencing showed broad conservation of transcriptional programmes. Genes associated with SARS-CoV-2 infection were not differentially expressed between children and adults, although individuals showed some variability in their expression of viral infection-associated genes. Our results chart important cell intrinsic differences in transcriptional profile and regenerative capacity between tracheal epithelial cells of children and adults.",
    "Comparison of airway RNAseq data between adults and children. Three related datasets are presented: fresh airway epithelial cells, fresh sorted basal cells and cultured basal cells.",
    "Adam Pennycuick",
    "Elizabeth F Maughan",
    "",
    ""
  )
)

geo.unsorted <- data.frame(
  "Sample name" = pheno.unsorted$sample,
  "title" = pheno.unsorted$sample,
  "source name" = "Airway epithelial cells",
  "organism" = "Homo sapiens",
  "characteristics: group" = pheno.unsorted$group,
  "characteristics: tag" = "",
  "characteristics: tag" = "",
  "molecule" = "total RNA",
  "description" = "",
  "processed data file" = "",
  "raw file" = "",
  "raw file" = "",
  "raw file" = "",
  "raw file" = "",
  "raw file" = ""
)

pheno.basal$r1_file <- sapply(pheno.basal$sample, function(x) {
  list.files("/Volumes/AdamData/RNAseq/lizzie_basal2/fastqs/", pattern = gsub("B_X", "", x))[1]
})
pheno.basal$r2_file <- sapply(pheno.basal$sample, function(x) {
  list.files("/Volumes/AdamData/RNAseq/lizzie_basal2/fastqs/", pattern = gsub("B_X", "", x))[2]
})

geo.basal <- data.frame(
  "Sample name" = pheno.basal$sample,
  "title" = pheno.basal$sample,
  "source name" = "Airway epithelial cells",
  "organism" = "Homo sapiens",
  "characteristics: group" = pheno.basal$group,
  "characteristics: tag" = "",
  "characteristics: tag" = "",
  "molecule" = "total RNA",
  "description" = "",
  "processed data file" = "",
  "raw file" = pheno.basal$r1_file,
  "raw file" = pheno.basal$r2_file,
  "raw file" = "",
  "raw file" = "",
  "raw file" = ""
)


pheno.cultured$r1_file <- sapply(pheno.cultured$sample, function(x) {
  list.files("/Volumes/AdamData/RNAseq/lizzie_cultured2/fastqs/", pattern = gsub("C_", "", x))[1]
})
pheno.cultured$r2_file <- sapply(pheno.cultured$sample, function(x) {
  list.files("/Volumes/AdamData/RNAseq/lizzie_cultured2/fastqs/", pattern = gsub("C_", "", x))[2]
})

geo.cultured <- data.frame(
  "Sample name" = pheno.cultured$sample,
  "title" = pheno.cultured$sample,
  "source name" = "Airway epithelial cells",
  "organism" = "Homo sapiens",
  "characteristics: group" = pheno.cultured$group,
  "characteristics: tag" = "",
  "characteristics: tag" = "",
  "molecule" = "total RNA",
  "description" = "",
  "processed data file" = "",
  "raw file" = pheno.cultured$r1_file,
  "raw file" = pheno.cultured$r2_file,
  "raw file" = "",
  "raw file" = "",
  "raw file" = ""
)

geo.dir <- paste0(opdir, "geo/")
dir.create(geo.dir, showWarnings = F)
WriteXLS(list(
    "shared" = geo.shared, 
    "unsorted" = geo.unsorted, 
    "basal" = geo.basal, 
    "cultured" = geo.cultured
  ), ExcelFileName = paste0(geo.dir, "geo_templates.xlsx"))

# Output count matrices for GEO as tab delimited text
write.table(counts.unsorted, file = paste0(geo.dir, "counts.unsorted.txt"), sep="\t", quote = F)
write.table(counts.basal, file = paste0(geo.dir, "counts.basal.txt"), sep="\t", quote = F)
write.table(counts.cultured, file = paste0(geo.dir, "counts.cultured.txt"), sep="\t", quote = F)




################################################################################
# Reviewer responses
################################################################################
# Compare non-smokers vs ex smokers in the unsorted group
# p <- pheno.unsorted[which(pheno.unsorted$group == 'Adult'),]
# d <- counts.unsorted[,which(pheno.unsorted$group == 'Adult')]
# p$smoking <- c("ex", "ex", "ex", "ex", "non", "non")
# dds <- DESeqDataSetFromMatrix(d, p, design = ~smoking)
# dds <- DESeq(dds)
# res <- results(dds, name = "smoking_non_vs_ex")
# res <- res[order(abs(res$stat), decreasing = T),]
# head(res)
# length(which(res$padj < 0.05))
