# Rank order plots
# Annotate where the expression of a given gene in a dataset sits within the distribution of all genes
# Expects RAW COUNTS as input
load('data/gtex.lung.cache.RData')
unexpressed.genes <- rownames(gtex.lung.counts)[which(rowSums(gtex.lung.counts) == 0)]
# rankOrderPlot(sample(unexpressed.genes, 100), counts.unsorted)

rankOrderPlot <- function(genes, counts, title=NA) {
  # For testing purposes:
  # genes <- c("GAPDH", genes)
  
  
  # Take mean CPM per sample as our aggregate expression measure
  c <- edgeR::cpm(counts)
  
  x <- data.frame(
    gene = rownames(c),
    val = rowSums(c) / dim(c)[2],
    label = NA,
    stringsAsFactors = F
  )
  sel <- which(x$gene %in% genes)
  x$label[sel] <- x$gene[sel]
  x <- x[order(x$val),]
  x$gene <- factor(x$gene, levels = x$gene)
  
  # Define a 'noise' cutoff
  # Take genes unexpressed in lung from GTEX, take their mean + 2sd as a cutoff
  # Knock off the top and bottom quartiles to avoid skew by outliers
  g2 <- intersect(rownames(counts), unexpressed.genes)
  sel <- which(x$gene %in% g2)
  vals <- x$val[sel]
  vals <- vals[which(vals >= quantile(vals)[2] & vals <= quantile(vals)[4])]
  cutoff <- mean(vals) + 2*sd(vals)
  
  # Print how many genes are above the cutoff:
  print(length(which(rowSums(counts) > cutoff)) / dim(counts)[1])
  
  fig <- ggplot(x, aes(x=gene, y=log(val+1e-5))) +
    geom_hline(yintercept = log(cutoff+1e-5), color='grey') +
    geom_point(aes(color = gene %in% genes), size = 0.1) +
    geom_label_repel(aes(label = label)) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank()) +
    guides(color=F) +
    ylab('Summed log normalised counts')
  
  
  if(!is.na(title)) {
    fig <- fig + ggtitle(title)
  }
  fig
}

