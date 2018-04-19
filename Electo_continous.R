###
# Author: Verena Friedl
# Created: April 27, 2017
# Last changed: April 27, 2017


#setwd('~/Google Drive/stuartlab/Treehouse/SCRIPTS')
source('./Electo_methods.R')

similarity_before_filename <- "../DATA/simMatrix.pancan.atlas.naAs0.nonaTis.noZeroVar.tsv"
sample_rankings <- readRDS(paste0(similarity_before_filename,"_sampleRankings.RDS"))

#expression <- read.delim("../DATA/TCGA_PANCANATLAS/mRNA/out.tsv",
#                       row.names = 1,
#                       header=T,
#                       stringsAsFactors=F,
#                       check.names=F)
#expression <- t(expression)
#saveRDS(expression,file = "../DATA/TCGA_PANCANATLAS/mRNA/out.RDS")

expression <- readRDS("../DATA/TCGA_PANCANATLAS/mRNA/out.RDS")
expression <- expression[rownames(sample_rankings),]

feature_correlation <- c()
#for(i in c(1:ncol(expression))){
genes <- unique(read.delim("../DATA/differentially_expressed_genes_in_cancer.txt",
                    header=T,
                    stringsAsFactors=F,
                    check.names=F)[,1])
genes_first <- sapply(genes,function(x) strsplit(x,"/",fixed=T)[[1]][1])
genes_second <- sapply(genes,function(x) strsplit(x,"/",fixed=T)[[1]][2])
genes_second <- genes_second[which(!(is.na(genes_second)))]
genes_all <- unique(c(genes_first,genes_second))
genes_intersect <- intersect(sapply(colnames(expression),function(x) strsplit(x,"|",fixed=T)[[1]][1]),genes_all)
for(gene in genes_intersect){
  ###
  gene <- "ABCD3"
  ###
  i <- which(sapply(colnames(expression),function(x) strsplit(x,"|",fixed=T)[[1]][1]) == gene)
  feature_name <- colnames(expression)[i]
  feature_values <- expression[,i]
  names(feature_values) <- rownames(expression)
  
  samples <- rownames(sample_rankings)
  correlations <- c()
  sample_feature_values <- c()
  for(sample in samples){
    ranks <- c(1:ncol(sample_rankings))
    names(ranks) <- sample_rankings[sample,]
    #ranks <- ranks[names(feature_values)]
    feature_values_ranked <- feature_values[sample_rankings[sample,]] 
    plot(ranks, feature_values_ranked, pch=16, col=rgb(0.5,0.5,0.5,0.5))
    correlation <- cor.test(ranks, feature_values_ranked, method="spearman")$estimate
    correlations <- c(correlations,correlation)
    sample_feature_values <- c(sample_feature_values,feature_values[sample])
  }
  png(paste0("../plots/expressionFeatures/",feature_name,".png"))
  plot(correlations,sample_feature_values, pch=16, col=rgb(0.5,0.5,0.5,0.5)
       ,main=feature_name
       ,xlab = "correlation(rank,attr value)"
       ,ylab = "query attribute value"
       ,xlim = c(-1,1))
  c <- cor.test(correlations,sample_feature_values, method="spearman")$estimate
  mtext(paste0("spearman correlation: ",round(c,digits=5)), side = 3)
  dev.off()
  feature_correlation <- c(feature_correlation,c)
}
saveRDS(feature_correlation, file = "../DATA/feature_correlation_values_continuousElecto_selectionDiffExpressedGenes.RDS")

