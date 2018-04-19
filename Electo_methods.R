###
# Author: Verena Friedl
# Created: April 19, 2017
# Last changed: May 11, 2017


# INPUTS:
  # similarity matrix: sample by sample
  # features: sample by feature

### METHODS
readSimilarityMatrix <- function(filename){
  similarity <- read.delim(filename,
                           header=T,
                           row.names=1,
                           stringsAsFactors=F,
                           check.names=F)
  saveRDS(similarity,file = paste0(filename,".RDS"))
  return(similarity)
}

# calculate ranking of samples by similarity to each sample
# input: similarity matrix
# return: sample by rank matrix, containing sample names ordered by similarity to row name
getSampleRankings <- function(similarity_matrix){
  samples <- rownames(similarity_matrix)
  sample_ranking_mat <- c()
  for(sample in samples){
    ranked_samples <- names(sort(similarity_matrix[sample,,drop=F],decreasing = T)[-c(1)])
    sample_ranking_mat <- rbind(sample_ranking_mat,ranked_samples)
  }
  rownames(sample_ranking_mat) <- samples
  colnames(sample_ranking_mat) <- c(1:(length(samples)-1))
  return(sample_ranking_mat)
}

# calculate KS test of ranks of positive samples to uniform distribution
# input: ranked samples, names of positive samples
# return: KS distance
KStoUniform <- function(ranked_samples, pos_samples){
  pos_ranks <- which(ranked_samples %in% pos_samples)
  ks <- suppressWarnings(ks.test(pos_ranks,"punif",1,length(ranked_samples),alternative="greater")) #### ToDo: change back to two.sided OR leave, but decide!!!!!!!!!!!!
  return(ks$statistic)
}

# get KS distance for all samples in a set
# input: sample set, the sample ranking matrix, names of positive samples
# return: vector of KS distances
getKSDistribution <- function(samples, sample_ranking_mat, pos_samples){
  distances <- c()
  for(sample in samples){
    ranked_samples <- sample_ranking_mat[sample,]
    ranked_samples <- ranked_samples[which(ranked_samples %in% rownames(sample_ranking_mat))]
    distances <- c(distances,KStoUniform(ranked_samples, pos_samples))
  }
  names(distances) <- samples
  return(distances)
}

testSeparationRaw <- function(feature.object){
  ks_greater <- suppressWarnings(ks.test(feature.object$pos_distribution_raw,feature.object$neg_distribution_raw,alternative="greater"))
  ks_less<- suppressWarnings(ks.test(feature.object$pos_distribution_raw,feature.object$neg_distribution_raw,alternative="less"))

  if(ks_less$statistic > ks_greater$statistic){
    direction <- 1
    ks_distance <- ks_less$statistic
    ks_pvalue <- ks_less$p.value
  } else{
    direction <- -1
    ks_distance <- ks_greater$statistic
    ks_pvalue <- ks_greater$p.value
  }
  return(list(ks_distance,ks_pvalue,direction))
}

testSeparation <- function(feature.object){
  pos_infered_data <- c()
  neg_infered_data <- c()
  number_pos_samples <- length(feature.object$pos_samples) + feature.object$pseudocounts
  number_neg_samples <- length(feature.object$neg_samples)
  
  for(i in c(1:length(feature.object$pos_distribution))){
    p_count <- round(feature.object$pos_distribution[i] * number_pos_samples)
    pos_infered_data <- c(pos_infered_data,
                          runif(p_count, min=feature.object$breaks[i], max=feature.object$breaks[i+1]))
    n_count <- round(feature.object$neg_distribution[i] * number_neg_samples)
    neg_infered_data <- c(neg_infered_data,
                          runif(n_count, min=feature.object$breaks[i], max=feature.object$breaks[i+1]))
  }
  
  ks_greater <- suppressWarnings(ks.test(pos_infered_data,neg_infered_data,alternative="greater"))
  ks_less<- suppressWarnings(ks.test(pos_infered_data,neg_infered_data,alternative="less"))

  if(ks_less$statistic > ks_greater$statistic){
    direction <- 1
    ks_distance <- ks_less$statistic
    ks_pvalue <- ks_less$p.value
  } else{
    direction <- -1
    ks_distance <- ks_greater$statistic
    ks_pvalue <- ks_greater$p.value
  }
  return(list(ks_distance,ks_pvalue,direction))
}

# smooth a histogram
# input: the histogram counts, smoothing parameter alpha (defaults to 1)
# return: smoothed histogram counts (normalized to sum up to 1)
smoothHistogram <- function(counts, breaks = c(0:20)/20, alpha = 1){

  counts_smoothed <- c()
  for(i in c(2:length(breaks))){
    #count = histogram$counts[i-1]
    bin <- breaks[i]
    
    bin_sum = 0
    for(j in c(2:length(breaks))){
      count <- counts[j-1]
      bin2 <- breaks[j]
      distance = abs(bin - bin2)*20

      bin_sum = bin_sum + (2^(-alpha * distance) * count)
    }
    counts_smoothed <- c(counts_smoothed,bin_sum)
  }
  counts_smoothed <- counts_smoothed/sum(counts_smoothed)
  return(counts_smoothed)
}

# add pseudocounts to histogram, based on the frequencies in another histogram
# input: the histogram, total number of pseudocounts to add (defaults to 10)
# return: histogram counts including pseudocounts
addPseudocounts <- function(histogram, frequencies, pseudocounts = 10){
  counts_pseudo <- c()
  for(i in c(1:length(frequencies))){
    base_frequency <- frequencies[i]
    new_count <- histogram$counts[i] + (base_frequency * pseudocounts)
    counts_pseudo <- c(counts_pseudo,new_count)
  }
  return(counts_pseudo)
}

calculateLogOddsRatio <- function(feature.object){
  log_ratios <- c()
  for(i in c(1:length(feature.object$neg_distribution))){
    p_neg <- feature.object$neg_distribution[i]
    p_pos <- feature.object$pos_distribution[i]
    log_ratios <- c(log_ratios,(log10(p_pos/p_neg)))
  }
  return(log_ratios)
}
calculateProbabilities <- function(feature.object,prior){
  probability_vector <- c()
  for(i in c(1:length(feature.object$neg_distribution))){
    p_neg <- feature.object$neg_distribution[i]
    p_pos <- feature.object$pos_distribution[i]
    prob <- (p_pos*prior) / (p_pos*prior + (p_neg*(1-prior)))
    probability_vector <- c(probability_vector,prob)
  }
  return(probability_vector)
}

getTissuePriors <- function(feature.object, disease_annotation){
  disease_types <- unique(disease_annotation)
  priors <- c()
  for(type in disease_types){
    type_samples <- names(disease_annotation[disease_annotation == type])
    pos_samples_type <- feature.object$pos_samples[which(feature.object$pos_samples %in% type_samples)]
    neg_samples_type <- feature.object$neg_samples[which(feature.object$neg_samples %in% type_samples)]
    priors <- c(priors,length(pos_samples_type)/(length(neg_samples_type)+length(pos_samples_type)))
  }
  names(priors) <- disease_types
  return(priors)
}

getTissueProbabilities <- function(feature.object){
  probabilities_list <- list()
  for(i in c(1:length(feature.object$tissue_priors))){
    prior <- feature.object$tissue_priors[i]
    if(!(is.na(prior)) & !(prior == 0)){
      probability_vector <- calculateProbabilities(feature.object,prior)
      tissue <- names(feature.object$tissue_priors)[i]
      probabilities_list[[tissue]] <- probability_vector
    }
  }
  return(probabilities_list)
}


# class to keep results for a binary feature
BinaryFeature <- function(name="",
                          feature_values=c(),
                          pos_samples=c(),
                          neg_samples=c(),
                          pos_distribution_raw=c(),
                          neg_distribution_raw=c(),
                          breaks = c(),
                          pos_distribution=c(),
                          neg_distribution=c(), 
                          pseudocounts = 0,
                          log_odds_ratio=c(), 
                          general_prior=0,
                          probability_general_prior=c(), 
                          tissue_priors=c(),
                          probabilities_tissue_priors=list(), 
                          separation_raw=0,
                          separation_raw_pvalue=0,
                          separation_raw_pvalue_fdrcorrected=0,
                          separation_raw_direction=0,
                          separation=0,
                          separation_pvalue=0,
                          separation_pvalue_fdrcorrected=0,
                          separation_direction=0) {
  me <- list(
    name = name,
    feature_values = feature_values,
    pos_samples = pos_samples,
    neg_samples = neg_samples,
    pos_distribution_raw = pos_distribution_raw,
    neg_distribution_raw = neg_distribution_raw,
    breaks = breaks,
    pos_distribution = pos_distribution,
    neg_distribution = neg_distribution,
    pseudocounts = pseudocounts,
    separation_raw = separation_raw,
    separation_raw_pvalue = separation_raw_pvalue,
    separation_raw_pvalue_fdrcorrected = separation_raw_pvalue_fdrcorrected,
    separation_raw_direction = separation_raw_direction,
    separation = separation,
    separation_pvalue = separation_pvalue,
    separation_pvalue_fdrcorrected = separation_pvalue_fdrcorrected,
    separation_direction = separation_direction,
    log_odds_ratio = log_odds_ratio,   
    general_prior = general_prior,
    probability_general_prior = probability_general_prior,  
    tissue_priors = tissue_priors,
    probabilities_tissue_priors = probabilities_tissue_priors
  )
  
  ## Set the name for the class
  class(me) <- append(class(me),"BinaryFeature")
  return(me)
}
set_pos_samples <- function(feature.object, newValue){
  feature.object$pos_samples <- newValue
  return(feature.object)
}
set_neg_samples <- function(feature.object, newValue){
  feature.object$neg_samples <- newValue
  return(feature.object)
}
set_pos_distribution_raw <- function(feature.object, newValue){
  feature.object$pos_distribution_raw <- newValue
  return(feature.object)
}
set_neg_distribution_raw <- function(feature.object, newValue){
  feature.object$neg_distribution_raw <- newValue
  return(feature.object)
}
set_separation_raw <- function(feature.object, newValueList){
  feature.object$separation_raw <- newValueList[[1]]
  feature.object$separation_raw_pvalue <- newValueList[[2]]
  feature.object$separation_raw_direction <- newValueList[[3]]
  return(feature.object)
}
set_separation <- function(feature.object, newValueList){
  feature.object$separation <- newValueList[[1]]
  feature.object$separation_pvalue <- newValueList[[2]]
  feature.object$separation_direction <- newValueList[[3]]
  return(feature.object)
}
set_separation_pvalue_fdrcorrected <- function(feature.object, newValue){
  feature.object$separation_pvalue_fdrcorrected <- newValue
  return(feature.object)
}
set_separation_raw_pvalue_fdrcorrected <- function(feature.object, newValue){
  feature.object$separation_raw_pvalue_fdrcorrected <- newValue
  return(feature.object)
}
set_breaks <- function(feature.object, newValue){
  feature.object$breaks <- newValue
  return(feature.object)
}
set_pseudocounts <- function(feature.object, newValue){
  feature.object$pseudocounts <- newValue
  return(feature.object)
}
set_pos_distribution <- function(feature.object, newValue){
  feature.object$pos_distribution <- newValue
  return(feature.object)
}
set_neg_distribution <- function(feature.object, newValue){
  feature.object$neg_distribution <- newValue
  return(feature.object)
}
set_log_odds_ratio <- function(feature.object, newValue){
  feature.object$log_odds_ratio <- newValue
  return(feature.object)
}
set_general_prior <- function(feature.object, newValue){
  feature.object$general_prior <- newValue
  return(feature.object)
}
set_tissue_priors <- function(feature.object, newValue){
  feature.object$tissue_priors <- newValue
  return(feature.object)
}
set_probability_general_prior <- function(feature.object, newValue){
  feature.object$probability_general_prior <- newValue
  return(feature.object)
}
set_probabilities_tissue_priors <- function(feature.object, newValue){
  feature.object$probabilities_tissue_priors <- newValue
  return(feature.object)
}

analyzeFeature <- function(feature.object,sample_ranking_mat, disease_annotation=c()){
  UseMethod("analyzeFeature",feature.object)
}
# Analyze one binary feature
# input: feature values as named vector, the sample ranking matrix
# return: 
analyzeFeature.BinaryFeature <- function(feature.object, sample_ranking_mat, disease_annotation=c()){
  feature_values <- feature.object$feature_values
  pos_samples <- names(feature_values[which(feature_values == 1)])
  feature.object <- set_pos_samples(feature.object,pos_samples)
  neg_samples <- names(feature_values[which(!(names(feature_values) %in% pos_samples))])
  feature.object <- set_neg_samples(feature.object,neg_samples)
  
  pos_KSdistribution <- getKSDistribution(pos_samples, sample_ranking_mat, pos_samples)
  feature.object <- set_pos_distribution_raw(feature.object,pos_KSdistribution)
  neg_KSdistribution <- getKSDistribution(neg_samples, sample_ranking_mat, pos_samples)
  feature.object <- set_neg_distribution_raw(feature.object,neg_KSdistribution)
  
  feature.object <- set_separation_raw(feature.object,testSeparationRaw(feature.object))
  
  breaks <- c(0:20)/20
  feature.object <- set_breaks(feature.object,breaks)
  neg_histogram <- hist(neg_KSdistribution,breaks = breaks,plot=F)
  pos_histogram <- hist(pos_KSdistribution,breaks = breaks,plot=F)
  
  neg_smoothed <- smoothHistogram(neg_histogram$counts)
  pseudocounts <- 10
  feature.object <- set_pseudocounts(feature.object, pseudocounts)
  pos_pseudo <- addPseudocounts(pos_histogram, neg_smoothed, pseudocounts)
  pos_pseudo_smoothed <- smoothHistogram(pos_pseudo)
  
  #pos_histogram$counts <- pos_pseudo_smoothed
  #plot(pos_histogram)
  #neg_histogram$counts <- neg_smoothed
  #plot(neg_histogram)
  
  feature.object <- set_pos_distribution(feature.object,pos_pseudo_smoothed)
  feature.object <- set_neg_distribution(feature.object,neg_smoothed)
  feature.object <- set_separation(feature.object,testSeparation(feature.object))
  
  
  feature.object <- set_log_odds_ratio(feature.object,calculateLogOddsRatio(feature.object))
  
  general_prior <- length(pos_samples)/(length(neg_samples)+length(pos_samples))
  feature.object <- set_probability_general_prior(feature.object,calculateProbabilities(feature.object,general_prior))

  feature.object <- set_tissue_priors(feature.object, getTissuePriors(feature.object,disease_annotation))
  feature.object <- set_probabilities_tissue_priors(feature.object, getTissueProbabilities(feature.object))
  
  return(feature.object)
}


# Extract KS-test distances and p-values from a list of feature objects
# Input: List of feature objects
# Output: Matrix with KS-test distance, p-value, FDR-corrected p-value as columns (sorted by FDR-corrected p-value) and features as rows
getSeparationValues <- function(feature_list){
  feature_name <- c()
  sep <- c()
  pval <- c()
  for(i in c(1:length(feature_list))){
    feature_name <- c(feature_name,strsplit(feature_list[[i]]$name,"_mut",fixed=T)[[1]][1])
    sep <- c(sep,feature_list[[i]]$separation * feature_list[[i]]$separation_direction)
    pval <- c(pval,feature_list[[i]]$separation_pvalue)
  }
  pval_fdrcorrected <- p.adjust(pval, method = "fdr")
  mat <- cbind(sep,pval,pval_fdrcorrected)
  rownames(mat) <- feature_name
  mat <- mat[order(pval_fdrcorrected),]
  return(mat)
}
