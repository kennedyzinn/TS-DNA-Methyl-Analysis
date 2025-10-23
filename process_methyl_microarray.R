process_methyl_microarray <- function(rgSet){
  
  # Detect p values
  p_values <- detectionP(rgSet, type = "m+u")
  cat("The number of probes with at least one sample having a detection p-value above 0.01 is:", nrow(p_values[rowSums(p_values>0.01)>0,]))
  
  mSet <- preprocessFunnorm(rgSet, sex = rep(0, nrow(sample_sheet)))
  
  # ensure probes are in the same order in the mSetSq and detP objects
  p_values <- p_values[match(featureNames(mSet),rownames(p_values)),]
  
  # remove any probes that have failed in one or more samples; this next line
  # checks for each row of p_values whether the number of values < 0.01 is equal
  # to the number of samples (TRUE) or not (FALSE)
  keep <- rowSums(p_values < 0.01) == ncol(mSet)
  table(keep)
  
  # Subset the GenomicRatioSet
  mSet <- mSet[keep,]
  
  # Remove probes with known SNPs
  gmSet <- gmSet <- dropLociWithSnps(mapToGenome(mSet))
  
  # Remove cross reactive probes
  xreactive_probes <- xreactive_probes(array_type="EPIC")
  keep <- !(featureNames(gmSet) %in% xreactive_probes)
  gmSet <- gmSet[keep,]
  
  # make gmSet available in the global environment so M values and Beta values can be saved
  return(gmSet)
}
