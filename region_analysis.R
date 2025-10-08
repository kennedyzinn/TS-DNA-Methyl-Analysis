region_analysis <- function(contrast_coeff, betacutoff, fdr){
  myAnnotation <- cpg.annotate(datatype = "array", object = Mval, what = "M", 
                               arraytype = c("EPICv1"), 
                               analysis.type = "differential", design = fullSv, contrasts = T, cont.matrix = contrast,
                               fdr = fdr, coef = contrast_coeff)
  DMRs <- suppressMessages(dmrcate(myAnnotation, C=2, betacutoff = betacutoff))
  results.ranges <- extractRanges(DMRs, genome = "hg38")
  return(results.ranges)
}
