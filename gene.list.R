gene.list <- function(methStatus){
  gene.list <- results$GeneSet[results$MethylationStatus==methStatus]
  gene.list <- mapIds(org.Hs.eg.db, keys = gene.list, column = "ENTREZID", keytype = "SYMBOL")
  return(gene.list)
}