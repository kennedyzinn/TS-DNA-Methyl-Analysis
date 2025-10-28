library(minfi)

# load the idat files
path_to_idat <- "/Users/kennedyzinn/Desktop/OneDrive - UTHealth Houston/GRA/TS/dna methyl/SPMTS2_Array data/"

rgSet <- read.metharray.exp(path_to_idat)

# Pre-processing

source("process_methyl_microarray.R")
gmSet <- process_methyl_microarray(rgSet)

saveRDS(gmSet, file = "gmSet.rds")