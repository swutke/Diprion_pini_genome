.libPaths(c("/path/to/R_packages_lib/", .libPaths()))
libpath <- .libPaths()[1]


# # # # # # # UpSet plot of OrthoFinder results  # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


### install (if necessary) and load packages
# install.packages("UpSetR")
library(UpSetR)

### set path to working directory
setwd("/path/to/wd/")

### get the data = OrthoFinder output
orthogroups <- read.table("Orthogroups.GeneCount.tsv", header=T, sep="\t", stringsAsFactors = F)
orthogroups[orthogroups > 0] <- 1

colnames(orthogroups) <- c("Orthogroup", "D. pini", "D. similis", "N. fabricii", "N. lecontei", "N. pinetum", "N. virginiana", "A. rosae", "O. abietus", "A. cephalotes", "P. dominula", "A. mellifera", "Total")

orthogroups <- orthogroups[,c(1,2,3,5,6,4,7,8,9,10,12,11,13)]
colnames(orthogroups)

selected_species <- colnames(orthogroups)[2:(ncol(orthogroups) -1)] 

upset(orthogroups, nsets = ncol(orthogroups), sets = rev(selected_species), keep.order = T, order.by = "freq", text.scale = 2, show.numbers = "no", set_size.scale_max = 12500)
