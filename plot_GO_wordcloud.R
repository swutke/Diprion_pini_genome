.libPaths(c("/path/to/R_packages_lib/", .libPaths()))
libpath <- .libPaths()[1]



# # # # # # # Word cloud of functional annotation  # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


### install (if necessary) and load packages
library(wordcloud2)
library(tidyverse)
library(webshot)
#webshot::install_phantomjs()
library("htmlwidgets")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GO.db")


### set path to working directory
setwd("/path/to/wd/")

### get the data = PANNZER output
go_table <- read.csv(file = "GO.out", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
counttable <- count(go_table, desc, sort = TRUE)

bp_table <- filter(go_table, ontology == "BP")
bp_counttable <- count(bp_table, desc, sort = TRUE)

mf_table <- filter(go_table, ontology == "MF")
mf_counttable <- count(mf_table, desc, sort = TRUE)

### Make the graph
bp_graph <- wordcloud2(slice_max(bp_counttable,order_by = n, n=30), size = 0.3, shape = 'circle', color='random-dark')
mf_graph <- wordcloud2(slice_max(mf_counttable,order_by = n, n=30), size = 0.3, shape = 'circle', color='random-dark')

### save it in html
saveWidget(bp_graph,"tmp.html",selfcontained = F)
### and in pdf
webshot("tmp.html","Wordcloud_BP.pdf", delay =5, vwidth = 1000, vheight=1000)

### save it in html
saveWidget(mf_graph,"tmp.html",selfcontained = F)
### and in pdf
webshot("tmp.html","Wordcloud_MF.pdf", delay =5, vwidth = 1000, vheight=1000)
