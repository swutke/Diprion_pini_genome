.libPaths(c("/path/to/R_packages_lib/", .libPaths()))
libpath <- .libPaths()[1]


# # # # # # # plot repeat landscape # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # #

### load packages
library(reshape)
library(ggplot2)
library(viridis)
library(hrbrthemes)
library(tidyverse)
library(gridExtra)

### set path to working directory
setwd("/path/to/wd/")

### get the data = RepeatMasker output
KimuraDistance_table <- read.csv("genome.kimuradist",sep="\t")

kd_melt = melt(KimuraDistance_table,id="Div")

kd_melt$norm = kd_melt$value/270000000 * 100   # genome_size=270000000

ggplot(kd_melt, aes(fill=variable, y=norm, x=Div)) +
  geom_bar(position="stack", stat="identity",color="black") +
  scale_fill_viridis(discrete = T) +
  theme_classic() +
  xlab("Kimura substitution level") +
  ylab("Percent of the genome") +
  labs(fill = "") +
  coord_cartesian(xlim = c(0, 70)) +
  theme(axis.text=element_text(size=11),axis.title =element_text(size=12)) +
  ggtitle("N. pinetum repeat landscape") +
  theme(plot.title=element_text(hjust=0.5, face='bold', size=20))
