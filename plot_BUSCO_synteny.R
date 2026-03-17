.libPaths(c("/path/to/R_packages_lib/", .libPaths()))
libpath <- .libPaths()[1]


# # # # # # # Synteny between two genomes # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # #


### load packages
require(RIdeogram)

### set path to working directory
setwd("/path/to/wd/")

# get the data / set paths to data = BUSCO output; create karyotype file
spec1_karyotype_file <- "karyotype_in/Dsimilis_karyotype.txt" # karyotype
spec1_coord_file <- "coordinates_in/Dsimilis_busco_coordinates.txt"       # busco coordinates

spec2_karyotype_file <- "karyotype_in/Dpini_karyotype.txt" # karyotype
spec2_coord_file <- "coordinates_in/Dpini_busco_coordinates.txt"   # busco coordinates

### load karyotypes ----
spec1_karyotype <- read.csv(file=spec1_karyotype_file, sep="\t", header = TRUE, stringsAsFactors = F)
spec2_karyotype <- read.csv(file=spec2_karyotype_file, sep="\t", header = TRUE, stringsAsFactors = F)

### load coordinates ----
spec1_coord <- read.csv(file=spec1_coord_file, sep="\t", header = FALSE, stringsAsFactors = F)
spec2_coord <- read.csv(file=spec2_coord_file, sep="\t", header = FALSE, stringsAsFactors = F)

### get chromosome names
chromosomes_spec1 <- spec1_karyotype$Chr
chromosomes_spec2 <- spec2_karyotype$Chr

spec1_coord <- spec1_coord[spec1_coord$V2=="Complete", ]
spec2_coord <- spec2_coord[spec2_coord$V2=="Complete", ]

# only keep results for selected sequences
spec1_coord <- spec1_coord[spec1_coord$V3 %in% spec1_karyotype$Chr,]
spec2_coord <- spec2_coord[spec2_coord$V3 %in% spec2_karyotype$Chr,]


# assign a name to species 1 and 2 
species1 <- "Dsimilis" # assign a name to species 1
species2 <- "Dpini" # assign a name to species 2

###merge kayotypes ----
karyotype_merge <- rbind(spec1_karyotype, spec2_karyotype, make.row.names=FALSE)                             # !!! CHANGE

karyotype_merge <- karyotype_merge[, c(1,2,3,7,4,5,6)]
karyotype_merge <- as.data.frame(karyotype_merge, stringsAsFactors = F )

synteny_coord <- merge(spec1_coord, spec2_coord, by="V1")                                                    # !!! CHANGE
synteny_coord <- synteny_coord[, c(-1, -2, -6)]
colnames(synteny_coord) <- c("Species_1", "Start_1", "End_1", "Species_2", "Start_2", "End_2")
synteny_coord$fill <- "cccccc"

# convert names into integer according to position in karyotype_merge for species 1
for(seq in chromosomes_spec1) {                                                                                # !!! CHANGE
  index <- which(grepl(seq, karyotype_merge$Chr))
  synteny_coord$Species_1 <- gsub(seq, index, synteny_coord$Species_1)
}

# replace chromosome name with chromosome number in synteny coordinate table for species 2
for(seq in chromosomes_spec2) {                                                                              # !!! CHANGE
  index <- which(grepl(seq, karyotype_merge$Chr))
  index_species_1 <- index - length(chromosomes_spec1)   # length of chromosome list of species 1             # !!! CHANGE
  synteny_coord$Species_2 <- gsub(seq, index_species_1, synteny_coord$Species_2)
}

# You can assign different colors to lines linking scaffolds
synteny_coord[synteny_coord$Species_1==1,]$fill <- "E3E418" # yellow
synteny_coord[synteny_coord$Species_1==2,]$fill <- "31688E" # petrol
synteny_coord[synteny_coord$Species_1==3,]$fill <- "FF3300" # orange-red
synteny_coord[synteny_coord$Species_1==4,]$fill <- "0000FF" # blue
synteny_coord[synteny_coord$Species_1==5,]$fill <- "FF33FF" # pink
synteny_coord[synteny_coord$Species_1==6,]$fill <- "669900" # green
synteny_coord[synteny_coord$Species_1==7,]$fill <- "9933CC" # purple
synteny_coord[synteny_coord$Species_1==8,]$fill <- "FF9966" # salmon
synteny_coord[synteny_coord$Species_1==9,]$fill <- "CCFF00" # green-yellow
synteny_coord[synteny_coord$Species_1==10,]$fill <- "993300" # brown
synteny_coord[synteny_coord$Species_1==11,]$fill <- "000066" # dark blue
synteny_coord[synteny_coord$Species_1==12,]$fill <- "FFCC00" # yellow-orange
synteny_coord[synteny_coord$Species_1==13,]$fill <- "CCFFCC" # mint-green
synteny_coord[synteny_coord$Species_1==14,]$fill <- "CC0000" # red

synteny_coord <- synteny_coord[with(synteny_coord, order(fill)),]


write.table(synteny_coord, file = paste0(species1, "_", species2, "_synteny.txt"), quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
write.table(karyotype_merge, file = paste0(species1, "_", species2, "_karyotype_merge.txt"), quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

synteny_coord <- read.table(paste0(species1, "_", species2, "_synteny.txt"), sep = "\t", header = T, stringsAsFactors = F)

# Plot with ideogram function
ideogram(karyotype = karyotype_merge, synteny = synteny_coord, output = paste0(species1, "_", species2, "_synteny.svg"), width = 1000)

svg2pdf(paste0(species1, "_", species2, "_synteny.svg"), file = paste0(species1, "_", species2, "_synteny.pdf"), width = 2, height = 2, dpi = 1200)
