

rm(list = ls())

library(GenomicRanges)

load("~/Desktop/Domain_manuscript/R_objects/PCA_species")
load("~/Desktop/Domain_manuscript/R_objects/Rep_info_species")
source(file="~/Desktop/Domain_manuscript/Domain_manuscript_scripts/functions.R")

### Our functions will only work in the correct environmnet
referenceDF <- data.frame(spec = "Human", genome = "hg19")
speciesDF <- data.frame(spec = c( "Dog", "Chimp", "Mouse", "Rhesus"), genome = c("canFam3","panTro4", "mm9", "rheMac3"))

# maybe write this as a function/functions? 

for(i in 1:nrow(speciesDF)){
  
  refSpecGenome <- referenceDF$genome[1] 
  queSpecGenome <- speciesDF$genome[i]
  refSpec <- referenceDF$spec[1]
  queSpec <- speciesDF$spec[i]
  
  align <- removeRepAlign(refSpecGenome=refSpecGenome, queSpecGenome=queSpecGenome, refSpec=refSpec, queSpec=queSpec)
  # now that i can remove the repeats all i have to do is get the bin map 
  align2 <- isolateBinAlign(align=align,refSpec=refSpec,queSpec=queSpec)
  binMap <- buildBinMap(align=align2,refSpec=refSpec, queSpec=queSpec)
  
  
  dat <- list(get(paste(refSpec, "PCA", sep ="")), get(paste(queSpec, "PCA", sep ="")), binMap)
  names(dat) <- c(paste(refSpec, "Ref", sep = ""), paste(queSpec, "Que", sep = ""), "binMap")
  assign(x=paste(refSpec, "Ref_", queSpec, "Que", sep = ""), dat)
}

saveList <- c(paste(referenceDF$spec[1], "Ref_", speciesDF$spec, "Que", sep = ""))
save(list=saveList, file="~/Desktop/Domain_manuscript/R_objects/binMaps")



# have a series of alignemnts thta can be used to map bins !!!
# the alignments shouldn't interfare with the bin boundries
# intersect and find the overlaps or just find Ols because boundries are resolved. 

# still Not sure if I should go back and remove the chromosomes that are unplaced across all species and start again




# now we roughly know which bins are next to each other and the bp number overlap get the bin number for each 
# calculate the fraction 
# then join them up 


# somthing is going on
# we are getting alignments ro a bin that are bigger than 1000000 bp



# this way we know whihc entries belong where


