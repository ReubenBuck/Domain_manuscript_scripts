# this can be the fig 1 pipeline 



library(GenomicRanges)
library(rtracklayer)





# set up the bins


# write a function that can pull out bin sizes we want

# we pull the bins we want to build the biplots

# we can pull in sorted repeat files 


rm(list = ls())

library(rtracklayer)

source("~/Desktop/Domain_manuscript/Domain_manuscript_scripts/rep_db.R")
source(file="~/Desktop/Domain_manuscript/Domain_manuscript_scripts/functions.R")
options(stringsAsFactors=FALSE)


# pull each species
bin.size = 1000000
speciesDF <- data.frame(spec = c("Human", "Dog", "Chimp", "Mouse", "Rhesus"), genome = c("hg19", "canFam3","panTro4", "mm9", "rheMac3"))


for(s in 1:(nrow(speciesDF))){
  rep <- rep_info(spec1=speciesDF$spec[s], genome=speciesDF$genome[s])
  # remove unplaced chromosomes from repeats
  TE.names <- names(rep)
  for(i in TE.names){
    if(length(grep("_", rep[[i]]$genoName)) > 0){
      rep[[i]] <- rep[[i]][-(grep("_", rep[[i]]$genoName)),]
    }
  }
  bins <- binned.genome.reader(genome=speciesDF$genome[s], bin.size=bin.size, keep.rate=.9)
  bins <- bins[[1]]
  # remove unplaced chromosomes from bins
  if(length(grep("_", bins$chr)) > 0){
    bins <- bins[-(grep("_", bins$chr)),]
  }
  bin.sort = binSort(repList = rep, bins=bins,TE.names = names(rep), repType = rep("repeats", length(rep)))
  assign(paste(speciesDF$spec[s], "covCount", sep = "_"), bin.sort$counts)
  assign(paste(speciesDF$spec[s], "covRate", sep = "_"), bin.sort$rates)
  assign(paste(speciesDF$spec[s], "repInfo", sep = "_"), rep)
}


# some rules


# read in the pca plots
# find the ancestral PC
# maybe also save proportion of variation explained. 


for(s in 1:nrow(speciesDF)){
  bin.rate <- get(paste(speciesDF$spec[s], "covRate", sep = "_"))
  pca <- prcomp(bin.rate[,5:ncol(bin.rate)], scale.=T)
  repStruct <- repFamStruct(speciesDF$spec[s])
  agg <- aggregate(pca$rotation[repStruct$TEname,1:2], by=list(repStruct$repType), FUN = mean)
  rownames(agg) <- agg[,1]
  agg <- agg[,c("PC1", "PC2")]
  sinePC <- colnames(agg)[max(abs(agg["new_SINE",])) == abs(agg["new_SINE",])]
  ancPC <- colnames(agg)[colnames(agg)!=sinePC]
#  ancPC <- colnames(agg)[max(abs(agg["ancient",])) == abs(agg["ancient",])]
 
#  if(ancPC == sinePC){
    # if they are poth the same we pick the anc PC that shows the biggest difference between ancestra elments
    # and new lines
#     ancPC <- colnames(agg)[max(abs(agg["new_LINE",] - agg["ancient",])) == abs(agg["new_LINE",] - agg["ancient",])]
#     sinePC <- colnames(agg)[colnames(agg)!=ancPC]
#  }
  sumPCA <- summary(pca)
  bin.ratePCA <- list(binInfo = data.frame(bin.rate[,1:4]),
                      x = data.frame(ancient_PC = pca$x[,ancPC] * (agg["ancient",ancPC] / abs(agg["ancient",ancPC])), 
                                     new_SINE_PC = pca$x[,sinePC] * (agg["new_SINE",sinePC] / abs(agg["new_SINE",sinePC]))
                      ),
                      rotation = data.frame(ancient_PC = pca$rotation[,ancPC] * (agg["ancient",ancPC] / abs(agg["ancient",ancPC])), 
                                            new_SINE_PC = pca$rotation[,sinePC] * (agg["new_SINE",sinePC] / abs(agg["new_SINE",sinePC]))
                      ),
                      importance = data.frame(ancient_PC = sumPCA$importance[2,ancPC], 
                                              new_SINE_PC = sumPCA$importance[2,sinePC]),
                      data = data.frame(bin.rate[,5:ncol(bin.rate)])
  )
  assign(paste(speciesDF$spec[s], "PCA", sep = ""), bin.ratePCA)
}



reuben.biplot(x=HumanPCA$x, y = HumanPCA$rotation, x.col = 8)
reuben.biplot(x=DogPCA$x, y = DogPCA$rotation, x.col = 8)
reuben.biplot(x=ChimpPCA$x, y = ChimpPCA$rotation, x.col = 8)
reuben.biplot(x=MousePCA$x, y = MousePCA$rotation, x.col = 8)
reuben.biplot(x=RhesusPCA$x, y = RhesusPCA$rotation, x.col = 8)


objectsPCA <- paste(speciesDF$spec, "PCA", sep = "")
objectsRep <- paste(speciesDF$spec, "_repInfo", sep = "")
save(list=objectsPCA, file="~/Desktop/Domain_manuscript/R_objects/PCA_species")
save(list=objectsRep, file="~/Desktop/Domain_manuscript/R_objects/Rep_info_species")



# now we are up to the stage where we have stable PCA plots
# and all the necesary information

# the next step is to read in the alignmnet data and convert the coordinates

# take out the repeats

# and work out which bin belongs to which bin. 

# we only have 2 dimensions now so it will be easier to do ks test. 


# work on each species and pull them all in and run each one through the rep_dp 


# if we get things to a stage where they can be plotted elsewhere 
# save R objects that can be read back in and we can play with the plotting features somewhere else
# to do PCA all we need is bin_info 


# also we should only remove repeats we are interested in from our pairwise stuff

















