
# here is somewhere to run through all the speceis and remodel them based on their coorinates
# no more long readign in and processing of data
# it can all be done on mannagable df



# also worth thinking about is the cutoff levels 

# and permution analysis on correlations. 


rm(list = ls())

load(file="~/Desktop/Domain_manuscript/R_objects/binMaps")
load("~/Desktop/Domain_manuscript/R_objects/PCA_species")


source(file="~/Desktop/element_curves/element_curves_scripts/functions.R")



referenceDF <- data.frame(spec = "Human", genome = "hg19")
speciesDF <- data.frame(spec = c( "Dog", "Chimp", "Mouse", "Rhesus"), genome = c("canFam3","panTro4", "mm9", "rheMac3"))






for(s in speciesDF$spec){
  assign(s,mapRemodeler(refSpec="Human", queSpec = s, cutoff = 150000))
}

 
pdf(file="~/Desktop/Domain_manuscript/plots/FIg1/PCA_wiggle/wiggle.pdf", width=20, height=10)
layout(matrix(c(1,2)))
par(mar=c(2,5,5,5))

plot(HumanPCA$binInfo$start[HumanPCA$binInfo$chr == "chr1"],
     as.numeric((HumanPCA$x$new_SINE_PC[HumanPCA$binInfo$chr == "chr1"])), 
     type="l",
     xlim = c(0,10^8),
     xlab = "chr1",
     ylab = "new_SINE_PC")
for(s in 1:nrow(speciesDF)){
  List <- get(as.character(speciesDF$spec[s]))
  coords = get(paste("HumanRef_", speciesDF$spec[s],"Que", sep = "" ))
  coords <- coords$HumanRef$binInfo[List$new_SINE_PC$remodeldPC$refNo,]
  datQue <- List$new_SINE_PC$remodeldPC$quePC[coords$chr == "chr1"]
  start = coords$start[coords$chr == "chr1"]
  lines(start, as.numeric((datQue)), col = s + 1)
}
legend("topright",c("Human", as.character(speciesDF$spec)), fill = 1:5)
abline(h = 0, lty = 2)

par(mar=c(5,5,2,5))
plot(HumanPCA$binInfo$start[HumanPCA$binInfo$chr == "chr1"],
     as.numeric((HumanPCA$x$ancient_PC[HumanPCA$binInfo$chr == "chr1"])), 
     type="l",
     xlim = c(0,10^8),
     xlab = "chr1",
     ylab = "ancient_PC")
for(s in 1:nrow(speciesDF)){
  List <- get(as.character(speciesDF$spec[s]))
  coords = get(paste("HumanRef_", speciesDF$spec[s],"Que", sep = "" ))
  coords <- coords$HumanRef$binInfo[List$ancient_PC$remodeldPC$refNo,]
  datQue <- List$ancient_PC$remodeldPC$quePC[coords$chr == "chr1"]
  start = coords$start[coords$chr == "chr1"]
  lines(start, as.numeric((datQue)), col = s + 1)
}
#legend("topright",c("Human", as.character(speciesDF$spec)), fill = 1:5)
abline(h = 0, lty = 2)

dev.off()

    # we need to get rid of all the human ref stuff
    
    
# so for now we will just do a cross section of chromosome1


# else where I'll plot the PCA with coloured arrows. 
    
    


hist(corDist)


probFunc <- ecdf(corDist)
probFunc(cor(remodeld.ANC[,2], remodeld.ANC[,3])) - 1

cor(remodeld.ANC[,2:3])













