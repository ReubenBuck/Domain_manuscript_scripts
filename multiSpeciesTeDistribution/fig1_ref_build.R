
# here is somewhere to run through all the speceis and remodel them based on their coorinates
# no more long readign in and processing of data
# it can all be done on mannagable df



# also worth thinking about is the cutoff levels 

# and permution analysis on correlations. 


rm(list = ls())

load(file="~/Desktop/Domain_manuscript/R_objects/binMaps")
load("~/Desktop/Domain_manuscript/R_objects/PCA_species")


source(file="~/Desktop/Domain_manuscript/Domain_manuscript_scripts/functions.R")



referenceDF <- data.frame(spec = "Human", genome = "hg19")
speciesDF <- data.frame(spec = c( "Dog", "Chimp", "Mouse", "Rhesus"), genome = c("canFam3","panTro4", "mm9", "rheMac3"))


for(s in speciesDF$spec){
  assign(s,mapRemodeler(refSpec="Human", queSpec = s, cutoff = .05))
}

 
pdf(file="~/Desktop/Domain_manuscript/plots/FIg1/PCA_wiggle/wiggle.pdf", width=20, height=10)
layout(matrix(c(1,2)))
par(mar=c(2,5,5,5))

plotdf <- smooth.spline(HumanPCA$binInfo$start[HumanPCA$binInfo$chr == "chr1"], as.numeric((HumanPCA$x$new_SINE_PC[HumanPCA$binInfo$chr == "chr1"])))

plot(plotdf$x,
     plotdf$y, 
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
    plotdf <- smooth.spline(start, datQue,all.knots = T)
  lines(plotdf$x, plotdf$y, col = s + 1)
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


# so lets read in replication domains


RDs <- read.table(file = "~/Desktop/Domain_manuscript/Data/DNN_HMM_repliDomains/GSE53984_GSM923452_Huvec_Rep1_segments.bed", 
                  col.names = c("chr", "start", "end", "value"), colClasses = c("character", "integer", "integer", "character"))
RDs$value[RDs$value == "ERD"] = 1
RDs$value[RDs$value == "LRD"] = -1
RDs$value[RDs$value == "UTZ"] = NA
RDs$value[RDs$value == "DTZ"] = NA
RDs <- RDs[complete.cases(RDs),]    
RDs$value <- as.numeric(RDs$value)
RDs <- RDs[RDs$end - RDs$start + 1 >= 2000000,]

library("circlize")

layout(1)
par(mar = c(1,1,1,1))



pdf(file = "~/Desktop/Domain_manuscript/plots/FIg1/ancPCAcircos.pdf",onefile = T)

circos.initializeWithIdeogram(species = "hg19", plotType = c( "labels"), sort.chr = T)
circos.par("track.height" = .05)
circos.genomicTrackPlotRegion(data = RDs, bg.border = "white", panel.fun = function(region,value, ...){
  circos.genomicRect(region,value,col = ifelse(value[[1]] > 0 , "red","white"), border = NA, ...)
})

for(s in nrow(speciesDF):1){
  List <- get(as.character(speciesDF$spec[s]))
  coords = get(paste("HumanRef_", speciesDF$spec[s],"Que", sep = "" ))
  coords <- coords$HumanRef$binInfo[List$ancient_PC$remodeldPC$refNo,]
  datQue <- List$ancient_PC$remodeldPC$quePC
  bedS <- data.frame(coords[,1:3], value = datQue)
  bedS = bedS[order(bedS$value),][c(1:as.integer(.2*nrow(bedS)), (nrow(bedS)-(as.integer(.2*nrow(bedS))-1)):nrow(bedS)),]
  circos.par("track.height" = .08)
  circos.genomicTrackPlotRegion(data = bedS, panel.fun = function(region,value, ...){
    circos.genomicRect(region,value,col = ifelse(value[[1]] > 0, "coral3", "aquamarine4"), border = NA, ...)
  })
  
}

bedH = data.frame(HumanPCA$binInfo[,1:3], value = (HumanPCA$x$ancient_PC))
bedH = bedH[order(bedH$value),][c(1:as.integer(.2*nrow(bedH)), (nrow(bedH)-(as.integer(.2*nrow(bedH))-1)):nrow(bedH)),]
circos.par("track.height" = .08)
circos.genomicTrackPlotRegion(data = bedH, panel.fun = function(region,value, ...){
  circos.genomicRect(region,value,col = ifelse(value[[1]] > 0, "coral3", "aquamarine4"), border = NA, ...)
})
dev.off()



pdf(file = "~/Desktop/Domain_manuscript/plots/FIg1/newPCAcircos.pdf",onefile = T)
circos.initializeWithIdeogram(species = "hg19", plotType = c( "labels"))
circos.par("track.height" = .05)
circos.genomicTrackPlotRegion(data = RDs, bg.border = "white", panel.fun = function(region,value, ...){
  circos.genomicRect(region,value,col = ifelse(value[[1]] > 0 , "red","white"), border = NA, ...)
})

for(s in nrow(speciesDF):1){
  List <- get(as.character(speciesDF$spec[s]))
  coords = get(paste("HumanRef_", speciesDF$spec[s],"Que", sep = "" ))
  coords <- coords$HumanRef$binInfo[List$new_SINE_PC$remodeldPC$refNo,]
  datQue <- List$new_SINE_PC$remodeldPC$quePC
  bedS <- data.frame(coords[,1:3], value = datQue)
  bedS = bedS[order(bedS$value),][c(1:as.integer(.2*nrow(bedS)), (nrow(bedS)-(as.integer(.2*nrow(bedS))-1)):nrow(bedS)),]
  circos.par("track.height" = .08)
  circos.genomicTrackPlotRegion(data = bedS, panel.fun = function(region,value, ...){
    circos.genomicRect(region,value,col = ifelse(value[[1]] > 0, "coral3", "aquamarine4"), border = NA, ...)
  })
  
}

bedH = data.frame(HumanPCA$binInfo[,1:3], value = (HumanPCA$x$new_SINE_PC))
bedH = bedH[order(bedH$value),][c(1:as.integer(.1*nrow(bedH)), (nrow(bedH)-(as.integer(.1*nrow(bedH))-1)):nrow(bedH)),]
circos.par("track.height" = .08)
circos.genomicTrackPlotRegion(data = bedH, panel.fun = function(region,value, ...){
  circos.genomicRect(region,value,col = ifelse(value[[1]] > 0, "coral3", "aquamarine4"), border = NA, ...)
})

plot.new()
legend("center", legend = c("> 80", "< 20"), fill = c("coral3", "aquamarine4"), title = "percentile",bty = "n")
legend("bottom", legend = c("ERD"), fill = c("red"), title = "domain", bty = "n")


dev.off()














