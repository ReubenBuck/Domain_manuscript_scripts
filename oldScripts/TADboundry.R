### the killobase scale of the genome at 50 kb bins? or at different sizes based on where we get our base approximations 

# so what do types matter?
# or should we look at levels of heterochromatin accross the boundries
# or shoudl we go with actuall coverage levels 
# maybe include the density of genes 
# or classes of genes if we bin it we can make the matrix smaller 
# get the 10kb level and look each side of the divide and compare to genome averge to under or over. 
# how far are we from getting domains?
# there is always the possibility of using lifted ones?
# for now lets take the human H1esc domain boundries 
# or go to regulatory domains?

#lets look at the cell type conserved domain boundries, look at accumulation of TEs either side of the line
# there are aboout 1000 conserved boundries and lets say we look at 10 kb windows, that mean we should have 50 windows wither side


# How to get the conserved boundries 
# maybe read in both data sets and pull out the gaps and look to see how they overlap
rm(list = ls())

library(GenomicRanges)
library(rtracklayer)

setwd("~/Desktop/Domain_manuscript/")

source("Domain_manuscript_scripts/functions.R")
source("Domain_manuscript_scripts/rep_db.R")
load("R_objects/chromStateCombined")


genome = "hg19"
spec1 = "Human"



rep <- rep_info(spec1 = spec1, genome = genome)


joinRep <- list(old_L1 = rbind(rep$L1ME, rep$L1MD, rep$L1MC, rep$L1MB),
                new_L1 = rbind(rep$L1MA, rep$L1PB, rep$L1PA, rep$L1HS),
                Alu = rbind(rep$AluJ, rep$AluS, rep$AluY),
                Ancient = rbind(rep$MIR, rep$L2)
                )
ChromoChoice = "R"
joinChromo <- c(list(`H1-hESC`[[ChromoChoice]]), list(HUVEC[[ChromoChoice]]), list(HepG2[[ChromoChoice]]), list(GM12878[[ChromoChoice]]), list(`HeLa-S3`[[ChromoChoice]]), list(K562[[ChromoChoice]]))
names(joinChromo) <- c("H1hESC", "HUVEC", "HepG2", "GM12878", "HelaS3", "K562")
names(joinChromo) <- paste(names(joinChromo), ChromoChoice, sep = "_")
joinRepChromatin <- c(joinRep, joinChromo)



if(spec1 == "Human"){
  escTad <- read.table("Data/TD_data_one/Human/hESC/combined/hglft_genome_42f6_7dfa40.bed")
  tisTad <- read.table("Data/TD_data_one/Human/IMR90/combined/hglft_genome_461b_7e0e40.bed")
}else if(spec1 == "Mouse"){
  escTad <- read.table("Data/TD_data_one/Mouse/mESC/HindIII_combined/total.HindIII.combined.domain")
  tisTad <- read.table("Data/TD_data_one/Mouse/cortex/combined/total.combined.domain")
}


esc.gr <- disjoin(
  GRanges(seqnames=Rle(escTad[,1]), 
                  ranges = IRanges(start = escTad[,2], end = escTad[,3])
  )
)
esc.gr <- esc.gr[width(esc.gr) > 1]
escBoundry <- gaps(esc.gr)
escBoundry <- escBoundry[width(escBoundry) < 100]


tis.gr <- disjoin(
  GRanges(seqnames = Rle(tisTad[,1]),
                  ranges = IRanges(start = tisTad[,2], end=tisTad[,3])
                  )
)
tis.gr <- tis.gr[width(tis.gr) > 1]
tisBoundry <- gaps(tis.gr)
tisBoundry <- tisBoundry[width(tisBoundry) < 100]


# so now we know the boudry regions maybe we can overlap them and get our conserved levels 
boundryOL <- as.matrix(findOverlaps(tisBoundry,escBoundry))

nrow(boundryOL)
length(unique(boundryOL[,1]))
length(unique(boundryOL[,2]))
length(tisBoundry)
length(escBoundry)

# we use a 40000 max gap and we set our line in the middle 

# for each one we take the min start and max stop and see where it centers

jointBoundry <- data.frame(chr = seqnames(tisBoundry[boundryOL[,1]]))
jointBoundry$start <- apply(X = cbind(start(tisBoundry[boundryOL[,1]]), start(escBoundry[boundryOL[,2]])),FUN = min,MARGIN = 1)
jointBoundry$end <- apply(X = cbind(end(tisBoundry[boundryOL[,1]]), end(escBoundry[boundryOL[,2]])),FUN = max,MARGIN = 1)
jointBoundry$line <- (jointBoundry$start + jointBoundry$end) /2

# turn it into bins and get the rate 
# we need an indicator of which row and column each section corrsponds to 
# strat form -500000 and g up every 10kb 

binSize <- 50000
regionSize <- 500000

# there is also the potential Ill have to get GC rates
boundryAnotation <- function(repList, binSize, regionSize, repTypes, boundryLine, boundryChr){
  TErateMatrixList <- NULL
  for(i in 1:length(repList)){
    TErateMatrixList <- c(TErateMatrixList, list(matrix(nrow = length(boundryLine), ncol = regionSize/binSize)))
  }
  names(TErateMatrixList) <- names(repList)
  
  for(i in 1:(regionSize/binSize)){
    print(i)
    bin <- data.frame(chr = boundryChr, 
                      start = (boundryLine - (regionSize/2)) + ((i - 1) * binSize), 
                      end = (boundryLine - (regionSize/2)) + ((i) * binSize))
    bin$Known <- binSize
    sorted <- binSort(rep = repList, bins = bin,TE.names = names(repList), repType = repTypes )
    for(te in 1:length(repList)){
      TErateMatrixList[[(names(repList)[te])]][,i] <- sorted$rates[,names(repList)[te]]
    }
    
  }
  return(TErateMatrixList)
}

matList <- boundryAnotation(repList = joinRepChromatin,binSize = 50000, regionSize = 1000000, 
                            boundryLine = jointBoundry$line, boundryChr = jointBoundry$chr,
                            repTypes = c(rep("repeats", length(joinRep)), rep("chromatin", length(joinChromo))))





# I think we need to smothen out the scores to actually look at the regions
# they tend to be a bit stochastic 
# pretty much everything except alu 

TE_choice <- "Ancient"

matChoice <- matList[[TE_choice]]
matChoice <- t(apply(X = matChoice, MARGIN = 1, smooth))
 flipDF <- data.frame(rowSums(matChoice[,(((ncol(matChoice)/2)+1) - 10):(ncol(matChoice)/2)]), 
                      rowSums(matChoice[,((ncol(matChoice)/2)+1):((ncol(matChoice)/2) + 10)]))
 fliper <- rep(0, nrow(matChoice))
 fliper[flipDF[[2]] > flipDF[[1]]] <- 1

 
 matChoice[fliper==1,] <-  rev(matChoice[fliper==1,])
  clust_pattern <- hclust(dist(matChoice))

 

layout(matrix(1:(4 *2), nrow = 2),heights = c(1,1.5))
par(mar=c(2,2,5,2))
for(i in 1:4){
  TEmat = matList[[i]]
 #TEmat <- t(apply(X = matList[[i]], MARGIN = 1, smooth))
  TEmatUnflip <- TEmat
  flipDF <- data.frame(rowSums(TEmat[,(((ncol(TEmat)/2)+1) - 10):(ncol(TEmat)/2)]), rowSums(TEmat[,((ncol(TEmat)/2)+1):((ncol(TEmat)/2) + 10)]))
   fliper <- rep(0, nrow(TEmat))
   fliper[flipDF[,2] > flipDF[,1]] <- 1
  TEmat[fliper==1,] <- rev(TEmat[fliper==1,])
  
  clust_pattern <- hclust(dist(TEmat))
  TEmat <- TEmat[clust_pattern$order,]
  par(mar=c(0,2,5,2))
  plot(colMeans(TEmat), type = "l", xaxt = "n", xlab = "", ylab= "", main = names(matList)[i],col = 2)
  lines(colMeans(TEmatUnflip), col = 1)
  abline(v=25.5)
  par(mar=c(2,2,0,2))
  image(t(TEmat), col = grey(seq(1,0,-.005)), xaxt = "n", yaxt = "n")
  abline(v = .5)
  axis(side = 1,at = c(0,.5,1), labels = c("-500 kb", "TAD boundry", "+500 kb"))
}




# shoild probably get some error bars on here
# based on how many I tend to see per 10000 bp or whatever
# so use the variance form the sampling or compute the actual rate?








# how does this correspond to heterochromatin 

# how TEs act according to genome compartments 
# maybe this is the wrong way to be looking at it
# maybe we should do coverage type plots again
# moving from one compartment to the next




# we may need to do a flip to look at the effects
# We can either flip each one as we did with the clustering 

# or we can flip according to a program, much like the single clustering 
# flipping can be made easy by counting the total on each side
# the higer value is on the right we revrese the order, of that row

# pretty easy we can get the rev function to do it
# first reverse then cluster. 


