### just need to put together L1s 
### get a boxplot going where we look at the distribution at various locations in intergenic space 
# this means overlapping regions with TEs 

## it is likly we have code to handle this already, its just being able to find it, I could use the boundary wave stuff
## it takes metadata scores and maybe makes a list of them.
## we want a list of elemement scores for each data point
## this way we can plot fluctuations 
## how can we write this so we are not overlapping other regions or does it nor matter given the stats. 

## another way to start will be to look at average age of elemetent correlated wuth GC or gene number





rm(list = ls())
setwd("~/Desktop/Domain_manuscript/")
library(GenomicRanges)
library(rtracklayer)


source("Domain_manuscript_scripts/functions.R")
source("Domain_manuscript_scripts/rep_db.R")
load("R_objects/chromStateCombined")


genome = "hg19"
spec1 = "Human"



rep <- rep_info(spec1 = spec1, genome = genome)





# L1s <- rbind(rep$L1ME, rep$L1MD, rep$L1MC, rep$L1MB, rep$L1MA,rep$L1PB, rep$L1PA, rep$L1HS)
# L1s$newClass = c(rep("L1ME", nrow(rep$L1ME)), rep("L1MD", nrow(rep$L1MD)), rep("L1MC", nrow(rep$L1MC)), rep("L1MB", nrow(rep$L1MB)), rep("L1MA", nrow(rep$L1MA)),
#                 rep("L1PB", nrow(rep$L1PB)), rep("L1PA", nrow(rep$L1PA)), rep("L1HS", nrow(rep$L1HS)))

L1s <- rbind(rep$AluJ, rep$AluS, rep$AluY)
#L1s$newClass = c(rep("L1ME", nrow(rep$L1ME)), rep("L1MD", nrow(rep$L1MD)), rep("L1MC", nrow(rep$L1MC)), rep("L1MB", nrow(rep$L1MB)), rep("L1MA", nrow(rep$L1MA)),
             #    rep("L1PB", nrow(rep$L1PB)), rep("L1PA", nrow(rep$L1PA)), rep("L1HS", nrow(rep$L1HS)))

# L1s <- rbind(rep$L1MA,rep$L1PB, rep$L1PA, rep$L1HS)
# L1s$newClass = c( rep("L1MA", nrow(rep$L1MA)),
# rep("L1PB", nrow(rep$L1PB)), rep("L1PA", nrow(rep$L1PA)), rep("L1HS", nrow(rep$L1HS)))


Rlist <- list(L1 = L1s)

SampGenomeMed <- localGCgenome(repList = Rlist,binSize = 10000,sampSize = 10000, genome = genome, repType = "repeats")

#SampGenomeMed <- SampGenomeMed[SampGenomeMed$L1>0,]

# GC reflects gene content 

L1.gr <- GRanges(seqnames = Rle(L1s$genoName), 
                 ranges = IRanges(start = L1s$genoStart, end = L1s$genoEnd), 
                 misMatch = rowSums(L1s[,c("milliDiv", "milliDel", "milliIns")]),
              #   newClass = L1s$newClass
                 )
sampGenomeMed.gr <- GRanges(seqnames = Rle(SampGenomeMed$chr), 
                            ranges = IRanges(start = SampGenomeMed$start, end = SampGenomeMed$end))


fOL <- as.matrix(findOverlaps(sampGenomeMed.gr, L1.gr))

# number of scores or number fo bases that belong to each score, then it becomes a weighted mean 

fOL.misMatch <- data.frame(fOL[,1], elementMetadata(L1.gr)$misMatch[fOL[,2]])
colnames(fOL.misMatch) = c("sampID", "misMatch")
fOL.misMatch$sampID = as.factor(fOL.misMatch$sampID)

fOL.misMatch$GC <- SampGenomeMed$GC[fOL[,1]]
fOL.misMatch$width <- width(L1.gr)[fOL[,2]]



# the average TE base belongs to an element of 


plot(fOL.misMatch$GC, fOL.misMatch$misMatch, cex = .2, pch = 16, ylim = c(0,500))
lo <- loess( fOL.misMatch$misMatch ~ fOL.misMatch$GC, weights = fOL.misMatch$width)
lines(seq(.3,.6,.001), predict(lo,newdata = seq(.3,.6,.001)), lwd = 3)

par(new = T)
plot(fOL.misMatch$GC, fOL.misMatch$misMatch, cex = .2, pch = 16, type = "n", ylim = c(0,1))
points(SampGenomeMed$GC, SampGenomeMed$L1, col = 2, pch = 16, cex = .2)
lo.L1 <- loess(SampGenomeMed$L1 ~ SampGenomeMed$GC)
lines(seq(.3,.6,.001), predict(lo.L1,newdata = seq(.3,.6,.001)), col = 2, lwd = 3)




fOL.weigthed <- NULL
for(i in 1:nrow(fOL.misMatch)){
  fOL.weigthed <- rbind(fOL.weigthed, matrix(fOL.misMatch[i,c(2,3)], nrow = as.integer(fOL.misMatch[i,"width"]/100), ncol = 2, byrow = T))
}

weightedDF <- data.frame(bin = cut(unlist(fOL.weigthed[,2]), breaks = 100), GC = unlist(fOL.weigthed[,2]), misMatch = unlist(fOL.weigthed[,1]))

bp <- boxplot(weightedDF$misMatch ~ weightedDF$bin, outline = F)

weightedDF = weightedDF[weightedDF$misMatch < 800,]

library(ggplot2)
gp <- ggplot(data = weightedDF, aes(x = bin, y = (misMatch)))
gp + geom_violin(trim = F  , scale = "width")


plot(density(weightedDF$misMatch[weightedDF$bin == "(0.333,0.349]"]))

