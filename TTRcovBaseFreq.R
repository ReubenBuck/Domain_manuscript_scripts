#### so we'll just properly fis the boundary analysis of the regions we are interested in



#### 




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


joinRep <- list(old_L1 = rbind(rep$L1ME, rep$L1MD, rep$L1MC, rep$L1MB),
                new_L1 = rbind(rep$L1PA, rep$L1HS),
                Alu = rbind(rep$AluJ, rep$AluS, rep$AluY),
                Ancient = rbind(rep$MIR, rep$L2)
)

#joinRep = rep
ChromoChoice = "R"
joinChromo <- c(list(`H1-hESC`[[ChromoChoice]]), list(HUVEC[[ChromoChoice]]), list(HepG2[[ChromoChoice]]), list(GM12878[[ChromoChoice]]), list(`HeLa-S3`[[ChromoChoice]]), list(K562[[ChromoChoice]]))
names(joinChromo) <- c("H1hESC", "HUVEC", "HepG2", "GM12878", "HelaS3", "K562")
names(joinChromo) <- paste(names(joinChromo), ChromoChoice, sep = "_")

chromoSample <- HUVEC
for(i in 1:length(chromoSample)){
  colnames(chromoSample[[i]]) <- c(colnames(chromoSample[[i]])[1],"genoName", "genoStart", "genoEnd", colnames(chromoSample[[i]])[5:10])
}
joinRepChromatin <- c(joinRep, chromoSample)


### lets genreate soem random regions and claculate our p

GenomeBins <- binned.genome.reader(genome = genome,bin.size = 1000000,keep.rate = .9)[[1]]
BinnedReps <- binSort(repList = joinRepChromatin, bins = GenomeBins,TE.names = names(joinRepChromatin),repType = c(rep("repeats", length(joinRep)), rep("repeats", length(chromoSample))))
GenomeP <- c(Alu = mean(BinnedReps$rates$Alu)/10000, new_L1 = mean(BinnedReps$rates$new_L1)/10000, old_L1 = mean(BinnedReps$rates$old_L1)/10000, Ancient = mean(BinnedReps$rates$Ancient)/10000, TSS = mean(BinnedReps$rates$TSS)/10000,
             T = mean(BinnedReps$rates$T)/10000, CTCF = mean(BinnedReps$rates$CTCF)/10000, E = mean(BinnedReps$rates$E)/10000, PF = mean(BinnedReps$rates$PF)/10000, R =mean(BinnedReps$rates$R)/10000, WE =mean(BinnedReps$rates$WE)/10000 )

hist(BinnedReps$rates$T, breaks = 100)
sd(BinnedReps$rates$new_L1)
plot(ecdf(BinnedReps$rates$new_L1))
A <- ecdf(BinnedReps$rates$new_L1)
#### this is where we get our regions 
# 

fileNames <- list.files("Data/DNN_HMM_repliDomains/")
cellName <- as.matrix(as.data.frame(strsplit(fileNames, "_"))[3,])[1:length(fileNames)]
RTdomain <- NULL
for(i in 1:length(fileNames)){
  cell <- read.table(paste("Data/DNN_HMM_repliDomains/", fileNames[i], sep = ""),header = F,
                     colClasses = c("character", "integer", "integer", "character"))
  colnames(cell) <- c("chr", "start", "end", "domain")
  RTdomain <- c(RTdomain,list(cell))
}
names(RTdomain) <- cellName

# really we want to pull out a boundry line?
# basicly the start of a down domain and the end of an up domaintwisted around
cellType <- "Huvec"

TTR <- RTdomain[[cellType]][RTdomain[[cellType]]$domain %in% c("UTZ","DTZ") ,]
Edomain <- RTdomain[[cellType]][RTdomain[[cellType]]$domain == "ERD" ,]
Ldomain <- RTdomain[[cellType]][RTdomain[[cellType]]$domain == "LRD" ,]

# RTdomains <- read.table("Data/ConsTimingDomains", header = T)
# Edomain <- RTdomains[RTdomains$domain == 1,]
# Ldomain <- RTdomains[RTdomains$domain == -1,]

InterE <- gaps(GRanges(seqnames = Rle(Edomain$chr),
                       ranges = IRanges(start = Edomain$start, end = Edomain$end)))

# InterE <- GRanges(seqnames = Rle(Ldomain$chr),
#                  ranges = IRanges(start = Ldomain$start, end = Ldomain$end))


for(te in 1:length(names(joinRepChromatin))){
  
  
  TEchoice <- names(joinRepChromatin)[te]
 # BIN <- data.frame(Ldomain[,1:3], Known = Ldomain$end- Ldomain$start + 1, Alu = 1, new_L1 = 1, old_L1 = 1, Ancient = 1, TSS = 1, T = 1, CTCF = 1 , E = 1, PF = 1, R = 1, WE = 1)
#  featCovL <- covCalcPlot(lenChoice = 1000000,repBins = BIN, repChoice = TEchoice, repList = joinRepChromatin)
  BIN <- data.frame(Edomain[,1:3], Known = Edomain$end- Edomain$start + 1, Alu = 1, new_L1 = 1, old_L1 = 1, Ancient = 1, TSS = 1, T = 1, CTCF = 1 , E = 1, PF = 1, R = 1, WE = 1)
#  BIN <- BIN[BIN$end - BIN$start > 1001,]
  featCovE <- covCalcPlot(lenChoice = 1000000,repBins = BIN, repChoice = TEchoice, repList = joinRepChromatin)
 # BIN <- data.frame(TTR[TTR$domain=="UTZ",1:3], Known = TTR[TTR$domain=="UTZ","end"]- TTR[TTR$domain=="UTZ","start"] + 1, Alu = 1, new_L1 = 1, old_L1 = 1, Ancient = 1, TSS = 1, T = 1, CTCF = 1 , E = 1, PF = 1, R = 1, WE = 1)
 # featCovUTZ <- covCalcPlot(lenChoice = 200000,repBins = BIN, repChoice = TEchoice, repList = joinRepChromatin)
  BIN <- data.frame(chr = seqnames(InterE), start = start(InterE), end = end(InterE), Known = width(InterE), Alu = 1, new_L1 = 1, old_L1 = 1, Ancient = 1, TSS = 1, T = 1, CTCF = 1 , E = 1, PF = 1, R = 1, WE = 1)
  featCovInterE <- covCalcPlot(lenChoice = 1000000,repBins = BIN, repChoice = TEchoice, repList = joinRepChromatin)
  
  
  
  #plot(featCovUTZ$rawRepCov5/featCovUTZ$baseFreq5, type = "l", ylim = c(0,.2))
  #plot(featCovUTZ$rawRepCov3/featCovUTZ$baseFreq3, type = "l", ylim = c(0,.2))
  png(filename = paste("plots/repliSup/", TEchoice, "_3.png", sep = ""), width = 40, height = 10, res = 300,units = "cm")
  upperY = max(c(featCovInterE$rawRepCov5/featCovInterE$baseFreq5, featCovE$rawRepCov3/featCovE$baseFreq3,
                 featCovInterE$rawRepCov3/featCovInterE$baseFreq3, featCovE$rawRepCov5/featCovE$baseFreq5))
  
  layout(matrix(c(1,2,3,4), nrow = 1))
  par(mar=c(5,5,5,0))
  plot(featCovInterE$rawRepCov5/featCovInterE$baseFreq5, type = "l", col = 1, ylim = c(0,upperY), main = paste(TEchoice,"late"))
  smLine <- smooth.spline(featCovInterE$rawRepCov5/featCovInterE$baseFreq5, nknots = 50)
  lines(x = smLine$x, y = smLine$y, col = 3)
  sd = sqrt(featCovInterE$baseFreq5 * GenomeP[TEchoice] * (1-GenomeP[TEchoice])) / featCovInterE$baseFreq5
  lines( GenomeP[TEchoice] - (2*sd), col =2)
  lines( GenomeP[TEchoice] + (2*sd), col =2)
  sd = (sqrt(featCovInterE$baseFreq5 * featCovInterE$p * (1-featCovInterE$p)) / featCovInterE$baseFreq5)
  lines(featCovInterE$p - (2*sd), col =4)
  lines(featCovInterE$p + (2*sd), col =4)
  
  par(mar=c(5,0,5,5))
  plot(featCovE$rawRepCov3/featCovE$baseFreq3, type = "l", col = 1,  ylim = c(0,upperY), main = paste(TEchoice,"early"))
  smLine <- smooth.spline(featCovE$rawRepCov3/featCovE$baseFreq3, nknots = 30)
  lines(x = smLine$x, y = smLine$y, col = 3)
  sd = sqrt(featCovE$baseFreq3 * GenomeP[TEchoice] * (1-GenomeP[TEchoice])) / featCovE$baseFreq3
  lines(GenomeP[TEchoice] - (2*sd), col =2)
  lines(GenomeP[TEchoice] + (2*sd), col =2)
  sd = sqrt(featCovE$baseFreq3 * featCovE$p * (1-featCovE$p)) / featCovE$baseFreq3
  lines(featCovE$p - (2*sd), col =4)
  lines(featCovE$p + (2*sd), col =4)
  
  par(mar=c(5,5,5,0))
  plot(featCovE$rawRepCov5/featCovE$baseFreq5, type = "l", col = 1,  ylim = c(0,upperY),  main = paste(TEchoice,"early"))
  smLine <- smooth.spline(featCovE$rawRepCov5/featCovE$baseFreq5, nknots = 30)
  lines(x = smLine$x, y = smLine$y, col = 3)
  sd = sqrt(featCovE$baseFreq5 * GenomeP[TEchoice] * (1-GenomeP[TEchoice])) / featCovE$baseFreq5
  lines(GenomeP[TEchoice] - (2*sd), col =2)
  lines(GenomeP[TEchoice] + (2*sd), col =2)
  sd = sqrt(featCovE$baseFreq5 * featCovE$p * (1-featCovE$p)) / featCovE$baseFreq5
  lines(featCovE$p - (2*sd), col =4)
  lines(featCovE$p + (2*sd), col =4)
  
  par(mar=c(5,0,5,5))
  plot(featCovInterE$rawRepCov3/featCovInterE$baseFreq3, type = "l", col = 1, ylim = c(0,upperY),  main = paste(TEchoice,"late"))
  smLine <- smooth.spline(featCovInterE$rawRepCov3/featCovInterE$baseFreq3, nknots = 50)
  lines(x = smLine$x, y = smLine$y, col = 3)
  sd = sqrt(featCovInterE$baseFreq3 * GenomeP[TEchoice] * (1-GenomeP[TEchoice])) / featCovInterE$baseFreq3
  lines( GenomeP[TEchoice] - (2*sd), col =2)
  lines( GenomeP[TEchoice] + (2*sd), col =2)
  sd = sqrt(featCovInterE$baseFreq3 * featCovInterE$p * (1-featCovInterE$p)) / featCovInterE$baseFreq3
  lines( featCovInterE$p - (2*sd), col =4)
  lines( featCovInterE$p + (2*sd), col =4)
  
  dev.off()
  
}

# so now we are looking at compartment increase and genome enrichment


# it is kind of interesting that nothing seems to be the same as what we anticipated based on our other work
# are some of our assumptions on that statistical model broken with L1s ?
# It might have something to do with the fact that L1s are bigger and there is less of them 
# the random model assumes that L1 bases are randomly distributed around the genome. 
# In reality L1 bases travel in packets 
# the probability of 1 base being an L1 is somewhat dependant on its neighbor being an L1


bins.early <- binSort(repList = joinRepChromatin,bins = Edomain,TE.names = names(joinRepChromatin), repType = rep("repeats", length(joinRepChromatin)))$counts
bins.late <- binSort(repList = joinRepChromatin,bins = Ldomain,TE.names = names(joinRepChromatin), repType = rep("repeats", length(joinRepChromatin)))$counts
bins.ttr <- binSort(repList = joinRepChromatin,bins = TTR,TE.names = names(joinRepChromatin), repType = rep("repeats", length(joinRepChromatin)))$counts



plot( log10(bins.early$end - bins.early$start), log10(bins.early$E))
points(log10(bins.late$end - bins.late$start), log10(bins.late$E), col = 2)

plot(((bins.early$R)/ (bins.early$end - bins.early$start)), ((bins.early$T) / (bins.early$end - bins.early$start)))
points(((bins.late$R) /(bins.late$end - bins.late$start)), ((bins.late$T)/ (bins.late$end - bins.late$start)), col = 2)
points(((bins.ttr$R) /(bins.ttr$end - bins.ttr$start)), ((bins.ttr$T)/ (bins.ttr$end - bins.ttr$start)), col = 3)





sd = sqrt(output$baseFreq5 * p * (1-p)) / output$baseFreq5

plot(output$rawRepCov5/output$baseFreq5, type = "l")
lines(p - (2*sd), col =2)
lines(p + (2*sd), col =2)


