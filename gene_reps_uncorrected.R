#### lets get the unadjusted regions working to get the figures looking better. 
### then we can finish up the writing
### concentrate on finshing 






setwd("~/Desktop/Domain_manuscript/")

rm(list = ls())


library(GenomicRanges)
library(rtracklayer)
library(zoo)

spec1 <- "Human"
genome = "hg19"


source(file="~/Desktop/Domain_manuscript/Domain_manuscript_scripts/functions.R")
source(file="~/Desktop/Domain_manuscript/Domain_manuscript_scripts/rep_db.R")
load(file = "~/Desktop/Domain_manuscript/R_objects/chromStateCombined")

# so lets read in the whole genome and find a way to build the neighbormatrix

# ref gene
web <- paste("http://hgdownload.soe.ucsc.edu/goldenPath/", genome, "/database/refGene.txt.gz", sep = "")
con <- gzcon(url(web))
txt <- readLines(con)
refgene <- read.delim(textConnection(txt), header = FALSE)

web <- paste("http://hgdownload.soe.ucsc.edu/goldenPath/", genome, "/database/chromInfo.txt.gz", sep = "")
con <- gzcon(url(web))
txt <- readLines(con)
chrom_info <- read.delim(textConnection(txt), header = FALSE)

# pull intergenic regions
refgene.gr <- GRanges(seqnames=Rle(refgene[,3]), ranges = IRanges(start = refgene[,5], end = refgene[,6]))

refgene_gap.gr <- filterIntergenic(refgene)

intronKeep.gr <- filterIntron(refgene)


# seq gaps
ses <- browserSession("UCSC")
genome(ses) <- genome
data.types <- c("gaps")
track.name <- c("gap")
table.name <- c("gap")
for(i in 1:length(data.types)){
  dat <- getTable(
    ucscTableQuery(
      ses, 
      track = track.name[i],
      table = table.name[i]
    )
  )
  assign(data.types[i], dat)        
}
gaps.gr <- GRanges(seqnames = Rle(gaps$chrom),
                   ranges = IRanges(start = gaps$chromStart, end = gaps$chromEnd)
)		 


# removing gaps intergenic
int <- intersect(refgene_gap.gr, gaps.gr)
OL <- as.matrix(findOverlaps(refgene_gap.gr, int))
agg <- aggregate(x=width(int)[OL[,2]], by=list(OL[,1]), FUN = sum)
bins_gene_gap <- data.frame(chr = seqnames(refgene_gap.gr), start = start(refgene_gap.gr), end = end(refgene_gap.gr), Known = width(refgene_gap.gr))
bins_gene_gap$Known[agg[,1]] = bins_gene_gap$Known[agg[,1]] - agg[,2]
bins_gene_gap <- (bins_gene_gap[bins_gene_gap$Known > bins_gene_gap$end - bins_gene_gap$start -101,])
bins_gene_gap <- (bins_gene_gap[(bins_gene_gap$end - bins_gene_gap$start) / bins_gene_gap$Known > .95,])

# remove gaps from intron
gap.int <- intersect(intronKeep.gr, gaps.gr)
gapOl <- as.matrix(findOverlaps(intronKeep.gr, gap.int))
bins_intron = as.data.frame(intronKeep.gr)[,1:4]
colnames(bins_intron)[c(1,4)] <- c("chr", "Known")
bins_intron$Knonw[gapOl[,1]] <- bins_intron$Knonw[gapOl[,1]] - width(gap.int)[gapOl[,2]]

# get repeat info
rep = rep_info(spec1=spec1,genome=genome)
# sort into intergenic bins and intronic bins
intergenic_reps <- binSort(rep=rep, bins=bins_gene_gap, TE.names=names(rep), repType = rep("repeats",length(rep)))
intron_reps <- binSort(rep=rep, bins=bins_intron, TE.names=names(rep), repType = rep("repeats",length(rep)))

intergenic_chromatin <- binSort(rep = `H1-hESC`, bins = bins_gene_gap, TE.names = names(`H1-hESC`), repType = rep("chromatin",length(`H1-hESC`)))
intron_chromatin <- binSort(rep = `H1-hESC`, bins = bins_intron, TE.names = names(`H1-hESC`), repType = rep("chromatin",length(`H1-hESC`)))



intergenicJoinedFam <- data.frame(intergenic_reps$counts[,c("chr", "start", "end","Known")],
                                  old_L1 = rowSums(intergenic_reps$counts[,c("L1ME", "L1MD", "L1MC", "L1MB")]),
                                  new_L1 = rowSums(intergenic_reps$counts[,c("L1MA", "L1PB",  "L1PA","L1HS")]),
                                  Alu = rowSums(intergenic_reps$counts[,c("AluS", "AluY", "AluJ")]),
                                  Ancient = rowSums(intergenic_reps$counts[,c("MIR", "L2")])
)

intronJoinedFam <- data.frame(intron_reps$counts[,c("chr", "start", "end","Known")],
                              old_L1 = rowSums(intron_reps$counts[,c("L1ME", "L1MD", "L1MC","L1MB")]),
                              new_L1 = rowSums(intron_reps$counts[,c("L1MA", "L1PB", "L1PA","L1HS")]),
                              Alu = rowSums(intron_reps$counts[,c("AluS", "AluY", "AluJ")]),
                              Ancient = rowSums(intron_reps$counts[,c("MIR", "L2")])
)




joinRep <- list(old_L1 = rbind(rep$L1ME, rep$L1MD, rep$L1MC, rep$L1MB),
                new_L1 = rbind(rep$L1MA, rep$L1PB  ,rep$L1PA, rep$L1HS),
                Alu = rbind(rep$AluJ, rep$AluS, rep$AluY),
                Ancient = rbind(rep$MIR, rep$L2))

joinRepChromatin <- c(joinRep, `H1-hESC`, HepG2, K562, HUVEC, `HeLa-S3`, GM12878)

#### lets get the average R chromatin profile and measure how that changes according to TE levels

###### intergenic section
#######
######


repTypes <- c(rep("repeats", length(joinRep)), rep("chromatin", length(`H1-hESC`)))

len = 40000

TEs_intergenic_Alu <- covCalcPlot5prime3prime(lenChoice = len,repChoice = "Alu",repBins = intergenicJoinedFam,repList = joinRep,refgene = refgene,type = "intergenic",repType = "repeats")
TEs_intergenic_old_L1 <- covCalcPlot5prime3prime(lenChoice = len,repChoice = "old_L1",repBins = intergenicJoinedFam,repList = joinRep,refgene = refgene,type = "intergenic",repType = "repeats")
TEs_intergenic_new_L1 <- covCalcPlot5prime3prime(lenChoice = len,repChoice = "new_L1",repBins = intergenicJoinedFam,repList = joinRep,refgene = refgene,type = "intergenic",repType = "repeats")
TEs_intergenic_Ancient <- covCalcPlot5prime3prime(lenChoice = len,repChoice = "Ancient",repBins = intergenicJoinedFam,repList = joinRep,refgene = refgene,type = "intergenic",repType = "repeats")

prime5rateR <- NULL
prime3rateR <- NULL
chromType = c("H1-hESC", "HepG2", "K562", "HUVEC", "HeLa-S3", "GM12878")
for(i in 1:length(chromType)){
  chromo <- get(chromType[i])
  intergenic_chromatin <- binSort(rep = chromo, bins = bins_gene_gap, TE.names = names(chromo), repType = rep("chromatin",length(chromo)))
  TEs_intergenic_R <- covCalcPlot5prime3prime(lenChoice = len,repChoice = "R",repBins = intergenic_chromatin$counts,repList = chromo,refgene = refgene,type = "intergenic",repType = "chromatin")
  assign(x = paste("TEs_intergenic_R", chromType[i], sep= "_"), TEs_intergenic_R)
  prime5rateR <- rbind(prime5rateR , TEs_intergenic_R$rawRepCov5/TEs_intergenic_R$baseFreq5prime)
  prime3rateR <- rbind(prime3rateR , TEs_intergenic_R$rawRepCov3/TEs_intergenic_R$baseFreq3prime)
}


### here we specify the colours 


# need to get colours working on 5 prime side

cutLine5 <- data.frame(minBP = rep(NA, 10),maxBP =rep(NA, 10), minBPsqrt = rep(NA, 10), 
                       maxBPsqrt = rep(NA, 10), 
                       freqBPmin = max(TEs_intergenic_Alu$baseFreq5prime)/10 * 0:9, 
                       freqBPmax = max(TEs_intergenic_Alu$baseFreq5prime)/10 * 1:10, 
                       cols = grey.colors(10,start = .9,end = 0))
for(i in 1:10){
  if(cutLine5$freqBPmax[i] >= min(TEs_intergenic_Alu$baseFreq5prime)){
      cutLine5$maxBP[i] = min((((len+1):1))[TEs_intergenic_Alu$baseFreq5prime <= cutLine5$freqBPmax[i]])
      cutLine5$maxBPsqrt[i] = sqrt(cutLine5$maxBP[i])
  }
  
  if(cutLine5$freqBPmax[i] >= min(TEs_intergenic_Alu$baseFreq5prime)){
    cutLine5$minBP[i] = max((((len+1):1))[TEs_intergenic_Alu$baseFreq5prime >= cutLine5$freqBPmin[i]])
    cutLine5$minBPsqrt[i] = sqrt(cutLine5$minBP[i])
  }
}
cutLine5$cols <- as.character(cutLine5$cols)


cutLine3 <- data.frame(minBP = rep(NA, 10),maxBP =rep(NA, 10), minBPsqrt = rep(NA, 10), 
                       maxBPsqrt = rep(NA, 10), 
                       freqBPmin = max(TEs_intergenic_Alu$baseFreq3prime)/10 * 0:9, 
                       freqBPmax = max(TEs_intergenic_Alu$baseFreq3prime)/10 * 1:10, 
                       cols = grey.colors(10,start = .9,end = 0))
for(i in 1:10){
  if(cutLine3$freqBPmax[i] >= min(TEs_intergenic_Alu$baseFreq3prime)){
    cutLine3$minBP[i] = min((1:(len+1))[TEs_intergenic_Alu$baseFreq3prime <= cutLine3$freqBPmax[i]])
    cutLine3$minBPsqrt[i] = sqrt(cutLine3$minBP[i])
  }
  
  if(cutLine3$freqBPmax[i] >= min(TEs_intergenic_Alu$baseFreq3prime)){
    cutLine3$maxBP[i] = max((1:(len+1))[TEs_intergenic_Alu$baseFreq3prime >= cutLine3$freqBPmin[i]])
    cutLine3$maxBPsqrt[i] = sqrt(cutLine3$maxBP[i])
  }
}
cutLine3$cols <- as.character(cutLine3$cols)



#### begin plotting 

pdf(file = "plots/geneRep/intergenic/chromo2.pdf", onefile = TRUE,width =  5,height = 3)
layout(matrix(c(1,2), nrow = 1))
par(mar = c(3,4,5,1))
het.smooth <- smooth.spline(sqrt((len+1):1)*-1, prime5rateR[1,])
plot(het.smooth$x,het.smooth$y, type = "n", ylim = c(0,1), ylab = "represive chromatin", xaxt = "n", xlab = "")
for(i in 1:nrow(prime5rateR)){
  het.smooth <- smooth.spline(sqrt((len+1):1)*-1, prime5rateR[i,])
  lines(het.smooth$x, het.smooth$y, col = 8, lwd = 1)
}
het.smooth <- smooth.spline(sqrt((len+1):1)*-1, colMeans(prime5rateR))
lines(het.smooth$x, het.smooth$y, col = 1, lwd = 3)
axis(side = 1,at = c(0, -50, -100, -150, -200), labels = c("", "", "", "", ""))


par(mar = c(3,1,5,4))
het.smooth <- smooth.spline(sqrt(1:(len+1)), prime3rateR[1,])
plot(het.smooth$x, het.smooth$y, type = "n", ylim = c(0,1), ylab = "", yaxt = "n", xaxt = "n", xlab = "")
for(i in 1:nrow(prime3rateR)){
  het.smooth <- smooth.spline(sqrt(1:(len+1)), prime3rateR[i,])
  lines(het.smooth$x, het.smooth$y, col = 8, lwd = 1)
}
het.smooth <- smooth.spline(sqrt(1:(len+1)), colMeans(prime3rateR))
lines(het.smooth$x, het.smooth$y, col = 1, lwd = 3)
axis(side = 1,at = c(0, 50, 100, 150, 200), labels = c("", "", "", "", ""))

dev.off()


for(i in 1:length(names(joinRep))){
  
  pdf(file = paste("plots/geneRep/intergenic/TE_2", i, sep  = "_"),width =  5,height = 3)
  layout(matrix(c(1,2), nrow = 1))
  par(mar = c(3,4,5,1))
  sdTE <- get(paste("TEs_intergenic_", names(joinRep)[i], sep = ""))
  sd <- (sqrt(sdTE$baseFreq5prime * sdTE$p * (1 - sdTE$p))/sdTE$baseFreq5prime)
  te.smooth5 <- smooth.spline(sqrt((len+1):1)*-1, sdTE$rawRepCov5/sdTE$baseFreq5prime)
  plot(te.smooth5$x, te.smooth5$y, col = 1,  type = "n", ylab = names(joinRep)[i],
       lwd = 3, xaxt = "n", xlab = "")
  for(c in 1:nrow(cutLine5)){
    if(!is.na(cutLine5$minBP[c])){
      lines(te.smooth5$x[te.smooth5$x >= -cutLine5$minBPsqrt[c] & te.smooth5$x <= -cutLine5$maxBPsqrt[c]],
            te.smooth5$y[te.smooth5$x >= -cutLine5$minBPsqrt[c] & te.smooth5$x <= -cutLine5$maxBPsqrt[c]],
            col = cutLine5$cols[c], lwd = 3)
    }
  }
  axis(side = 1,at = c(0, -50, -100, -150, -200), labels = c("", "", "", "", ""))
  lines(sqrt((len+1):1)*-1, sdTE$p + 2*sd, col = 2)
  lines(sqrt((len+1):1)*-1, sdTE$p - 2*sd, col = 2)
  
  par(mar = c(3,1,5,4))
  sdTE <- get(paste("TEs_intergenic_", names(joinRep)[i], sep = ""))
  sd <- (sqrt(sdTE$baseFreq3prime * sdTE$p * (1 - sdTE$p))/sdTE$baseFreq3prime)
  te.smooth3 <- smooth.spline(sqrt(1:(len+1)), sdTE$rawRepCov3/sdTE$baseFreq3prime)
  plot(te.smooth5$x*-1, te.smooth5$y, col = 1, type = "n",
       yaxt = "n", lwd = 3, xaxt = "n", xlab = "")
  for(c in 1:nrow(cutLine3)){
    if(!is.na(cutLine3$minBP[c])){
      lines(te.smooth3$x[te.smooth3$x >= cutLine3$minBPsqrt[c] & te.smooth3$x <= cutLine3$maxBPsqrt[c]],
            te.smooth3$y[te.smooth3$x >= cutLine3$minBPsqrt[c] & te.smooth3$x <= cutLine3$maxBPsqrt[c]],
            col = cutLine5$cols[c], lwd = 3)
    }
  }
  
  
  axis(side = 1,at = c(0, 50, 100, 150, 200), labels = c("", "", "", "", ""))
  lines(sqrt(1:(len+1)), sdTE$p + 2*sd, col = 2)
  lines(sqrt(1:(len+1)), sdTE$p - 2*sd, col = 2)
  
  dev.off()
}



pdf(file = "plots/geneRep/intergenic/TEextra2", onefile = T,width =  5,height = 3)
layout(matrix(c(1,2), nrow = 1))
par(mar = c(3,4,5,1))
te.smooth <- smooth.spline(sqrt((len+1):1)*-1, sdTE$rawRepCov5/sdTE$baseFreq5prime)
plot(te.smooth$x, te.smooth$y, col = 1, ylim = c(0,.2), type = "n", ylab = "",
     lwd = 3, xaxt = "n", xlab = "", yaxt = "n")
axis(side = 1,at = c(0, -50, -100, -150, -200), labels = (c(0, -50, -100, -150, -200)^2)/1000 )

par(mar = c(3,1,5,4))
te.smooth <- smooth.spline(sqrt(1:(len+1)), sdTE$rawRepCov3/sdTE$baseFreq3prime)
plot(te.smooth$x, te.smooth$y, col = 1, ylim = c(0,.2), type = "n",
     yaxt = "n", lwd = 3, xaxt = "n", xlab = "")
axis(side = 1,at = c(0, 50, 100, 150, 200), labels = (c(0, -50, -100, -150, -200)^2)/1000 )

par(mar = c(6,5,6,5))
image(matrix(c(0,.25,.5,.75,1)), col = grey.colors(n=5,1,0), xaxt = "n", yaxt = "n", main = "position frequency (kb)")
axis(side = 1,at = c(0,.5,1), labels = c(0, 10, 20))
dev.off()

#### before we do introns we should make a legend for bp frequency




###### lets have a go at intronic and get that up to scratch 








repTypes <- c(rep("repeats", length(joinRep)), rep("chromatin", length(`H1-hESC`)))

len = 20000

TEs_intron_Alu <- covCalcPlot5prime3prime(lenChoice = len,repChoice = "Alu",repBins = intronJoinedFam,repList = joinRep,refgene = refgene,type = "intron",repType = "repeats")
TEs_intron_old_L1 <- covCalcPlot5prime3prime(lenChoice = len,repChoice = "old_L1",repBins = intronJoinedFam,repList = joinRep,refgene = refgene,type = "intron",repType = "repeats")
TEs_intron_new_L1 <- covCalcPlot5prime3prime(lenChoice = len,repChoice = "new_L1",repBins = intronJoinedFam,repList = joinRep,refgene = refgene,type = "intron",repType = "repeats")
TEs_intron_Ancient <- covCalcPlot5prime3prime(lenChoice = len,repChoice = "Ancient",repBins = intronJoinedFam,repList = joinRep,refgene = refgene,type = "intron",repType = "repeats")

prime5rateR <- NULL
prime3rateR <- NULL
chromType = c("H1-hESC", "HepG2", "K562", "HUVEC", "HeLa-S3", "GM12878")
for(i in 1:length(chromType)){
  chromo <- get(chromType[i])
  intron_chromatin <- binSort(rep = chromo, bins = bins_intron, TE.names = names(chromo), repType = rep("chromatin",length(chromo)))
  TEs_intron_R <- covCalcPlot5prime3prime(lenChoice = len,repChoice = "R",repBins = intron_chromatin$counts,repList = chromo,refgene = refgene,type = "intron",repType = "chromatin")
  assign(x = paste("TEs_intron_R", chromType[i], sep= "_"), TEs_intron_R)
  prime5rateR <- rbind(prime5rateR , TEs_intron_R$rawRepCov5/TEs_intron_R$baseFreq5prime)
  prime3rateR <- rbind(prime3rateR , TEs_intron_R$rawRepCov3/TEs_intron_R$baseFreq3prime)
}


### here we specify the colours 


# need to get colours working on 5 prime side

cutLine5 <- data.frame(minBP = rep(NA, 10),maxBP =rep(NA, 10), minBPsqrt = rep(NA, 10), 
                       maxBPsqrt = rep(NA, 10), 
                       freqBPmin = max(TEs_intron_Alu$baseFreq5prime)/10 * 0:9, 
                       freqBPmax = max(TEs_intron_Alu$baseFreq5prime)/10 * 1:10, 
                       cols = grey.colors(10,start = .9,end = 0))
for(i in 1:10){
  if(cutLine5$freqBPmax[i] >= min(TEs_intron_Alu$baseFreq5prime)){
    cutLine5$maxBP[i] = min((((len+1):1))[TEs_intron_Alu$baseFreq5prime <= cutLine5$freqBPmax[i]])
    cutLine5$maxBPsqrt[i] = sqrt(cutLine5$maxBP[i])
  }
  
  if(cutLine5$freqBPmax[i] >= min(TEs_intron_Alu$baseFreq5prime)){
    cutLine5$minBP[i] = max((((len+1):1))[TEs_intron_Alu$baseFreq5prime >= cutLine5$freqBPmin[i]])
    cutLine5$minBPsqrt[i] = sqrt(cutLine5$minBP[i])
  }
}
cutLine5$cols <- as.character(cutLine5$cols)


cutLine3 <- data.frame(minBP = rep(NA, 10),maxBP =rep(NA, 10), minBPsqrt = rep(NA, 10), 
                       maxBPsqrt = rep(NA, 10), 
                       freqBPmin = max(TEs_intron_Alu$baseFreq3prime)/10 * 0:9, 
                       freqBPmax = max(TEs_intron_Alu$baseFreq3prime)/10 * 1:10, 
                       cols = grey.colors(10,start = .9,end = 0))
for(i in 1:10){
  if(cutLine3$freqBPmax[i] >= min(TEs_intron_Alu$baseFreq3prime)){
    cutLine3$minBP[i] = min((1:(len+1))[TEs_intron_Alu$baseFreq3prime <= cutLine3$freqBPmax[i]])
    cutLine3$minBPsqrt[i] = sqrt(cutLine3$minBP[i])
  }
  
  if(cutLine3$freqBPmax[i] >= min(TEs_intron_Alu$baseFreq3prime)){
    cutLine3$maxBP[i] = max((1:(len+1))[TEs_intron_Alu$baseFreq3prime >= cutLine3$freqBPmin[i]])
    cutLine3$maxBPsqrt[i] = sqrt(cutLine3$maxBP[i])
  }
}
cutLine3$cols <- as.character(cutLine3$cols)



#### begin plotting 

pdf(file = "plots/geneRep/intton/chromo.pdf", onefile = TRUE, width = 5 ,height = 3)
layout(matrix(c(1,2), nrow = 1))
par(mar = c(3,4,5,1))
het.smooth <- smooth.spline(sqrt((len+1):1)*-1, prime5rateR[1,])
plot(het.smooth$x,het.smooth$y, type = "n", ylim = c(0,1), ylab = "represive chromatin", xlab = "", xaxt = "n")
for(i in 1:nrow(prime5rateR)){
  het.smooth <- smooth.spline(sqrt((len+1):1)*-1, prime5rateR[i,])
  lines(het.smooth$x, het.smooth$y, col = 8, lwd = 1)
}
het.smooth <- smooth.spline(sqrt((len+1):1)*-1, colMeans(prime5rateR))
lines(het.smooth$x, het.smooth$y, col = 1, lwd = 3)
axis(side = 1,at = c(0, -35, -70, -105, -140), labels = c("", "", "", "", ""))


par(mar = c(3,1,5,4))
het.smooth <- smooth.spline(sqrt(1:(len+1)), prime3rateR[1,])
plot(het.smooth$x, het.smooth$y, type = "n", ylim = c(0,1), ylab = "", yaxt = "n", xaxt = "n", xlab = "")
for(i in 1:nrow(prime3rateR)){
  het.smooth <- smooth.spline(sqrt(1:(len+1)), prime3rateR[i,])
  lines(het.smooth$x, het.smooth$y, col = 8, lwd = 1)
}
het.smooth <- smooth.spline(sqrt(1:(len+1)), colMeans(prime3rateR))
lines(het.smooth$x, het.smooth$y, col = 1, lwd = 3)
axis(side = 1,at = c(0, 35, 70, 105, 140), labels = c("", "", "", "", ""))
dev.off()


for(i in 1:length(names(joinRep))){
  pdf(file = paste("plots/geneRep/intton/TE",i, ".pdf",sep = ""), onefile = TRUE, width = 5 ,height = 3)
  layout(matrix(c(1,2), nrow = 1))
  par(mar = c(3,4,5,1))
  sdTE <- get(paste("TEs_intron_", names(joinRep)[i], sep = ""))
  sd <- (sqrt(sdTE$baseFreq5prime * sdTE$p * (1 - sdTE$p))/sdTE$baseFreq5prime)
  te.smooth5 <- smooth.spline(sqrt((len+1):1)*-1, sdTE$rawRepCov5/sdTE$baseFreq5prime)
  plot(te.smooth5$x, te.smooth5$y, col = 1, type = "n", ylab = names(joinRep)[i],
       lwd = 3, xaxt = "n", xlab = "")
  for(c in 1:nrow(cutLine5)){
    if(!is.na(cutLine5$minBP[c])){
      lines(te.smooth5$x[te.smooth5$x >= -cutLine5$minBPsqrt[c] & te.smooth5$x <= -cutLine5$maxBPsqrt[c]],
            te.smooth5$y[te.smooth5$x >= -cutLine5$minBPsqrt[c] & te.smooth5$x <= -cutLine5$maxBPsqrt[c]],
            col = cutLine5$cols[c], lwd = 3)
    }
  }
  axis(side = 1,at = c(0, -35, -70, -105, -140), labels = c("", "", "", "", ""))
  lines(sqrt((len+1):1)*-1, sdTE$p + 2*sd, col = 2)
  lines(sqrt((len+1):1)*-1, sdTE$p - 2*sd, col = 2)
  
  par(mar = c(3,1,5,4))
  sdTE <- get(paste("TEs_intron_", names(joinRep)[i], sep = ""))
  sd <- (sqrt(sdTE$baseFreq3prime * sdTE$p * (1 - sdTE$p))/sdTE$baseFreq3prime)
  te.smooth3 <- smooth.spline(sqrt(1:(len+1)), sdTE$rawRepCov3/sdTE$baseFreq3prime)
  plot(te.smooth5$x*-1, te.smooth5$y, col = 1, type = "n",
       yaxt = "n", lwd = 3, xaxt = "n", xlab = "")
  for(c in 1:nrow(cutLine3)){
    if(!is.na(cutLine3$minBP[c])){
      lines(te.smooth3$x[te.smooth3$x >= cutLine3$minBPsqrt[c] & te.smooth3$x <= cutLine3$maxBPsqrt[c]],
            te.smooth3$y[te.smooth3$x >= cutLine3$minBPsqrt[c] & te.smooth3$x <= cutLine3$maxBPsqrt[c]],
            col = cutLine5$cols[c], lwd = 3)
    }
  }
  
  
  axis(side = 1,at = c(0, 35, 70, 105, 140), labels = c("", "", "", "", ""))
  lines(sqrt(1:(len+1)), sdTE$p + 2*sd, col = 2)
  lines(sqrt(1:(len+1)), sdTE$p - 2*sd, col = 2)
  dev.off()
}

pdf(file = "plots/geneRep/intton/extra.pdf", onefile = TRUE, width = 5 ,height = 3)
layout(matrix(c(1,2), nrow = 1))

par(mar = c(3,4,5,1))
te.smooth <- smooth.spline(sqrt((len+1):1)*-1, sdTE$rawRepCov5/sdTE$baseFreq5prime)
plot(te.smooth$x, te.smooth$y, col = 1, ylim = c(0,.2), type = "n", ylab = "",
     lwd = 3, xaxt = "n", xlab = "", yaxt = "n")
axis(side = 1,at = c(0, -35, -70, -105, -140), labels = (c(0, -35, -70, -105, -140)^2)/1000 )

par(mar = c(3,1,5,4))
te.smooth <- smooth.spline(sqrt(1:(len+1)), sdTE$rawRepCov3/sdTE$baseFreq3prime)
plot(te.smooth$x, te.smooth$y, col = 1, ylim = c(0,.2), type = "n",
     yaxt = "n", lwd = 3, xaxt = "n", xlab = "")
axis(side = 1,at = c(0, 35, 70, 105, 140), labels = (c(0, 35, 70, 105, 140)^2)/1000 )

par(mar = c(6,5,6,5))
image(matrix(c(0,.25,.5,.75,1)), col = grey.colors(n=5,1,0), xaxt = "n", yaxt = "n", main = "position frequency (kb)")
axis(side = 1,at = c(0,.5,1), labels = c(0, 70, 140))
dev.off()


















# maybe if we get all the R from each cell line and get an average over the region 

# how could we find these insertion channels 

# if we look at the mean age profile of elements do we see an increse in % identity 







