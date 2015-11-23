# now we can begin to combine the anlysis of both kinds of regions
# this will be our fig 2

# where we look at the relationship between TE accumulation and genes

# the idea is that TEs accumulate according to our model around genes

# on both fractions theres constraint pertaining to the amount but the constraint on position is different





setwd("~/Desktop/topological_domains/")

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
intergenic_reps <- binSort(rep=rep, bins=bins_gene_gap, TE.names=names(rep))

intron_reps <- binSort(rep=rep, bins=bins_intron, TE.names=names(rep))

intergenic_chromatin <- binSort(rep = `H1-hESC`, bins = bins_gene_gap, TE.names = names(`H1-hESC`), repType = "chromatin")

intron_chromatin <- binSort(rep = `H1-hESC`, bins = bins_intron, TE.names = names(`H1-hESC`), repType = "chromatin")

#interesting that intervals with repeats in them tend to be over 1000 in length

intergenicSample <- sample(x = 1:nrow(intergenic_reps$counts), size = 10000, replace = FALSE)
intronSample <- sample(x = 1:nrow(intron_reps$counts), size = 10000, replace = FALSE)


TEs_intergenic <- covCalcPlot5prime3primeGC(lenChoice=100000,repChoice="L1PA",repBins=intergenic_reps$counts[intergenicSample,],refgene=refgene,type="intergenic",genome=genome,repList=rep)

TEs_intron <- covCalcPlot5prime3primeGC(lenChoice=100000,repChoice="L1PA",repBins=intron_reps$counts[intronSample,],refgene=refgene,type="intron",genome=genome,repList=rep)

Intergenic_break <- LocalGClevel(maxLen = 20000,repBins = intergenic_reps$counts[intergenicSample,], repChoice = "L1PA", genome = genome,repList = rep)
Intron_break <- LocalGClevel(maxLen = 20000,repBins = intron_reps$counts[intronSample,], repChoice = "L1PA", genome = genome,repList = rep)



# 
# Intron_breakSamp <- Intron_break[sample(1:nrow(Intron_break),size = 1000, replace=FALSE),]
# plot((Intron_breakSamp$GC), Intron_breakSamp$cov/(Intron_breakSamp$end - Intron_breakSamp$start + 1), pch = 16, cex =.2, type = "n")
# 
# 
# for(i in 1:10){
# Intron_breakSamp <- Intron_break[Intron_break$width>10000,][sample(1:nrow(Intron_break[Intron_break$width>10000,]),size = 1000, replace=FALSE),]
# points(Intron_breakSamp$GC, Intron_breakSamp$cov/Intron_breakSamp$width, pch = 16, cex =.1)
# IntronModLO <- loess((Intron_breakSamp$cov/Intron_breakSamp$width) ~ Intron_breakSamp$GC)
# lines(seq(0,1,.01), predict(IntronModLO, newdata = seq(0,1,.01)), lty = 2, lwd = 2, col =2)
# 
# Intergenic_breakSamp <- Intergenic_break[Intergenic_break$width>10000,][sample(1:nrow(Intergenic_break[Intergenic_break$width>10000,]),size = 1000, replace=FALSE),]
# points(Intergenic_breakSamp$GC, Intergenic_breakSamp$cov/Intergenic_breakSamp$width, pch = 16, cex =.1, col = 3)
# IntergenicModLO <- loess((Intergenic_breakSamp$cov/Intergenic_breakSamp$width) ~ Intergenic_breakSamp$GC)
# lines(seq(0,1,.01), predict(IntergenicModLO, newdata = seq(0,1,.01)), lty = 2, lwd = 2, col =4)
# 
# }


###### 
####
##
#
#     Get the sampel Regions 
#   get the TE groups
#  run the groups through

###   use our old GC calculations for normalization
#



intergenicJoinedFam <- data.frame(intergenic_reps$counts[,c("chr", "start", "end","Known")],
                                  old_L1 = rowSums(intergenic_reps$counts[,c("L1ME", "L1MD", "L1MC","L1MB")]),
                                  new_L1 = rowSums(intergenic_reps$counts[,c("L1MA", "L1PB", "L1PA","L1HS")]),
                                  Alu = rowSums(intergenic_reps$counts[,c("AluS", "AluY", "AluJ")]),
                                  Ancient = rowSums(intergenic_reps$counts[,c("MIR", "L2")])
)

intronJoinedFam <- data.frame(intron_reps$counts[,c("chr", "start", "end","Known")],
                              old_L1 = rowSums(intron_reps$counts[,c("L1ME", "L1MD", "L1MC","L1MB")]),
                              new_L1 = rowSums(intron_reps$counts[,c("L1MA", "L1PB", "L1PA","L1HS")]),
                              Alu = rowSums(intron_reps$counts[,c("AluS", "AluY", "AluJ")]),
                              Ancient = rowSums(intron_reps$counts[,c("MIR", "L2")])
)


intergenicJoinedFam <- intergenicJoinedFam[intergenicSample,]
intronJoinedFam <- intronJoinedFam[intronSample,]


joinRep <- list(old_L1 = rbind(rep$L1ME, rep$L1MD, rep$L1MC, rep$L1MB),
                new_L1 = rbind(rep$L1MA, rep$L1PB, rep$L1PA, rep$L1HS),
                Alu = rbind(rep$AluJ, rep$AluS, rep$AluY),
                Ancient = rbind(rep$MIR, rep$L2))



joinSampGenomeBig <- localGCgenome(repList = joinRep,binSize = 20000,sampSize = 5000, genome = genome)
joinSampGenomeSmall <- localGCgenome(repList = joinRep,binSize = 1000,sampSize = 5000, genome = genome)



layout(matrix(c(1,2,3,4), nrow = 2, ncol = 2))
par(mar=c(4,4,4,4))

plot(joinSampGenomeBig$GC, joinSampGenomeBig$old_L1, main = "old_L1", pch = 16, cex = .2, xlab = "GC content", ylab = "TE fraction", ylim = c(0,1))
old_L1Mod <- loess(joinSampGenomeBig$old_L1 ~ joinSampGenomeBig$GC)
lines(seq(0,1,0.01), predict(old_L1Mod, newdata = seq(0,1,0.01)), col = 2, lwd = 2)

plot(joinSampGenomeBig$GC, joinSampGenomeBig$new_L1, main = "new_L1", pch = 16, cex = .2, xlab = "GC content", ylab = "TE fraction", ylim = c(0,1))
new_L1Mod <- loess(joinSampGenomeBig$new_L1 ~ joinSampGenomeBig$GC)
lines(seq(0,1,0.01), predict(new_L1Mod, newdata = seq(0,1,0.01)), col = 2, lwd = 2)

plot(joinSampGenomeSmall$GC, joinSampGenomeSmall$Alu, main = "Alu", pch = 16, cex = .2, xlab = "GC content", ylab = "TE fraction", ylim = c(0,1))
AluMod <- loess(joinSampGenomeSmall$Alu ~ joinSampGenomeSmall$GC)
lines(seq(0,1,0.01), predict(AluMod, newdata = seq(0,1,0.01)), col = 2, lwd = 2)

plot(joinSampGenomeSmall$GC, joinSampGenomeSmall$Ancient, main = "Ancient", pch = 16, cex = .2, xlab = "GC content", ylab = "TE fraction", ylim = c(0,1))
AncientMod <- loess(joinSampGenomeSmall$Ancient ~ joinSampGenomeSmall$GC)
lines(seq(0,1,0.01), predict(AncientMod, newdata = seq(0,1,0.01)), col = 2, lwd = 2)




minS <- NULL
maxS <- NULL
minRepCov <- NULL
maxRepCov <- NULL
lenChoice = 100000
genome_type <- "intergenic"
binnedGenome <- get(paste(genome_type , "JoinedFam", sep = ""))



repChoice <- "old_L1"
old_L1 <- covCalcPlot5prime3prime(lenChoice=lenChoice,repChoice=repChoice, repBins=binnedGenome,repList=joinRep, refgene=refgene, type= genome_type,minBinSize=minS, repType = "repeats",maxBinSize=maxS,minRepCov=NULL, maxRepCov=NULL)
repChoice <- "new_L1"
new_L1 <- covCalcPlot5prime3prime(lenChoice=lenChoice,repChoice=repChoice, repBins=binnedGenome,repList=joinRep, refgene=refgene, type= genome_type,minBinSize=minS, repType = "repeats",maxBinSize=maxS,minRepCov=NULL, maxRepCov=NULL)
repChoice <- "Alu"
Alu <- covCalcPlot5prime3prime(lenChoice=lenChoice,repChoice=repChoice, repBins=binnedGenome,repList=joinRep, refgene=refgene, type= genome_type,minBinSize=minS, repType = "repeats",maxBinSize=maxS,minRepCov=NULL, maxRepCov=NULL)
repChoice <- "Ancient"
Ancient <- covCalcPlot5prime3prime(lenChoice=lenChoice,repChoice=repChoice, repBins=binnedGenome,repList=joinRep, refgene=refgene, type= genome_type,minBinSize=minS, repType = "repeats",maxBinSize=maxS,minRepCov=NULL, maxRepCov=NULL)

sqrtCoordinates <- -(sqrt(1:(lenChoice + 1))[(lenChoice + 1):1])
Coordinates <- -((lenChoice + 1))


TEset <- get(paste("TEs_", genome_type,sep=""))

k = 1000
GCsum3 <- rollsum(TEset$prime3gc,k = k, fill="extend", align = "center")
baseSum3 <- rollsum(TEset$baseFreq3prime,k = k, fill="extend", align = "center")
gcRate3 <- GCsum3/baseSum3
n3 = TEs_intergenic$baseFreq3prime

GCsum5 <- rollsum(TEset$prime5gc,k = k, fill="extend", align = "center")
baseSum5 <- rollsum(TEset$baseFreq5prime,k = k, fill="extend", align = "center")
gcRate5 <- GCsum5/baseSum5
n5 = TEs_intergenic$baseFreq5prime

ylims <- c(0,.2)



greys <- grey(level = seq(.9,0,-.1))
prime5BaseFreq <- data.frame(baseFreq =TEset$baseFreq5prime, 
                             group = cut(TEset$baseFreq5prime, breaks = length(greys), labels = as.character(1:length(greys))),
                             position = sqrtCoordinates)
prime3BaseFreq <- data.frame(baseFreq =TEset$baseFreq3prime, 
                             group = cut(TEset$baseFreq3prime, breaks = length(greys), labels = as.character(1:length(greys))),
                             position = sqrtCoordinates)

layout(matrix(c(1,3,5,7,2,4,6,8), nrow = 4, ncol = 2), heights = c(10,10,10,15,10,10,10,15))

par(mar=c(0,5,1,0))
plot(sqrtCoordinates,(old_L1$rawRepCov5/TEset$baseFreq5prime), type = "n", ylab = "old_L1", xlab = "", ylim = ylims, xaxt = "n" )
for(c in 1:length(greys)){
  lines(sqrtCoordinates[prime5BaseFreq$group == c], (old_L1$rawRepCov5/TEset$baseFreq5prime)[prime5BaseFreq$group == c], col = greys[c])
}
lines(sqrtCoordinates,predict(old_L1Mod, gcRate5) + (2 * sqrt(n5 * predict(old_L1Mod, gcRate5) *  (1 - predict(old_L1Mod, gcRate5))))/n5, col = 2)
lines(sqrtCoordinates,predict(old_L1Mod, gcRate5) - (2 * sqrt(n5 * predict(old_L1Mod, gcRate5) *  (1 - predict(old_L1Mod, gcRate5))))/n5, col = 2)
grid()

par(mar=c(0,0,1,5))
plot(sqrt(1:100001),(old_L1$rawRepCov3/TEset$baseFreq3prime), type = "l", ylab = "old_L1", xlab = "" , ylim = ylims, xaxt = "n", yaxt = "n")
 lines(sqrt(1:100001),predict(old_L1Mod, gcRate3) + (2 * sqrt(n3 * predict(old_L1Mod, gcRate3) *  (1 - predict(old_L1Mod, gcRate3))))/n3, col = 2)
 lines(sqrt(1:100001),predict(old_L1Mod, gcRate3) - (2 * sqrt(n3 * predict(old_L1Mod, gcRate3) *  (1 - predict(old_L1Mod, gcRate3))))/n3, col = 2)
grid()

par(mar=c(0,5,1,0))
plot(sqrtCoordinates,(new_L1$rawRepCov5/TEset$baseFreq5prime), type = "l", ylab = "new_L1", xlab = "", ylim = ylims, xaxt = "n"  )
 lines(sqrtCoordinates,predict(new_L1Mod, gcRate5) + (2 * sqrt(n5 * predict(new_L1Mod, gcRate5) *  (1 - predict(new_L1Mod, gcRate5))))/n5, col = 2)
 lines(sqrtCoordinates,predict(new_L1Mod, gcRate5) - (2 * sqrt(n5 * predict(new_L1Mod, gcRate5) *  (1 - predict(new_L1Mod, gcRate5))))/n5, col = 2)
grid()

par(mar=c(0,0,1,5))
plot(sqrt(1:100001),(new_L1$rawRepCov3/TEset$baseFreq3prime), type = "l", ylab = "new_L1", xlab = "" , ylim = ylims, xaxt = "n", yaxt = "n")
 lines(sqrt(1:100001),predict(new_L1Mod, gcRate3) + (2 * sqrt(n3 * predict(new_L1Mod, gcRate3) *  (1 - predict(new_L1Mod, gcRate3))))/n3, col = 2)
 lines(sqrt(1:100001),predict(new_L1Mod, gcRate3) - (2 * sqrt(n3 * predict(new_L1Mod, gcRate3) *  (1 - predict(new_L1Mod, gcRate3))))/n3, col = 2)
grid()

par(mar=c(0,5,1,0))
plot(sqrtCoordinates,(Alu$rawRepCov5/TEset$baseFreq5prime), type = "l", ylab = "Alu", xlab = "" , ylim = ylims, xaxt = "n" )
 lines(sqrtCoordinates,predict(AluMod, gcRate5) + (2 * sqrt(n5 * predict(AluMod, gcRate5) *  (1 - predict(AluMod, gcRate5))))/n5, col = 2)
 lines(sqrtCoordinates,predict(AluMod, gcRate5) - (2 * sqrt(n5 * predict(AluMod, gcRate5) *  (1 - predict(AluMod, gcRate5))))/n5, col = 2)
grid()

par(mar=c(0,0,1,5))
plot(sqrt(1:100001),(Alu$rawRepCov3/TEset$baseFreq3prime), type = "l", ylab = "Alu", xlab = "" , ylim = ylims, xaxt = "n", yaxt = "n")
 lines(sqrt(1:100001),predict(AluMod, gcRate3) + (2 * sqrt(n3 * predict(AluMod, gcRate3) *  (1 - predict(AluMod, gcRate3))))/n3, col = 2)
 lines(sqrt(1:100001),predict(AluMod, gcRate3) - (2 * sqrt(n3 * predict(AluMod, gcRate3) *  (1 - predict(AluMod, gcRate3))))/n3, col = 2)
grid()

par(mar=c(6,5,1,0))
plot(sqrtCoordinates,(Ancient$rawRepCov5/TEset$baseFreq5prime), type = "l", ylab = "Ancient", xlab = "" , ylim = ylims, xaxt = "n" )
 lines(sqrtCoordinates,predict(AncientMod, gcRate5) + (2 * sqrt(n5 * predict(AncientMod, gcRate5) *  (1 - predict(AncientMod, gcRate5))))/n5, col = 2)
 lines(sqrtCoordinates,predict(AncientMod, gcRate5) - (2 * sqrt(n5 * predict(AncientMod, gcRate5) *  (1 - predict(AncientMod, gcRate5))))/n5, col = 2)
grid()

axis(1,at = seq(0,-350, by = -50), labels = seq(0,350, by = 50)^2)


par(mar=c(6,0,1,5))
plot(sqrt(1:100001),(Ancient$rawRepCov3/TEset$baseFreq3prime), type = "l", ylab = "Ancient", xlab = "", ylim = ylims, xaxt = "n", yaxt = "n" )
 lines(sqrt(1:100001),predict(AncientMod, gcRate3) + (2 * sqrt(n3 * predict(AncientMod, gcRate3) *  (1 - predict(AncientMod, gcRate3))))/n3, col = 2)
 lines(sqrt(1:100001),predict(AncientMod, gcRate3) - (2 * sqrt(n3 * predict(AncientMod, gcRate3) *  (1 - predict(AncientMod, gcRate3))))/n3, col = 2)
grid()

axis(1,at = seq(0,350, by = 50), labels = seq(0,350, by = 50)^2)

### are we controlling for smaller elements ?
### it may be better to use a 1 kb SINE estimate



##### lets get the chormatin states in here

###### its probably worth thinking about GC content
###### so the ideas is that insertions are a result of chromatin at locations near a gene 

heteroBig <- localGCgenome(repList = `H1-hESC`,binSize = 20000,sampSize = 5000, genome = genome, repType = "chromatin")

layout(1)
par(mar=c(5,5,5,5))
plot(heteroBig$GC, heteroBig$T, pch = 16, cex = .3)
heterolo <- loess(heteroBig$T~heteroBig$GC)
lines(seq(.1,1,.01), predict(heterolo, seq(.1,1,.01)))


Rlo  <- loess(heteroBig$R~heteroBig$GC)
Elo  <- loess(heteroBig$E~heteroBig$GC)
PFlo <- loess(heteroBig$PF~heteroBig$GC)
CTCFlo <- loess(heteroBig$CTCF~heteroBig$GC)
Ttlo <- loess(heteroBig$T~heteroBig$GC)
TSSlo <- loess(heteroBig$TSS~heteroBig$GC)
WElo <- loess(heteroBig$WE~heteroBig$GC)



# we need to get bins that have chromatin status included 


colnames(intergenic_chromatin$counts)[5:ncol(intergenic_chromatin$counts)]

R <- covCalcPlot5prime3prime(lenChoice=lenChoice,repChoice="R", repBins=intergenic_chromatin$counts,repList=`H1-hESC`, refgene=refgene, type= genome_type, repType = "chromatin")
E <- covCalcPlot5prime3prime(lenChoice=lenChoice,repChoice="E", repBins=intergenic_chromatin$counts,repList=`H1-hESC`, refgene=refgene, type= genome_type, repType = "chromatin")
PF <- covCalcPlot5prime3prime(lenChoice=lenChoice,repChoice="PF", repBins=intergenic_chromatin$counts,repList=`H1-hESC`, refgene=refgene, type= genome_type, repType = "chromatin")
CTCF <- covCalcPlot5prime3prime(lenChoice=lenChoice,repChoice="CTCF", repBins=intergenic_chromatin$counts,repList=`H1-hESC`, refgene=refgene, type= genome_type, repType = "chromatin")
Tt <- covCalcPlot5prime3prime(lenChoice=lenChoice,repChoice="T", repBins=intergenic_chromatin$counts,repList=`H1-hESC`, refgene=refgene, type= genome_type, repType = "chromatin")
TSS <- covCalcPlot5prime3prime(lenChoice=lenChoice,repChoice="TSS", repBins=intergenic_chromatin$counts,repList=`H1-hESC`, refgene=refgene, type= genome_type, repType = "chromatin")
WE <- covCalcPlot5prime3prime(lenChoice=lenChoice,repChoice="WE", repBins=intergenic_chromatin$counts,repList=`H1-hESC`, refgene=refgene, type= genome_type, repType = "chromatin")




layout(matrix(c(1,2), nrow=1))

par(mar=c(5,5,5,0))
plot(-(sqrt(100001:1)),((R$rawRepCov5/R$baseFreq5prime) - predict(Rlo, gcRate5)) + R$p, type = "l", ylim = c(0,1), ylab = "",xaxt="n", xlab = "")
lines(-(sqrt(100001:1)),((E$rawRepCov5/E$baseFreq5prime)- predict(Elo, gcRate5)) + E$p, col = 2)
lines(-(sqrt(100001:1)),((PF$rawRepCov5/PF$baseFreq5prime)- predict(PFlo, gcRate5)) + PF$p, col = 3)
lines(-(sqrt(100001:1)),((CTCF$rawRepCov5/CTCF$baseFreq5prime)- predict(CTCFlo, gcRate5)) + CTCF$p, col = 4)
lines(-(sqrt(100001:1)),((Tt$rawRepCov5/Tt$baseFreq5prime)- predict(Ttlo, gcRate5)) + Tt$p, col = 5)
lines(-(sqrt(100001:1)),((TSS$rawRepCov5/TSS$baseFreq5prime)- predict(TSSlo, gcRate5)) + TSS$p, col = 6)
lines(-(sqrt(100001:1)),((WE$rawRepCov5/WE$baseFreq5prime)- predict(WElo, gcRate5)) + WE$p, col = 7)
par(new=TRUE)
new_L1_Z <- (new_L1$rawRepCov5/new_L1$baseFreq5prime - predict(new_L1Mod, gcRate5))   / ((1 * sqrt(n5 * predict(new_L1Mod, gcRate5) *  (1 - predict(new_L1Mod, gcRate5))))/n5)
plot(-(sqrt(100001:1)), new_L1_Z, ylim = c(-10,10) , col = 8, type = "l", axes = FALSE, ylab = "", xlab = "")
abline(h=2, lty=2)
abline(h=-2, lty=2)

axis(1,at = seq(0,-350, by = -50), labels = seq(0,350, by = 50)^2)



par(mar=c(5,0,5,5))

plot((sqrt(1:100001)),((R$rawRepCov3/R$baseFreq3prime)- predict(Rlo, gcRate3)) + R$p, type = "l", ylim = c(0,1), xaxt = "n", yaxt="n", ylab = "n", xlab = "")
lines((sqrt(1:100001)),((E$rawRepCov3/E$baseFreq3prime)- predict(Elo, gcRate3)) + E$p, col = 2)
lines((sqrt(1:100001)),((PF$rawRepCov3/PF$baseFreq3prime)- predict(PFlo, gcRate3)) + PF$p, col = 3)
lines((sqrt(1:100001)),((CTCF$rawRepCov3/CTCF$baseFreq3prime)- predict(CTCFlo, gcRate3)) + CTCF$p, col = 4)
lines((sqrt(1:100001)),((Tt$rawRepCov3/Tt$baseFreq3prime)- predict(Ttlo, gcRate3)) + Tt$p, col = 5)
lines((sqrt(1:100001)),((TSS$rawRepCov3/TSS$baseFreq3prime)- predict(TSSlo, gcRate3)) + TSS$p, col = 6)
lines((sqrt(1:100001)),((WE$rawRepCov3/WE$baseFreq3prime)- predict(WElo, gcRate3)) + WE$p, col = 7)
par(new=TRUE)
new_L1_Z <- (new_L1$rawRepCov3/new_L1$baseFreq3prime - predict(new_L1Mod, gcRate3))   / ((1 * sqrt(n3 * predict(new_L1Mod, gcRate3) *  (1 - predict(new_L1Mod, gcRate3))))/n3)
plot((sqrt(1:100001)),(new_L1_Z), ylim = c(-10,10)  , col = 8, type = "l",xlab = "", axes=FALSE, ylab = "")
abline(h=2, lty=2)
abline(h=-2, lty=2)

axis(1,at = seq(0,350, by = 50), labels = seq(0,350, by = 50)^2)
axis(4,at=seq(-10,10), labels=seq(-10,10))





# how can we have TEs with less than the expected TE rate the whole time 
## so need to proberly normalize for GC content
## lets get the rolling windo wto work properly










# region comparison

minS <- NULL
maxS <- NULL
minRepCov <- 5000
maxRepCov <- 7000
lenChoice = 100000
repChoice <- "L1PA"
TEintergenic <- covCalcPlot5prime3prime(repType = "repeats",lenChoice=lenChoice,repChoice=repChoice, repBins=intergenic_reps$counts,repList=rep, refgene=refgene, type= "intergenic",minBinSize=minS, maxBinSize=maxS,minRepCov=minRepCov, maxRepCov=maxRepCov)
TEintron <- covCalcPlot5prime3prime(repType = "repeats",lenChoice=lenChoice,repChoice=repChoice, repBins=intron_reps$counts,repList=rep, refgene=refgene, type= "intron", minBinSize=minS, maxBinSize=maxS,minRepCov=minRepCov, maxRepCov=maxRepCov)

logCoordinates <- c(-(log10(1:length(TEintergenic$zRepCov5))[length(TEintergenic$zRepCov5):1]), 0, log10(1:length(TEintergenic$zRepCov3)))
Coordinates <- -(length(TEintergenic$zRepCov5)):length(TEintergenic$zRepCov3)

plot(logCoordinates, c(TEintergenic$zRepCov5,NA ,TEintergenic$zRepCov3), type = "l", 
     ylim = c(min(c(TEintergenic$zRepCov5,TEintergenic$zRepCov3, 
                    TEintron$zRepCov5, TEintron$zRepCov3)),
              max(c(TEintergenic$zRepCov5,TEintergenic$zRepCov3, 
                    TEintron$zRepCov5, TEintron$zRepCov3))
              )
     )
lines(logCoordinates, c(TEintron$zRepCov5,NA ,TEintron$zRepCov3), col = 2)

abline(h=2,lty=2)
abline(h=-2,lty=2)



# accumulation comparison intergenic
# we look at families fo rthe convayerbelt effct
# we could look for sites with both quite easily 


## ALU

minS <- 50000
maxS <- NULL
minRepCov <- NULL
maxRepCov <- NULL
lenChoice = 50000
genome_type <- "intergenic"
binnedGenome <- get(paste(genome_type, "_reps", sep = ""))$counts


repChoice <- "AluJ"
AluJ <- covCalcPlot5prime3prime(lenChoice=lenChoice,repChoice=repChoice, repBins=binnedGenome,repList=rep, refgene=refgene, type= genome_type,minBinSize=minS, maxBinSize=maxS,minRepCov=minRepCov, maxRepCov=maxRepCov)
repChoice <- "AluS"
AluS <- covCalcPlot5prime3prime(lenChoice=lenChoice,repChoice=repChoice, repBins=binnedGenome,repList=rep, refgene=refgene, type= genome_type,minBinSize=minS, maxBinSize=maxS,minRepCov=minRepCov, maxRepCov=maxRepCov)
repChoice <- "AluY"
AluY <- covCalcPlot5prime3prime(lenChoice=lenChoice,repChoice=repChoice, repBins=binnedGenome,repList=rep, refgene=refgene, type= genome_type,minBinSize=minS, maxBinSize=maxS,minRepCov=minRepCov, maxRepCov=maxRepCov)


logCoordinates <- c(-(log10(1:length(AluJ$zRepCov5))[length(AluJ$zRepCov5):1]), 0, log10(1:length(AluJ$zRepCov3)))
Coordinates <- -(length(AluJ$zRepCov5)):length(AluJ$zRepCov3)

plot(logCoordinates, c(AluJ$zRepCov5,NA ,AluJ$zRepCov3), type = "l", 
     ylim = c(min(c(AluJ$zRepCov5,AluJ$zRepCov3,
                    AluS$zRepCov5,AluS$zRepCov3,
                    AluY$zRepCov5,AluY$zRepCov3)),
              max(c(AluJ$zRepCov5,AluJ$zRepCov3,
                   AluS$zRepCov5,AluS$zRepCov3,
                   AluY$zRepCov5,AluY$zRepCov3))
     ),
     main = genome_type, ylab = "TE accumulation (Z score)"
)
lines(logCoordinates, c(AluS$zRepCov5,NA ,AluS$zRepCov3), col = 2)
lines(logCoordinates, c(AluY$zRepCov5,NA ,AluY$zRepCov3), col = 3)
legend("bottomright", c("AluJ", "AluS", "AluY"), fill = c(1,2,3))

abline(h=2,lty=2)
abline(h=-2,lty=2)





#### L1 Families


minS <- NULL
maxS <- NULL
minRepCov <- 5000
maxRepCov <- 7000
lenChoice = 50000
genome_type <- "intergenic"
binnedGenome <- get(paste(genome_type, "_reps", sep = ""))$counts

repChoice <- "L1ME"
L1ME <- covCalcPlot5prime3prime(lenChoice=lenChoice,repChoice=repChoice, repBins=binnedGenome,repList=rep, refgene=refgene, type= genome_type,minBinSize=minS, maxBinSize=maxS,minRepCov=minRepCov, maxRepCov=maxRepCov)
repChoice <- "L1MD"
L1MD <- covCalcPlot5prime3prime(lenChoice=lenChoice,repChoice=repChoice, repBins=binnedGenome,repList=rep, refgene=refgene, type= genome_type,minBinSize=minS, maxBinSize=maxS,minRepCov=minRepCov, maxRepCov=maxRepCov)
repChoice <- "L1MC"
L1MC <- covCalcPlot5prime3prime(lenChoice=lenChoice,repChoice=repChoice, repBins=binnedGenome,repList=rep, refgene=refgene, type= genome_type,minBinSize=minS, maxBinSize=maxS,minRepCov=minRepCov, maxRepCov=maxRepCov)
repChoice <- "L1MB"
L1MB <- covCalcPlot5prime3prime(lenChoice=lenChoice,repChoice=repChoice, repBins=binnedGenome,repList=rep, refgene=refgene, type= genome_type,minBinSize=minS, maxBinSize=maxS,minRepCov=minRepCov, maxRepCov=maxRepCov)
repChoice <- "L1MA"
L1MA <- covCalcPlot5prime3prime(lenChoice=lenChoice,repChoice=repChoice, repBins=binnedGenome,repList=rep, refgene=refgene, type= genome_type,minBinSize=minS, maxBinSize=maxS,minRepCov=minRepCov, maxRepCov=maxRepCov)
repChoice <- "L1PB"
L1PB <- covCalcPlot5prime3prime(lenChoice=lenChoice,repChoice=repChoice, repBins=binnedGenome,repList=rep, refgene=refgene, type= genome_type,minBinSize=minS, maxBinSize=maxS,minRepCov=minRepCov, maxRepCov=maxRepCov)
repChoice <- "L1PA"
L1PA <- covCalcPlot5prime3prime(lenChoice=lenChoice,repChoice=repChoice, repBins=binnedGenome,repList=rep, refgene=refgene, type= genome_type,minBinSize=minS, maxBinSize=maxS,minRepCov=minRepCov, maxRepCov=maxRepCov)

logCoordinates <- c(-(log10(1:length(L1ME$zRepCov5))[length(L1ME$zRepCov5):1]), 0, log10(1:length(L1ME$zRepCov3)))
Coordinates <- -(length(L1ME$zRepCov5)):length(L1ME$zRepCov3)

plot(Coordinates, c(L1MA$zRepCov5,NA ,L1MA$zRepCov3), type = "n", 
     ylim = c(min(c(L1ME$zRepCov5,L1ME$zRepCov3,
                    L1MD$zRepCov5,L1MD$zRepCov3,
                    L1MC$zRepCov5,L1MC$zRepCov3,
                    L1MB$zRepCov5,L1MB$zRepCov3,
                    L1MA$zRepCov5,L1MA$zRepCov3,
                    L1PB$zRepCov5,L1PB$zRepCov3,
                    L1PA$zRepCov5,L1PA$zRepCov3
                    )),
              max(c(L1ME$zRepCov5,L1ME$zRepCov3,
                    L1MD$zRepCov5,L1MD$zRepCov3,
                    L1MC$zRepCov5,L1MC$zRepCov3,
                    L1MB$zRepCov5,L1MB$zRepCov3,
                    L1MA$zRepCov5,L1MA$zRepCov3,
                    L1PB$zRepCov5,L1PB$zRepCov3,
                    L1PA$zRepCov5,L1PA$zRepCov3
              ))
     ),
     main = genome_type, ylab = "TE accumulation (Z score)"
)
legend("bottomright", c("L1ME", "L1MD", "L1MC", "L1MB", "L1MA", "L1PB", "L1PA"), fill = c(1,2,3,4,5,6,8))

abline(h=2,lty=2)
abline(h=-2,lty=2)

lines(Coordinates, c(L1ME$zRepCov5,NA ,L1ME$zRepCov3), col = 1)
lines(Coordinates, c(L1MD$zRepCov5,NA ,L1MD$zRepCov3), col = 2)
lines(Coordinates, c(L1MC$zRepCov5,NA ,L1MC$zRepCov3), col = 3)
lines(Coordinates, c(L1MB$zRepCov5,NA ,L1MB$zRepCov3), col = 4)
lines(Coordinates, c(L1MA$zRepCov5,NA ,L1MA$zRepCov3), col = 5)
lines(Coordinates, c(L1PB$zRepCov5,NA ,L1PB$zRepCov3), col = 6)
lines(Coordinates, c(L1PA$zRepCov5,NA ,L1PA$zRepCov3), col = 8)









### lets try joiing families

intergenicJoinedFam <- data.frame(intergenic_reps$counts[,c("chr", "start", "end","Known")],
                        old_L1 = rowSums(intergenic_reps$counts[,c("L1ME", "L1MD", "L1MC","L1MB")]),
                        new_L1 = rowSums(intergenic_reps$counts[,c("L1MA", "L1PB", "L1PA","L1HS")]),
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
                new_L1 = rbind(rep$L1MA, rep$L1PB, rep$L1PA, rep$L1HS),
                Alu = rbind(rep$AluJ, rep$AluS, rep$AluY),
                Ancient = rbind(rep$MIR, rep$L2))


minS <- NULL
maxS <- NULL
minRepCov <- NULL
maxRepCov <- NULL
lenChoice = 50000
genome_type <- "intron"
binnedGenome <- get(paste(genome_type , "JoinedFam", sep = ""))



repChoice <- "old_L1"
old_L1 <- covCalcPlot5prime3prime(lenChoice=lenChoice,repChoice=repChoice, repBins=binnedGenome,repList=joinRep, refgene=refgene, type= genome_type,minBinSize=minS, maxBinSize=maxS,minRepCov=NULL, maxRepCov=NULL)
repChoice <- "new_L1"
new_L1 <- covCalcPlot5prime3prime(lenChoice=lenChoice,repChoice=repChoice, repBins=binnedGenome,repList=joinRep, refgene=refgene, type= genome_type,minBinSize=minS, maxBinSize=maxS,minRepCov=NULL, maxRepCov=NULL)
repChoice <- "Alu"
Alu <- covCalcPlot5prime3prime(lenChoice=lenChoice,repChoice=repChoice, repBins=binnedGenome,repList=joinRep, refgene=refgene, type= genome_type,minBinSize=minS, maxBinSize=maxS,minRepCov=NULL, maxRepCov=NULL)
repChoice <- "Ancient"
Ancient <- covCalcPlot5prime3prime(lenChoice=lenChoice,repChoice=repChoice, repBins=binnedGenome,repList=joinRep, refgene=refgene, type= genome_type,minBinSize=minS, maxBinSize=maxS,minRepCov=NULL, maxRepCov=NULL)

logCoordinates <- c(-(log10(1:(lenChoice + 1))[(lenChoice + 1):1]), 0, log10(1:(lenChoice + 1)))
Coordinates <- -((lenChoice + 1)):(lenChoice + 1)


plot(logCoordinates, c(old_L1$zRepCov5,NA ,old_L1$zRepCov3), type = "l", 
     ylim = c(min(c(old_L1$zRepCov5,old_L1$zRepCov3,
                    new_L1$zRepCov5,new_L1$zRepCov3,
                    Alu$zRepCov5,Alu$zRepCov3,
                    Ancient$zRepCov5,Ancient$zRepCov3)),
              max(c(old_L1$zRepCov5,old_L1$zRepCov3,
                    new_L1$zRepCov5,new_L1$zRepCov3,
                    Alu$zRepCov5,Alu$zRepCov3,
                    Ancient$zRepCov5,Ancient$zRepCov3))
     ),
     main = genome_type, ylab = "TE accumulation (Z score)"
)

lines(logCoordinates, c(new_L1$zRepCov5,NA ,new_L1$zRepCov3), col = 2)
lines(logCoordinates, c(Alu$zRepCov5,NA ,Alu$zRepCov3), col = 3)
lines(logCoordinates, c(Ancient$zRepCov5,NA ,Ancient$zRepCov3), col = 4)


legend("bottomright", names(joinRep), fill = c(1,2,3,4))

abline(h=2,lty=2)
abline(h=-2,lty=2)








plot((1:length(ALuSintron$zRepCov)), ALuSintron$zRepCov,col = 1, type = "l")
lines((1:length(L1MEintron$zRepCov)), L1MEintron$zRepCov,col = 2, type = "l")



abline(h = 2)


# maybe new insertions are accoicated with a specific type of histone mark 

plot((1:length(ALuSintron$zRepCov)), ALuSintron$rawRepCov, type = "l")
lines((1:length(ALuSintron$zRepCov)), ALuS$rawRepCov, col = 2)
par(new = T)
plot((1:length(ALuSintron$zRepCov)), ALuSintron$rawRepCov/ALuSintron$baseFreq, type = "l")
lines((1:length(ALuSintron$zRepCov)), ALuS$rawRepCov/ALuSintergenic$baseFreq, col= 2)


par(new=T)
plot((1:length(ALuSintron$zRepCov)), ALuSintron$zRepCov, type = "l", col = 1)
lines((1:length(ALuSintron$zRepCov)), ALuSintergenic$zRepCov, type = "l", col = 2)

# as far as i can tell, enrichment tends to be at these smaller sizes. 


plot((1:length(ALuSintron$zRepCov)), ALuSintron$baseFreq/sum(ALuSintron$baseFreq), type = "l")
lines((1:length(ALuSintergenic$zRepCov)), ALuSintergenic$baseFreq/sum(ALuSintergenic$baseFreq), type = "l", col = 2)


plot((1:length(ALuSintron$zRepCov)), ALuSintron$baseFreq, type = "l")
lines((1:length(ALuSintergenic$zRepCov)), ALuSintergenic$baseFreq, type = "l", col = 2)









# There may be an effect where shorter intorns tend to be from shorter genes that 
# tend to get expressed more. This means they are more chromaticly open 
# leading to accumulation of TEs and a confounding factor in our results
# One straturdy is to stratify on Intron size. 
# I should cut out small introns and then aother patterns may begin to immerge

# if we are going to do the bidirectional plot we need to incorperate strand information. 
# Essentially know if something is up or downstream of a gene 
# also this will change the number of up and downstream regions we are looking at





TEchioce = "AluY"
datType = "counts"

plot(log10(intron_reps$rates$Known),  log10(intron_reps[[datType]][,TEchioce]) , pch = 16, cex = .1, col = 2, main = TEchioce,
     ylab = paste("log10", datType), xlab = "log10 interval size")

points(log10(intergenic_reps$rates$Known),log10(intergenic_reps[[datType]][,TEchioce]),  pch = 16, cex = .1)
legend("topleft", legend = c("intergenic", "intron"), fill = c(1,2))
a = (1:100)
b = (1:100)

fit <- lm(b~a)
abline(lm(b~a))



plot(log10(intron_reps$rates$Known -(intron_reps[[datType]][,TEchioce]) )  ,  log10(intron_reps[[datType]][,TEchioce]) , pch = 16, cex = .1, col = 2, main = TEchioce,
     ylab = paste("log10", datType), xlab = "log10 interval size")
points(log10( (intergenic_reps$rates$Known) -intergenic_reps[[datType]][,TEchioce] ) ,log10(intergenic_reps[[datType]][,TEchioce]),  pch = 16, cex = .1)



# can these differences just be a result of more smaller spaces

# should we stratify on region size

# also we are thinking about similar outcomes for both regions but by different mechanisms. 


plot(density(log10(intron_reps$counts$Known)), col = 2, xlab = "log10 interval length")
lines(density(log10(intergenic_reps$counts$Known)), col = 1)
legend("topleft", legend=c("intergenic", "intron"), fill=c(1,2))


bin_strat <- c(seq(100, 900, by = 100), seq(1000, 9000, by = 1000), seq(10000, 100000, by = 10000))
result_tabInter <- data.frame(intervalSize = bin_strat[1:(length(bin_strat)-1)], mean = NA, sd = NA)
result_tabIntron <- data.frame(intervalSize = bin_strat[1:(length(bin_strat)-1)], mean = NA, sd = NA)


repChoice = "AluS"
InterDat = NULL
for(i in 1:nrow(result_tabInter)){
  selection = intergenic_reps$rates[,repChoice][intergenic_reps$rates$Known < bin_strat[i+1] & intergenic_reps$rates$Known > bin_strat[i]]
  selection = selection[selection > 0]
  result_tabInter$mean[i] = mean(selection)
  result_tabInter$sd[i] = sd(selection)
  InterDat <- c(InterDat, list(selection))
}

IntronDat = NULL
for(i in 1:nrow(result_tabIntron)){
  selection = intron_reps$rates[,repChoice][intron_reps$rates$Known < bin_strat[i+1] & intron_reps$rates$Known > bin_strat[i]]
  selection = selection[selection > 0]
  result_tabIntron$mean[i] = mean(selection)
  result_tabIntron$sd[i] = sd(selection)
  IntronDat <- c(IntronDat, list(selection))
  
}
names(InterDat) <- result_tabInter$intervalSize
names(IntronDat) <- result_tabIntron$intervalSize



plot(log10(result_tabInter$intervalSize), result_tabInter$mean, type = "b", ylim=c(1,10000))
lines(log10(result_tabInter$intervalSize), result_tabInter$mean + (result_tabInter$sd), lty = 2)
lines(log10(result_tabInter$intervalSize), result_tabInter$mean-(result_tabInter$sd), lty = 2)

lines(log10(result_tabIntron$intervalSize), result_tabIntron$mean, type = "b", col = 2)
lines(log10(result_tabIntron$intervalSize), result_tabIntron$mean + (result_tabIntron$sd), lty = 2, col = 2)
lines(log10(result_tabIntron$intervalSize),  result_tabIntron$mean - (result_tabIntron$sd), lty = 2, col = 2)


boxplot(IntronDat,outline=T, ylim = c(1,10000), col = 2)
par(new = T)
boxplot(InterDat,outline=T, ylim = c(1,10000))


boxplot(InterDat,outline=T, ylim = c(1,10000))
par(new = T)
boxplot(IntronDat,outline=T, ylim = c(1,10000), col = 2)


boxplot(IntronDat,outline=T, ylim = c(1,10000), col = 2)
boxplot(InterDat,outline=T, ylim = c(1,10000))




# other experiments to try will be to look at heterochromatin coverage
# and DNase1 coverage
# also coorect for GC content

# not too hard, there is a bioconductor packakage we can use to get that info. 
# we could potentially do this today 
# we know which ranges we are grabbing for both intron and intergenic 
# all we need is to factor it in
### lets try pull out GC content

Seq=Hsapiens[["chr1"]]


chrChoice <- "chr21"

Seq <- Hsapiens[["chr21"]]

# we can do this for each region via chromosome. 
Seq.set=DNAStringSet(Hsapiens[[chrChoice]], start=start(refgene_gap.gr[seqnames(refgene_gap.gr) == chrChoice]), end=end(refgene_gap.gr[seqnames(refgene_gap.gr) == chrChoice]))
bases=alphabetFrequency(Seq.set,baseOnly=TRUE)
CGcontent <- rowSums(bases[,c("G", "C")]) / rowSums(bases[,c("A", "C", "T", "G")])



# find 

chromos <- as.character(unique(intergenic_reps$counts$chr))
IntergenicGC <- rep(NA, nrow(intergenic_reps$counts))
for(i in 1:length(chromos)){
  Seq.set=DNAStringSet(Hsapiens[[chromos[i]]], 
                       start=intergenic_reps$counts$start[intergenic_reps$counts$chr == chromos[i]], 
                       end=intergenic_reps$counts$end[intergenic_reps$counts$chr == chromos[i]])
  bases=data.frame(alphabetFrequency(Seq.set,baseOnly=TRUE))
  CGcontent <- rowSums(bases[,c("G", "C")]) / rowSums(bases[,c("A", "C", "T", "G")])
  IntergenicGC[intergenic_reps$counts$chr == chromos[i]] <- CGcontent
  
}

IntronGC <- rep(NA, nrow(intron_reps$counts))
for(i in 1:length(chromos)){
  Seq.set=DNAStringSet(Hsapiens[[chromos[i]]], 
                       start=intron_reps$counts$start[intron_reps$counts$chr == chromos[i]], 
                       end=intron_reps$counts$end[intron_reps$counts$chr == chromos[i]])
  bases=data.frame(alphabetFrequency(Seq.set,baseOnly=TRUE))
  CGcontent <- rowSums(bases[,c("G", "C")]) / rowSums(bases[,c("A", "C", "T", "G")])
  IntronGC[intron_reps$counts$chr == chromos[i]] <- CGcontent
  
}


repChoice <- "AluS"
intronGCset <- IntronGC[intron_reps$rates[,repChoice] > 0]
intronTEset <- intron_reps$rates[intron_reps$rates[,repChoice] > 0,repChoice]
intergenicGCset <- IntergenicGC[intergenic_reps$rates[,repChoice] > 0]
intergenicTEset <- intergenic_reps$rates[intergenic_reps$rates[,repChoice] > 0,repChoice]

plot((intronGCset),
      (intronTEset),
      pch = 16, cex = .05, col = 2 , main = repChoice,
     ylab = "TE coverage rate (log10 per 10000bp)",
     xlab = "GC content"
      )
points((intergenicGCset),
     (intergenicTEset), 
     pch = 16, cex = .05)
interMod <- lm((intergenicTEset) ~ (intergenicGCset))
intronMod <- lm((intronTEset) ~ (intronGCset))


intronMod2 <- loess((intronTEset) ~ (intronGCset))
interMod2 <- loess((intergenicTEset) ~ (intergenicGCset))
lines(seq(.1,.65,.01), predict(object=intronMod2, newdata=seq(.1,.65,.01)), col = 2, lwd = 2)
lines(seq(.1,.65,.01), predict(object=interMod2, newdata=seq(.1,.65,.01)), col = 1, lwd = 2)



abline(intronMod, col = 2, lwd = 2, lty = 2)
abline(interMod, col = 1, lwd = 2, , lty = 2)
legend("bottomright", legend=c("intergenic", "intron"), fill = c(1,2))

legend("topright", bty="n", legend=c(paste("intergenic r =", format(sqrt(summary(interMod)$r.squared), digits=3)), 
                                     paste("intron r =", format(sqrt(summary(intronMod)$r.squared), digits=3))
                                     )
       )




### if we can extract the Cs and Gs and put them in a rep style list opject
### then we can plug them straight into our annalysis 
### the Cs and Gs can be in any order



# this is essentially another way to work it out 
# it is probably quiclker then calculating each one individually and then sending it through our script
# alternativly 
# we could just add the start position of the interval 
# this function gcRepObject


# lets see if we can vectorise when we calculate the mean 
 chromos <- as.character(unique(repBins$chr))
  repBinsGC <- rep(NA, nrow(repBins))
  for(i in 1:length(chromos)){
    Seq.set=DNAStringSet(Hsapiens[[chromos[i]]], 
                         start=repBins$start[repBins$chr == chromos[i]], 
                         end=repBins$end[repBins$chr == chromos[i]])
    bases=data.frame(alphabetFrequency(Seq.set,baseOnly=TRUE))
    CGcontent <- rowSums(bases[,c("G", "C")]) / rowSums(bases[,c("A", "C", "T", "G")])
    repBinsGC[repBins$chr == chromos[i]] <- CGcontent
    
  }
  
  

  chromos <- as.character(unique(repBins$chr))
  repBinsGC <- rep(NA, nrow(repBins))
  
  
  for(i in 1:length(chromos)){
    Seq.set=DNAStringSet(Hsapiens[[chromos[i]]], 
                         start=repBins$start[repBins$chr == chromos[i]], 
                         end=repBins$end[repBins$chr == chromos[i]])
    bases=data.frame(alphabetFrequency(Seq.set,baseOnly=TRUE))
    CGcontent <- rowSums(bases[,c("G", "C")]) / rowSums(bases[,c("A", "C", "T", "G")])
    repBinsGC[repBins$chr == chromos[i]] <- CGcontent
    
  }
  
  
# calculate residual 


# maybe just get bins 10,000 bp 
# do the regression analysis for each repeat 



library(BSgenome)
bsGenome <- available.genomes()[grep(genome,available.genomes())]
bsGenome <- bsGenome[-(grep("masked", bsGenome))]
library(bsGenome, character.only=TRUE)
bsSpec <- get(strsplit(bsGenome,split="\\.")[[1]][2])

bins <- binned.genome.reader(genome=genome, bin.size=20000, keep.rate=.9)
bins <- bins[[1]]
# remove unplaced chromosomes from bins
if(length(grep("_", bins$chr)) > 0){
  bins <- bins[-(grep("_", bins$chr)),]
}


bin.samp <- bins[sample(1:nrow(bins),10000),]

bin.sort = binSort(rep=rep, bins=bin.samp,TE.names=names(rep))
bin.samp = bin.sort$rates
bin.samp[,names(rep)] = bin.sort$rates[,names(rep)]/10000

repChoice = "L1PA"
rep.gr <- GRanges(seqnames = Rle(rep[[repChoice]]$genoName),
                  ranges = IRanges(start = rep[[repChoice]]$genoStart, end = rep[[repChoice]]$genoEnd))
bin.gr <- GRanges(seqnames=Rle(bin.samp$chr), 
                  ranges = IRanges(start=bin.samp$start, end=bin.samp$end-1))

diff.gr <- setdiff(bin.gr,rep.gr)
bin.diff <- as.matrix(findOverlaps(bin.gr, diff.gr))
if(length((1:nrow(bin.diff))[duplicated(bin.diff[,2])])!=0){
  print("some bins merged")
}
# bin diff provides our map for going back 



# we will get the GC percent of these regions and weight them by their length and get the mean

chromos <- unique(as.character(seqnames(diff.gr)))
bin.sortGC <- rep(NA, length(diff.gr))

for(i in 1:length(chromos)){
  Seq.set=DNAStringSet(bsSpec[[chromos[i]]], 
                       start=start(diff.gr)[as.character(seqnames(diff.gr)) == chromos[i]], 
                       end=end(diff.gr)[as.character(seqnames(diff.gr)) == chromos[i]]
                       )
  bases=data.frame(alphabetFrequency(Seq.set,baseOnly=TRUE))
  CGcontent <- rowSums(bases[,c("G", "C")]) / rowSums(bases[,c("A", "C", "T", "G")])
  bin.sortGC[as.character(seqnames(diff.gr)) == chromos[i]] <- CGcontent
  print(chromos[i])
}

getContent <- data.frame(binID = unique(bin.diff[,1]), GCcontent = NA)
for( i in 1:nrow(getContent)){
  w <- width(diff.gr)[bin.diff[getContent$binID[i] == bin.diff[,1],2]]
  getContent[i,"GCcontent"] <- weighted.mean(bin.sortGC[bin.diff[getContent$binID[i] == bin.diff[,1],2]], w=w)
}



### old way

chromos <- as.character(unique(bin.samp$chr))
bin.sortGC <- rep(NA, nrow(bin.samp))

for(i in 1:length(chromos)){
  Seq.set=DNAStringSet(bsSpec[[chromos[i]]], 
                       start=bin.samp$start[bin.samp$chr == chromos[i]], 
                       end=bin.samp$end[bin.samp$chr == chromos[i]])
  bases=data.frame(alphabetFrequency(Seq.set,baseOnly=TRUE))
  CGcontent <- rowSums(bases[,c("G", "C")]) / rowSums(bases[,c("A", "C", "T", "G")])
  bin.sortGC[bin.samp$chr == chromos[i]] <- CGcontent
  print(chromos[i])
}



# this way we can look for changes





plot(bin.sortGC[bin.samp[,repChoice] > 0],bin.samp[,repChoice][bin.samp[,repChoice] > 0], pch = 16, cex = .3, ylim = c(0,1), xlim = c(0,1))



modLM <- lm(bin.samp[,repChoice][bin.samp[,repChoice] > 0] ~ bin.sortGC[bin.samp[,repChoice] > 0])
abline(modLM, lty = 2, lwd = 2, col = 2)

modLO <- loess(bin.samp[,repChoice][bin.samp[,repChoice] > 0] ~ bin.sortGC[bin.samp[,repChoice] > 0])
lines(seq(0,1,.01), predict(modLO, newdata=seq(0,1,.01)), col = 2, lwd = 2)






plot(getContent$GCcontent,bin.samp[,repChoice][getContent$binID], pch = 16, cex = .3, ylim = c(0,1), xlim = c(0,1))


modLM <- lm(bin.samp[,repChoice][bin.samp[,repChoice] > 0] ~ bin.sortGC[bin.samp[,repChoice] > 0])
abline(modLM, lty = 2, lwd = 2, col = 2)

modLO <- loess(bin.samp[,repChoice][bin.samp[,repChoice] > 0] ~ bin.sortGC[bin.samp[,repChoice] > 0])
lines(seq(0,1,.01), predict(modLO, newdata=seq(0,1,.01)), col = 2, lwd = 2)


# GC content is being interesting, I wonder if getting the content of the regions then taking the mean over all the bases the right approach

# we seem to be taking the avarge rate
# rather than the rate over a position. 

# also it might be worth calculating the non TE GC rate. 
# it would be quite easy # we just take the setdiff
# and use those coordinates and map it back to the bin



# we could look at the difference between the two approaches 
# diff in GC contnet 
# against GC content

plot(bin.samp[getContent$binID,repChoice],getContent$GCcontent / bin.sortGC[getContent$binID])

plot(bin.samp[getContent$binID,repChoice],getContent$GCcontent - bin.sortGC[getContent$binID])
plot(getContent$GCcontent,bin.sortGC[getContent$binID])


# maybe we need to think about a different function for capturing cg adjustments



# thinking about how to control some potential sources of confunding 
# There will always be the problem of the repeat sequence contributing to the GC content 


# it is probably best to find at what level is the system robust to this 

# the probability a repeat picked at random has surrounding gc content of X




###### OLD APPROACH


chromos <- as.character(unique(bin.samp$chr))
bin.sortGC <- rep(NA, nrow(bin.samp))

for(i in 1:length(chromos)){
  Seq.set=DNAStringSet(bsSpec[[chromos[i]]], 
                       start=bin.samp$start[bin.samp$chr == chromos[i]], 
                       end=bin.samp$end[bin.samp$chr == chromos[i]])
  bases=data.frame(alphabetFrequency(Seq.set,baseOnly=TRUE))
  CGcontent <- rowSums(bases[,c("G", "C")]) / rowSums(bases[,c("A", "C", "T", "G")])
  bin.sortGC[bin.samp$chr == chromos[i]] <- CGcontent
  print(chromos[i])
}


# maybe just sample 1000 or something 




plot(bin.sortGC[bin.samp[,repChoice] > 0],bin.samp[,repChoice][bin.samp[,repChoice] > 0], pch = 16, cex = .3, ylim = c(0,1), xlim = c(0,1))

modLM <- lm(bin.samp[,repChoice][bin.samp[,repChoice] > 0] ~ bin.sortGC[bin.samp[,repChoice] > 0])
abline(modLM, lty = 2, lwd = 2, col = 2)

modLO <- loess(bin.samp[,repChoice] ~ bin.sortGC)
lines(seq(0,1,.01), predict(modLO, newdata=seq(0,1,.01)), col = 2, lwd = 2)


predict(modLO, .8)


TEs <- covCalcPlot5prime3primeGC(lenChoice=100000,repChoice="L1PA",repBins=intergenic_reps$counts,refgene=refgene,type="intergenic",genome=genome,repList=rep)


# so now we have the GC regression slope

plot(TEs_intergenic$prime3gc/TEs_intergenic$baseFreq3prime, type = "l", col = 2)
par(new=TRUE)
plot(TEs_intergenic$rawRepCov3/TEs_intergenic$baseFreq3prime, type = "l", yaxt ="n")

axis(4,at=seq(0,1,by=.01), labels=seq(0,1,by=.01))

# how do we get the expected TE level based on GC content 
plot(predict(modLO, TEs_intergenic$prime3gc/TEs_intergenic$baseFreq3prime), type = "l")

lines(TEs_intergenic$rawRepCov3/TEs_intergenic$baseFreq3prime, col=2)




#TEs_intorn <- covCalcPlot5prime3primeGC(lenChoice=100000,repChoice="L1PA",repBins=intron_reps$counts,refgene=refgene,type="intron",genome=genome,repList=rep)



hist(TEs_intorn$prime3gc/TEs_intorn$baseFreq3prime)

plot(TEs_intorn$prime3gc/TEs_intorn$baseFreq3prime, type = "l", col = 2)
par(new=TRUE)
plot(TEs_intorn$rawRepCov3/TEs_intorn$baseFreq3prime, type = "l", yaxt ="n")

axis(4,at=seq(0,1,by=.01), labels=seq(0,1,by=.01))

# how do we get the expected TE level based on GC content 
plot((1:100001),predict(modLO, TEs_intorn$prime3gc/TEs_intorn$baseFreq3prime), type = "l", ylim = c(0, .1))

lines((1:100001),TEs_intorn$rawRepCov3/TEs_intorn$baseFreq3prime, col=2)


# why do the values never go above a certain number



plot(TEs_intorn$prime3gc/TEs_intorn$baseFreq3prime, TEs$rawRepCov3/TEs$baseFreq3prime)

plot(TEs$rawRepCov3/TEs$baseFreq3prime , predict(modLO, TEs_intorn$prime3gc/TEs_intorn$baseFreq3prime))



tail(predict(modLO, TEs_intorn$prime3gc/TEs_intorn$baseFreq3prime))

# why is it the expected value does not go above a certain limit
# is it becuase the GC range is out of the scope of the regression function

plot(TEs_intorn$prime3gc/TEs_intorn$baseFreq3prime, type = "l")


gc.no <- TEs_intorn$prime3gc
# lets start doing the rolling mean appraoch of the gc content about 10000bp

# maybe it should be a rolling weighted mean
# one easy way would be to set up the coordinates 
# fill in the weights at each stage 
# fill in the values
# at each point we have the GC number or rate and the base frequency

library(zoo)

CG_coorect <- rollsum()


# we have a weighted mean problem 



z <- zoo(c(NA, 2, NA, 1, 4, 5, 2, NA))
na.fill(z)
na.fill(z, "extend")
na.fill(z, c("extend", NA))
na.fill(z, -(1:3))
na.fill(z, list(NA, NULL, NA))


x.Date <- as.Date(paste(2004, rep(1:4, 4:1), sample(1:28, 10), sep = "-"))
x <- zoo(rnorm(12), x.Date)

x
rollsum(x, 3,align="left" ,fill="extend")
rollsum(x, 3,align="right" ,fill="extend")

k = 20000
lenChoice = 100000


names(gc.no) <- 1:(lenChoice+1)

k = 1000
GCsum <- rollsum(TEs_intorn$prime3gc,k = k, fill="extend", align = "left")
baseSum <- rollsum(TEs_intorn$baseFreq3prime,k = k, fill="extend", align = "left")
TEsum <- rollsum(TEs_intorn$rawRepCov3,k = k, fill="extend", align = "left")
gcRate <- GCsum/baseSum


plot(TEsum/baseSum, type = "l")
lines(predict(modLO, gcRate), col = 2)

p = TEs_intorn$totalRep/TEs_intorn$totalCov
n = TEs_intorn$baseFreq3prime


# how would we adjust by GC content 

#lines((p*n)/n - predict(modLO, gcRate) + ((p*n)/n))
plot((TEsum/baseSum) - predict(modLO, gcRate) + ((p*n)/n) , type = "l")
lines( ((p*n + 2*((p*n) / sqrt( p*n*(1-p) )))/n)   , col = 2  )
lines( ((p*n - 2*((p*n) / sqrt( p*n*(1-p) )))/n)  , col = 2  )

### we need to sample and heatmap this thing to see the dynamics



lines( ((p*n + 2*((p*n) / sqrt( p*n*(1-p) )))/n)  - predict(modLO, gcRate) + ((p*n)/n) , col = 2  )
lines( ((p*n - 2*((p*n) / sqrt( p*n*(1-p) )))/n)  - predict(modLO, gcRate) + ((p*n)/n) , col = 2  )

par(new=TRUE)
plot(n, type = "l")



lines( (p*n - 2*((p*n) / sqrt( p*n*(1-p) )))/n   )


##### the next thing is to normalise aginst all bins 
#### that is the best way to adhust the regression line 
#### furthermore we know this will change the CG dynamics



# why are we still getting these increases right near the front.

# we want p + the resdidual 




k = 1000


GCsum <- rollsum(TEs$prime3gc,k = k, fill="extend", align = "left")
baseSum <- rollsum(TEs$baseFreq3prime,k = k, fill="extend", align = "left")
TEsum <- rollsum(TEs$rawRepCov3,k = k, fill="extend", align = "left")
gcRate <- GCsum/baseSum


plot(TEsum/baseSum, type = "l", ylim=c(0, max(TEsum/baseSum)))
lines(predict(modLO, gcRate), col = 2)
lines(n*p/n)


p = TEs$totalRep/TEs$totalCov
n = TEs$baseFreq3prime


# how would we adjust by GC content 

#lines((p*n)/n - predict(modLO, gcRate) + ((p*n)/n))
plot((TEs$rawRepCov3/TEs$baseFreq3prime) - predict(modLO, gcRate) + ((p*n)/n) , type = "l")
lines( ((p*n + 2*((p*n) / sqrt( p*n*(1-p) )))/n)   , col = 2  )
lines( ((p*n - 2*((p*n) / sqrt( p*n*(1-p) )))/n), col = 2  )





lines( ((p*n + 2*((p*n) / sqrt( p*n*(1-p) )))/n)  - predict(modLO, gcRate) + ((p*n)/n) , col = 2  )
lines( ((p*n - 2*((p*n) / sqrt( p*n*(1-p) )))/n)  - predict(modLO, gcRate) + ((p*n)/n) , col = 2  )



# work on correct sampling stratergies. 









plot(TEs_intorn$rawRepCov3/TEs_intorn$baseFreq3prime, type = "l")
lines(predict(modLO, gcRate),type = "l", col = 2)
lines(predict(modLO, TEs_intorn$prime3gc/TEs_intorn$baseFreq3prime), col = 3)

coords <- data.frame(start = seq(0-k, lenChoice, by = 1), end = seq(2, lenChoice + 1 +k +1, by = 1))
coords$start[coords$start < 1] <- 1
coords$end[coords$end > lenChoice + 1] <- lenChoice +1




A <- gc.no[coords$start : coords$end]




plot(coords$start,type = "l")
lines(coords$end)
