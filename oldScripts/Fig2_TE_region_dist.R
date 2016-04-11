# now we can begin to combine the anlysis of both kinds of regions
# this will be our fig 2

# where we look at the relationship between TE accumulation and genes

# the idea is that TEs accumulate according to our model around genes

# on both fractions theres constraint pertaining to the amount but the constraint on position is different


#### maybe tomorrow i can compare both methods
### now that I've found a way of normalising for this
# still havnt been able to get the 




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
intergenic_reps <- binSort(rep=rep, bins=bins_gene_gap, TE.names=names(rep), repType = rep("repeats",length(rep)))
intron_reps <- binSort(rep=rep, bins=bins_intron, TE.names=names(rep), repType = rep("repeats",length(rep)))

intergenic_chromatin <- binSort(rep = `H1-hESC`, bins = bins_gene_gap, TE.names = names(`H1-hESC`), repType = rep("chromatin",length(`H1-hESC`)))
intron_chromatin <- binSort(rep = `H1-hESC`, bins = bins_intron, TE.names = names(`H1-hESC`), repType = rep("chromatin",length(`H1-hESC`)))

#interesting that intervals with repeats in them tend to be over 1000 in length

intergenicSample <- sample(x = 1:nrow(intergenic_reps$counts), size = 10000, replace = FALSE)
intronSample <- sample(x = 1:nrow(intron_reps$counts), size = 10000, replace = FALSE)


intergenic_chromatinSamp <- intergenic_chromatin$counts[intergenicSample,]
intron_chromatinSamp <- intron_chromatin$counts[intronSample,]

TEs_intergenic <- covCalcPlot5prime3primeGC(lenChoice=100000,repBins=intergenic_reps$counts[intergenicSample,],refgene=refgene,type="intergenic",genome=genome)
TEs_intron <- covCalcPlot5prime3primeGC(lenChoice=40000,repBins=intron_reps$counts[intronSample,],refgene=refgene,type="intron",genome=genome)

save(TEs_intergenic, file = "~/Desktop/Domain_manuscript/R_objects/TEs_intergenicGC")
save(TEs_intron, file = "~/Desktop/Domain_manuscript/R_objects/TEs_intronGC")

save(intergenicSample, file = "~/Desktop/Domain_manuscript/R_objects/intergenicSample")
save(intronSample, file = "~/Desktop/Domain_manuscript/R_objects/intronSample")


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

joinRepChromatin <- c(joinRep, `H1-hESC`)

repTypes <- c(rep("repeats", length(joinRep)), rep("chromatin", length(`H1-hESC`)))
#joinSampGenomeBig <- localGCgenome(repList = joinRepChromatin,binSize = 20000,sampSize = 5000, genome = genome, repType = repTypes)
#joinSampGenomeSmall <- localGCgenome(repList = joinRepChromatin,binSize = 1000,sampSize = 5000, genome = genome, repType = repTypes)
joinSampGenomeMed <- localGCgenome(repList = joinRepChromatin,binSize = 10000,sampSize = 10000, genome = genome, repType = repTypes)
#joinSampGenomeLarge <- localGCgenome(repList = joinRepChromatin,binSize = 30000,sampSize = 10000, genome = genome, repType = repTypes)

# increased the sampling size of our random genome for our adjustments 


# maybe do it all at the 1000bp level. 
 
# we can put the heterochromatin in here 

### adjust GC contnent 
layout(matrix(c(1)))
par(mar=c(5,5,5,5))

plot(joinSampGenomeMed$R, joinSampGenomeMed$GC, pch = 16, cex = .2)
RgcLoMed <- loess(joinSampGenomeMed$GC ~ joinSampGenomeMed$R)
lines(seq(0,1,0.01), predict(RgcLoMed, newdata = seq(0,1,0.01)), col = 2, lwd = 2)
GC.adjMed <- (joinSampGenomeMed$GC - predict(RgcLoMed, joinSampGenomeMed$R) ) + mean(joinSampGenomeMed$GC)

plot(joinSampGenomeMed$GC, joinSampGenomeMed$R, pch = 16, cex = .2)
rMod <- loess(joinSampGenomeMed$R ~ joinSampGenomeMed$GC)
lines(seq(0,1,0.01), predict(rMod, newdata = seq(0,1,0.01)), col = 2, lwd = 2)


layout(matrix(c(1,2,3,4), nrow = 2, ncol = 2))
par(mar=c(4,4,4,4))

plot(GC.adjMed, joinSampGenomeMed$old_L1, main = "old_L1", pch = 16, cex = .1, xlab = "GC content", ylab = "TE fraction", ylim = c(0,1), col = "darkred")
points(joinSampGenomeMed$GC, joinSampGenomeMed$old_L1, main = "old_L1", pch = 16, cex = .1, ylim = c(0,1), col = "darkblue")
old_L1Mod <- loess(joinSampGenomeMed$old_L1 ~ joinSampGenomeMed$GC)
lines(seq(0,1,0.01), predict(old_L1Mod, newdata = seq(0,1,0.01)), lwd = 3, col = "blue")
old_L1ModAdj <- loess(joinSampGenomeMed$old_L1 ~ GC.adjMed)
lines(seq(0,1,0.01), predict(old_L1ModAdj, newdata = seq(0,1,0.01)), lwd = 3, col = "red")
legend("topright", legend = c("Adjusted", "Raw"), fill = c(2,4), title = "GC data")

plot(GC.adjMed, joinSampGenomeMed$new_L1, main = "new_L1", pch = 16, cex = .1, xlab = "GC content", ylab = "TE fraction", ylim = c(0,1), col = "darkred")
points(joinSampGenomeMed$GC, joinSampGenomeMed$new_L1, main = "new_L1", pch = 16, cex = .1, ylim = c(0,1), col = "darkblue")
new_L1Mod <- loess(joinSampGenomeMed$new_L1 ~ joinSampGenomeMed$GC)
lines(seq(0,1,0.01), predict(new_L1Mod, newdata = seq(0,1,0.01)), lwd = 3, col = "blue")
new_L1ModAdj <- loess(joinSampGenomeMed$new_L1 ~ GC.adjMed)
lines(seq(0,1,0.01), predict(new_L1ModAdj, newdata = seq(0,1,0.01)), lwd = 3, col = "red")
legend("topright", legend = c("Adjusted", "Raw"), fill = c(2,4), title = "GC data")


plot(GC.adjMed, joinSampGenomeMed$Alu, main = "Alu", pch = 16, cex = .1, xlab = "GC content", ylab = "TE fraction", ylim = c(0,1), col = "darkred")
points(joinSampGenomeMed$GC, joinSampGenomeMed$Alu, main = "Alu", pch = 16, cex = .1, ylim = c(0,1), col = "darkblue")
AluMod <- loess(joinSampGenomeMed$Alu ~ joinSampGenomeMed$GC)
lines(seq(0,1,0.01), predict(AluMod, newdata = seq(0,1,0.01)), lwd = 3, col = "blue")
AluModAdj <- loess(joinSampGenomeMed$Alu ~ GC.adjMed)
lines(seq(0,1,0.01), predict(AluModAdj, newdata = seq(0,1,0.01)), lwd = 3, col = "red")
legend("topright", legend = c("Adjusted", "Raw"), fill = c(2,4), title = "GC data")


plot(GC.adjMed, joinSampGenomeMed$Ancient, main = "Ancient", pch = 16, cex = .1, xlab = "GC content", ylab = "TE fraction", ylim = c(0,1), col = "darkred")
points(joinSampGenomeMed$GC, joinSampGenomeMed$Ancient, main = "Ancient", pch = 16, cex = .1, ylim = c(0,1), col = "darkblue")
AncientMod <- loess(joinSampGenomeMed$Ancient ~ joinSampGenomeMed$GC)
lines(seq(0,1,0.01), predict(AncientMod, newdata = seq(0,1,0.01)), lwd = 3, col = "blue")
AncientModAdj <- loess(joinSampGenomeMed$Ancient ~ GC.adjMed)
lines(seq(0,1,0.01), predict(AncientModAdj, newdata = seq(0,1,0.01)), lwd = 3, col = "red")
legend("topright", legend = c("Adjusted", "Raw"), fill = c(2,4), title = "GC data")




#### 



#### here we actually do our plotting 
### before I get onto this I need to correct 






genome_type <- "intron"
if(genome_type == "intergenic"){
  lenChoice = 100000
}else if (genome_type == "intron"){
  lenChoice = 40000
}

TEset <- get(paste("TEs_", genome_type,sep=""))
binnedGenomeRep <- get(paste(genome_type , "JoinedFam", sep = ""))
binnedGenomeChr <- get(paste(genome_type , "_chromatinSamp", sep = ""))

rChromatin <- covCalcPlot5prime3prime(lenChoice=lenChoice,
                             repChoice="R", 
                             repBins=binnedGenomeChr,
                             repList=`H1-hESC`, 
                             refgene=refgene, 
                             type= genome_type, 
                             repType = "chromatin")


repChoice <- names(joinRep)
for(i in repChoice){
  assign(i, covCalcPlot5prime3prime(lenChoice=lenChoice,
                                  repChoice=i, 
                                  repBins=binnedGenomeRep,
                                  repList=joinRep, 
                                  refgene=refgene, 
                                  type= genome_type,
                                  repType = "repeats")
  )
}

sqrtCoordinates <- -(sqrt(1:(lenChoice + 1))[(lenChoice + 1):1])
Coordinates <- -((1:(lenChoice + 1))[(lenChoice + 1):1])

# our normalization facotrs 
meanRes <- meanRoller(k = 10001, primeGC3 = TEset$prime3gc, primeGC5 = TEset$prime5gc, primeBF3 = TEset$baseFreq3prime, primeBF5 = TEset$baseFreq5prime)
meanRes.pre_3 <- predict(RgcLoMed,rChromatin$rawRepCov3/rChromatin$baseFreq3prime)
meanRes.adj_3 <- (meanRes$prime3gc - meanRes.pre_3) + mean(joinSampGenomeMed$GC)
meanRes.pre_5 <- predict(RgcLoMed,rChromatin$rawRepCov5/rChromatin$baseFreq5prime)
meanRes.adj_5 <- (meanRes$prime5gc - meanRes.pre_5) + mean(joinSampGenomeMed$GC)

layout(matrix(c(1,2), nrow = 1))
par(mar=c(5,5,5,0))
plot(sqrtCoordinates,meanRes$prime5gc, type = "l", ylab = "GC fraction", ylim = c(.3,.6), main = "upstream", xaxt = "n", xlab = "")
lines(sqrtCoordinates,meanRes.pre_5, col = 2)
lines(sqrtCoordinates,meanRes.adj_5, col = 4)
par(new=TRUE)
plot(sqrtCoordinates,rChromatin$rawRepCov5/TEset$baseFreq5prime, type = "l", ylab = "", xlab = "", yaxt = "n",xaxt ="n" ,col= 3, ylim = c(.3,1))
axis(1,at = seq(0,-sqrt(lenChoice), by = -40), labels = seq(0,sqrt(lenChoice), by = 40)^2)
legend("bottomleft", legend = c("raw GC","predicted GC" ,"adjusted GC", "heterochromatin"), fill = c(1,2,4,3), title = "GC data")

par(mar=c(5,0,5,5))
plot(sqrt(1:(lenChoice+1)),meanRes$prime3gc, type = "l",ylim = c(.3,.6) , main = "3 prime", xaxt = "n", xlab = "", yaxt = "n")
lines(sqrt(1:(lenChoice+1)),meanRes.pre_3, col = 2)
lines(sqrt(1:(lenChoice+1)),meanRes.adj_3, col = 4)
par(new=TRUE)
plot(sqrt(1:(lenChoice+1)),rChromatin$rawRepCov3/TEset$baseFreq3prime, type = "l", ylab = "", xlab = "", yaxt = "n", xaxt = "n", col= 3, ylim = c(.3,1))
axis(4,at = seq(0,1,.1), labels = seq(0,1,.1))
mtext(text = "repressive chromatin fraction",side = 4, outer = F,line = 3)

axis(1,at = seq(0,sqrt(lenChoice), by = 40), labels = seq(0,sqrt(lenChoice), by = 40)^2)

layout(1)

#### lets start plotting and adjusting
greys <- grey(level = seq(.9,0,-.1))
prime5BaseFreq <- data.frame(baseFreq =TEset$baseFreq5prime, 
                             group = cut(TEset$baseFreq5prime, breaks = length(greys), labels = as.character(1:length(greys))),
                             position = sqrtCoordinates)
prime3BaseFreq <- data.frame(baseFreq =TEset$baseFreq3prime, 
                             group = cut(TEset$baseFreq3prime, breaks = length(greys), labels = as.character(1:length(greys))),
                             position = sqrtCoordinates) 



ylims <- c(0,.2)

#### the loop will start around here
repChoice <- names(joinRep)

  pdf(file = "~/Desktop/Domain_manuscript/plots/geneRep/intergenic.pdf",width = 9,height = 5,onefile = T)
    layout(matrix(nrow = 4, c(1,3,5,7,2,4,6,8)))

for(te in 1:length(repChoice)){
  
  
  # i can probably change the MOD and gc rate depending on wheather i want to adjust chromatin a certain way
  element <- get(repChoice[te])
  Mod <- get(paste(repChoice[te],"ModAdj", sep =""))
  gcRate5 <- meanRes.adj_5
  gcRate3 <- meanRes.adj_3
  
  
  par(mar=c(0,5,0,1))
  plot(sqrtCoordinates,(element$rawRepCov5/TEset$baseFreq5prime), type = "n", ylab = paste(repChoice[te], "fraction"), xlab = "", ylim = ylims, xaxt = "n" , main = "", yaxt = "n")
  grid()
  axis(2,at = seq(.05,.15,by=.05), labels = seq(.05,.15,by=.05),col =  1, cex.axis = 1.5)

  for(c in 1:length(greys)){
    lines(sqrtCoordinates[prime5BaseFreq$group == c], (element$rawRepCov5/TEset$baseFreq5prime)[prime5BaseFreq$group == c], col = greys[c])
  }
  lines(sqrtCoordinates,predict(Mod, gcRate5) + (2 * sqrt(TEset$baseFreq5prime * predict(Mod, gcRate5) *  (1 - predict(Mod, gcRate5))))/TEset$baseFreq5prime, col = 2)
  lines(sqrtCoordinates,predict(Mod, gcRate5) - (2 * sqrt(TEset$baseFreq5prime * predict(Mod, gcRate5) *  (1 - predict(Mod, gcRate5))))/TEset$baseFreq5prime, col = 2)
  par(new=TRUE)
  plot(sqrtCoordinates,rChromatin$rawRepCov5/TEset$baseFreq5prime, type = "l", ylab = "", xlab = "", yaxt = "n", xaxt = "n", col= "darkblue", ylim = c(.3,1))
  
#  axis(1,at = seq(0,-sqrt(lenChoice), by = -40), labels = seq(0,sqrt(lenChoice), by = 40)^2)
  
  
  par(mar=c(0,1,0,5))
  plot(sqrt(1:(lenChoice+1)),(element$rawRepCov3/TEset$baseFreq3prime), type = "n", ylab = "", xlab = "", ylim = ylims, xaxt = "n" , yaxt = "n", main = "")
  grid()
  for(c in 1:length(greys)){
    lines(sqrt(1:(lenChoice+1))[prime3BaseFreq$group == c], (element$rawRepCov3/TEset$baseFreq3prime)[prime3BaseFreq$group == c], col = greys[c])
  }
  lines(sqrt(1:(lenChoice+1)),predict(Mod, gcRate3) + (2 * sqrt(TEset$baseFreq3prime * predict(Mod, gcRate3) *  (1 - predict(Mod, gcRate3))))/TEset$baseFreq3prime, col = 2)
  lines(sqrt(1:(lenChoice+1)),predict(Mod, gcRate3) - (2 * sqrt(TEset$baseFreq3prime * predict(Mod, gcRate3) *  (1 - predict(Mod, gcRate3))))/TEset$baseFreq3prime, col = 2)
  par(new=TRUE)
  plot(sqrt(1:(lenChoice+1)),rChromatin$rawRepCov3/TEset$baseFreq3prime, type = "l", ylab = "", xlab = "", yaxt = "n", xaxt = "n", col=  "darkblue", ylim = c(.3,1))
  axis(4,at = seq(.475,1,by=.35), labels = seq(.475,1,by=.35),col =  "darkblue", cex.axis = 1.5, col.axis="darkblue")
 # mtext(text = "repressive chromatin fraction",side = 4, outer = F,line = 3)
  
#  axis(1,at = seq(0,sqrt(lenChoice), by = 40), labels = seq(0,sqrt(lenChoice), by = 40)^2)
  
  #### two different effects

}
    par(mar=c(5,5,0,1))
    plot(sqrtCoordinates,(element$rawRepCov5/TEset$baseFreq5prime), type = "n", ylab = paste(repChoice[te], "fraction"), xlab = "", ylim = ylims, xaxt = "n" , main = "", yaxt = "n")
    axis(1,at = seq(0,-sqrt(lenChoice), by = -40), labels = seq(0,sqrt(lenChoice), by = 40)^2, cex.axis=1.5)
  par(mar=c(5,1,0,5))
  plot(sqrt(1:(lenChoice+1)),(element$rawRepCov3/TEset$baseFreq3prime), type = "n", ylab = "", xlab = "", ylim = ylims, xaxt = "n" , yaxt = "n", main = "")
  axis(1,at = seq(0,sqrt(lenChoice), by = 40), labels = seq(0,sqrt(lenChoice), by = 40)^2, cex.axis=1.5)
  dev.off()

## Trying to get the range thing to work
  
  greys
  max(c(max(TEset$baseFreq5prime),max(TEset$baseFreq3prime)))
  min(c(min(TEset$baseFreq5prime),min(TEset$baseFreq3prime)))
  mean(c(TEset$baseFreq5prime, TEset$baseFreq3prime))
  pdf("~/Desktop/Domain_manuscript/plots/geneRep/intergenic_legend.pdf")
  layout(1)
  par(mar=c(5,5,5,27))
  image(x = t(matrix(1:10)),col = greys, xaxt = "n", yaxt = "n")
  axis(side = 4, at = seq(from = 0, to = 1,length.out = 10), seq(1000, 10000, 1000))
  dev.off()
  # this is with adjusted heterochromatin

  
  
  #### lets work out which thinks we should save
  
  
  
  
  
  
  
  

for(te in 1:length(repChoice)){

element <- get(repChoice[te])
Mod <- get(paste(repChoice[te],"Mod", sep =""))
gcRate5 <- meanRes$prime5gc
gcRate3 <- meanRes$prime3gc


layout(matrix(nrow = 1, c(1,2)))
par(mar=c(5,5,5,0))
plot(sqrtCoordinates,(element$rawRepCov5/TEset$baseFreq5prime), type = "n", ylab = paste(repChoice[te], "fraction"), xlab = "", ylim = ylims, xaxt = "n" , main = "upsteam")
grid()
for(c in 1:length(greys)){
  lines(sqrtCoordinates[prime5BaseFreq$group == c], (element$rawRepCov5/TEset$baseFreq5prime)[prime5BaseFreq$group == c], col = greys[c])
}
lines(sqrtCoordinates,predict(Mod, gcRate5) + (2 * sqrt(TEset$baseFreq5prime * predict(Mod, gcRate5) *  (1 - predict(Mod, gcRate5))))/TEset$baseFreq5prime, col = 2)
lines(sqrtCoordinates,predict(Mod, gcRate5) - (2 * sqrt(TEset$baseFreq5prime * predict(Mod, gcRate5) *  (1 - predict(Mod, gcRate5))))/TEset$baseFreq5prime, col = 2)
par(new=TRUE)
#plot(sqrtCoordinates,((rChromatin$rawRepCov5/TEset$baseFreq5prime) - predict(rMod, gcRate5)) + mean(joinSampGenomeMed$R), type = "l", ylab = "", xlab = "", yaxt = "n", xaxt = "n", col= "darkblue", ylim = c(.3,1))
plot(sqrtCoordinates,(rChromatin$rawRepCov5/TEset$baseFreq5prime), type = "l", ylab = "", xlab = "", yaxt = "n", xaxt = "n", col= "darkblue", ylim = c(.3,1))

axis(1,at = seq(0,-350, by = -50), labels = seq(0,350, by = 50)^2)


par(mar=c(5,0,5,5))
plot(sqrt(1:(lenChoice+1)),(element$rawRepCov3/TEset$baseFreq3prime), type = "n", ylab = "", xlab = "", ylim = ylims, xaxt = "n" , yaxt = "n", main = "downsteam")
grid()
for(c in 1:length(greys)){
  lines(sqrt(1:(lenChoice+1))[prime3BaseFreq$group == c], (element$rawRepCov3/TEset$baseFreq3prime)[prime3BaseFreq$group == c], col = greys[c])
}
lines(sqrt(1:(lenChoice+1)),predict(Mod, gcRate3) + (2 * sqrt(TEset$baseFreq3prime * predict(Mod, gcRate3) *  (1 - predict(Mod, gcRate3))))/TEset$baseFreq3prime, col = 2)
lines(sqrt(1:(lenChoice+1)),predict(Mod, gcRate3) - (2 * sqrt(TEset$baseFreq3prime * predict(Mod, gcRate3) *  (1 - predict(Mod, gcRate3))))/TEset$baseFreq3prime, col = 2)
par(new=TRUE)
#plot(sqrt(1:(lenChoice+1)),((rChromatin$rawRepCov3/TEset$baseFreq3prime) - predict(rMod, gcRate3)) + mean(joinSampGenomeMed$R), type = "l", ylab = "", xlab = "", yaxt = "n", xaxt = "n", col=  "darkblue", ylim = c(.3,1))
plot(sqrt(1:(lenChoice+1)),((rChromatin$rawRepCov3/TEset$baseFreq3prime)), type = "l", ylab = "", xlab = "", yaxt = "n", xaxt = "n", col=  "darkblue", ylim = c(.3,1))

axis(4,at = seq(.3,1,length.out = 5), labels = seq(.3,1,length.out = 5),col =  "darkblue")
mtext(text = "repressive chromatin fraction",side = 4, outer = F,line = 3)

axis(1,at = seq(0,350, by = 50), labels = seq(0,350, by = 50)^2)

}

### we use the GC rate to get the expected level of things 
### we redict the mean rate based on GC, 




plot(intron_chromatinSamp$Known,(intron_chromatinSamp$Known- intron_chromatinSamp$R)/ intron_chromatinSamp$Known)

plot((intergenic_chromatinSamp$Known),(intergenic_chromatinSamp$Known- intergenic_chromatinSamp$R)/ intergenic_chromatinSamp$Known, pch = 16, cex = .2)
points((intron_chromatinSamp$Known),(intron_chromatinSamp$Known- intron_chromatinSamp$R)/ intron_chromatinSamp$Known, col = 2, pch = 16, cex = .2)

mod.gap.intorn <- loess(((intron_chromatinSamp$Known- intron_chromatinSamp$R)/ intron_chromatinSamp$Known) ~ (intron_chromatinSamp$Known))
lines(seq(1,200000,.5),predict(mod.gap.intorn,seq(1,200000,.5) ), col = 2)
mod.gap.intergenic <- loess(((intergenic_chromatinSamp$Known- intergenic_chromatinSamp$R)/ intergenic_chromatinSamp$Known) ~ (intergenic_chromatinSamp$Known))
lines(seq(6,14,.5),predict(mod.gap.intergenic,seq(6,14,.5) ), col = 1)











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




