# now we can begin to combine the anlysis of both kinds of regions
# this will be our fig 2

# where we look at the relationship between TE accumulation and genes

# the idea is that TEs accumulate according to our model around genes

# on both fractions theres constraint pertaining to the amount but the constraint on position is different


setwd("~/Desktop/topological_domains/")

rm(list = ls())


library(GenomicRanges)
library(rtracklayer)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)

spec1 <- "Human"
genome = "hg19"


source(file="~/Desktop/element_curves/element_curves_scripts/functions.R")
source(file="~/Desktop/Domain_manuscript/Domain_manuscript_scripts/rep_db.R")


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



plot(density(log10(intron_reps$counts$Known)), col = 2)
lines(density(log10(intergenic_reps$counts$Known)), col = 1)

length(intron_reps$counts$Known[intron_reps$counts$Known > 60000])
length(intergenic_reps$counts$Known[intergenic_reps$counts$Known > 60000])



log10(60000)




# region comparison

minS <- 30000
maxS <- NULL
minRepCov <- 5000
maxRepCov <- 7000
lenChoice = 100000
repChoice <- "L1PA"
TEintergenic <- covCalcPlot5prime3prime(lenChoice=lenChoice,repChoice=repChoice, repBins=intergenic_reps$counts,repList=rep, refgene=refgene, type= "intergenic",minBinSize=minS, maxBinSize=maxS,minRepCov=minRepCov, maxRepCov=maxRepCov)
TEintron <- covCalcPlot5prime3prime(lenChoice=lenChoice,repChoice=repChoice, repBins=intron_reps$counts,repList=rep, refgene=refgene, type= "intron", minBinSize=minS, maxBinSize=maxS,minRepCov=minRepCov, maxRepCov=maxRepCov)

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


repChoice = "L1PB"
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


repChoice <- "Known"
intronGCset <- IntronGC[intron_reps$rates[,repChoice] > 0]
intronTEset <- intron_reps$rates[intron_reps$rates[,repChoice] > 0,repChoice]
intergenicGCset <- IntergenicGC[intergenic_reps$rates[,repChoice] > 0]
intergenicTEset <- intergenic_reps$rates[intergenic_reps$rates[,repChoice] > 0,repChoice]

plot((intronGCset),
      log10(intronTEset),
      pch = 16, cex = .2, col = 2 , main = repChoice,
     ylab = "TE coverage rate (log10 per 10000bp)",
     xlab = "GC content"
      )
points((intergenicGCset),
     log10(intergenicTEset), 
     pch = 16, cex = .2)
interMod <- lm(log10(intergenicTEset) ~ (intergenicGCset))
intronMod <- lm(log10(intronTEset) ~ (intronGCset))

abline(intronMod, col = 2, lwd = 2)
abline(interMod, col = 1, lwd = 2)
legend("bottomright", legend=c("intergenic", "intron"), fill = c(1,2))

legend("topright", bty="n", legend=c(paste("intergenic r =", format(sqrt(summary(interMod)$r.squared), digits=3)), 
                                     paste("intron r =", format(sqrt(summary(intronMod)$r.squared), digits=3))
                                     )
       )








