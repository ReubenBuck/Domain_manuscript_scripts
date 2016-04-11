##### geneAnalysis in domains 


rm(list = ls())
setwd("~/Desktop/Domain_manuscript/")
library(GenomicRanges)
library(rtracklayer)


source("Domain_manuscript_scripts/functions.R")
source("Domain_manuscript_scripts/rep_db.R")
load("R_objects/chromStateCombined")


genome = "hg19"
spec1 = "Human"

domainRanges <- read.table("Data/ConsTimingDomains", header = T)




dd <- as.character(domainRanges$domain)
ER <- (1:length(dd))[dd == "ERD"]
UT = ER-1
DT = ER+1

df <- data.frame(UT = dd[UT],ER = dd[ER], DT = dd[DT])
triplet <- ER[df$UT == "UTZ" & df$ER == "ERD" & df$DT == "DTZ"]
ranges <- unique(c(triplet[domainRanges[triplet+1,"start"] - domainRanges[triplet,"end"] == 1],
                   triplet[domainRanges[triplet,"start"] - domainRanges[triplet-1,"end"] == 1]))

cleanE <- domainRanges[ranges,]

LR <- (1:length(dd))[dd == "LRD"]
UT = LR-1
DT = LR+1

df <- data.frame(UT = dd[UT],LR = dd[LR], DT = dd[DT])
triplet <- LR[df$UT == "DTZ" & df$LR == "LRD" & df$DT == "UTZ"]
ranges <- unique(c(triplet[domainRanges[triplet+1,"start"] - domainRanges[triplet,"end"] == 1],
                   triplet[domainRanges[triplet,"start"] - domainRanges[triplet-1,"end"] == 1]))

cleanL <- domainRanges[ranges,]


cleanE.gr <- GRanges(seqnames = Rle(cleanE$chr), 
                     ranges = IRanges(start = cleanE$start, end = cleanE$end))
cleanE.gr

cleanL.gr <- GRanges(seqnames = Rle(cleanL$chr), 
                     ranges = IRanges(start = cleanL$start, end = cleanL$end))
cleanL.gr

findOverlaps(cleanL.gr,cleanE.gr)


### now we have our two genome sections
### we can start to look at the TE content within
rep <- rep_info(spec1 = spec1, genome = genome)
joinRep <- list(old_L1 = rbind(rep$L1ME, rep$L1MD, rep$L1MC, rep$L1MB),
                new_L1 = rbind(rep$L1MA, rep$L1PB, rep$L1PA, rep$L1HS),
                Alu = rbind(rep$AluJ, rep$AluS, rep$AluY),
                Ancient = rbind(rep$MIR, rep$L2)
)

web <- paste("http://hgdownload.soe.ucsc.edu/goldenPath/", genome, "/database/refGene.txt.gz", sep = "")
con <- gzcon(url(web))
txt <- readLines(con)
refgene <- read.delim(textConnection(txt), header = FALSE)
refgene.gr <- GRanges(seqnames=Rle(refgene[,3]), ranges = IRanges(start = refgene[,5], end = refgene[,6]))
refgene_gap.gr <- filterIntergenic(refgene)
intronKeep.gr <- filterIntron(refgene)


# Intergenic 

OLintergenicE <- as.matrix(findOverlaps(refgene_gap.gr, cleanE.gr, type = "within"))
InterGenicE <- refgene_gap.gr[OLintergenicE[,1]]

OLintergenicL <- as.matrix(findOverlaps(refgene_gap.gr, cleanL.gr, type = "within"))
InterGenicL <- refgene_gap.gr[OLintergenicL[,1]]

plot(density(log10(width(InterGenicE))))

findOverlaps(InterGenicL,InterGenicE)


plot(density(log10(width(InterGenicL))))
lines(density(log10(width(InterGenicE))), col = 2)


histL <- hist(log10(width(InterGenicL)), breaks = 100)
loL <- loess(histL$counts ~ histL$mids)
histE <- hist(log10(width(InterGenicE)), breaks = histL$breaks)
loE <- loess(histE$counts ~ histE$mids)
plot(seq(2,7,.01), predict(loL, seq(2,7,.01)), type = "l", ylim = c(0,200))
lines(seq(2,7,.01), predict(loE, seq(2,7,.01)), col = 2)

plot((sort(width(InterGenicE))), (length(InterGenicE):1), ylab = "number of ranges > range width",
     xlab = "range width (log10 bp)", type = "l", main = "Intergenic position frequency")

lines((sort(width(InterGenicL))), (length(InterGenicL):1), col = 2)
legend("topright", legend = c("Early", "Late"), fill = c(1,2), title = "Replication Domain")






# Intron

OLintronE <- as.matrix(findOverlaps(intronKeep.gr, cleanE.gr, type = "within"))
IntronE <- intronKeep.gr[OLintronE[,1]]

OLintronL <- as.matrix(findOverlaps(intronKeep.gr, cleanL.gr, type = "within"))
IntronL <- intronKeep.gr[OLintronL[,1]]

plot(density(log10(width(IntronE))))

findOverlaps(IntronL,IntronE)


plot(density(log10(width(IntronL))), col = 2)
lines(density(log10(width(IntronE))), col = 1)
legend("topright", legend = c("Early", "Late"), fill = c(1,2), title = "Replication Domain")


histL <- hist(log10(width(IntronL)), breaks = 100)
loL <- loess(histL$counts ~ histL$mids)
histE <- hist(log10(width(IntronE)), breaks = histL$breaks)
loE <- loess(histE$counts ~ histE$mids)
plot(seq(2,7,.01), predict(loL, seq(2,7,.01)), type = "l", ylim = c(0,200))
lines(seq(2,7,.01), predict(loE, seq(2,7,.01)), col = 2)

plot(log10(sort(width(IntronE))), (length(IntronE):1), ylab = "number of ranges > range width",
     xlab = "range width (log10 bp)", type = "l", main = "Intergenic position frequency")

lines(log10(sort(width(IntronL))), (length(IntronL):1), col = 2)
legend("topright", legend = c("Early", "Late"), fill = c(1,2), title = "Replication Domain")


# so looking at the 5' and 3' positions of the genome we can get expected region estimates. 
# the fold decrease in bases called as 




# to explain the distribution of L1s in late replicating regions, 
# we find that the majority of bases are too close to genes. 

# our expected number of bases will be this curve multiplied by the relationship with position summed 

# so we get this curve that says how many spots are available this far from the edge
# so we get that curve 

# this is strange, it seems that the LRD has the same size distribution as the ERD in terms of intergenic regions.
# the idea is that there are less genes and that the reduced amount of genes allows for the insertion of TEs.






# we'd also expect that the numbers of genes should correlate with the size of the region in E domains, not so much in L domains




plot(log10(width(InterGenicE)[OLintergenicE[,2]]), (width(refgene_gap.gr)[OLintergenicE[,1]]))


points(log10(width(InterGenicL)[OLintergenicL[,2]]), (width(refgene_gap.gr)[OLintergenicL[,1]]), col = 2)


Lsize <- loess((width(refgene_gap.gr)[OLintergenicL[,1]]) ~ log10(width(InterGenicL)[OLintergenicL[,2]]))
Lpred <- predict(Lsize,newdata = seq(2,7,.01), se = T)
Lyline <- c(Lpred$fit+(2*Lpred$se.fit) , rev(Lpred$fit-(2*Lpred$se.fit)))


Esize <- loess((width(refgene_gap.gr)[OLintergenicE[,1]]) ~ log10(width(InterGenicL)[OLintergenicE[,2]]))
Epred <- predict(Esize, newdata = seq(2,7,.01), se = T)
Eyline <- c(Epred$fit+(2*Epred$se.fit) , rev(Epred$fit-(2*Epred$se.fit)))


plot(seq(2,7,.01), Epred$fit, col = 2, type = "l", lwd = 3, ylim = c(0,300000))

polygon(x = c(seq(2,7,.01),seq(7,2,-.01))[!is.na(Eyline)], y = Eyline[!is.na(Eyline)], density = 20, col = 2,border = 0)

lines(seq(2,7,.01), Lpred$fit, col = 1, type = "l", lwd = 3)
polygon(x = c(seq(2,7,.01),seq(7,2,-.01))[!is.na(Lyline)], y = Lyline[!is.na(Lyline)], density = 20, col = 1,border = 0)




Ecut <- cut(x = log10(width(InterGenicE)[OLintergenicE[,2]]),breaks = seq(2,7,.3))
boxplot((width(refgene_gap.gr)[OLintergenicE[,1]]) ~ Ecut, outline = F, ylim = c(0,1500000))


Lcut <- cut(x = log10(width(InterGenicL)[OLintergenicL[,2]]), breaks = seq(2,7,.3))
boxplot((width(refgene_gap.gr)[OLintergenicL[,1]]) ~ Lcut, outline = F,add = T)



Ecut <- cut(x = log10(width(InterGenicE)[OLintergenicE[,2]]),breaks = seq(2,7,.1))
Lcut <- cut(x = log10(width(InterGenicL)[OLintergenicL[,2]]), breaks = seq(2,7,.1))

dfE <- data.frame(IntervalSize = (width(refgene_gap.gr)[OLintergenicE[,1]]), cut = Ecut, group = "Early")
dfL <- data.frame(IntervalSize = (width(refgene_gap.gr)[OLintergenicL[,1]]), cut = Lcut, group = "Late")


both <- rbind(dfE, dfL)


boxplot((both$IntervalSize) ~ both$group + both$cut, col = c(0,2), outline = F, varwidth = F, range = 1)



head(both)

summary(both$cut[both$group == "Early"])
# how do we get a random model?
# the Ranges, they are going to vary 

# we look at the number of bases in inntergenic regions of different sizes 
# the number of bases 

# is a random model the right thing?
# the idea is that 




plot(log10(width(InterGenicE)[OLintergenicE[,2]]), (width(refgene_gap.gr)[OLintergenicE[,1]]))


points(log10(width(InterGenicL)[OLintergenicL[,2]]), (width(refgene_gap.gr)[OLintergenicL[,1]]), col = 2)



interE <- summary(as.factor(OLintergenicE[,2]), maxsum = length(unique(OLintergenicE[,2])))

plot(sqrt(width(cleanE.gr)[as.integer(names(interE))]), log10((interE) / (width(cleanE.gr)[as.integer(names(interE))])), ylim = c(-7,-4))


interL <- summary(as.factor(OLintergenicL[,2]), maxsum = length(unique(OLintergenicL[,2])))

points(sqrt(width(cleanL.gr)[as.integer(names(interL))]), log10((interL)/(width(cleanL.gr)[as.integer(names(interL))])), col = 2)


# we choose a gene or exon and see what happens 

# this will go to show that the size around these genes is somewhat regulated 
# also the genes in these domains are regulated and generealy more conserved across species. 

LrefOL <- as.matrix(findOverlaps(reduce(refgene_gap.gr), cleanL.gr, type = "within"))
interL <- summary(as.factor(LrefOL[,2]), maxsum = length(unique(LrefOL[,2])))

ErefOL <- as.matrix(findOverlaps(reduce(refgene_gap.gr), cleanE.gr, type = "within"))
interE <- summary(as.factor(ErefOL[,2]), maxsum = length(unique(ErefOL[,2])))


x <- width(cleanE.gr)[as.integer(names(interE))]
y = interE
y = y[order(x)]
x = sort(x)

n = sum(y)
p = (x/sum(x))

m <- n*p
SD = sqrt(n*p*(1-p))


plot((x), (y/x))
lines((x), (m/x))
lines((x), ((m +(3*SD) )/x ))
lines((x), ((m - (3*SD) )/x ))

df2 <- data.frame((width(reduce(refgene.gr))[ErefOL[,1]]),as.factor(ErefOL[,2])) 
agg <- aggregate(x = as.integer(df2[,1]), by = list(df2[,2]), FUN = sum)

x <- width(cleanL.gr)[as.integer(names(interL))] - 
  

y = interL
y = y[order(x)]
x = sort(x)

n = sum(y)
p = (x/sum(x))

m <- n*p
SD = sqrt(n*p*(1-p))


plot((x), (y/x))
lines((x), (m/x))
lines((x), ((m +(3*SD) )/x ))
lines((x), ((m - (3*SD) )/x ))


# are we better off binning 
# so we are not getting accurate measurments becasue we are not considering the ranges with no genes
# should we also be looking at the 


hist(log10(width(cleanE.gr)), breaks = 20)

#### not looking for weather the number of genes correlates to the size of the region 
### trying to test for constraint over the gene density 








