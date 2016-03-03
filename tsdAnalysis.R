#### 


## what are we going to do with TSD analysis 

## havnt really joined the repeats yet


### that is ok I just need to make a pipeline i can later update to feed data into 

# get some sort of gff output or somehting 

####### Important questions this will answear. 
##### newest repeats, GC bias , proximit to chromatin boudaries. 

#### what does a chromatin boundary look like?
### ttr might be interesting too. 

## gc will take a while to figuer out
## lets start bu looking at gene proximity


### so where do we get non tsd repeats from? 

### just get the rm gff file or move everyhting over to the big computers



rm(list = ls())
setwd("~/Desktop/Domain_manuscript/")
library(GenomicRanges)
library(rtracklayer)


reps <- read.table(file = "~/Desktop/quilt/repfile/hg19.rm.gff",
                   colClasses = c("character", "character","character", 
                                  "integer", "integer", "integer", 
                                  "character", "character", "character", 
                                  "character", "character", "integer", 
                                  "integer", "integer" ))
colnames(reps) <- c("chr", "program", "interval", "start", "end", "length", "strand", "dot", "featType", "repName", "repType", "repStart", "repEnd", "repLeft")

tsd <- read.table(file = "~/Desktop/quilt/tsdFinder/hg19_chrKnown.tsd.txt", skip = 1)
colnames(tsd) <- c("chr", "strat", "end", "leftStart", 
                   "rightStart", "tsdLength", "repFam", 
                   "repType", "repStart", "repEnd", "repLeft",
                   "repFam2", "repType2" , "tsdSeq", "strand")

unique(tsd$chr)


L1hist <- hist(tsd$V6[tsd$V8 == "LINE/L1"])
Aluhist <- hist(tsd$V6[tsd$V8 == "SINE/Alu"])
plot(L1hist$mids, L1hist$counts, type = "l")
lines(Aluhist$mids, Aluhist$counts, col = 2)

hist(tsd$V6)

plot(tsd[tsd$repType == "LINE/L1", "tsdLength"], log10(tsd[tsd$repType == "LINE/L1", "repEnd"]), cex = .3, pch = 16)

stripchart((tsd[tsd$repType == "LINE/L1", "repEnd"]) ~ as.factor(tsd[tsd$repType == "LINE/L1", "tsdLength"]), 
           method = "jitter", pch = 16, vertical = TRUE, jitter = .3)
abline(h=(6000), col = 2)
abline(h=(8000), col = 2)


colours <- rep(0, nrow(tsd[tsd$repType == "LINE/L1",]))
colours[grep("L1ME", tsd$repFam[tsd$repType == "LINE/L1"])] =1
colours[grep("L1MD", tsd$repFam[tsd$repType == "LINE/L1"])] =2
colours[grep("L1MC", tsd$repFam[tsd$repType == "LINE/L1"])] =3
colours[grep("L1MB", tsd$repFam[tsd$repType == "LINE/L1"])] =4
colours[grep("L1MA", tsd$repFam[tsd$repType == "LINE/L1"])] =5
colours[grep("HAL1", tsd$repFam[tsd$repType == "LINE/L1"])] =0
colours[grep("L1PB", tsd$repFam[tsd$repType == "LINE/L1"])] =6
colours[grep("L1PA", tsd$repFam[tsd$repType == "LINE/L1"])] =7
colours[grep("L1HS", tsd$repFam[tsd$repType == "LINE/L1"])] =8



cexL1 <- rep(.3, length(colours))
cexL1[colours == 8] = 1

plot(tsd[tsd$repType == "LINE/L1", "tsdLength"][colours > 0] + runif(n = length(colours[colours > 0]), min = -.4, max = .4), 
     (tsd[tsd$repType == "LINE/L1", "repEnd"])[colours > 0], 
     cex = cexL1[colours > 0], pch = 16, col = colours[colours > 0],
     xlab = "TSD size",
     ylab = "3' coordinate")
legend("bottomright", legend = c("L1ME", "L1MD", "L1MC", "L1MB", "L1MA", "L1PB", "L1PA","L1HS"), fill = c(1,2,3,4,5,6,7,8))







# the average distance to a gene for different groups we could probably put these on a contour plot 

# maybe in 3D we could look at each distribution 

# or we could colur the points that we are already looking at


#### gnene tables 

spec1 <- "Human"
genome = "hg19"


web <- paste("http://hgdownload.soe.ucsc.edu/goldenPath/", genome, "/database/refGene.txt.gz", sep = "")
con <- gzcon(url(web))
txt <- readLines(con)
refgene <- read.delim(textConnection(txt), header = FALSE)
colnames(refgene) <- c("bin", "ID", "chr", "strand", "start", "end", "txnStart", "txnEnd", "exonCount",
                       "exonStarts", "exonEnds", "score", "name2", "cdsStartStat", "cdsEndStat", "exonFrames")




# maybe remove intronic elements 
# 
tsdL1 <- tsd[tsd$repType == "LINE/L1",]
tsdL1.gr <- GRanges(seqnames = Rle(tsdL1$chr),
                    ranges = IRanges(start = tsdL1$strat,end=tsdL1$end))

refgene.gr <- GRanges(seqnames = Rle(refgene$chr), 
                      ranges = IRanges(start = refgene$start, end = refgene$end))





# we can get a curve that shows the cumulative distribution of elements within x space of gnenes

sizes <- c(0,1,5,10,50,100,500,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,20000,30000,40000,50000,100000, 200000,300000,400000,500000, 600000,1000000)
OLcount<- matrix(0,ncol = length(sizes), nrow = length(colours))
OLmat<- matrix(0,ncol = length(sizes), nrow = length(colours))

for(i in 1:length(sizes)){
  geneGrow <- reduce(refgene.gr)
  start(geneGrow) = start(geneGrow) - sizes[i]
  end(geneGrow) = end(geneGrow) + sizes[i]
  OLcount[,i] <- countOverlaps(tsdL1.gr, geneGrow)
  OLmat[,i] <- OLcount[,i] > 0
}


tsdL = 13

plot(log10(sizes), (colSums(OLmat)), type = "n", ylim = c(0,1), xlim= c(3,6))
for(i in c(5:8)){
  lines(log10(sizes), 
        (  (colSums(OLmat[colours == i & tsdL1$tsdLength >=tsdL & tsdL1$repEnd > 6000,]) - 
              sum(OLmat[colours == i & tsdL1$tsdLength >=tsdL & tsdL1$repEnd > 6000,1])) /
           (length(colours[colours == i & tsdL1$tsdLength >=tsdL & tsdL1$repEnd > 6000])- 
              sum(OLmat[colours == i & tsdL1$tsdLength >=tsdL & tsdL1$repEnd > 6000,1])) ), 
         col=i, lwd = 3)
}




for(i in 5:8){
  lines(log10(sizes), 
        (  (colSums(OLmat[colours == i & tsdL1$tsdLength <tsdL & tsdL1$repEnd > 6000,]) - 
              sum(OLmat[colours == i & tsdL1$tsdLength <tsdL & tsdL1$repEnd > 6000,1])) /
             (length(colours[colours == i & tsdL1$tsdLength <tsdL & tsdL1$repEnd > 6000])- 
                sum(OLmat[colours == i & tsdL1$tsdLength <tsdL & tsdL1$repEnd > 6000,1])) ), 
        col=i, lty = 2, lwd = 3)
}


legend("bottomright", legend = c("L1ME", "L1MD", "L1MC", "L1MB", "L1MA", "L1PB", "L1PA","L1HS"), fill = c(1,2,3,4,5,6,7,8))



# looking at tsd repeats and their distance from genes. 

# it would be interesting to compare tsd with nontsd. 


# when i get back pull these things right






head(reps)

##reps_IDs <- paste(reps$chr, reps$start, reps$end)

repSampL1 <- reps[reps$repType == "LINE/L1", ]

#sample(size = 100,x = 1:nrow(repSampL1))
repSampL1 <- repSampL1[sample(size = 100000,x = 1:nrow(repSampL1)),]

repSampL1.gr <- GRanges(seqnames = Rle(repSampL1$chr),
                        ranges = IRanges(start = repSampL1$start, end = repSampL1$end))



coloursGenome <- rep(0, nrow(repSampL1))
coloursGenome[grep("L1ME", repSampL1$repName)] =1
coloursGenome[grep("L1MD", repSampL1$repName)] =2
coloursGenome[grep("L1MC", repSampL1$repName)] =3
coloursGenome[grep("L1MB", repSampL1$repName)] =4
coloursGenome[grep("L1MA", repSampL1$repName)] =5
coloursGenome[grep("HAL1", repSampL1$repName)] =0
coloursGenome[grep("L1PB", repSampL1$repName)] =6
coloursGenome[grep("L1PA", repSampL1$repName)] =7
coloursGenome[grep("L1HS", repSampL1$repName)] =8




OLcountGenome <- matrix(0,ncol = length(sizes), nrow = nrow(repSampL1))
OLmatGenome <- matrix(0,ncol = length(sizes), nrow = nrow(repSampL1))

for(i in 1:length(sizes)){
  geneGrow <- reduce(refgene.gr)
  start(geneGrow) = start(geneGrow) - sizes[i]
  end(geneGrow) = end(geneGrow) + sizes[i]
  OLcountGenome[,i] <- countOverlaps(repSampL1.gr, geneGrow)
  OLmatGenome[,i] <- OLcountGenome[,i] > 0
}


plot(c(3,6), c(0,1), type = "n", 
     xlab = "distance from gene (log 10)", 
     ylab = "proprtion of L1s")

for(i in 5:8){
  lines(log10(sizes), 
        (  (colSums(OLmatGenome[coloursGenome == i ,]) - 
              sum(OLmatGenome[coloursGenome == i ,1])) /
             (length(coloursGenome[coloursGenome == i ])- 
                sum(OLmatGenome[coloursGenome == i ,1])) ), 
        col=i-4, lty = 2, lwd = 3)
}

for(i in c(5:8)){
  lines(log10(sizes), 
        (  (colSums(OLmat[colours == i & tsdL1$tsdLength >=tsdL & tsdL1$repEnd > 6000,]) - 
              sum(OLmat[colours == i & tsdL1$tsdLength >=tsdL & tsdL1$repEnd > 6000,1])) /
             (length(colours[colours == i & tsdL1$tsdLength >=tsdL & tsdL1$repEnd > 6000])- 
                sum(OLmat[colours == i & tsdL1$tsdLength >=tsdL & tsdL1$repEnd > 6000,1])) ), 
        col=i-4, lwd = 3)
}

legend("bottomright", legend = c( "L1MA", "L1PB", "L1PA","L1HS"), fill = c(1,2,3,4))
legend("right", legend = c( "TSD", "Random Sample"), lty = c(1,2), lwd = c(3,3))


