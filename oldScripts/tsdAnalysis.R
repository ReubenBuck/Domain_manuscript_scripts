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
source("Domain_manuscript_scripts/functions.R")

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
colnames(tsd) <- c("chr", "start", "end", "leftStart", 
                   "rightStart", "tsdLength", "repFam", 
                   "repType", "repStart", "repEnd", "repLeft",
                   "repFam2", "repType2" , "tsdSeq", "strand")

unique(tsd$chr)


L1hist <- hist(tsd$V6[tsd$V8 == "LINE/L1"])
Aluhist <- hist(tsd$V6[tsd$V8 == "SINE/Alu"])
plot(L1hist$mids, L1hist$counts, type = "l")
lines(Aluhist$mids, Aluhist$counts, col = 2)

hist(tsd$V6)

Rtype = "LINE/L1"

plot(tsd[tsd$repType == Rtype, "tsdLength"], log10(tsd[tsd$repType == Rtype, "repEnd"]), cex = .3, pch = 16)

stripchart((tsd[tsd$repType == Rtype, "repEnd"]) ~ as.factor(tsd[tsd$repType == Rtype, "tsdLength"]), 
           method = "jitter", pch = 16, vertical = TRUE, jitter = .3)
abline(h=(6000), col = 2)
abline(h=(8000), col = 2)

if(Rtype == "LINE/L1"){
  colours <- rep(0, nrow(tsd[tsd$repType == Rtype,]))
  colours[grep("L1ME", tsd$repFam[tsd$repType == Rtype])] =1
  colours[grep("L1MD", tsd$repFam[tsd$repType == Rtype])] =2
  colours[grep("L1MC", tsd$repFam[tsd$repType == Rtype])] =3
  colours[grep("L1MB", tsd$repFam[tsd$repType == Rtype])] =4
  colours[grep("L1MA", tsd$repFam[tsd$repType == Rtype])] =5
  colours[grep("HAL1", tsd$repFam[tsd$repType == Rtype])] =0
  colours[grep("L1PB", tsd$repFam[tsd$repType == Rtype])] =6
  colours[grep("L1PA", tsd$repFam[tsd$repType == Rtype])] =7
  colours[grep("L1HS", tsd$repFam[tsd$repType == Rtype])] =8
}


if(Rtype == "SINE/Alu"){
  colours <- rep(0, nrow(tsd[tsd$repType == Rtype,]))
  colours[grep("AluJ", tsd$repFam[tsd$repType == Rtype])] =1
  colours[grep("AluS", tsd$repFam[tsd$repType == Rtype])] =2
  colours[grep("AluY", tsd$repFam[tsd$repType == Rtype])] =3
}


cexL1 <- rep(.3, length(colours))
cexL1[colours == 8] = 1

plot(tsd[tsd$repType == Rtype, "repEnd"][colours > 0] + runif(n = length(colours[colours > 0]), min = -.4, max = .4), 
     (tsd[tsd$repType == Rtype, "repEnd"])[colours > 0], 
     cex = cexL1[colours > 0], pch = 16, col = colours[colours > 0],
     xlab = "TSD size",
     ylab = "3' coordinate")


if(Rtype == "LINE/L1"){
  legend("bottomright", legend = c("L1ME", "L1MD", "L1MC", "L1MB", "L1MA", "L1PB", "L1PA","L1HS"), fill = c(1,2,3,4,5,6,7,8))
}
if(Rtype == "SINE/Alu"){
  legend("bottomright", legend = c("AluJ", "AluS", "AluY"), fill = c(1,2,3))
}






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
tsdL1 <- tsd[tsd$repType == Rtype,]
tsdL1.gr <- GRanges(seqnames = Rle(tsdL1$chr),
                    ranges = IRanges(start = tsdL1$start,end=tsdL1$end))

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
prime3endMin = 290


plot(log10(sizes), (colSums(OLmat)), type = "n", ylim = c(0,1), xlim= c(3,6))
for(i in c(1:max(colours))){
  lines(log10(sizes), 
        (  (colSums(OLmat[colours == i & tsdL1$tsdLength >=tsdL & tsdL1$repEnd > prime3endMin,]) - 
              sum(OLmat[colours == i & tsdL1$tsdLength >=tsdL & tsdL1$repEnd > prime3endMin,1])) /
           (length(colours[colours == i & tsdL1$tsdLength >=tsdL & tsdL1$repEnd > prime3endMin])- 
              sum(OLmat[colours == i & tsdL1$tsdLength >=tsdL & tsdL1$repEnd > prime3endMin,1])) ), 
         col=i, lwd = 3)
}




for(i in 1:max(colours)){
  lines(log10(sizes), 
        (  (colSums(OLmat[colours == i & tsdL1$tsdLength <tsdL & tsdL1$repEnd > prime3endMin,]) - 
              sum(OLmat[colours == i & tsdL1$tsdLength <tsdL & tsdL1$repEnd > prime3endMin,1])) /
             (length(colours[colours == i & tsdL1$tsdLength <tsdL & tsdL1$repEnd > prime3endMin])- 
                sum(OLmat[colours == i & tsdL1$tsdLength <tsdL & tsdL1$repEnd > prime3endMin,1])) ), 
        col=i, lty = 2, lwd = 3)
}



if(Rtype == "LINE/L1"){
  legend("bottomright", legend = c("L1ME", "L1MD", "L1MC", "L1MB", "L1MA", "L1PB", "L1PA","L1HS"), fill = c(1,2,3,4,5,6,7,8))
}
if(Rtype == "SINE/Alu"){
  legend("bottomright", legend = c("AluJ", "AluS", "AluY"), fill = c(1,2,3))
}


# looking at tsd repeats and their distance from genes. 

# it would be interesting to compare tsd with nontsd. 


# when i get back pull these things right






head(reps)

##reps_IDs <- paste(reps$chr, reps$start, reps$end)

repSampL1 <- reps[reps$repType == Rtype, ]

#sample(size = 100,x = 1:nrow(repSampL1))
repSampL1 <- repSampL1[sample(size = 100000,x = 1:nrow(repSampL1)),]

repSampL1.gr <- GRanges(seqnames = Rle(repSampL1$chr),
                        ranges = IRanges(start = repSampL1$start, end = repSampL1$end))


if(Rtype == "LINE/L1"){
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
}

if(Rtype == "SINE/Alu"){
  coloursGenome <- rep(0, nrow(repSampL1))
  coloursGenome[grep("AluJ", repSampL1$repName)] =1
  coloursGenome[grep("AluS", repSampL1$repName)] =2
  coloursGenome[grep("AluY", repSampL1$repName)] =3
}
  



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
     ylab = "cumulative proportion of TEs")

for(i in 1:max(coloursGenome)){
  lines(log10(sizes), 
        (  (colSums(OLmatGenome[coloursGenome == i ,]) - 
              sum(OLmatGenome[coloursGenome == i ,1])) /
             (length(coloursGenome[coloursGenome == i ])- 
                sum(OLmatGenome[coloursGenome == i ,1])) ), 
        col=i, lty = 2, lwd = 3)
}

for(i in 1:max(coloursGenome)){
  lines(log10(sizes), 
        (  (colSums(OLmat[colours == i & tsdL1$tsdLength >=tsdL & tsdL1$repEnd > prime3endMin,]) - 
              sum(OLmat[colours == i & tsdL1$tsdLength >=tsdL & tsdL1$repEnd > prime3endMin,1])) /
             (length(colours[colours == i & tsdL1$tsdLength >=tsdL & tsdL1$repEnd > prime3endMin])- 
                sum(OLmat[colours == i & tsdL1$tsdLength >=tsdL & tsdL1$repEnd > prime3endMin,1])) ), 
        col=i, lwd = 3)
}

if(Rtype == "LINE/L1"){
  legend("bottomright", legend = c( "L1MA", "L1PB", "L1PA","L1HS"), fill = c(1,2,3,4))
}
if(Rtype == "SINE/Alu"){
  legend("bottomright", legend = c( "AluJ", "AluS", "AluY"), fill = c(1,2,3,4))
}

legend("right", legend = c( "TSD", "Random Sample"), lty = c(1,2), lwd = c(3,3))


#### maybe some summary stats on the layout of TSD families 

layout(matrix(c(1,2), nrow = 1))
pie(summary(tsd$repType), main = "TSD")
pie(summary(as.factor(reps$repType[reps$repType == "LINE/L1" | reps$repType == "SINE/MIR" | reps$repType == "SINE/Alu" | reps$repType == "LINE/L2" ])),  main = "Genome")

#pie(summary(as.factor(reps$repType[c("SINE/Alu", "SINE/MIR",  "LINE/L1", "LINE/L2") %in% reps$repType, ])))


#### My intron dataset too, it would be interesting to see what is happening there 


# the distance to boundaries aswell as the size of the thing. 

intronKeep.gr <- filterIntron(refgene)


# introns we see over representation for one of the compartments 
# one way we could look at it is different kinds of L1s and the sizes of introns they are in 
# another way is we could look at each intron and how many L1s does it have

layout(matrix(1))



tsdListL1 <- NULL
repListL1 <- NULL
L1names <- c("L1ME", "L1MD", "L1MC", "L1MB", "L1MA", "L1PB", "L1PA", "L1HS")


L1names <- c("AluJ", "AluS", "AluY")

for(i in 1:length(L1names)){
  tsdListL1 <- c(tsdListL1, list(tsd[grep(L1names[i], tsd$repFam),]))
  repListL1 <- c(repListL1, list(repSampL1[grep(L1names[i], repSampL1$repName),]))
}
names(tsdListL1) <- names(repListL1) <- L1names


# 
# l <- density(log10(width(intronKeep.gr)))
# sp <- smooth.spline(l$x,l$y,all.knots = F)
# 
# sp <- loess(l$y~l$x)
# plot(sp$x,(sp$y*(10^sp$x) ), type = "l")
# par(new=T)
# plot(sp$x,(sp$y*(10^sp$x) ), type = "n", ylim=c(0,.7))
# 

den <- density((width(intronKeep.gr)), weights = width(intronKeep.gr)/sum(width(intronKeep.gr)))


IntronSamp <- sample(width(intronKeep.gr), 10000,prob = width(intronKeep.gr), replace = F)
IntronHist <- hist((IntronSamp), breaks = 100)


plot(den, xaxt = "n")

abline(v = weighted.mean(log2(width(intronKeep.gr)), width(intronKeep.gr)/sum(width(intronKeep.gr))  ), lwd = 2)

axis(side = 1, at = seq(from = 0, to = 1000,by = 100), labels = seq(from = 0, to = 1000,by = 100)^2)

axis(side = 1, at = seq(from = 0, to = 100,by = .5), labels = round(2^seq(from = 0, to = 100,by = .5)))

plot(den, xlim = c(1,200000))
#plot(c(0,1500000),c(0,25), type = "n")

for(i in 1:8){
  L1FAM <- L1names[i]
  
#  L1groupTSD <- tsdListL1[[L1FAM]][tsdListL1[[L1FAM]][,"tsdLength"] >= 13 & tsdListL1[[L1FAM]][,"repEnd"] > 6000,]
  L1groupTSD <- tsdListL1[[L1FAM]][tsdListL1[[L1FAM]][,"tsdLength"] >= 13,]
  
  tsdL1.gr <- GRanges(seqnames = Rle(L1groupTSD$chr),
                      ranges = IRanges(start = L1groupTSD$start, end = L1groupTSD$end)
  )
  
  L1groupREP <- repListL1[[L1FAM]]
  repL1.gr <- GRanges(seqnames = Rle(L1groupREP$chr),
                      ranges = IRanges(start = L1groupREP$start, end = L1groupREP$end)
  )
  
  fOLtsd <- as.matrix(findOverlaps(tsdL1.gr, intronKeep.gr))
  fOLrep <- as.matrix(findOverlaps(repL1.gr, intronKeep.gr))
  
    denTSD <- density((width(intronKeep.gr)[fOLtsd[,2]]))
     lines((denTSD$x),denTSD$y, type = "l", col = i, lwd = 3)
     abline(v = mean(sqrt(width(intronKeep.gr)[fOLtsd[,2]])), col = i)
  #   
#   denREP <- density(log2(width(intronKeep.gr)[fOLrep[,2]]))
#   lines((denREP$x),denREP$y, type = "l", col = i, lwd = 3, lty = 3)
#   abline(v = mean(log2(width(intronKeep.gr)[fOLrep[,2]])), col = i)
}

#### 


# likly our results are due to a size effect rather than a boundary effect. 

myhist <- hist((width(intronKeep.gr)[fOLtsd[,2]]), breaks = IntronHist$breaks)


myhist$counts/IntronHist$counts *1000

plot(sqrt(myhist$mids),(myhist$counts/IntronHist$counts *1000), xaxt = "n")
axis(side = 1,at = seq(0,20,2),labels = 2^seq(0,20,2))

# so now we are able to plot how many per site

# how can we get the expected number on a random distribution model
# we have x elements acroos b bases per size level 
#par(new = T)




x = (myhist$mids)
y = (myhist$counts/IntronHist$counts *1000)

df <- data.frame(x = x, y = y)
df = df[complete.cases(df),]

df$y = df$y/max(df$y)
df$x = df$x/max(df$x)
plot(df$y ~ df$x)
fit= lm(df$y~df$x)
abline(fit)
par(new = T)
plot((df$y - fit$coefficients[1])~df$x, col = 2, xlim = c(0,1), ylim = c(0,1))
abline(a=0,b=1,col =2)


plot(myhist$mids, y/sum(y,na.rm = T), xlim = c(0,20), type = "p", pch = 16)




sampDF <- data.frame(x=seq(0,1,by = .0001) * max(2^myhist$mids), y = seq(0,1,by = .0001) * max(y,na.rm = T))
lines(log2(sampDF$x),(sampDF$y) + 100)



plot((myhist$mids), y/sum(y,na.rm = T),xaxt = "n")
axis(side = 1,at = seq(0,1000,100), labels = seq(0,1000,100)^2)
lines(myhist$mids, myhist$mids/sum(myhist$mids))
# how do we get our variance

n = sum(y, na.rm = T)
p = myhist$mids/sum(2^myhist$mids)
tsdMean = n*p

plot(myhist$mids, y)
lines(myhist$mids, tsdMean + 300)

tsdSD <- sqrt(n*p*(1-p))
lines(myhist$mids, tsdMean + 300 + (2*tsdSD) , col = 2)
lines(myhist$mids, tsdMean  + 300 - (2*tsdSD) , col = 2)


plot(2^myhist$mids , y - (tsdMean + 300))
lines(2^myhist$mids, (2*tsdSD) , col = 2)
lines(2^myhist$mids, -(2*tsdSD) , col = 2)
abline(h=0)
# so with this model we have the size of the boundary depletion zone
# and also a way to measure the accumulation
# doenst seem like we're capturing all the variation
# some varation may be related to 












par(new = T)
plot(IntronHist$mids, myhist$counts, xaxt = "n", yaxt = "n", col = 3, pch = 16)


plot(10^(denTSD$x),(denTSD$y/(predict(sp,newdata = denTSD$x)$y * predict(sp,newdata = denTSD$x)$x)), type = "l", col = i, lwd = 3)
plot((denTSD$x),(predict(sp,newdata = denTSD$x)$y * predict(sp,newdata = denTSD$x)$x) / denTSD$y, type = "l", col = i, lwd = 3)

plot(denREP$x, denREP$y)
# can we get some sort of rate, like a per bp thing, so turn it into counts

denTSD <- hist(log10(width(intronKeep.gr)[fOLtsd[,2]]), breaks = 100)
plot((denTSD$mids),denTSD$counts / (predict(sp,newdata = denTSD$mids) * (10^denTSD$mids) ), type = "l", col = i, lwd = 3)

sp$y*(10^sp$x)

plot(density(tsdListL1$L1ME$tsdLength[tsdListL1$L1ME$repEnd>=6000]))
lines(density(tsdListL1$L1MA$tsdLength[tsdListL1$L1MA$repEnd>=6000]), col = 2)
lines(density(tsdListL1$L1PA$tsdLength[tsdListL1$L1PA$repEnd>=6000]), col = 3)
lines(density(tsdListL1$L1HS$tsdLength[tsdListL1$L1HS$repEnd>=6000]), col = 4)





# we peobably need to control for how much a certain size is represented 

# we just need to work out what expected is
# sample a bunch of genomic regions and find out the sizes of introns they're in. 
# so line up the whole genome and pick from it



# for the intron anlysis 
# maybe then our small differences will become amplified. 

# so probably the best way to see whats happening with L1 insertions in introns is to take a random sample
# if we take multiple samples from the genome then we can see the probability of getting such a distribhution



## when it comes to introns, what are we looking at? 

# the number of bases in an inton of size x == how many bases are within intron size x

# how many insertions are within an intron of size x

# looking at the denisty, what we find is a slight shift , we were asking about wheather or not there was enrichment 

plot(predict(sp))
sp <- smooth.spline(l$x,l$y,all.knots = T)
predict(low,c(5,6))










