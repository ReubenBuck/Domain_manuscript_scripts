
library(IRanges)
rm(list = ls())

setwd("./Domain_manuscript/")

# so we make a series of repeats 

set.seed(89)

xl <- rep(0,20)
xr <- NULL
xr[seq(2,20,2)] <- sort(c(runif(n = 6,min = .05,max = .4), runif(n = 1,min = .4,max = .6), runif(n = 2,min = .7,max = 1),1))
xr[seq(1,19,2)] <- xr[seq(2,20,2)]
yb <- seq(.95,0,length.out = 20)
yt <- seq(.96,.01,length.out = 20)


# lets make our genome
sampInt <- sample(10, size = 10,replace = FALSE)
intervalLength <- (xr[seq(1,19,2)]/sum(xr[seq(1,19,2)]))[sampInt]
labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K")
intervalDF <- data.frame(length = intervalLength, start = labels[1:10], end = labels[2:11], startNo = 1:10, endNo = 2:11 )
RTNs <- sort(sample(1:100,size = 20, replace = FALSE))/100


featurePosition <- 0
for(i in 1:10){featurePosition <- c(featurePosition, sum(intervalLength[1:i]))}


intervalDF[order(intervalDF$length),"lineStart"] <- seq(20,2,-2)
intervalDF[order(intervalDF$length),"lineEnd"] <- seq(19,1,-2)







#pdf(file = "Desktop/Domain_manuscript/writing/round2_20160503/draftsTex/supmaterial/TexFigs/supFig/boundaryAnalysis/boundaryExample.pdf",
#    height = 7, width = 6)

layout(matrix(c(1,2,3,4,5), nrow = 5), heights = c(.5,2,.7,.7,.7))

par(mar = c(0,0,0,0), oma = c(5,5,5,5))
plot(1, xlim = c(0,1), ylim = c(0,1), axes = F, type = "n")
#abline(h=0)
for(i in 1:(length(featurePosition)-1)){
#  lines(c(featurePosition[i], featurePosition[i] + (featurePosition[i+1] - featurePosition[i])/2), c(0.1,0.1), col = "grey80", lwd = 3)
  lines(c(featurePosition[i], featurePosition[i] + (featurePosition[i+1] - featurePosition[i])/2) ,c(.1,.1), col = "grey70",  lwd = 5)
}
for(i in 1:(length(featurePosition)-1)){
#  lines(c(featurePosition[i] + (featurePosition[i+1] - featurePosition[i])/2, featurePosition[i+1]), c(0.1,0.1), col = "grey60", lwd = 3)
  lines(c(featurePosition[i] + (featurePosition[i+1] - featurePosition[i])/2,  featurePosition[i+1]),c(.1,.1), col = "grey30",  lwd = 5)
  
}
text(x = featurePosition, y = rep(0.1, 10), labels = labels, pos = 3, cex = 1)
#for(i in 1:length(RTNs)){
 # points(x= RTNs[i], y = rep(0.1), pch = 15, cex = 1.4, col = 1)
  #points(x= RTNs[i], y = rep(0.1), pch = 15, cex = 1.2, col = topo.colors(length(RTNs))[i])
#  print(i)
#}
rect(xleft = RTNs - .01, xright = RTNs, ybottom = 0,ytop = .2,density = -1, col = topo.colors(length(RTNs)), lwd = .5)

points(x = featurePosition, y = rep(0.1, 11), cex = 1.6, pch = 15, col = "grey40")

#legend("topleft", legend = c("feature", "retrotransposon"), col = c("grey40", topo.colors(1)), pch = c(15,16), bty = "n")
# split intervals and align by boundary    

par(mar = c(1,5,.5,5))
plot.new()
for(i in 1:20){
 lines(c(xl[i],xr[i]), c(yb[i],yb[i]), col = c("grey70","grey30")[(i %% 2 == 0) + 1], lwd = 5)
#  rect(xleft = xl[i], xright = xr[i], ybottom = yb[i]-.01, ytop = yb[i]+.01, col = c("grey70","grey30")[(i %% 2 == 0) + 1], angle = c(135,45)[(i %% 2 == 0) + 1], density = 20, border = c("grey70","grey30")[(i %% 2 == 0) + 1],lwd = 2)
}
mtext( at = intervalDF$lineStart * .05 - .05, text = intervalDF$start, side = 2, las = 2, cex = .7)
mtext( at = intervalDF$lineEnd * .05 - .05, text = intervalDF$end, side = 2, las = 2, cex = .7)

insertionLoci <- NULL
cols <- topo.colors(length(RTNs))
for(i in 1:length(RTNs)){
  distance <- (RTNs[i] - featurePosition)[order(abs(RTNs[i] - featurePosition))[1]]
  closest <- order(abs(RTNs[i] - featurePosition))[1]
  if(distance > 0){
#    points(x = abs(distance) * sum(xr[seq(1,19,2)]) * 2,  y = intervalDF$lineStart[intervalDF$startNo == closest]  *.05 -.05, pch = 15, cex = 1.7, col = 1)
#    points(x = abs(distance) * sum(xr[seq(1,19,2)]) * 2,  y = intervalDF$lineStart[intervalDF$startNo == closest]  *.05 -.05, pch = 15, cex = 1.5, col = cols[i])
    rect(xleft = abs(distance) * sum(xr[seq(1,19,2)]) * 2-(.01 * sum(xr)), xright = (abs(distance) * sum(xr[seq(1,19,2)]) * 2),ybottom= intervalDF$lineStart[intervalDF$startNo == closest]  *.05 -.07, ytop= intervalDF$lineStart[intervalDF$startNo == closest]  *.05 -.03, density = -1, col = cols[i])
  }
  if(distance < 0){
#    points(x = abs(distance) * sum(xr[seq(1,19,2)]) * 2, y =  intervalDF$lineEnd[intervalDF$endNo == closest]  *.05 - .05, pch = 15, cex = 1.7, col = 1)
#    points(x = abs(distance) * sum(xr[seq(1,19,2)]) * 2, y =  intervalDF$lineEnd[intervalDF$endNo == closest]  *.05 - .05, pch = 15, cex = 1.5, col = cols[i])
    rect(xleft = abs(distance) * sum(xr[seq(1,19,2)]) * 2-(.01 * sum(xr)), xright = (abs(distance) * sum(xr[seq(1,19,2)]) * 2) + .005,ybottom= intervalDF$lineEnd[intervalDF$endNo == closest]  *.05 -.07, ytop= intervalDF$lineEnd[intervalDF$endNo == closest]  *.05 -.03, density = -1, col = cols[i])
  }
  insertionLoci <- c(insertionLoci, abs(distance) * sum(xr[seq(1,19,2)]) * 2)
}
points(rep(0,20), seq(0,.95,.05), pch = 15 ,col = "grey40", cex =1.5)

range <- IRanges(start = round(insertionLoci*10000)-100, end = round(insertionLoci*10000) + 100)


# need to sort out the problems of numbering the intervals approproatly first
par(mar=c(1,5,.5,5))

insertionLociH <- hist(insertionLoci, breaks = seq(0,1,.05), plot = F)
insertionLociH <- c(t(matrix(rep(smooth(insertionLociH$counts),5), byrow = FALSE, nrow = 20))[1:100],0)



# so how do we get the distribution
bpFreq <- NULL
for(i in 1:length(xr)){
  bpFreq <- c(bpFreq, seq(0, xr[i], .01))
}
bpFreq <- table(as.factor(bpFreq))

plot(names(bpFreq), bpFreq[1:101], type = "l", xlim = c(0,1), xaxt = "n", las = 2, ylim = c(0,20), ylab = "position depth", lwd = 3)

plot(seq(0,1,.01), insertionLociH, type = "l", xlim = c(0,1), xaxt = "n", las = 2 , ylab = "retrotransposon\ncoverage", lwd = 3)

plot(seq(0,1,.01),insertionLociH/bpFreq, type = "l", ylim = c(0,.6), xaxt = "n", yaxt = "n", ylab = "retrotransposon\ndensity", xlab = "distance from boundary (kb)", lwd = 3)
axis(side = 2, at = seq(0,1,.1), las = 2)
axis(side = 1, at = seq(0,1,.1), label = seq(0,1,.1)*10000)

mtext("distance from boundary (bp)", outer = TRUE, side = 1, line = 2, cex = .8)

#dev.off()





# dark needs to be darker 
# explain smoothing 
# text needs to explain the transitioon from genome to half intervals 
# maybe give the RTNs some length so that 
# change bp frequency to depth
# rectangles with hatching
range <- IRanges(start = round(RTNs*10000)- 50, end = round(RTNs*10000) + 50)
plot(coverage(range), type = "l")

  intervalLength

  intervalDF[order(intervalDF$length),]

#IntervalDF <- data.frame(xl = xl, xr = xr, yb = yb, yt = yt)
labels[sampInt]

densitiesBoundary <- rep(.1,20)
# if there is someway to put stuff 