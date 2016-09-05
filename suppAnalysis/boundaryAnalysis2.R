
library(IRanges)

set.seed(100)

features <- c(0,sort(sample(x = 1:98,9,replace = FALSE)),99)/100
RTNs <- sort(sample((0:99)[-(features*100 +1)], 18, replace = FALSE))/100



pdf(file = "Desktop/Domain_manuscript/writing/round2_20160503/draftsTex/supmaterial/TexFigs/supFig/boundaryAnalysis/boundaryExample.pdf",
        height = 7, width = 6, onefile = TRUE)

layout(matrix(c(1,2,3,4,5), nrow = 5), heights = c(.5,2,.7,.7,.7))

letters <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K")



par(mar = c(0,5,0,5), oma = c(5,5,5,5))
plot(x = features, y = rep(.5, length(features)), axes = FALSE, pch = 15, ylim = c(0,1), type = "n", ylab = "", xlab = "")
rect(xleft = features[2:11] +.01, xright = features[2:11] - (features[2:11] - features[1:10])/2, ybottom = .47,ytop = .53, col = "grey30")
rect(xleft = features[1:10], xright = features[1:10] + (features[2:11] - features[1:10])/2 + .005, ybottom = .47,ytop = .53, col = "grey70")
rect(xleft = features , xright = features + .01, ybottom = .45,ytop = .55, col= "white")
rect(xleft = RTNs,xright = RTNs + .01, ybottom = .45,ytop = .55, col = topo.colors(length(RTNs)), density = -1)
text(x = features + .005,y = .5, labels = letters,pos = 3)

featIntLeft <- diag(as.matrix(dist(features))[2:11,1:10])/2 -.005
featIntRight <- diag(as.matrix(dist(features))[2:11,1:10])/2 -.005



par(mar = c(1,5,.5,5))
plot(0,0,xlim=c(0,1), ylim = c(0,1), axes = FALSE, type = "n", ylab = "", xlab ="")
rect(xleft = 0,xright = sort(featIntLeft/max(featIntLeft)),ybottom = seq(.91,0.01,by = -.1), ytop = seq(.89,-.01,by = -.1), col = "grey30", density = -1)
rect(xleft = 0,xright = sort(featIntRight/max(featIntRight)),ybottom = seq(.96,.06,by = -.1), ytop = seq(.94,.04,by = -.1), col= "grey60", density = -1)
rect(xleft = -1, xright = 0, ybottom = seq(-.02,.93,.05), ytop = seq(.02,.97,.05), col = "white")
# how do we get t


mtext(text = letters[1:10][order(featIntLeft)],side = 2,at = seq(.95,.05,-.1), las = 2,adj = 0, line = 1, cex = .6)
mtext(text = letters[2:11][order(featIntLeft)],side = 2,at = seq(.90,0,-.1), las = 2,adj = 0, line = 1, cex = .6)

 

# we need to take into account featue size

insertionLoci <- list(start = NULL, end = NULL)
cols <- topo.colors(length(RTNs))
for(i in 1:length(RTNs)){
  distance <- (RTNs[i] - features)[order(abs(RTNs[i] - features))[1]] /max(featIntRight)
  closest <- order(abs(RTNs[i] - features))[1]
  if(distance > 0){
    rect(xleft = abs(distance) - (.01 / max(featIntLeft)) , xright = abs(distance) ,ybottom= seq(.93,-.02,-.1)[order(featIntLeft) == closest], ytop= seq(.97,.02,-.1)[order(featIntLeft) == closest], density = -1, col = cols[i])
    insertionLoci$start = c(insertionLoci$start, abs(distance) - (.01 / max(featIntLeft)) )
    insertionLoci$end = c(insertionLoci$end, abs(distance))
  }
  if(distance < 0){
    rect(xleft = abs(distance) - (.01 / max(featIntRight)) , xright = abs(distance) ,ybottom= seq(.88,-.02,-.1)[order(featIntRight) == closest-1], ytop= seq(.92,.02,-.1)[order(featIntLeft) == closest-1], density = -1, col = cols[i])
    insertionLoci$start = c(insertionLoci$start, abs(distance) - (.01 / max(featIntRight)))
    insertionLoci$end = c(insertionLoci$end, abs(distance))
  }
}


rangeRTN <- IRanges(start = round(insertionLoci$start*10000), end = round(insertionLoci$end*10000) -1)

# so do a coverage plot for what we see on the other plot
rangeInt <- IRanges(start = rep(1,20), end = round(c(featIntLeft/max(featIntLeft), featIntLeft/max(featIntLeft)) * 10000) - 1)


plot(c(coverage(rangeRTN),rep(0,10000))[1:9999], type = "l", xlim = c(0,10000), xaxt = "n", las = 2, ylab = "retrotransposon\ncoverage", lwd = 3)
grid()

plot(coverage(rangeInt), type = "l", ylim = c(0,20), xaxt = "n", las = 2, ylab = "position depth", lwd = 3)
grid()

plot(c(coverage(rangeRTN),rep(0,5000))[1:9999]/coverage(rangeInt), type = "l", ylim = c(0,.6), yaxt = "n", ylab = "retrotransposon\ndensity", xlab = "distance from boundary (kb)", lwd = 3)
axis(side = 2, at = seq(0,1,.1), las = 2)
mtext("distance from boundary (bp)", outer = TRUE, side = 1, line = 2, cex = .8)
grid()


layout(1)
plot.new()
legend("topleft", legend = c("feature", "retrotransposon"),fill = c("white", "green"), bty = "n", title = "DNA element")
legend("right", legend = c("left", "right"),fill = c("grey70", "grey30"), bty = "n", title = "inter-feature interval")

axis(side = 1,at = c(0,.2,.4), labels = c(0,NA,"2 kb"))
axis(side = 1,at = c(.6,.8,1), labels = c(0,NA,"20 kb"))

dev.off()

max()

plot.new()
axis(side = 1,at = c(,.2,.4), labels = c(0,NA,"200 kb"))
axis(side = 1,at = c(.6,.8,1), labels = c(0,NA,"20 kb"))
#text(.7,-.05,"kb")


