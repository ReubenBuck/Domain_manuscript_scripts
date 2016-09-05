


rm(list = ls())

setwd("./Domain_manuscript/")

# so we make a series of repeats 

set.seed(89)

xl <- rep(0,20)
xr <- sort(c(runif(n = 13,min = .05,max = .4), runif(n = 2,min = .4,max = .6), runif(n = 4,min = .7,max = 1),1))
yb <- seq(.95,0,length.out = 20)
yt <- seq(.96,.01,length.out = 20)

densitiesInterval <- c(rep(.2,10),rep(.05,10))
densitiesBoundary <- c(rep(.1,10),rep(.1,10))


for(p in 0:1){
 pdf(file = paste("writing/round2_20160503/draftsTex/supmaterial/TexFigs/supFig/biasDemonstrate/boundary",p,".pdf", sep =""))
  layout(matrix(c(4,1,2,3),nrow = 2),heights = c(2,1), widths = c(1,2))
  plot.new()
  par(mar = c(0,0,5,5))
  AllSampled <- NULL
  plot.new()
  Boundary = p
  for(i in 1:20){
    lines(c(xl[i],xr[i]), c(yb[i],yb[i]))
    if (Boundary) {
      samplingNumber <- as.integer(densitiesBoundary[i] * (xr[i]/.02))
    }else if (!Boundary) {
      samplingNumber <- as.integer(densitiesInterval[i] * (xr[i]/.02))
    }
    
    if(samplingNumber == 0){next()}
    
    if (!Boundary) {
      sampled <- sample(x = seq(0,xr[i]-.019,by=.02),size = samplingNumber, replace = FALSE)
    } else if(Boundary){
      sampled <- sample(x = seq(0,xr[i] - runif(n = 1,min = 0,max = xr[i] - .019),by=.02),size = samplingNumber, replace = FALSE)
    }
    
    AllSampled <- c(AllSampled, list(sampled))
    names(AllSampled)[length(AllSampled)] <- i
    for(j in sampled){
      rect(xleft = j, xright = j + .02, ybottom = yb[i] , ytop = yt[i], col = 1)
    }
  }
  
  
  
  counts <- hist(unlist(AllSampled), breaks = seq(0,1,length.out = 21), plot = FALSE)$counts
  bpFreq <- hist(sample(x = seq(0,max(xr),length.out = 20),size = 10000,prob = rev(xr), replace = TRUE), 
                 breaks = seq(0,1,length.out = 21), plot = FALSE)$counts
  bpFreq <- (bpFreq/max(bpFreq)) * 20
  
  par(mar = c(5,0,1,5))
  
  plot(smooth(counts/bpFreq),type = "l", ylab = "retrotransposon\ndensity", xlab = "distance from feature boundary (bp)", 
       xaxt = "n",xlim = c(0,20), yaxt = "n", lwd = 3)
  axis(side = 1,at = seq(0,20,5),seq(0,1,.25)*10000)
  axis(side = 2, labels = c("low", "mid", "high"), at = c(0, max(smooth(counts/bpFreq)/2),max(smooth(counts/bpFreq))))
  mtext(text = "retrotransposon\ndensity",side = 2, line = 2)
  
  
  intervalDense <- NULL
  for(i in 1:20){
    intervalDense <- c(intervalDense, length(AllSampled[[as.character(i)]]))
  }
  
  par(mar = c(0,5,5,0))
  plot(x = rev(smooth(intervalDense/xr)), y = 1:20,
       xlim = c(max(smooth(intervalDense/xr))*(1+p),0), type = "l", xlab = "", ylab = "",
       xaxt="n", yaxt = "n", lwd = 3)
  axis(side = 3,labels = c("low", "mid", "high"), at=c(0,(max(smooth(intervalDense/xr))*(1+p))/2,max(smooth(intervalDense/xr))*(1+p) ))
  mtext(text = "retrotransposon\ndensity per interval",side = 3, line = 2)
  
  dev.off()
}




