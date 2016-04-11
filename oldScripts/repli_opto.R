x <- c(hist(rnorm(10000), breaks = 100)$counts, -hist(rnorm(10000,sd=2), breaks = 100)$counts)


rand.point <- function(){
  min = runif(min = -10, max = 0, n = 1)
  return(runif(100, min = min, max = min+10))
}


x = replicate(n = 10,expr = rand.point(),simplify = TRUE)
x <- x[1:length(x)]
plot(x, ylim = c(-10,10))

lineSP <- smooth.spline(1:1000,y=x,all.knots = TRUE)
driv <- predict(lineSP,deriv = 1)

lines(lineSP)
par(new=TRUE)
plot(driv$x,abs(driv$y), col = 3, type = "l", ylim = c(0,.25))

hist(abs(driv$y))

# just use an spline and clacualte the dreivitive and find a cut off
# not too hard

plot(lineSP$y)
lineSP <- smooth.spline(lineSP$x,y=lineSP$y,all.knots = TRUE)
lines(lineSP$y)


RT <- read.table("Downloads/RT_H9_Neural Crest_All.txt",header = TRUE,sep = "\t", skip = 13, nrows = 100000)
RT <- RT[RT$Chromosome=="chr1",]
x = (RT$End_Position- RT$Start_Position)/2 + RT$Start_Position
y = RT$Data_Value
xlim = c(7000000,132000000)
plot(x,y, type = "l", xlim = xlim, pch = 16, cex = .1)
# data points per bp, to how many knots per bp
spLine <- smooth.spline(x,y, nknots = as.integer(length(y)/300))
lines(spLine,col = 3)
deriv <- predict(spLine,deriv = 1)
par(new=TRUE)
plot(x = deriv$x,y=(deriv$y), type = "l", col = 2, yaxt = "n", xlim = xlim)
abline(h = 2.75 * (10^(-6)), lty = 1)
abline(h = -2.75 * (10^(-6)), lty = 1)
abline(h = .8 * (10^(-6)), lty = 2)
abline(h = -.8 * (10^(-6)), lty = 2)

# alternativly we could use a HMM to get three states

plot(deriv$x,(abs(deriv$y)), type = "l", xlim = c(10000000,50000000))
sm.spline <- smooth.spline(x = deriv$x,y=abs(deriv$y),nknots = as.integer((max(x)- min(x))/ 1000000))
lines(sm.spline$x, sm.spline$y, col = 2)

library(depmixS4)

nStates = 4
mod <- depmix(list((y)~1),nstates=nStates,
              family=list(gaussian()), ntimes = length(y))
fit <- fit(mod)
pos <- posterior(fit)

df = data.frame(firstOrder = y, state = as.factor(pos$state))

layout(matrix(c(1)))

plot(x, y,type = "n",xlim = c(50000000,100000000))
for(i in 1:nStates){
  state <- rep(NA, length(y))
  state[df$state == i] = y[df$state==i]
  lines(x,state,col = i)
}



#boxplot((df$firstOrder) ~ df$state)
stripchart(abs(df$firstOrder) ~ df$state, method = "jitter", vertical = TRUE, jitter = .5, pch = 16, cex = .1)


plot(c(min(spLine$x), max(spLine$x)), c(log(1),log(2)), type = "n", xlim = c(10000000,50000000))
lines(spLine$x, log(pos$S1+1), col = 2)
lines(spLine$x, log(pos$S2+1), col = 3)
lines(spLine$x, log(pos$S3+1), col = 4)



plot(deriv$x,deriv$y, type = "l")
#lines(spLine$x,spLine$y, col = 3)
par(new=TRUE)
plot(spLine$x,pos$state, type = "p", col = 2, pch = 16, cex = .5)


plot(spLine$x, spLine$y, type = "l", xlim = c(50000000,100000000))
par(new=TRUE)
plot(sm.spline$x,pos$state, type = "p", col = 2, pch = 16, cex = .5, xlim = c(50000000,100000000))




##### not sure how to find these regions, maybe i shoudl sample random segments of the genome
##### cluster them 



# so so far the more accurate approach is to use the slope classification
# so doing a hmm on the first order derivitive we can expand the measures a bit further
# the process become more a priori 
# then it becomes easier to select things that cross a certain range 

# don't like the way its cutting straight lines
# one thing we can do is smooth over the derviative a second time 

library(DNAcopy)

x = (RT$End_Position- RT$Start_Position)/2 + RT$Start_Position

y = RT$Data_Value
xlim = c(1000000,52000000)
plot(x,y, type = "p", xlim = xlim, pch = 16, cex = .1)
CNA.object <- CNA(genomdat = RT$Data_Value, maploc = RT$Start_Position, chrom = RT$Chromosome, data.type="logratio")
#CNA.sm <- smooth.CNA(CNA.object,smooth.region = 1000)
#lines(CNA.sm$maploc,CNA.sm$Sample.1, col = 2)

CNA.seg <- segment(CNA.object, nperm=10000, alpha=1e-10, undo.splits="sdundo", undo.SD=1.5, verbose=2)

plot(CNA.seg, plot.type = "w")


hist(y, breaks = 1000)
plot(x,y, type = "l", xlim = c(200000000,240000000))
hist((abs(y[1:1000]- y[2:1001])), breaks = 100)




###### Lets try with spat stats


RT <- read.table("~/Downloads/RT_GM12801_Lymphoblastoid_All.txt",header = TRUE,sep = "\t", skip = 13, nrows = 100000)
RT <- RT[RT$Chromosome=="chr1",]
x = (RT$End_Position- RT$Start_Position)/2 + RT$Start_Position
y = RT$Data_Value


source("~/Desktop/Domain_manuscript/Domain_manuscript_scripts/functions.R")


library(spdep)
NB <- cell2nb(nrow = length(x), ncol = 1)
wList <- nb2listw(NB)
LM <- localmoran(y,wList)

xlim = c(50000000,70000000)
layout(matrix(c(1,2), nrow = 2))
plot(x = x, y = y, type = "l", xlim = xlim)
plot(x,LM[,1], type = "l", xlim = xlim)
lines(x,sqrt(LM[,3]) )

# so here we can just use it as a meassure and decided at what level we can calassify regions 
### this method has potential cause it is able to capture a think we want which is how similar are adjacent points
### However it will be dificult to properly prepare the weight matrix


### alternativly there are already pre identified regions based on the university of washington data where we cna observe interesting cahnges 
### that migh be enough to show at the meeting on friday






