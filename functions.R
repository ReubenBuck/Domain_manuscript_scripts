# write out each function so it can handle the genome 
# this means work off the alternative matricies




# how to use th enew matrix to do this

# we have to worry about zero values telling us there are no neighbors
# when we calculate the w matrix

neighbor.matrix = function(binned_genome, X){
  B <- binned_genome
  C <- B
  C$start.end <- NA
  C$end.start <- NA
  
  C[1:(nrow(C)-1), "start.end"] = C[2:(nrow(C)),"start"]
  C[2:(nrow(C)), "end.start"] = C[1:(nrow(C)-1),"end"]
  D <- data.frame(C$start.end - C$end, C$start - C$end.start)
  D[is.na(D)] = 0
  
  coord.1 <- coord.2 <- matrix(0, nrow = nrow(B), ncol = 2)
  coord.1[1:nrow(B),1] = 1:nrow(B)
  coord.1[1:(nrow(B)-1),2] = 2:nrow(B)
  coord.1[D[,1] != 1,2] = coord.1[D[,1] != 1,1]
  coord.2[1:nrow(B),1] = 1:nrow(B)
  coord.2[2:(nrow(B)),2] = 1:(nrow(B)-1)
  coord.2[D[,2] != 1,2] = coord.2[D[,2] != 1,1]
  
  ### Removing all lonly bins from B and rewrite the above
  keep = coord.1[,2] != coord.2[,2]
  B <- B[keep,]
  
  C <- B
  C$start.end <- NA
  C$end.start <- NA
  
  C[1:(nrow(C)-1), "start.end"] = C[2:(nrow(C)),"start"]
  C[2:(nrow(C)), "end.start"] = C[1:(nrow(C)-1),"end"]
  D <- data.frame(C$start.end - C$end, C$start - C$end.start)
  D[is.na(D)] = 0
  
  coord.1 <- coord.2 <- matrix(0, nrow = nrow(B), ncol = 2)
  
  coord.1[1:nrow(B),1] = 1:nrow(B)
  coord.1[1:(nrow(B)-1),2] = 2:nrow(B)
  coord.1[D[,1] != 1,2] = coord.1[D[,1] != 1,1]
  
  
  coord.2[1:nrow(B),1] = 1:nrow(B)
  coord.2[2:(nrow(B)),2] = 1:(nrow(B)-1)
  coord.2[D[,2] != 1,2] = coord.2[D[,2] != 1,1]
  
  keepers = coord.1[,1]
  
  Neighbors = matrix(ncol = 3, nrow = nrow(D))
  
  Neighbors[,1] = X[coord.2[,2]]
  Neighbors[,2] = X[coord.2[,1]]
  Neighbors[,3] = X[coord.1[,2]]
  Neighbors[Neighbors[,1] == Neighbors[,2],1] = 0
  Neighbors[Neighbors[,2] == Neighbors[,3],3] = 0
  
  w.mat = matrix(ncol = 3, nrow = nrow(D))
  w.mat[,2] = coord.2[,1]
  w.mat[,1] = coord.2[,2]
  w.mat[,3] = coord.1[,2]
  w.mat[w.mat[,1] == w.mat[,2],1] = 0
  w.mat[w.mat[,2] == w.mat[,3],3] = 0
  w.mat[w.mat > 0] = 1
  
  mats <- c(list(Neighbors), list(w.mat), list(keepers))
  names(mats) = c("neighbors", "w.mat", "keepers")
  return(mats)
  
}






hotspot.G.WG <- function(X, neighbors, w.mat){
  Weight.matrix = w.mat
  n = length(X)
  X.bar = mean(X)
  S = sqrt((sum(X^2) / n) - X.bar^2)
  
  NUM = (rowSums(neighbors) - (X.bar * rowSums(Weight.matrix)))
  DEN = S * sqrt(((n * rowSums(Weight.matrix^2)) - (rowSums(Weight.matrix)^2))/(n-1))
  Z_score = NUM/DEN
  
  
  P_value = 2*pnorm(-abs(Z_score))
  FDR = p.adjust(p=P_value, method="fdr")
  return(data.frame(G.Z_socre = Z_score, P_value = P_value, FDR = FDR))
}



hotspot.G <- function(X, w){
  n = length(X)
  X.bar = mean(X)
  S = sqrt((sum(X^2) / n) - X.bar^2)
  Z_score = (rowSums(t(t(w) * X)) - (X.bar * rowSums(w)))/(S * sqrt(((n * rowSums(w^2))- rowSums(w)^2)/(n-1)))
  P_value = 2*pnorm(-abs(Z_score))
  FDR = p.adjust(p=P_value, method="fdr")
  
  return(data.frame(G.Z_socre = Z_score, P_value = P_value, FDR = FDR))
}




Z.moran.WG = function(X, neighbors, w.mat){
  
  neighbors[,2] = 0
  w.mat[,2] = 0
  w = w.mat
  
  X.bar = mean(X)
  n = length(X)
  S0 <- sum(w) 
  # calculate I first
  # n is the number of elements 
  I = (n/S0) * (sum((X-X.bar) * (neighbors-X.bar) * w)/sum((X - X.bar)^2))
  # so i can get the I value 
  # next we get  expected I
  E.I <- -1 / (n - 1)
  # Now to clculate the variance 
  #    S1
  S1 = sum((w + (w))^2)/2
  #    S2
  S2 = sum((rowSums(w) + rowSums(w))^2)
  ### D
  D = sum((X - X.bar)^4)/(sum((X - X.bar)^2))^2
  # A, B, C no longer require any loops, so pretty easy 
  A = n * ( ((n^2 - (3*n) + 3 ) * S1) - (n * S2) + (3 * S0^2) )
  B = D * ( ((n^2 - n) * S1) - ((2 * n) * S2) + (6 * S0^2) )
  C = (n - 1) * (n - 2) * (n -3) * S0^2
  E.I2 <- (A - B) / C
  # calculate the variance
  V.I = E.I2 - E.I^2
  # calculate z.score
  Z.I <- (I - E.I)/sqrt(V.I)
  P_value = 2*pnorm(-abs(Z.I))
  return(c(Moran_I = I, Z_score = Z.I, Var = V.I, P_value = P_value))
}

# the last one to do will be Global_G
# annoyingly general G won't seem to show a negative Z score
# our peak example shows that there is definatly positve levels of clustering but our 
# trough example does not,
# for some reason it is too hard to get below the mean 
# also significance seems really fragile.
# it might be the way we calculate variance 

# as we increase each X by the same number the varaince gets smaller
# an interesting result considering we are really looking at this as a ratio


General.G <- function(X,w){
  n = length(X)
  X.bar = mean(X)
  w1 <- matrix(1,ncol = n, nrow = n)
  diag(w1) = 0
  General_G = sum(outer(X,X)*w)/sum(outer(X,X)*w1)
  # Expected
  W = sum(w)
  E.G <- W/(n*(n-1))
  
  #expected G square
  n.r4 = n * (n-1) * (n-2) * (n-4)
  
  # i think
  S1 = sum((w + t(w))^2)/2
  # there may be a better way to do this for the genome scale analysis. 
  # essentially becasue we know we are only doing one dimension
  S2 = sum((rowSums(w) + colSums(w))^2)
  
  
  # maybe we could adjust the Bs up here
  
  B0 = (((n^2) - (3*n) + 3)*S1) - (n * S2) + (3*(W^2))
  B1 = -( (((n^2) - n) *S1) - ((2*n) * S2) + (3*(W^2))  )
  B2 = -(  ( 2*n*S1) - ((n + 3) * S2) + (6*(W^2))  )
  B3 = ((4*(n - 1))*S1) - ((2 * (n +1))* S2) + (8*(W^2))
  B4 = S1 - S2 + (W^2)
  
  B0 = B0*(M.j(X,2)^2)
  B1 = B1*(M.j(X,4))
  B2 = B2*(M.j(X,1)^2)*(M.j(X,2))
  B3 = B3*(M.j(X,1))*(M.j(X,3))
  B4 = B4*(M.j(X,1)^4)
  
  
  
  EG.2 <- (1/( (((M.j(X,1)^2)- (M.j(X,2)))^2) * n.r4  )) *(B0 + B1 + B2 + B3 + B4) 
  
  V.G <- EG.2 - (E.G^2)
  
  
  General_G.Z = (General_G - E.G)/sqrt(V.G)
  P_value = 2*pnorm(-abs(General_G.Z))
  
  return(c(General_G = General_G, Z_score = General_G.Z, Var = V.G, P_value = P_value))
  
  
}


Z.moran = function(X, w){
  diag(w) <- 0
  X.bar = mean(X)
  n = length(X)
  S0 <- sum(w) 
  # calculate I first
  # n is the number of elements 
  I = (n/S0) * (sum(outer(X - X.bar, X - X.bar) * w)/sum((X - X.bar)^2))
  # so i can get the I value 
  # next we get  expected I
  E.I <- -1 / (n - 1)
  # Now to clculate the variance 
  #    S1
  S1 = sum((w + t(w))^2)/2
  #    S2
  S2 = sum((colSums(w) + rowSums(w))^2)
  ### D
  D = sum((X - X.bar)^4)/((sum((X - X.bar)^2))^2)
  # A, B, C no longer require any loops, so pretty easy 
  A = n * ( (( (n^2) - (3*n) + 3 ) * S1) - (n * S2) + (3 * (S0^2)) )
  B = D * ( ((n^2 - n) * S1) - ((2 * n) * S2) + (6 * (S0^2)) )
  C = (n - 1) * (n - 2) * (n -3) * (S0^2)
  E.I2 <- (A - B) / C
  # calculate the variance
  V.I = E.I2 - E.I^2
  # calculate z.score
  Z.I <- (I - E.I)/sqrt(V.I)
  
  return(c(Moran_I = I, Z_score = Z.I))
}


M.j <- function(X,j){
  return(sum(X^j))
}





# biggest challange right now is to get the outer of the X vaector

sample.space.Gen_G <- function(X){
  X <- as.numeric(X)
  A = NULL
  for(i in 1:(length(X) -1)){
    A <- c(A,X[i]*(sum(X[(i+1):length(X)])))
  }
  return(sum(A,na.rm=T)*2)
}


General.G.WG <- function(X,neighbors,w.mat){
  
  neighbors[,2] = 0
  w.mat[,2] = 0
  w = w.mat

  
  n = length(X)
  X.bar = mean(X)
    
  General_G = sum((X) * (neighbors) * w)/sample.space.Gen_G(X)
  # Expected
  W = sum(w)
  E.G <- W/(n*(n-1))
  #expected G square
  n.r4 = n * (n-1) * (n-2) * (n-4)
  S1 = sum((w + (w))^2)/2
  S2 = sum((rowSums(w) + rowSums(w))^2)
  # maybe we could adjust the Bs up here
  B0 = (((n^2) - (3*n) + 3)*S1) - (n * S2) + (3*(W^2))
  B1 = -( (((n^2) - n) *S1) - ((2*n) * S2) + (3*(W^2))  )
  B2 = -(  ( 2*n*S1) - ((n + 3) * S2) + (6*(W^2))  )
  B3 = ((4*(n - 1))*S1) - ((2 * (n +1))* S2) + (8*(W^2))
  B4 = S1 - S2 + (W^2)
  B0 = B0*(M.j(X,2)^2)
  B1 = B1*(M.j(X,4))
  B2 = B2*(M.j(X,1)^2)*(M.j(X,2))
  B3 = B3*(M.j(X,1))*(M.j(X,3))
  B4 = B4*(M.j(X,1)^4)
  EG.2 <- (B0 + B1 + B2 + B3 + B4)/ ( (((M.j(X,1)^2)- (M.j(X,2)))^2) * n.r4  )
  V.G <- EG.2 - (E.G^2)
  General_G.Z = (General_G - E.G)/sqrt(V.G)
  P_value = 2*pnorm(-abs(General_G.Z))
  return(c(General_G = General_G, Z_score = General_G.Z, Var = V.G, P_value = P_value))  
}




# (X - mean) / sd

#seems i found the bug
#we were getting the score for the whole group




### this function isn't working
# we are getting diffent answears if we feed it different bin sizes
binned.genome.reader <- function(genome = "hg19", bin.size = c(1500000), keep.rate = .9){
  
  # getting chrom info 
  web <- paste("http://hgdownload.soe.ucsc.edu/goldenPath/", genome, "/database/chromInfo.txt.gz", sep = "")
  con <- gzcon(url(web))
  txt <- readLines(con)
  chrom_info <- read.delim(textConnection(txt), header = FALSE)
  
  
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
  
  Final <- NULL
  for(z in 1:length(bin.size)){
    chrom_info1 <- chrom_info[chrom_info[,2] > bin.size[z],]
    
    bins <- NULL
    for( i in 1:dim(chrom_info1)[1]){
      start <- seq(from = 1, to = chrom_info1[i,2], by = bin.size[z])
      end <- c(seq(from = bin.size[z], to = chrom_info1[i,2], by = bin.size[z]), chrom_info1[i,2])
      pre.bin <- data.frame(chr = chrom_info1[i,1], start = start, end = end)
      bins <- rbind(bins, pre.bin)
    }
    
    bins$Known <- bins$end - bins$start + 1
    bins.gr <- GRanges(seqnames = Rle(bins$chr), 
                       ranges = IRanges(start = bins$start, end = bins$end - 1)
    )
    
    gap.int.gr <- intersect(gaps.gr, bins.gr)	   
    gap.bin.ol <- as.matrix(
      findOverlaps(bins.gr, gap.int.gr)
    )	   
    gap.size <- data.frame(bin = gap.bin.ol[,1], width = width(gap.int.gr[gap.bin.ol[,2]]))
    gap.size.bin <- aggregate(gap.size[,2],list(bin = gap.size$bin), sum)
    bins$Known[gap.size.bin$bin] <- bins$Known[gap.size.bin$bin] - gap.size.bin$x + 1
    # so now we have all the known bp for bins for a genome of any bin size
    bins <- bins[bins$Known >= keep.rate*bin.size[z],1:4]
    rownames(bins) <- 1:nrow(bins)
    Final <- c(Final, list(bins))
  }
  names(Final) <- paste(genome,"bin.size", bin.size, sep = "_")
  return(Final)
  
  
}





SquareNBqueen <- function(M, cbound, rbound, NBdist) {
  
  diag <- c(1:ncol(M) * ncol(M)) - ncol(M) + 1:nrow(M)
  diag.all <- diag
  for(i in 1:NBdist){
    diag.all <- c(diag.all, (diag + i)[1:(length(diag) - i)], (diag - i)[(1 + i):length(diag)])
  }
  
  zeros = NULL
  for(cut in 1:length(cbound)){
    out = NULL
    for(i in 1 : (nrow(M) - cbound[cut])){
      out <- c(out, c( (1:nrow(M)* nrow(M)) - nrow(M) + cbound[cut] + i , ((cbound[cut] + i) * nrow(M)) - nrow(M) + 1:nrow(M)))
    }
    out <- out[!( out %in% out[duplicated(out)] )]
    zeros <- c(zeros, out)
  }
  zeros.col <- unique(zeros)
  
  min.module <- diag.all[!(diag.all %in% zeros.col)]
  
  # get our max module immage and plot it out here 
  zeros = NULL
  for(cut in 1:length(rbound)){
    out = NULL
    for(i in 1 : (nrow(M) - rbound[cut])){
      out <- c(out, c( (1:nrow(M) * nrow(M)) - nrow(M) + rbound[cut] + i , ((rbound[cut] + i) * nrow(M)) - nrow(M) + 1:nrow(M)))
    }
    out <- out[!( out %in% out[duplicated(out)] )]
    zeros <- c(zeros, out)
  }
  zeros.row <- unique(zeros)
  max.module <- diag.all[!(diag.all %in% zeros.row)]
  
  # we just need to caculte the factors to add on to place them in the big matrix
  # we need the row column info to place them out 
  
  rows = matrix(data=rep(1:nrow(M), nrow(M)), nrow = nrow(M), ncol = ncol(M))
  cols = matrix(data=rep(1:nrow(M), nrow(M)), nrow = nrow(M), ncol = ncol(M), byrow=T)
  
  rows.max <- rows[max.module]
  cols.max <- cols[max.module]
  rows.min <- rows[min.module]
  cols.min <- cols[min.module]
  
  w.mat <- matrix(data = 0, nrow=length(M), ncol = length(M))
  
  for(i in 1:length(max.module)){
    w.mat[(cols.min * length(M) - length(M) + rows.min)  +   ((cols.max[i] * (length(M) * nrow(M))) - (length(M) * nrow(M)) + rows.max[i] * nrow(M) - nrow(M))] <- 1
  }
  
  return(w.mat)
}



geneModelplot <- function(chr.choice = NULL, pos.coor = NULL, neg.coor = NULL, intron.width = NULL,
                          refgene. = NULL){
  
  refgene. = refgene.[refgene.[,3] == chr.choice,]
  
  Estart <- strsplit(x=as.character(refgene.[,10]),split=",")
  Eend <- strsplit(x=as.character(refgene.[,10]),split=",")
  Istart <- Eend
  Iend <- Estart
  for(i in 1:nrow(refgene.)){
    if(refgene[i,9] > 1){
      Istart[[i]] = as.numeric(Istart[[i]][1:(length(Istart[[i]])-1)])
      Iend[[i]] = as.numeric(Iend[[i]][2:(length(Iend[[i]]))])
    }else{
      Istart[[i]] = as.numeric(Estart[[i]])
      Iend[[i]] = as.numeric(Eend[[i]])
    }
    Estart[[i]] = as.numeric(Estart[[i]])
    Eend[[i]] = as.numeric(Eend[[i]])
  }
  
  
  Exons.start.pos <- Estart[refgene.[,4] == "+"]
  Exons.start.neg <- Estart[refgene.[,4] == "-"]
  Exons.end.pos <- Eend[refgene.[,4] == "+"]
  Exons.end.neg <- Eend[refgene.[,4] == "-"]
  Introns.start.pos <- Istart[refgene.[,4] == "+"]
  Introns.start.neg <- Istart[refgene.[,4] == "-"]
  Introns.end.pos <- Iend[refgene.[,4] == "+"]
  Introns.end.neg <- Iend[refgene.[,4] == "-"]
  
  
  # exon rect
  for(i in 1:length(Exons.start.pos)){
    for(j in 1:length(Exons.start.pos[[i]])){
      rect(xleft=Exons.start.pos[[i]][j], xright=Exons.start.pos[[i]][j],
           ytop=pos.coor+intron.width, ybottom=pos.coor-intron.width,
           density=-1,col = 1)
    }
  }
  for(i in 1:length(Exons.start.neg)){
    for(j in 1:length(Exons.start.neg[[i]])){
      rect(xleft=Exons.start.neg[[i]][j], xright=Exons.start.neg[[i]][j],
           ytop=neg.coor+intron.width, ybottom=neg.coor-intron.width,
           density=-1,col = 1)
    }
  }
  
  ### intron
  for(i in 1:length(Introns.start.pos)){
    for(j in 1:length(Introns.start.pos[[i]])){
      rect(xleft=Introns.start.pos[[i]][j], xright=Introns.start.pos[[i]][j],
           ytop=pos.coor+(intron.width/2), ybottom=pos.coor-(intron.width/2),
           density=-1,col = 1)
    }
  }
  for(i in 1:length(Introns.start.neg)){
    for(j in 1:length(Introns.start.neg[[i]])){
      rect(xleft=Introns.start.neg[[i]][j], xright=Introns.start.neg[[i]][j],
           ytop=neg.coor+(intron.width/2), ybottom=neg.coor-(intron.width/2),
           density=-1,col = 1)
    }
  }
}



reuben.biplot <- function (x, y, var.axes = TRUE, col, cex = rep(par("cex"), 2), 
                           xlabs = NULL, ylabs = NULL, expand = 1, xlim = NULL, ylim = NULL, 
                           arrow.len = 0.1, main = NULL, sub = NULL, xlab = NULL, ylab = NULL, x.col = 1, y.col = 2,text.col = 1,
                           text.cex = 1,arrow.lwd = 1, ratio = NULL, ...) {
  n <- nrow(x)
  p <- nrow(y)
  if (missing(xlabs)) {
    xlabs <- dimnames(x)[[1L]]
    if (is.null(xlabs)) 
      xlabs <- 1L:n
  }
  xlabs <- as.character(xlabs)
  dimnames(x) <- list(xlabs, dimnames(x)[[2L]])
  if (missing(ylabs)) {
    ylabs <- dimnames(y)[[1L]]
    if (is.null(ylabs)) 
      ylabs <- paste("Var", 1L:p)
  }
  ylabs <- as.character(ylabs)
  dimnames(y) <- list(ylabs, dimnames(y)[[2L]])
  if (length(cex) == 1L) 
    cex <- c(cex, cex)
  if (missing(col)) {
    col <- par("col")
    if (!is.numeric(col)) 
      col <- match(col, palette(), nomatch = 1L)
    col <- c(col, col + 1L)
  }
  else if (length(col) == 1L) 
    col <- c(col, col)
  unsigned.range <- function(x) c(-abs(min(x, na.rm = TRUE)), 
                                  abs(max(x, na.rm = TRUE)))
  rangx1 <- unsigned.range(x[, 1L])
  rangx2 <- unsigned.range(x[, 2L])
  rangy1 <- unsigned.range(y[, 1L])
  rangy2 <- unsigned.range(y[, 2L])
  if (missing(xlim) && missing(ylim)) 
    xlim <- ylim <- rangx1 <- rangx2 <- range(rangx1, rangx2)
  else if (missing(xlim)) 
    xlim <- rangx1
  else if (missing(ylim)) 
    ylim <- rangx2
  if (is.null(ratio)) {
    ratio <- max(rangy1/rangx1, rangy2/rangx2)/expand
  }
  on.exit(par(op))
  op <- par(pty = "s")
  if (!is.null(main)) 
    op <- c(op, par(mar = par("mar") + c(0, 0, 1, 0)))
  plot(
    x, 
    type = "p", 
    xlim = xlim, 
    ylim = ylim, 
    col = x.col, 
    xlab = xlab, 
    ylab = ylab, 
    sub = sub, 
    main = main,
    pch = 16,
    cex = cex, 
    ...)
  #    text(
  #      x, 
  #    	xlabs, 
  #    	cex = cex[1L], 
  #    	col = x.col, 
  #    	...)
  par(new = TRUE)
  dev.hold()
  on.exit(dev.flush(), add = TRUE)
  plot(
    y, 
    axes = FALSE, 
    type = "n", 
    xlim = xlim * ratio, 
    ylim = ylim * ratio, 
    xlab = "", 
    ylab = "", 
    col = y.col,
    ...)
  # axis(3, col = col[2L], ...)
  #  axis(4, col = col[2L], ...)
  box(col = col[1L])
  text(y, 
       labels = ylabs, 
       cex = text.cex, 
       col = text.col, 
       ...)
  if (var.axes) 
    arrows(0, 0, y[, 1L] * 0.8, y[, 2L] * 0.8, col = y.col, 
           length = arrow.len,
           lwd = arrow.lwd)
  invisible()
}



resetPar <- function() {
  dev.new()
  op <- par(no.readonly = TRUE)
  dev.off()
  op
}


# defult bin sort is repeats
binSort <- function(repList, bins, TE.names, repType){
  bin.gr <- GRanges(seqnames=Rle(bins$chr),
                    ranges = IRanges(start=bins$start, end = bins$end - 1))
  
  for(te in 1:length(TE.names)){
    print(TE.names[te])
    if(repType[te] == "repeats"){
      te.gr <- repList[[TE.names[te]]]
      GR <- GRanges(seqnames = Rle(te.gr$genoName),
                  ranges = IRanges(start = te.gr$genoStart, end = te.gr$genoEnd))
    }else if(repType[te] == "chromatin"){
      te.gr <- repList[[TE.names[te]]]
      GR <- GRanges(seqnames = Rle(te.gr$chrom),
                    ranges = IRanges(start = te.gr$chromStart, end = te.gr$chromEnd))
    }else if(repType[te] == "GRange"){
      GR <- repList[te]
    }else{
      stop("no repType")
    }
    
    
    # here we can do repeat coverage 
    OL <- intersect(x=bin.gr, y=GR)
    OL.f <- as.matrix(findOverlaps(bin.gr,OL))
    if(length(OL.f) > 0){
      OL.agg <- aggregate( x= width(OL[OL.f[,2]]) , by= list(OL.f[,1]), FUN=sum)
      bins[OL.agg[,1],TE.names[te]] <- OL.agg[,2]
      bins[is.na(bins[,TE.names[te]]),TE.names[te]] <- 0
    }else{
      bins[,TE.names[te]] <- 0
    }
    
  }
  
  # it might be a good idea to make sure to sort out the seqnames issue first
  bin.counts <- bins
  bin.rates <- bins
  for(i in 1:length(TE.names)){
    bin.rates[,TE.names[i]] <- (bin.rates[,TE.names[i]] * 10000) / bin.rates$Known
  }
  
  return(list(counts = bin.counts, rates = bin.rates))
}




removeRepAlign <- function(refSpecGenome, queSpecGenome, refSpec, queSpec){
  
  web <- paste("http://hgdownload.soe.ucsc.edu/goldenPath/", queSpecGenome, "/database/chromInfo.txt.gz", sep = "")
  con <- gzcon(url(web))
  txt <- readLines(con)
  chrom_info <- read.delim(textConnection(txt), header = FALSE)
  chrom_infoQue <- chrom_info[,2]
  names(chrom_infoQue) <- chrom_info[,1]
  
  
  align <- read.table(paste("~/Desktop/Domain_manuscript/Data/usable_alignment/",refSpecGenome,".",queSpecGenome,".axt.txt", sep = ""), 
                      colClasses= c("character", "character", "numeric", "numeric", "character", "numeric", "numeric", "character", "numeric"))
  colnames(align) <- c("Alignment_number", "chrRef", "startRef", "endRef", "chrQue", "startQue", "endQue", "strand", "blastzScore")
  
  # remove unplaced chromosomes
  if(length(grep("_", align$chrQue)) > 0){
    align <- align[-(grep("_", align$chrQue)),]
  }
  
  if(length(grep("_", align$chrRef)) > 0){
    align <- align[-(grep("_", align$chrRef)),]
  }
  # neg strand converiosn 
  
  align[align$strand == "-","startQue"] <-  chrom_infoQue[align[align$strand == "-","chrQue"]] - align[align$strand == "-","startQue"]
  tmp <-  chrom_infoQue[align[align$strand == "-","chrQue"]] - align[align$strand == "-","endQue"]
  align[align$strand == "-","endQue"] <- align[align$strand == "-","startQue"]
  align[align$strand == "-","startQue"] <- tmp
  
  # remove all multi mapping regions
  Que.mm.gr <- GRanges(seqnames=Rle(align$chrQue), 
                       ranges=IRanges(start=align$startQue, end = align$endQue)
  )
  mm.ol <- as.matrix(findOverlaps(Que.mm.gr))
  # perform self overlaps to identify if regions overlap themselves
  multiQue <- as.matrix(findOverlaps(Que.mm.gr, Que.mm.gr))
  # however species 2 does
  # therefore regions in species 1 are mapping to the same place in species two
  # non multis
  singleQue <- multiQue[!(multiQue[,1] %in% unique(multiQue[duplicated(multiQue[,1]),1])),1]
  align <- align[singleQue,]
  rownames(align) <- 1:nrow(align)
  
  refRep <- NULL
  refRep_list <- get(paste(refSpec, "_repInfo", sep = ""))
  for(i in 1:length(refRep_list)){
    refRep <- rbind(refRep, refRep_list[[i]])
  }
  refRep.gr <- GRanges(seqnames=Rle(refRep$genoName),
                       ranges = IRanges(start=refRep$genoStart, end = refRep$genoEnd)
  )
  refRep.gr <- reduce(refRep.gr)
  
  aliRef.gr <- GRanges(seqnames=Rle(align$chrRef),
                       ranges = IRanges(start=align$startRef, end=align$endRef- 1)
  )
  
  if(length(aliRef.gr) == length(reduce(aliRef.gr))){
    print("no overlapping alignment fragments")
  }
  
  # how do we translate coordinates 
  sd <- setdiff(aliRef.gr, refRep.gr)
  ol.sd <- as.matrix(findOverlaps(aliRef.gr, sd))
  comp <- data.frame(align[ol.sd[,1],], as.data.frame(sd)[ol.sd[,2],])
  
  # make sure all our broken fragments fit within the correct area
  if(!all(comp$start >= comp$startRef)){
    stop("New ref start alignments are outside of old ones")
  }
  if(!all(comp$end <= comp$endRef)){
    stop("New end start alignments are outside of old ones")
  }
  # the difference between new start and old start as a portion of the length 
  comp$newStartQue <- round(((comp$start - comp$startRef)/(comp$endRef - comp$startRef + 1)) * (comp$endQue - comp$startQue + 1)) + comp$startQue
  comp$newEndQue <- round(((comp$end - comp$endRef)/(comp$endRef - comp$startRef + 1)) * (comp$endQue - comp$startQue + 1)) + comp$endQue
  
  
  
  # there is an inequality between the lengths of the aligning fragments
  # we are remving parts of the aligning fragments, we awant to make sure we do not disturb these inequalitites. 
  
  # the fraction that has been removed by repeats should be equall across both datasets
  print("proportion ovrerlapping repeats in Ref")
  print((sum(as.numeric(comp$endRef-comp$startRef)) - 
           sum(as.numeric(comp$end-comp$start))) / 
          sum(as.numeric(comp$endRef-comp$startRef)))
  
  print("proportion overlapping repeats in Que")
  print((sum(as.numeric(comp$endQue-comp$startQue)) - 
           sum(as.numeric(comp$newEndQue-comp$newStartQue))) / 
          sum(as.numeric(comp$endQue-comp$startQue)))
  
  print("cor between repeat potions removed from alignmnets")
  print(cor(((comp$endRef-comp$startRef) - (comp$end-comp$start)) / (comp$endRef-comp$startRef), 
            ((comp$endQue-comp$startQue) - (comp$newEndQue-comp$newStartQue)) / (comp$endQue-comp$startQue)))
  
  ## fractions maintained
  new.align <- comp[,c("Alignment_number", "chrRef", "start", "end", "chrQue", "newStartQue", "newEndQue", "strand", "blastzScore")]
  colnames(new.align) <- colnames(align)
  
  new.align <- new.align[new.align$endRef - new.align$startRef > 0,]
  new.align <- new.align[new.align$endQue - new.align$startQue > 0,]
  rownames(new.align) <- 1:nrow(new.align)
  
  queRep <- NULL
  queRep_list <- get(paste(queSpec, "_repInfo", sep = ""))
  for(i in 1:length(queRep_list)){
    queRep <- rbind(queRep, queRep_list[[i]])
  }
  queRep.gr <- GRanges(seqnames=Rle(queRep$genoName),
                       ranges = IRanges(start=queRep$genoStart, end = queRep$genoEnd)
  )
  queRep.gr <- reduce(queRep.gr)
  aliQue.gr <- GRanges(seqnames=Rle(new.align$chrQue),
                       ranges = IRanges(start=new.align$startQue, end=new.align$endQue - 1)
  )
  
  if(length(reduce(aliQue.gr)) == length(aliQue.gr)){
    print("no overlapping alignmnets")
  }
  # how do we translate coordinates 
  sd <- setdiff(aliQue.gr, queRep.gr)
  ol.sd <- as.matrix(findOverlaps(aliQue.gr, sd))
  comp <- data.frame(new.align[ol.sd[,1],], as.data.frame(sd)[ol.sd[,2],])
  if(!all(comp$start >= comp$startQue)){
    stop("New que start alignments are outside of old ones")
  }
  if(!all(comp$end <= comp$endQue)){
    stop("New que end alignments are outside of old ones")
  }
  # the difference between new start and old start as a portion of the length 
  comp$newStartRef <- round(((comp$start - comp$startQue)/(comp$endQue - comp$startQue + 1)) * (comp$endRef - comp$startRef + 1)) + comp$startRef
  comp$newEndRef <- round(((comp$end - comp$endQue)/(comp$endQue - comp$startQue + 1)) * (comp$endRef - comp$startRef + 1)) + comp$endRef
  
  print("fraction of repeats removed from Que")
  print(
    (sum(as.numeric(comp$endQue-comp$startQue)) - sum(as.numeric(comp$end-comp$start))) / sum(as.numeric(comp$endQue-comp$startQue))
  )
  
  print("fraction of repeats removed from Ref")
  print(
    (sum(as.numeric(comp$endRef-comp$startRef)) - sum(as.numeric(comp$newEndRef-comp$newStartRef))) / sum(as.numeric(comp$endRef-comp$startRef))
  )
  
  print("corelation of fractions of repeats removed")
  print(
    cor((((comp$endRef-comp$startRef) - (comp$newEndRef-comp$newStartRef)) / (comp$endRef-comp$startRef)), (((comp$endQue-comp$startQue) - (comp$end-comp$start)) / (comp$endQue-comp$startQue)))
  )
  # amybe need to do some analysis on the residual fractions left over
  
  #res_frac <- (((comp$endRef-comp$startRef) - (comp$newEndRef-comp$newStartRef)) / (comp$endRef-comp$startRef)) -  (((comp$endQue-comp$startQue) - (comp$end-comp$start)) / (comp$endQue-comp$startQue))
  
  # it tends to be small things that have the biggest fractional differences
  # maybe if we remove all our zeroed alignmnets we will be done 
  
  new.new.align <- comp[,c("Alignment_number", "chrRef", "newStartRef", "newEndRef", "chrQue", "start", "end", "strand", "blastzScore")]
  colnames(new.new.align) <- colnames(align)
  
  new.new.align <- new.new.align[new.new.align$endRef - new.new.align$startRef > 0,]
  new.new.align <- new.new.align[new.new.align$endQue - new.new.align$startQue > 0,]
  rownames(new.new.align) <- 1:nrow(new.new.align)
  
  return(new.new.align)
}



isolateBinAlign <- function(align, refSpec, queSpec){
  refPCA <- get(paste(refSpec, "PCA", sep = ""))
  refBin <- refPCA$binInfo
  refBin.gr <- GRanges(seqnames=Rle(refBin$chr), 
                       ranges=IRanges(start=refBin$start, end = refBin$end -1 )
  )
  
  refAli.gr <- GRanges(seqnames=Rle(align$chrRef),
                       ranges = IRanges(start=align$startRef, end = align$endRef)
  )
  
  intRef <- intersect(refBin.gr, refAli.gr)
  olRef <- as.matrix(findOverlaps(refAli.gr, intRef)) 
  comp <- data.frame(align[olRef[,1], ], as.data.frame(intRef)[olRef[,2],])
  
  # make sure all our broken fragments fit within the correct area
  if(!all(comp$start >= comp$startRef)){
    stop("New ref start alignments are outside of old ones")
  }
  if(!all(comp$end <= comp$endRef)){
    stop("New ref end alignments are outside of old ones")
  }
  # the difference between new start and old start as a portion of the length 
  comp$newStartQue <- round(((comp$start - comp$startRef)/(comp$endRef - comp$startRef + 1)) * (comp$endQue - comp$startQue + 1)) + comp$startQue
  comp$newEndQue <- round(((comp$end - comp$endRef)/(comp$endRef - comp$startRef + 1)) * (comp$endQue - comp$startQue + 1)) + comp$endQue
  
  # it looks 
  
  print("proportion of Ref bins outside alignmnets")
  print((sum(as.numeric(comp$endRef-comp$startRef)) - 
           sum(as.numeric(comp$end-comp$start))) / 
          sum(as.numeric(comp$endRef-comp$startRef)))
  
  print("proportion of Que bins outside alignmnets")
  print((sum(as.numeric(comp$endQue-comp$startQue)) - 
           sum(as.numeric(comp$newEndQue-comp$newStartQue))) / 
          sum(as.numeric(comp$endQue-comp$startQue)))
  
  print("cor between repeat potions removed from alignmnets")
  print(cor(((comp$endRef-comp$startRef) - (comp$end-comp$start)) / (comp$endRef-comp$startRef), 
            ((comp$endQue-comp$startQue) - (comp$newEndQue-comp$newStartQue)) / (comp$endQue-comp$startQue)))
  
  ## fractions maintained
  new.align <- comp[,c("Alignment_number", "chrRef", "start", "end", "chrQue", "newStartQue", "newEndQue", "strand", "blastzScore")]
  colnames(new.align) <- colnames(align)
  
  ####
  #
  #
  #
  
  
  quePCA <- get(paste(queSpec, "PCA", sep = ""))
  queBin <- quePCA$binInfo
  queBin.gr <- GRanges(seqnames=Rle(queBin$chr), 
                       ranges=IRanges(start=queBin$start, end = queBin$end -1 )
  )
  
  queAli.gr <- GRanges(seqnames=Rle(new.align$chrQue),
                       ranges = IRanges(start=new.align$startQue, end = new.align$endQue - 1)
  )
  
  intQue <- intersect(queBin.gr, queAli.gr)
  olQue <- as.matrix(findOverlaps(queAli.gr, intQue)) 
  comp <- data.frame(new.align[olQue[,1], ], as.data.frame(intQue)[olQue[,2],])
  
  if(!all(comp$start >= comp$startQue)){
    stop("New que start alignments are outside of old ones")
  }
  if(!all(comp$end <= comp$endQue)){
    stop("New que end alignments are outside of old ones")
  }
  
  # the difference between new start and old start as a portion of the length 
  comp$newStartRef <- round(((comp$start - comp$startQue)/(comp$endQue - comp$startQue + 1)) * (comp$endRef - comp$startRef + 1)) + comp$startRef
  comp$newEndRef <- round(((comp$end - comp$endQue)/(comp$endQue - comp$startQue + 1)) * (comp$endRef - comp$startRef + 1)) + comp$endRef
  
  print("fraction of repeats removed from Que")
  print(
    (sum(as.numeric(comp$endQue-comp$startQue)) - sum(as.numeric(comp$end-comp$start))) / sum(as.numeric(comp$endQue-comp$startQue))
  )
  
  print("fraction of repeats removed from Ref")
  print(
    (sum(as.numeric(comp$endRef-comp$startRef)) - sum(as.numeric(comp$newEndRef-comp$newStartRef))) / sum(as.numeric(comp$endRef-comp$startRef))
  )
  
  print("corelation of fractions of repeats removed")
  print(
    cor((((comp$endRef-comp$startRef) - (comp$newEndRef-comp$newStartRef)) / (comp$endRef-comp$startRef)), (((comp$endQue-comp$startQue) - (comp$end-comp$start)) / (comp$endQue-comp$startQue)))
  )
  
  
  # so this will find perfectly overlapping fragments between the bins and the alignmnets
  
  new.new.align <- comp[,c("Alignment_number", "chrRef", "newStartRef", "newEndRef", "chrQue", "start", "end", "strand", "blastzScore")]
  colnames(new.new.align) <- colnames(align)
  
  new.new.align <- new.new.align[new.new.align$endRef - new.new.align$startRef > 0,]
  new.new.align <- new.new.align[new.new.align$endQue - new.new.align$startQue > 0,]
  rownames(new.new.align) <- 1:nrow(new.new.align)
  
  return(new.new.align)
  
}


buildBinMap <- function(align, refSpec, queSpec){
  refPCA <- get(paste(refSpec, "PCA", sep = ""))
  refBin <- refPCA$binInfo
  refBin.gr <- GRanges(seqnames=Rle(refBin$chr), 
                       ranges=IRanges(start=refBin$start, end = refBin$end -1 )
  )
  refAli.gr <- GRanges(seqnames=Rle(align$chrRef),
                       ranges = IRanges(start=align$startRef, end = align$endRef -1)
  )
  refOL <- as.matrix(findOverlaps(refBin.gr, refAli.gr))
  if(nrow(refOL) == nrow(align)){
    print("no overlapping sequence")
  }
  
  quePCA <- get(paste(queSpec, "PCA", sep = ""))
  queBin <- quePCA$binInfo
  queBin.gr <- GRanges(seqnames=Rle(queBin$chr), 
                       ranges=IRanges(start=queBin$start, end = queBin$end -1 )
  )
  queAli.gr <- GRanges(seqnames=Rle(align$chrQue),
                       ranges = IRanges(start=align$startQue, end = align$endQue - 1)
  )
  queOL <- as.matrix(findOverlaps(queBin.gr, queAli.gr))
  if(nrow(queOL) == nrow(align)){
    print("no overlapping sequence")
  }
  
  mergeAlignBins <- merge(refOL, queOL, by = 2)
  colnames(mergeAlignBins) <- c("alignNo", "refNo", "queNo")
  mergeAlignBins$alignGroup <- as.factor(paste(mergeAlignBins$refNo, mergeAlignBins$queNo, sep = "_"))
  nrow(mergeAlignBins) == nrow(align)
  
  # all the information is there, just need to sum it up 
  aggAlignRef <- aggregate(x=align$endRef[mergeAlignBins$alignNo] - align$startRef[mergeAlignBins$alignNo] + 1, by=list(mergeAlignBins$alignGroup), FUN=sum)
  colnames(aggAlignRef) <- c("alignGroup", "refWidth")
  aggAlignRef$refNo <- mergeAlignBins$refNo[match(aggAlignRef$alignGroup, table=mergeAlignBins$alignGroup)]
  aggAlignRef$refFrac <- aggAlignRef$refWidth/refBin$Known[aggAlignRef$refNo]
  
  aggAlignQue <- aggregate(x=align$endQue[mergeAlignBins$alignNo] - align$startQue[mergeAlignBins$alignNo] + 1, by=list(mergeAlignBins$alignGroup), FUN=sum)
  colnames(aggAlignQue) <- c("alignGroup", "queWidth")
  aggAlignQue$queNo <- mergeAlignBins$queNo[match(aggAlignQue$alignGroup, table=mergeAlignBins$alignGroup)]
  aggAlignQue$queFrac <- aggAlignQue$queWidth/queBin$Known[aggAlignQue$queNo]
  
  BinMapping <- merge(aggAlignRef, aggAlignQue, by = "alignGroup")
  BinMapping <- BinMapping[,c("refNo", "refFrac","refWidth", "queNo", "queFrac", "queWidth")]
  
  return(BinMapping)
}


### Takes a bin mapping of two species genome and remodels the query species according to the reference 
### plot also takes a cutoff value that decided the minnimum representation of a bin 


mapRemodeler <- function(refSpec, queSpec, cutoff, PCs = c("ancient_PC", "new_SINE_PC")){
  
  analysis <- get(paste(refSpec, "Ref_", queSpec,"Que", sep = ""))
  BinMap <- analysis$binMap
  
  aggRefFrac <- aggregate(x = BinMap$refFrac, by = list(BinMap$refNo), FUN = sum)
  aggQueFrac <- aggregate(x = BinMap$queFrac, by = list(BinMap$queNo), FUN = sum)
  rmRbins <- aggRefFrac$Group.1[aggRefFrac$x < cutoff]
  rmQbins <- aggQueFrac$Group.1[aggQueFrac$x < cutoff]
  
  BinMap <- BinMap[!(BinMap$refNo %in% rmRbins),]
  BinMap <- BinMap[!(BinMap$queNo %in% rmQbins),]
  
  
  # so we need to look at the binMap differntly 
  # remove bins form the analysis in which only a small fraction of them map, instead of removing small fractions of bins
  
  metaList = NULL
  
  for(pc in PCs){
    quePC <- analysis[[paste(queSpec, "Que", sep = "")]]$x[,pc]
    refPC <- analysis[[paste(refSpec, "Ref", sep = "")]]$x[,pc]
    
    
    permutePC <- matrix(data=NA, nrow=length(quePC), ncol=1000)
    for(i in 1:1000){
      permutePC[,i] <- sample(quePC,replace=F, size=length(quePC))
    }
    
    quePC <- data.frame(real = quePC, permutePC)
    
    
    agg.quePC <- aggregate(quePC[BinMap$queNo,] * BinMap$queFrac, by=list(BinMap$refNo), FUN = sum)
    agg.refFrac <- aggregate(x=BinMap$refFrac, by=list(BinMap$refNo), FUN = sum)
    
    all(agg.quePC$Group.1 == agg.refFrac$Group.1)
    
    remodeldPC <- data.frame(refNo = agg.quePC$Group.1, 
                             refPC = refPC[agg.quePC$Group.1],
                             quePC = agg.quePC$real/agg.refFrac$x
    )
    remodeldPermutePC <- data.frame(refNo = agg.quePC$Group.1,
                                    refPC = refPC[agg.quePC$Group.1], 
                                    agg.quePC[,3:ncol(agg.quePC)]/agg.refFrac$x
    )
    
    corDist <- NULL
    for(i in 1:1000){
      corDist <- c(corDist, cor(remodeldPermutePC[,2], remodeldPermutePC[,2 + i]))
    }
    
    ksRes <- data.frame(ref_cutoff = ks.test(refPC,y=refPC[unique(BinMap$refNo)])$p.value,
                        que_cutoff = ks.test(remodeldPC$quePC,y= quePC$real[unique(BinMap$queNo)])$p.value,
                        que_initial = ks.test(remodeldPC$quePC,y= quePC$real)$p.value
    )
    resultList <- list(remodeldPC = remodeldPC, ksRes = ksRes, corDist = corDist)
    metaList <- c(metaList,list(resultList))
  }
  
  specs <- data.frame(ref = refSpec, que = queSpec)
  metaList <- c(metaList, list(unique(BinMap$queNo)), list(cutoff), list(specs))
  names(metaList) <- c(PCs, "queNo", "cutoff", "specs")
  
  return(metaList)
  
  
  
}


covCalcPlot <- function(lenChoice, repChoice, repBins , repList ,minRepCov = NULL, maxRepCov = NULL, minRepSize = NULL, maxRepSize = NULL, minBinSize = NULL, maxBinSize = NULL){
  
  
  repBins <- repBins[repBins[,repChoice] > 0, c("chr", "start", "end", "Known", repChoice)]
  
  repBins.gr <- GRanges(seqnames = Rle(repBins$chr), 
                        ranges = IRanges(start = repBins$start, end = repBins$end))
  
  if(!is.null(maxBinSize)){
    repBins <- repBins[repBins$end - repBins$start + 1 < maxBinSize,]
  }
  if(!is.null(minBinSize)){
    repBins <- repBins[repBins$end - repBins$start + 1 > minBinSize,]
  }
  repBins2 <- repBins
  
  colnames(repBins)[5] = "repCov"
  if(!is.null(maxRepCov)){
    repBins <- repBins[repBins$repCov < maxRepCov,]
  }
  if(!is.null(minRepCov)){
    repBins <- repBins[repBins$repCov > minRepCov,]
  }
  
  repGR <- GRanges(seqnames=Rle(repList[[repChoice]]$genoName),
                   ranges = IRanges(start=repList[[repChoice]]$genoStart, end=repList[[repChoice]]$genoEnd)
  )
  
  if(!is.null(maxRepSize)){
    repGR <- repGR[width(repGR) < maxRepSize,]
  }
  if(!is.null(minRepSize)){
    repGR <- repGR[width(repGR) > minRepSize,]
  }
  seqLen <- lenChoice
  names(seqLen) <- "seq"
  
  us.ends <- repBins$start + ((repBins$end - repBins$start)/2)
  us.ends[us.ends - repBins$start+1 > lenChoice] <- repBins$start[us.ends - repBins$start+1 > lenChoice] + lenChoice
  us.bins <- data.frame(chr = repBins$chr, start = repBins$start, end = us.ends)
  us.bins.gr <- GRanges(seqnames=Rle(us.bins$chr), 
                        ranges = IRanges(start=us.bins$start, end = us.bins$end))
  us.rep.int <- intersect(repGR, us.bins.gr)
  us.Ol <- as.matrix(findOverlaps(us.bins.gr, us.rep.int))
  us.cov <- data.frame(start = start(us.rep.int[us.Ol[,2]]) - start(us.bins.gr[us.Ol[,1]]) + 1, 
                       end =end(us.rep.int[us.Ol[,2]]) - start(us.bins.gr[us.Ol[,1]]) + 1)
  
  us.cov.r <- GRanges(seqnames=Rle("seq"), ranges = IRanges(start=us.cov$start, end = us.cov$end), seqlengths = seqLen + 1)
  
  
  ds.starts <- as.integer(repBins$end - ((repBins$end - repBins$start)/2))
  ds.starts[repBins$end - ds.starts + 1> lenChoice] <- repBins$end[repBins$end - ds.starts + 1> lenChoice] - lenChoice
  ds.bins <- data.frame(chr = repBins$chr, start = ds.starts, end = repBins$end)
  ds.bins.gr <- GRanges(seqnames=Rle(ds.bins$chr), 
                        ranges = IRanges(start=ds.bins$start, end = ds.bins$end))
  ds.rep.int <- intersect(repGR, ds.bins.gr)  
  ds.Ol <- as.matrix(findOverlaps(ds.bins.gr, ds.rep.int))
  ds.cov <- data.frame(start =   (lenChoice - (end(ds.bins.gr[ds.Ol[,1]]) - start(ds.bins.gr[ds.Ol[,1]]))) + (start(ds.rep.int[ds.Ol[,2]]) - start(ds.bins.gr[ds.Ol[,1]])) + 1 ,
                       end =  (lenChoice - (end(ds.bins.gr[ds.Ol[,1]]) - start(ds.bins.gr[ds.Ol[,1]]))) + (end(ds.rep.int[ds.Ol[,2]]) - start(ds.bins.gr[ds.Ol[,1]])) + 1)
  ds.cov.r <- GRanges(seqnames=Rle("seq"),ranges = IRanges(start=ds.cov$start, end=ds.cov$end),seqlengths = seqLen + 1)
  
  
  # we can get our base frequency another way
  
  us.bins.bf <- as.integer(coverage(IRanges(start = 1, width = width(us.bins.gr))))
  ds.bins.bf <- as.integer(coverage(IRanges(end = lenChoice +1, width = width(ds.bins.gr))))
  
  
  
  # totCov <- as.numeric(coverage(ds.cov.r)$seq[(lenChoice+1):1]) + as.numeric(coverage(us.cov.r)$seq)
  int <- intersect(repBins.gr, repGR)
  p = sum(width(int))/sum(as.numeric(width(repBins.gr)))
  
  output = list(rawRepCov5 = as.integer(coverage(ds.cov.r)$seq), rawRepCov3 = as.integer(coverage(us.cov.r)$seq), p = p, baseFreq5 = ds.bins.bf, baseFreq3 = us.bins.bf)
  return(output)
  
}


# about getting GC content 
# can now also take chromatin
# a chromatin replist is the chromatin of one cell type
# each element in the list corresponds to a particular chromatin type 

covCalcPlot5prime3prime <- function(lenChoice, repChoice, repBins , repList , 
                                    refgene , type ,repType ,minRepCov = NULL, maxRepCov = NULL, 
                                    minRepSize = NULL, maxRepSize = NULL, minBinSize = NULL, 
                                    maxBinSize = NULL){
  
  
  
  refgene.gr <- GRanges(seqnames=Rle(refgene[,3]), ranges = IRanges(start = refgene[,5], end = refgene[,6]))
  

  repBins <- repBins[, c("chr", "start", "end", "Known", repChoice)]
  
  if(!is.null(maxBinSize)){
    repBins <- repBins[repBins$end - repBins$start + 1 < maxBinSize,]
  }
  if(!is.null(minBinSize)){
    repBins <- repBins[repBins$end - repBins$start + 1 > minBinSize,]
  }
  colnames(repBins)[5] = "repCov"
  if(!is.null(maxRepCov)){
    repBins <- repBins[repBins$repCov < maxRepCov,]
  }
  if(!is.null(minRepCov)){
    repBins <- repBins[repBins$repCov > minRepCov,]
  }
  
  
  
  if(repType == "repeats"){
    repGR <- GRanges(seqnames=Rle(repList[[repChoice]]$genoName),
                     ranges = IRanges(start=repList[[repChoice]]$genoStart, end=repList[[repChoice]]$genoEnd)
    )
  }else if(repType == "chromatin"){
    repGR <- GRanges(seqnames=Rle(repList[[repChoice]]$chrom),
                     ranges = IRanges(start=repList[[repChoice]]$chromStart, end=repList[[repChoice]]$chromEnd)
    )
  }else{
    stop("no repType")
  }
  
  if(!is.null(maxRepSize)){
    repGR <- repGR[width(repGR) < maxRepSize,]
  }
  if(!is.null(minRepSize)){
    repGR <- repGR[width(repGR) > minRepSize,]
  }
  seqLen <- lenChoice + 1
  names(seqLen) <- "seq"
  
  # it might be easier if we break down the four groups at the start

  us.ends <- repBins$start + ((repBins$end - repBins$start)/2)
  us.ends[us.ends - repBins$start+1 > lenChoice] <- repBins$start[us.ends - repBins$start+1 > lenChoice] + lenChoice
  us.bins <- data.frame(chr = repBins$chr, start = repBins$start, end = us.ends)
  us.bins.gr <- GRanges(seqnames=Rle(us.bins$chr), 
                        ranges = IRanges(start=us.bins$start, end = us.bins$end))
 
  ds.starts <- as.integer(repBins$end - ((repBins$end - repBins$start)/2))
  ds.starts[repBins$end - ds.starts + 1> lenChoice] <- repBins$end[repBins$end - ds.starts + 1> lenChoice] - lenChoice
  ds.bins <- data.frame(chr = repBins$chr, start = ds.starts, end = repBins$end)
  ds.bins.gr <- GRanges(seqnames=Rle(ds.bins$chr), 
                        ranges = IRanges(start=ds.bins$start, end = ds.bins$end))
  
  if(type == "intergenic"){
    us.s.ol <- as.matrix(findOverlaps(us.bins.gr, refgene.gr, maxgap=1))
    if(!(length(unique(us.s.ol[,1])) == length(us.bins.gr))){
      stop("some intergenic regions have no upstream gene")
    }
    us.strand <- data.frame(intergenicID = us.s.ol[,1],strand = refgene[us.s.ol[,2], 4])
    us.plus <- unique(us.strand$intergenicID[us.strand$strand == "+"])
    us.minus <- unique(us.strand$intergenicID[us.strand$strand == "-"])
    
    if(length(c(us.plus[us.plus %in% us.minus],us.minus[us.minus %in% us.plus])) != 0){
      stop("strand confusion in upstream intergenic regions")
    }
    
    ds.s.ol <- as.matrix(findOverlaps(ds.bins.gr, refgene.gr, maxgap=1))
    if(!(length(unique(ds.s.ol[,1])) == length(ds.bins.gr))){
      stop("some intergenic regions have no downstream gene")
    }
    ds.strand <- data.frame(intergenicID = ds.s.ol[,1],strand = refgene[ds.s.ol[,2], 4])
    ds.plus <- unique(ds.strand$intergenicID[ds.strand$strand == "+"])
    ds.minus <- unique(ds.strand$intergenicID[ds.strand$strand == "-"])
    
    if(length(c(ds.plus[ds.plus %in% ds.minus],ds.minus[ds.minus %in% ds.plus])) != 0){
      stop("strand confusion in downstream intergenic regions")
    }
    
  }else if(type == "intron"){
    
    us.s.ol <- as.matrix(findOverlaps(us.bins.gr, refgene.gr, minoverlap=2))
    if(!(length(unique(us.s.ol[,1])) == length(us.bins.gr))){
      stop("some upstream intron regions don't belong to a gene")
    }
    us.strand <- data.frame(intronID = us.s.ol[,1],strand = refgene[us.s.ol[,2], 4])
    us.plus <- unique(us.strand$intronID[us.strand$strand == "+"])
    us.minus <- unique(us.strand$intronID[us.strand$strand == "-"])
    
    if(length(c(us.plus[us.plus %in% us.minus],us.minus[us.minus %in% us.plus])) != 0){
      stop("strand confusion in upstream intron region")
    }
    
    ds.s.ol <- as.matrix(findOverlaps(ds.bins.gr, refgene.gr, minoverlap=2))
    if(!(length(unique(ds.s.ol[,1])) == length(ds.bins.gr))){
      stop("some downstream intron regions don't belong to a gene")
    }
    ds.strand <- data.frame(intergenicID = ds.s.ol[,1],strand = refgene[ds.s.ol[,2], 4])
    ds.plus <- unique(ds.strand$intergenicID[ds.strand$strand == "+"])
    ds.minus <- unique(ds.strand$intergenicID[ds.strand$strand == "-"])
    
    if(length(c(ds.plus[ds.plus %in% ds.minus],ds.minus[ds.minus %in% ds.plus])) != 0){
      stop("strand confusion in downstream intron region")
    }
  }else{
    stop("type needs to be intron or intergenic")
  }
  ### everything that is minus we have to reverse the coordinates 
  # upstream region of intergenic in plus strand is the 3prime region of a gene
  #us.Ol[,1] tells us which bin the repeats belong to
  
  us.bins.gr_5 <- us.bins.gr[us.minus]
  us.rep.int_5 <- intersect(repGR, us.bins.gr_5)
  us.Ol_5 <- as.matrix(findOverlaps(us.bins.gr_5, us.rep.int_5))
  us.cov_5 <- data.frame(start = start(us.rep.int_5[us.Ol_5[,2]]) - start(us.bins.gr_5[us.Ol_5[,1]]) + 1, 
                         end =end(us.rep.int_5[us.Ol_5[,2]]) - start(us.bins.gr_5[us.Ol_5[,1]]) + 1)
  us.cov_5GR <- GRanges(seqnames = Rle("seq"), 
                        ranges = IRanges(start = us.cov_5$start, end = us.cov_5$end),
                        seqlengths = seqLen)
  us.binRange5 <- IRanges(start = 1 , width = width(us.bins.gr_5)) 
  
  us.bins.gr_3 <- us.bins.gr[us.plus]
  us.rep.int_3 <- intersect(repGR, us.bins.gr_3)
  us.Ol_3 <- as.matrix(findOverlaps(us.bins.gr_3, us.rep.int_3))
  us.cov_3 <- data.frame(start = start(us.rep.int_3[us.Ol_3[,2]]) - start(us.bins.gr_3[us.Ol_3[,1]]) + 1, 
                         end =end(us.rep.int_3[us.Ol_3[,2]]) - start(us.bins.gr_3[us.Ol_3[,1]]) + 1)
  us.cov_3GR <- GRanges(seqnames = Rle("seq"), 
                        ranges = IRanges(start = us.cov_3$start, end = us.cov_3$end),
                        seqlengths = seqLen)
  
  us.binRange3 <- IRanges(start = 1 , width = width(us.bins.gr_3)) 
  
  
  
  ds.bins.gr_5 <- ds.bins.gr[ds.plus]
  ds.rep.int_5 <- intersect(repGR, ds.bins.gr_5)
  ds.Ol_5 <- as.matrix(findOverlaps(ds.bins.gr_5, ds.rep.int_5))
  ds.cov_5 <- data.frame(start =   (lenChoice - (end(ds.bins.gr_5[ds.Ol_5[,1]]) - start(ds.bins.gr_5[ds.Ol_5[,1]]))) + (start(ds.rep.int_5[ds.Ol_5[,2]]) - start(ds.bins.gr_5[ds.Ol_5[,1]])) + 1 ,
                         end =  (lenChoice - (end(ds.bins.gr_5[ds.Ol_5[,1]]) - start(ds.bins.gr_5[ds.Ol_5[,1]]))) + (end(ds.rep.int_5[ds.Ol_5[,2]]) - start(ds.bins.gr_5[ds.Ol_5[,1]])) + 1)
  
  
  ds.cov_5GR <- GRanges(seqnames = Rle("seq"), 
                        ranges = IRanges(start = ds.cov_5$start, end = ds.cov_5$end),
                        seqlengths = seqLen)
  ds.binRange5 <- IRanges(end = lenChoice+1 , width = width(ds.bins.gr_5)) 
  
  ds.bins.gr_3 <- ds.bins.gr[ds.minus]
  ds.rep.int_3 <- intersect(repGR, ds.bins.gr_3)
  ds.Ol_3 <- as.matrix(findOverlaps(ds.bins.gr_3, ds.rep.int_3))
  ds.cov_3 <- data.frame(start =   (lenChoice - (end(ds.bins.gr_3[ds.Ol_3[,1]]) - start(ds.bins.gr_3[ds.Ol_3[,1]]))) + (start(ds.rep.int_3[ds.Ol_3[,2]]) - start(ds.bins.gr_3[ds.Ol_3[,1]])) + 1 ,
                       end =  (lenChoice - (end(ds.bins.gr_3[ds.Ol_3[,1]]) - start(ds.bins.gr_3[ds.Ol_3[,1]]))) + (end(ds.rep.int_3[ds.Ol_3[,2]]) - start(ds.bins.gr_3[ds.Ol_3[,1]])) + 1)
  
  ds.cov_3GR <- GRanges(seqnames = Rle("seq"), 
                        ranges = IRanges(start = ds.cov_3$start, end = ds.cov_3$end),
                        seqlengths = seqLen)
  
  ds.binRange3 <- IRanges(end = lenChoice+1 , width = width(ds.bins.gr_3)) 
  
  us.cov_5GRnum <- as.numeric(coverage(us.cov_5GR)$seq)[(lenChoice+1):1]
  ds.cov_5GRnum <- as.numeric(coverage(ds.cov_5GR)$seq)
  
  us.cov_3GRnum <- as.numeric(coverage(us.cov_3GR)$seq)
  ds.cov_3GRnum <- as.numeric(coverage(ds.cov_3GR)$seq)[(lenChoice+1):1]
  
  us.5bf <- as.numeric(coverage(us.binRange5))[(lenChoice+1):1]
  ds.5bf <- as.numeric(coverage(ds.binRange5))[1:(lenChoice+1)]
  us.3bf <- as.numeric(coverage(us.binRange3))[1:(lenChoice+1)]
  ds.3bf <- as.numeric(coverage(ds.binRange3))[(lenChoice+1):1]
  
  p = sum(repBins$repCov)/sum(repBins$end - repBins$start + 1)
  output = list(rawRepCov5 = ds.cov_5GRnum + us.cov_5GRnum,
                rawRepCov3 = ds.cov_3GRnum + us.cov_3GRnum,
                p = p, 
                baseFreq5prime = us.5bf + ds.5bf, 
                baseFreq3prime = us.3bf + ds.3bf
  )
  return(output)
  
}  




### takes out regions that don't have up and downstream or are both an up and downstream coordinate
filterIntergenic <- function(refgene){
  refgene.gr <- GRanges(seqnames=Rle(refgene[,3]), ranges = IRanges(start = refgene[,5], end = refgene[,6]))
  refgene_gap.gr <- gaps(refgene.gr)
  
  refgene_gap_startCH.gr <- refgene_gap.gr 
  end(refgene_gap_startCH.gr) <- end(refgene_gap.gr) - 1
  Sol <- as.matrix(findOverlaps(refgene_gap_startCH.gr, refgene.gr, maxgap=1))
  
  us.strand <- data.frame(intergenicID = Sol[,1],strand = refgene[Sol[,2], 4])
  plus <- unique(us.strand$intergenicID[us.strand$strand == "+"])
  minus <- unique(us.strand$intergenicID[us.strand$strand == "-"])
  strandConflict <- unique(c(plus[plus %in% minus] , minus[minus %in% plus]))
  
  refgene_gap_endCH.gr <- refgene_gap.gr 
  start(refgene_gap_endCH.gr) <- start(refgene_gap.gr) + 1
  Eol <- as.matrix(findOverlaps(refgene_gap_endCH.gr, refgene.gr, maxgap=1))
  
  ds.strand <- data.frame(intergenicID = Eol[,1],strand = refgene[Eol[,2], 4])
  plus <- unique(ds.strand$intergenicID[ds.strand$strand == "+"])
  minus <- unique(ds.strand$intergenicID[ds.strand$strand == "-"])
  strandConflict<- c(strandConflict, unique(c(plus[plus %in% minus] , minus[minus %in% plus])))
  
  gapKeep <- unique(c(Eol[Eol[,1] %in% Sol[,1],1],   Sol[Sol[,1] %in% Eol[,1],1]))
  gapKeep <- gapKeep[!(gapKeep %in% strandConflict)]
  
  refgene_gap.gr <- refgene_gap.gr[gapKeep]
  refgene_gap.gr <- refgene_gap.gr[-(grep("_", seqnames(refgene_gap.gr)))]
  return(refgene_gap.gr[width(refgene_gap.gr)>200])
}



# pulls out a set of both non overlapping and perfectly overlapping intronic regions
# strandedness is consistent between introns from multiple transcripts. 
# Introns overlapping any exons are also removed

filterIntron <- function(refgene){
  
  
  refgene.gr <- GRanges(seqnames=Rle(refgene[,3]), ranges = IRanges(start = refgene[,5], end = refgene[,6]))
  
  Echr <- NULL
  Ichr <- NULL
  Estart <- strsplit(x=as.character(refgene[,10]),split=",")
  Eend <- strsplit(x=as.character(refgene[,11]),split=",")
  Istart <- Eend
  Iend <- Estart
  for(i in 1:nrow(refgene)){
    if(refgene[i,9] > 1){
      Istart[[i]] = as.numeric(Istart[[i]][1:(length(Istart[[i]])-1)])
      Iend[[i]] = as.numeric(Iend[[i]][2:(length(Iend[[i]]))])
    }else{
      Istart[[i]] = NA
      Iend[[i]] = NA
    }
    Estart[[i]] = as.numeric(Estart[[i]])
    Eend[[i]] = as.numeric(Eend[[i]])
    Echr <- c(Echr, list(rep(as.character(refgene[i,3]), refgene[i,9])))
    Ichr <- c(Ichr, list(rep(as.character(refgene[i,3]), length(Istart[[i]]))))
  }
  
  Exons <- data.frame(chr = do.call(c, Echr), start = do.call(c,Estart), end = do.call(c,Eend))
  Introns <- data.frame(chr = do.call(c, Ichr), start = do.call(c,Istart), end = do.call(c,Iend))
  Introns <- Introns[!(is.na(Introns$start)),]                                        
  
  
  # intron filtering, probably could be a function
  
  intron.gr <- GRanges(seqnames = Rle(Introns$chr), 
                       ranges = IRanges(start=Introns$start, end=Introns$end)
  )
  I.OL <- as.matrix(findOverlaps(intron.gr))
  intronPull <- (1:length(intron.gr))[-unique(I.OL[duplicated(I.OL[,1]),1])]
  intronKeep.gr <- intron.gr[intronPull]
  intronReamin.gr <- intron.gr[unique(I.OL[duplicated(I.OL[,1]),1])]
  I.OL <- as.matrix(findOverlaps(intronReamin.gr, type = "equal"))
  eqOl.gr <- intronReamin.gr[unique(I.OL[duplicated(I.OL[,1]),1])]
  red <- reduce(eqOl.gr)
  ol2 <- as.matrix(findOverlaps(red, intronReamin.gr[-unique(I.OL[duplicated(I.OL[,1]),1])]))
  intronKeep.gr = c(intronKeep.gr, red[-unique(ol2[,1])])
  
  if(length(intronKeep.gr) == length(reduce(intronKeep.gr))){
    print("Single layer of intronic regions")
  }
  
  exon.gr <- GRanges(seqnames=Rle(Exons$chr),
                     ranges = IRanges(start=Exons$start, end = Exons$end))
  
  IolE <- as.matrix(findOverlaps(intronKeep.gr,exon.gr,minoverlap=2))
  
  intronKeep.gr <- intronKeep.gr[-unique(IolE[,1])]
  # strand consistancy 
  
  IolS <- as.matrix(findOverlaps(intronKeep.gr,refgene.gr,minoverlap=2))
  strandCon <- data.frame(intronID = IolS[,1],strand = refgene[IolS[,2], 4])
  plus <- unique(strandCon$intronID[strandCon$strand == "+"])
  minus <- unique(strandCon$intronID[strandCon$strand == "-"])
  strandConflict <- unique(c(plus[plus %in% minus] , minus[minus %in% plus]))
  if(length(strandConflict) > 0 ){
      intronKeep.gr <- intronKeep.gr[-strandConflict]
  }
  if(length(grep("_", seqnames(intronKeep.gr))) > 0){
    intronKeep.gr <- intronKeep.gr[-(grep("_", seqnames(intronKeep.gr)))]
  }
  return(intronKeep.gr[width(intronKeep.gr)>200])
  
}




#### maybe we don't have to remove our non zero bins. 
### Also we may be able to get a gnereal coverage amount
### This way we can account for the lack of TEs in regions relative to the genome. 

### This will be showing us where TEs accumulate in the genome 
### Rather than where do TEs accumulate in this region. 

# cahnge this function so we don't have to worry about the repeats 


covCalcPlot5prime3primeGC <- function(lenChoice, repBins  , refgene , type ,genome=genome ,
                                      minBinSize = NULL, maxBinSize = NULL){
  #### CG stats
  
  library(BSgenome)
  bsGenome <- available.genomes()[grep(genome,available.genomes())]
  bsGenome <- bsGenome[-(grep("masked", bsGenome))]
  library(bsGenome, character.only=TRUE)
  bsSpec <- get(strsplit(bsGenome,split="\\.")[[1]][2])
  
  
  refgene.gr <- GRanges(seqnames=Rle(refgene[,3]), ranges = IRanges(start = refgene[,5], end = refgene[,6]))
  
#  repBins <- repBins[repBins[,repChoice] > 0, c("chr", "start", "end", "Known", repChoice)]
  repBins <- repBins[, c("chr", "start", "end", "Known")]
  
  if(!is.null(maxBinSize)){
    repBins <- repBins[repBins$end - repBins$start + 1 < maxBinSize,]
  }
  if(!is.null(minBinSize)){
    repBins <- repBins[repBins$end - repBins$start + 1 > minBinSize,]
  }
  
  seqLen <- lenChoice
  names(seqLen) <- "seq"
  
  us.ends <- repBins$start + ((repBins$end - repBins$start)/2)
  us.ends[us.ends - repBins$start+1 > lenChoice] <- repBins$start[us.ends - repBins$start+1 > lenChoice] + lenChoice
  us.bins <- data.frame(chr = repBins$chr, start = repBins$start, end = us.ends)
  us.bins.gr <- GRanges(seqnames=Rle(us.bins$chr), 
                        ranges = IRanges(start=us.bins$start, end = us.bins$end))

  ds.starts <- as.integer(repBins$end - ((repBins$end - repBins$start)/2))
  ds.starts[repBins$end - ds.starts + 1> lenChoice] <- repBins$end[repBins$end - ds.starts + 1> lenChoice] - lenChoice
  ds.bins <- data.frame(chr = repBins$chr, start = ds.starts, end = repBins$end)
  ds.bins.gr <- GRanges(seqnames=Rle(ds.bins$chr), 
                        ranges = IRanges(start=ds.bins$start, end = ds.bins$end))
  
  if(type == "intergenic"){
    us.s.ol <- as.matrix(findOverlaps(us.bins.gr, refgene.gr, maxgap=1))
    if(!(length(unique(us.s.ol[,1])) == length(us.bins.gr))){
      stop("some intergenic regions have no upstream gene")
    }
    us.strand <- data.frame(intergenicID = us.s.ol[,1],strand = refgene[us.s.ol[,2], 4])
    us.plus <- unique(us.strand$intergenicID[us.strand$strand == "+"])
    us.minus <- unique(us.strand$intergenicID[us.strand$strand == "-"])
    
    if(length(c(us.plus[us.plus %in% us.minus],us.minus[us.minus %in% us.plus])) != 0){
      stop("strand confusion in upstream intergenic regions")
    }
    
    ds.s.ol <- as.matrix(findOverlaps(ds.bins.gr, refgene.gr, maxgap=1))
    if(!(length(unique(ds.s.ol[,1])) == length(ds.bins.gr))){
      stop("some intergenic regions have no downstream gene")
    }
    ds.strand <- data.frame(intergenicID = ds.s.ol[,1],strand = refgene[ds.s.ol[,2], 4])
    ds.plus <- unique(ds.strand$intergenicID[ds.strand$strand == "+"])
    ds.minus <- unique(ds.strand$intergenicID[ds.strand$strand == "-"])
    
    if(length(c(ds.plus[ds.plus %in% ds.minus],ds.minus[ds.minus %in% ds.plus])) != 0){
      stop("strand confusion in downstream intergenic regions")
    }
    
  }else if(type == "intron"){
    
    us.s.ol <- as.matrix(findOverlaps(us.bins.gr, refgene.gr, minoverlap=2))
    if(!(length(unique(us.s.ol[,1])) == length(us.bins.gr))){
      stop("some upstream intron regions don't belong to a gene")
    }
    us.strand <- data.frame(intronID = us.s.ol[,1],strand = refgene[us.s.ol[,2], 4])
    us.plus <- unique(us.strand$intronID[us.strand$strand == "+"])
    us.minus <- unique(us.strand$intronID[us.strand$strand == "-"])
    
    if(length(c(us.plus[us.plus %in% us.minus],us.minus[us.minus %in% us.plus])) != 0){
      stop("strand confusion in upstream intron region")
    }
    
    ds.s.ol <- as.matrix(findOverlaps(ds.bins.gr, refgene.gr, minoverlap=2))
    if(!(length(unique(ds.s.ol[,1])) == length(ds.bins.gr))){
      stop("some downstream intron regions don't belong to a gene")
    }
    ds.strand <- data.frame(intergenicID = ds.s.ol[,1],strand = refgene[ds.s.ol[,2], 4])
    ds.plus <- unique(ds.strand$intergenicID[ds.strand$strand == "+"])
    ds.minus <- unique(ds.strand$intergenicID[ds.strand$strand == "-"])
    
    if(length(c(ds.plus[ds.plus %in% ds.minus],ds.minus[ds.minus %in% ds.plus])) != 0){
      stop("strand confusion in downstream intron region")
    }
  }else{
    stop("type needs to be intron or intergenic")
  }
  # GC coverage calc
  All.seq<- DNAStringSet()
  prime5grDS <- ds.bins[ds.plus,]
  prime5grUS <- us.bins[us.minus,]
  
  chromos <- as.character(unique(repBins$chr))
  for(i in 1:length(chromos)){
    Seq.setDS=DNAStringSet(Hsapiens[[chromos[i]]], 
                           start=prime5grDS$start[prime5grDS$chr == chromos[i]], 
                           end=prime5grDS$end[prime5grDS$chr == chromos[i]])
    
    
    Seq.setUS=DNAStringSet(Hsapiens[[chromos[i]]], 
                           start=prime5grUS$start[prime5grUS$chr == chromos[i]], 
                           end=as.integer(prime5grUS$end[prime5grUS$chr == chromos[i]]))
    
    All.seq <- c(All.seq, reverse(Seq.setDS), Seq.setUS)
    print(paste("5 prime", chromos[i]))
  }
  
  gcPos <- c(unlist(vmatchPattern("C", All.seq)),unlist(vmatchPattern("G", All.seq)))
  prime5gcCov <- as.numeric(coverage(gcPos))
  prime5gcCov <- prime5gcCov[(lenChoice+1):1]
  
  
  All.seq<- DNAStringSet()
  prime3grDS <- ds.bins[ds.minus,]
  prime3grUS <- us.bins[us.plus,]
  
  chromos <- as.character(unique(repBins$chr))
  for(i in 1:length(chromos)){
    Seq.setDS=DNAStringSet(Hsapiens[[chromos[i]]], 
                           start=prime3grDS$start[prime3grDS$chr == chromos[i]], 
                           end=prime3grDS$end[prime3grDS$chr == chromos[i]])
    
    
    Seq.setUS=DNAStringSet(Hsapiens[[chromos[i]]], 
                           start=prime3grUS$start[prime3grUS$chr == chromos[i]], 
                           end=as.integer(prime3grUS$end[prime3grUS$chr == chromos[i]]))
    
    All.seq <- c(All.seq, reverse(Seq.setDS), Seq.setUS)
    print(paste("3 prime", chromos[i]))
  }
  
  gcPos <- c(unlist(vmatchPattern("C", All.seq)),unlist(vmatchPattern("G", All.seq)))
  prime3gcCov <- as.numeric(coverage(gcPos))
  
  
  ## gc pos might be wrong because it hasnt been end adjusted 
  
  
  
  #base frequency
  
  prime5Lengths <- c(width(us.bins.gr[us.minus]), width(ds.bins.gr[ds.plus]))
  prime3Lengths <- c(width(us.bins.gr[us.plus]), width(ds.bins.gr[ds.minus]))
  
  len.s <- sort(prime5Lengths)
  lenSum <- summary(as.factor(len.s), maxsum=length(unique(len.s)))
  CS <- cumsum(as.integer(lenSum)[length(lenSum):1])
  CS <- CS[length(CS):1]
  step <- stepfun(as.numeric(names(lenSum)[2:(length(lenSum))]), c(CS))
  prime5stepRes <- step(1:(lenChoice+1))[(lenChoice+1):1]
  
  len.s <- sort(prime3Lengths)
  lenSum <- summary(as.factor(len.s), maxsum=length(unique(len.s)))
  CS <- cumsum(as.integer(lenSum)[length(lenSum):1])
  CS <- CS[length(CS):1]
  step <- stepfun(as.numeric(names(lenSum)[2:(length(lenSum))]), c(CS))
  prime3stepRes <- step(1:(lenChoice+1))
  
  
  
  output = list(
                baseFreq5prime = prime5stepRes, 
                baseFreq3prime = prime3stepRes,
                prime5gc = prime5gcCov,
                prime3gc = prime3gcCov
  )
  return(output)
  
}  



LocalGClevel <- function(maxLen, repChoice,repBins,genome,repList){
  
  #Region_break <- repBins[repBins[, repChoice] > 0, c("chr", "start", "end")] 
  Region_break <- repBins[, c("chr", "start", "end")]
  print(c("Over max size", "max length", "mean length"))
  
  while(nrow(Region_break[Region_break$end - Region_break$start + 1 > maxLen,]) > 0){
    # pull long ones break them 
    Region_long <- Region_break[Region_break$end - Region_break$start > maxLen,]
    Region_break <- Region_break[!(Region_break$end - Region_break$start > maxLen),]
    Region_long2 <- Region_long1 <- Region_long
    longLengths <- (Region_long$end - Region_long$start + 1)
    cuts <- sapply(longLengths, FUN = sample, size = 1)
    if(length(cuts > 0)){
      Region_long1$end <- round(Region_long$start + cuts )
      Region_long2$start <- Region_long1$end + 1
      Region_break <- rbind(Region_break, Region_long1, Region_long2)
      print(c(nrow(Region_break[Region_break$end - Region_break$start > maxLen,]), 
              max(Region_break$end - Region_break$start),
              mean(Region_break$end - Region_break$start)
      )
      )
    }else{break()}
    
  }
  
  ### this is where we do the repeat overlaps and get the GC content 
  Region_break <- Region_break[Region_break$end - Region_break$start + 1 >1 & Region_break$end - Region_break$start + 1 < maxLen,]
  
  # now we can try and capture some local effects
  
  Intorn_break.gr <- GRanges(seqnames = Rle(Region_break$chr),
                             ranges = IRanges(start = Region_break$start, end = Region_break$end - 1))
  repsGR <- GRanges(seqnames = Rle(rep[[repChoice]]$genoName),
                    ranges = IRanges(start = rep[[repChoice]]$genoStart, end = rep[[repChoice]]$genoEnd))
  breakInt <- intersect(Intorn_break.gr, repsGR)
  breakOL <- as.matrix(findOverlaps(Intorn_break.gr, breakInt))
  breakCov <- aggregate(width(breakInt[breakOL[,2]]), by = list(breakOL[,1]), FUN = sum)
  Region_break$cov <- 0
  Region_break$cov[breakCov$Group.1] <- breakCov$x
  
  # I've got repeat coverage 
  # next is to collect GC content 
  
  
  library(BSgenome)
  bsGenome <- available.genomes()[grep(genome,available.genomes())]
  bsGenome <- bsGenome[-(grep("masked", bsGenome))]
  library(bsGenome, character.only=TRUE)
  bsSpec <- get(strsplit(bsGenome,split="\\.")[[1]][2])
  
  
  chromos <- as.character(unique(Region_break$chr))
  Region_break$GC <- rep(NA, nrow(Region_break))
  
  for(i in 1:length(chromos)){
    Seq.set=DNAStringSet(bsSpec[[chromos[i]]], 
                         start=Region_break$start[Region_break$chr == chromos[i]], 
                         end=Region_break$end[Region_break$chr == chromos[i]])
    bases=data.frame(alphabetFrequency(Seq.set,baseOnly=TRUE))
    CGcontent <- rowSums(bases[,c("G", "C")]) / rowSums(bases[,c("A", "C", "T", "G")])
    Region_break$GC[Region_break$chr == chromos[i]] <- CGcontent
    print(chromos[i])
  }
  
  
  Region_break$width <- Region_break$end - Region_break$start + 1
  
  #### here we get the final piece of the puzzle 
  return(Region_break)
  
}



localGCgenome <- function(repList, binSize, sampSize, genome, repType){
  
  library(BSgenome)
  bsGenome <- available.genomes()[grep(genome,available.genomes())]
  bsGenome <- bsGenome[-(grep("masked", bsGenome))]
  library(bsGenome, character.only=TRUE)
  bsSpec <- get(strsplit(bsGenome,split="\\.")[[1]][2])
  
  bins <- binned.genome.reader(genome=genome, bin.size=binSize, keep.rate=.95)
  bins <- bins[[1]]
  # remove unplaced chromosomes from bins
  if(length(grep("_", bins$chr)) > 0){
    bins <- bins[-(grep("_", bins$chr)),]
  }
  
  
  bin.samp <- bins[sample(1:nrow(bins),sampSize),]
  bin.sort = binSort(rep=repList, bins=bin.samp,TE.names=names(repList), repType = repType)
  bin.samp = bin.sort$rates
  bin.samp[,names(repList)] = bin.sort$rates[,names(repList)]/10000
  
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
  
  bin.samp$GC <- bin.sortGC
  return(bin.samp)
  
}



# motif analysis

LocalMotifLevel <- function(repBins,genome,motif){
  
  
  library(BSgenome)
  bsGenome <- available.genomes()[grep(genome,available.genomes())]
  bsGenome <- bsGenome[-(grep("masked", bsGenome))]
  library(bsGenome, character.only=TRUE)
  bsSpec <- get(strsplit(bsGenome,split="\\.")[[1]][2])
  
  Region = repBins[,c("chr", "start", "end")]
  chromos <- as.character(unique(Region$chr))
  Region$motif <- rep(NA, nrow(Region))
  motif <- DNAString(motif)
  

  
  for(i in 1:length(chromos)){
    Seq.set=DNAStringSet(bsSpec[[chromos[i]]], 
                         start=Region$start[Region$chr == chromos[i]], 
                         end=Region$end[Region$chr == chromos[i]])
    patternContent <- unlist(lapply(vmatchPattern(pattern = motif,subject = Seq.set), length))
    patternContentRC <- unlist(lapply(vmatchPattern(pattern = reverseComplement(motif),subject = Seq.set), length))
    
    
    Region$motif[Region$chr == chromos[i]] <- patternContent + patternContentRC
    print(chromos[i])
  }
  
  
  Region$width <- Region$end - Region$start + 1
  #### here we get the final piece of the puzzle 
  return(Region)
}



rateRoller <- function(primeGC3, primeGC5, primeBF3, primeBF5, k){
  
  if(k%%2 != 1){
    stop("k must be odd")
  }
  
  GCsum3 <- rollsum(primeGC3,k = k, align = "center")
  baseSum3 <- rollsum(primeBF3,k = k, align = "center")
  
  GCsum5 <- rollsum(primeGC5,k = k, align = "center")
  baseSum5 <- rollsum(primeBF5,k = k, align = "center")
  
  roll3start <- roll3end <- roll5start <- roll5end <- roll3startGC <- roll3endGC <- roll5startGC <- roll5startGC <- roll5endGC <- rep(NA,k/2)
  
  for(i in 1:(k/2)){
    roll3start[i] <- sum(primeBF3[1:(i+(i-1))])
    roll3end[length(roll3end) - i + 1] <- sum(primeBF3[length(primeBF3):1][1:(i+(i-1))])
    
    roll5start[i] <- sum(primeBF5[1:(i+(i-1))])
    roll5end[length(roll5end) - i +1 ] <- sum(primeBF5[length(primeBF5):1][1:(i+(i-1))])
    
    roll3startGC[i] <- sum(primeGC3[1:(i+(i-1))])
    roll3endGC[length(roll3endGC) - i + 1] <- sum(primeGC3[length(primeGC3):1][1:(i+(i-1))])
    
    roll5startGC[i] <- sum(primeGC5[1:(i+(i-1))])
    roll5endGC[length(roll5endGC) - i + 1] <- sum(primeGC5[length(primeGC5):1][1:(i+(i-1))])
  }
  
  return(list(prime3gc = c(roll3startGC, GCsum3, roll3endGC)/c(roll3start, baseSum3, roll3end), 
              prime5gc = c(roll5startGC, GCsum5, roll5endGC)/c(roll5start, baseSum5, roll5end)
  ) 
  )
  
}



meanRoller35 <- function(primeGC3, primeGC5, primeBF3, primeBF5, k){
  
  if(k%%2 != 1){
    stop("k must be odd")
  }
  
  rate3 <- primeGC3/primeBF3
  rate5 <- primeGC5/primeBF5
  
  Mrate3 <- rollmean(rate3,k = k, align = "center")
  Mrate5 <- rollmean(rate5,k = k, align = "center")
  rate3start <- rate3end <- rate5start <- rate5end <- rep(NA, k/2)
  for(i in 1:(k/2)){
    rate3start[i] <- mean(rate3[1:(i+(i-1))])
    rate3end[length(rate3end) - i + 1] <- mean(rate3[length(rate3):1][1:(i+(i-1))])
    rate5start[i] <- mean(rate5[1:(i+(i-1))])
    rate5end[length(rate5end) - i + 1] <- mean(rate5[length(rate5):1][1:(i+(i-1))])
  }
  return(list(prime3gc = c(rate3start, Mrate3, rate3end), prime5gc=c(rate5start, Mrate5, rate5end)))
}

meanSdRoller <- function(repCov, bpFreq,k){
  
  if(k%%2 != 1){
    stop("k must be odd")
  }
  
  rate <- repCov/bpFreq
  
  Mrate <- rollmean(rate,k = k, align = "center")
  SDrate <- rollapply(data = zoo(rate), width = k, FUN = sd)
  
  rateStartM <- rateEndM <- rateEndSD <- rateStartSD <- rep(NA, k/2)
  for(i in 1:(k/2)){
    rateStartM[i] <- mean(rate[1:(i+(i-1))])
    rateEndM[length(rateEndM) - i + 1] <- mean(rate[length(rate):1][1:(i+(i-1))])
    rateStartSD[i] <- sd(rate[1:(i+(i-1))])
    rateEndSD[length(rateEndSD) - i + 1] <- sd(rate[length(rate):1][1:(i+(i-1))])
  }
  return(list(mean = c(rateStartM, Mrate, rateEndM), SD = c(rateStartSD, Mrate, rateEndSD)))
}


covImage <- function(lenChoice, repBins , repList ,chromList, 
                     refgene , type ,polyL1 = NULL , minRepCov = NULL, maxRepCov = NULL, 
                     minRepSize = NULL, maxRepSize = NULL, minBinSize = NULL, 
                     maxBinSize = NULL){
  
  chromR <- chromList$R
  colnames(chromR)[2:4] <- c("genoName", "genoStart", "genoEnd")
  repList = c(repList, list(represedChromatin = chromR))
  if(length(polyL1) != 0){
    repList = c(repList, list(polyL1 = polyL1))
  }
  
  refgene.gr <- GRanges(seqnames=Rle(refgene[,3]), ranges = IRanges(start = refgene[,5], end = refgene[,6]))
  
  
  repBins <- repBins[, c("chr", "start", "end", "Known", repChoice)]
  
  if(!is.null(maxBinSize)){
    repBins <- repBins[repBins$end - repBins$start + 1 < maxBinSize,]
  }
  if(!is.null(minBinSize)){
    repBins <- repBins[repBins$end - repBins$start + 1 > minBinSize,]
  }
  colnames(repBins)[5] = "repCov"
  if(!is.null(maxRepCov)){
    repBins <- repBins[repBins$repCov < maxRepCov,]
  }
  if(!is.null(minRepCov)){
    repBins <- repBins[repBins$repCov > minRepCov,]
  }
  
  if(!is.null(maxRepSize)){
    repGR <- repGR[width(repGR) < maxRepSize,]
  }
  if(!is.null(minRepSize)){
    repGR <- repGR[width(repGR) > minRepSize,]
  }
  seqLen <- lenChoice + 1
  names(seqLen) <- "seq"
  
  # it might be easier if we break down the four groups at the start
  
  us.ends <- repBins$start + ((repBins$end - repBins$start)/2)
  us.ends[us.ends - repBins$start+1 > lenChoice] <- repBins$start[us.ends - repBins$start+1 > lenChoice] + lenChoice
  us.bins <- data.frame(chr = repBins$chr, start = repBins$start, end = us.ends)
  us.bins.gr <- GRanges(seqnames=Rle(us.bins$chr), 
                        ranges = IRanges(start=us.bins$start, end = us.bins$end))
  
  ds.starts <- as.integer(repBins$end - ((repBins$end - repBins$start)/2))
  ds.starts[repBins$end - ds.starts + 1> lenChoice] <- repBins$end[repBins$end - ds.starts + 1> lenChoice] - lenChoice
  ds.bins <- data.frame(chr = repBins$chr, start = ds.starts, end = repBins$end)
  ds.bins.gr <- GRanges(seqnames=Rle(ds.bins$chr), 
                        ranges = IRanges(start=ds.bins$start, end = ds.bins$end))
  
  if(type == "intergenic"){
    us.s.ol <- as.matrix(findOverlaps(us.bins.gr, refgene.gr, maxgap=1))
    if(!(length(unique(us.s.ol[,1])) == length(us.bins.gr))){
      stop("some intergenic regions have no upstream gene")
    }
    us.strand <- data.frame(intergenicID = us.s.ol[,1],strand = refgene[us.s.ol[,2], 4])
    us.plus <- unique(us.strand$intergenicID[us.strand$strand == "+"])
    us.minus <- unique(us.strand$intergenicID[us.strand$strand == "-"])
    
    if(length(c(us.plus[us.plus %in% us.minus],us.minus[us.minus %in% us.plus])) != 0){
      stop("strand confusion in upstream intergenic regions")
    }
    
    ds.s.ol <- as.matrix(findOverlaps(ds.bins.gr, refgene.gr, maxgap=1))
    if(!(length(unique(ds.s.ol[,1])) == length(ds.bins.gr))){
      stop("some intergenic regions have no downstream gene")
    }
    ds.strand <- data.frame(intergenicID = ds.s.ol[,1],strand = refgene[ds.s.ol[,2], 4])
    ds.plus <- unique(ds.strand$intergenicID[ds.strand$strand == "+"])
    ds.minus <- unique(ds.strand$intergenicID[ds.strand$strand == "-"])
    
    if(length(c(ds.plus[ds.plus %in% ds.minus],ds.minus[ds.minus %in% ds.plus])) != 0){
      stop("strand confusion in downstream intergenic regions")
    }
    
  }else if(type == "intron"){
    
    us.s.ol <- as.matrix(findOverlaps(us.bins.gr, refgene.gr, minoverlap=2))
    if(!(length(unique(us.s.ol[,1])) == length(us.bins.gr))){
      stop("some upstream intron regions don't belong to a gene")
    }
    us.strand <- data.frame(intronID = us.s.ol[,1],strand = refgene[us.s.ol[,2], 4])
    us.plus <- unique(us.strand$intronID[us.strand$strand == "+"])
    us.minus <- unique(us.strand$intronID[us.strand$strand == "-"])
    
    if(length(c(us.plus[us.plus %in% us.minus],us.minus[us.minus %in% us.plus])) != 0){
      stop("strand confusion in upstream intron region")
    }
    
    ds.s.ol <- as.matrix(findOverlaps(ds.bins.gr, refgene.gr, minoverlap=2))
    if(!(length(unique(ds.s.ol[,1])) == length(ds.bins.gr))){
      stop("some downstream intron regions don't belong to a gene")
    }
    ds.strand <- data.frame(intergenicID = ds.s.ol[,1],strand = refgene[ds.s.ol[,2], 4])
    ds.plus <- unique(ds.strand$intergenicID[ds.strand$strand == "+"])
    ds.minus <- unique(ds.strand$intergenicID[ds.strand$strand == "-"])
    
    if(length(c(ds.plus[ds.plus %in% ds.minus],ds.minus[ds.minus %in% ds.plus])) != 0){
      stop("strand confusion in downstream intron region")
    }
  }else{
    stop("type needs to be intron or intergenic")
  }
  
  
  us.bins.gr_5 <- us.bins.gr[us.minus]
  binStartsUs5 <- as.data.frame(us.bins.gr_5)
  binStartsUs5 <- data.frame( binStartsUs5, matrix(NA, nrow = length(us.bins.gr_5), ncol = length(repList),dimnames = list(1:length(us.bins.gr_5), names(repList))))
  covListUs5 <- NULL
  
  
  us.bins.gr_3 <- us.bins.gr[us.plus]
  binStartsUs3 <- as.data.frame(us.bins.gr_3)
  binStartsUs3 <- data.frame( binStartsUs3, matrix(NA, nrow = length(us.bins.gr_3), ncol = length(repList),dimnames = list(1:length(us.bins.gr_3), names(repList))))
  covListUs3 <- NULL
  
  
  ds.bins.gr_5 <- ds.bins.gr[ds.plus]
  binStartsDs5 <- as.data.frame(ds.bins.gr_5)
  binStartsDs5 <- data.frame( binStartsDs5, matrix(NA, nrow = length(ds.bins.gr_5), ncol = length(repList),dimnames = list(1:length( ds.bins.gr_5), names(repList))))
  covListDs5 <- NULL
  
  ds.bins.gr_3 <- ds.bins.gr[ds.minus]
  binStartsDs3 <- as.data.frame(ds.bins.gr_3)
  binStartsDs3 <- data.frame( binStartsDs3, matrix(NA, nrow = length(ds.bins.gr_3), ncol = length(repList),dimnames = list(1:length( ds.bins.gr_3), names(repList))))
  covListDs3 <- NULL
  
  
  # set the loop up here
  for(r in 1:length(repList)){
    repGR <- GRanges(seqnames=Rle(repList[[r]]$genoName),
                     ranges = IRanges(start = repList[[r]]$genoStart, end = repList[[r]]$genoEnd -1))
    us.rep.int_5 <- intersect(repGR, us.bins.gr_5)
    us.Ol_5 <- as.matrix(GenomicRanges::findOverlaps(us.bins.gr_5, us.rep.int_5, select = "first"))
    us.Ol_5 <- data.frame(1:nrow(us.Ol_5), us.Ol_5)
    us.Ol_5 <- us.Ol_5[complete.cases(us.Ol_5),]
    binStartsUs5[us.Ol_5[,1],names(repList)[r]] <- start(us.rep.int_5[us.Ol_5[,2]]) - start(us.bins.gr_5[us.Ol_5[,1]]) + 1
    
    us.Ol_5 <- as.matrix(GenomicRanges::findOverlaps(us.bins.gr_5, us.rep.int_5, select = "all"))
    us.cov_5 <- data.frame(start = start(us.rep.int_5[us.Ol_5[,2]]) - start(us.bins.gr_5[us.Ol_5[,1]]) + 1, 
                           end =end(us.rep.int_5[us.Ol_5[,2]]) - start(us.bins.gr_5[us.Ol_5[,1]]) + 1,
                           binID = us.Ol_5[,1], width.rank = (rank(width(us.bins.gr_5),ties.method = "random"))[us.Ol_5[,1]], width= width(us.bins.gr_5)[us.Ol_5[,1]])
    
    covListUs5 <- c(covListUs5, list(us.cov_5))
    
    
    us.rep.int_3 <- intersect(repGR, us.bins.gr_3)
    us.Ol_3 <- as.matrix(findOverlaps(us.bins.gr_3, us.rep.int_3, select = "first"))
    us.Ol_3 <- data.frame(1:nrow(us.Ol_3), us.Ol_3)
    us.Ol_3 <- us.Ol_3[complete.cases(us.Ol_3),]
    binStartsUs3[us.Ol_3[,1],names(repList)[r]] <- start(us.rep.int_3[us.Ol_3[,2]]) - start(us.bins.gr_3[us.Ol_3[,1]]) + 1
    
    
    us.Ol_3 <- as.matrix(findOverlaps(us.bins.gr_3, us.rep.int_3))
    us.cov_3 <- data.frame(start = start(us.rep.int_3[us.Ol_3[,2]]) - start(us.bins.gr_3[us.Ol_3[,1]]) + 1, 
                           end =end(us.rep.int_3[us.Ol_3[,2]]) - start(us.bins.gr_3[us.Ol_3[,1]]) + 1,
                           binID = us.Ol_3[,1], width.rank = (rank(width(us.bins.gr_3),ties.method = "random"))[us.Ol_3[,1]], width= width(us.bins.gr_3)[us.Ol_3[,1]])
    
    covListUs3 <- c(covListUs3, list(us.cov_3))
    
    
    
    ds.rep.int_5 <- intersect(repGR, ds.bins.gr_5)
    ds.Ol_5 <- as.matrix(findOverlaps(ds.bins.gr_5, ds.rep.int_5, select = "first"))
    ds.Ol_5 <- data.frame(1:nrow(ds.Ol_5), ds.Ol_5)
    ds.Ol_5 <- ds.Ol_5[complete.cases(ds.Ol_5),]
    binStartsUs3[ds.Ol_5[,1],names(repList)[r]] <- (lenChoice - (end(ds.bins.gr_5[ds.Ol_5[,1]]) - start(ds.bins.gr_5[ds.Ol_5[,1]]))) + (start(ds.rep.int_5[ds.Ol_5[,2]]) - start(ds.bins.gr_5[ds.Ol_5[,1]])) + 1
    
    ds.Ol_5 <- as.matrix(findOverlaps(ds.bins.gr_5, ds.rep.int_5))
    ds.cov_5 <- data.frame(start =   (lenChoice - (end(ds.bins.gr_5[ds.Ol_5[,1]]) - start(ds.bins.gr_5[ds.Ol_5[,1]]))) + (start(ds.rep.int_5[ds.Ol_5[,2]]) - start(ds.bins.gr_5[ds.Ol_5[,1]])) + 1 ,
                           end =  (lenChoice - (end(ds.bins.gr_5[ds.Ol_5[,1]]) - start(ds.bins.gr_5[ds.Ol_5[,1]]))) + (end(ds.rep.int_5[ds.Ol_5[,2]]) - start(ds.bins.gr_5[ds.Ol_5[,1]])) + 1,
                           binID = ds.Ol_5[,1], width.rank = (rank(width(ds.bins.gr_5),ties.method = "random"))[ds.Ol_5[,1]], width= width(ds.bins.gr_5)[ds.Ol_5[,1]])
    
    covListDs5 <- c(covListDs5, list(ds.cov_5))
    
    ds.rep.int_3 <- intersect(repGR, ds.bins.gr_3)
    ds.Ol_3 <- as.matrix(findOverlaps(ds.bins.gr_3, ds.rep.int_3, select = "first"))
    ds.Ol_3 <- data.frame(1:nrow(ds.Ol_3), ds.Ol_3)
    ds.Ol_3 <- ds.Ol_3[complete.cases(ds.Ol_3),]
    binStartsUs3[ds.Ol_3[,1],names(repList)[r]] <- (lenChoice - (end(ds.bins.gr_3[ds.Ol_3[,1]]) - start(ds.bins.gr_3[ds.Ol_3[,1]]))) + (start(ds.rep.int_3[ds.Ol_3[,2]]) - start(ds.bins.gr_3[ds.Ol_3[,1]])) + 1
    
    ds.Ol_3 <- as.matrix(findOverlaps(ds.bins.gr_3, ds.rep.int_3))
    ds.cov_3 <- data.frame(start =   (lenChoice - (end(ds.bins.gr_3[ds.Ol_3[,1]]) - start(ds.bins.gr_3[ds.Ol_3[,1]]))) + (start(ds.rep.int_3[ds.Ol_3[,2]]) - start(ds.bins.gr_3[ds.Ol_3[,1]])) + 1 ,
                           end =  (lenChoice - (end(ds.bins.gr_3[ds.Ol_3[,1]]) - start(ds.bins.gr_3[ds.Ol_3[,1]]))) + (end(ds.rep.int_3[ds.Ol_3[,2]]) - start(ds.bins.gr_3[ds.Ol_3[,1]])) + 1,
                           binID = ds.Ol_3[,1], width.rank = (rank(width(ds.bins.gr_3),ties.method = "random"))[ds.Ol_3[,1]], width= width(ds.bins.gr_3)[ds.Ol_3[,1]])
    
    covListDs3 <- c(covListDs3, list(ds.cov_3))
    
    
  }
  names(covListUs5) <- names(repList)
  names(covListUs3) <- names(repList)
  names(covListDs5) <- names(repList)
  names(covListDs3) <- names(repList)
  
  
  # we need to trun some of these around to match their relative gene position
  
  
  output = list(covListUs5 = covListUs5, 
                covListDs5 = covListDs5,
                covListUs3 = covListUs3,
                covListDs3 = covListDs3,
                startsUs5 = binStartsUs5,
                startsDs5 = binStartsDs5,
                startsUs3 = binStartsUs3,
                startsDs3 = binStartsDs3
  )
  return(output)
  
}  


### this will get us a matrix with repeat info for equally sized regions around a boundary
### there will be a matrix for each repeat type looked at, and all matricies will be stored in a list

# something going wrong with consesnsus set


boundaryAnotation <- function(repList, binSize, regionSize, repTypes, boundaryLine, boundaryChr){
  TErateMatrixList <- NULL
  for(i in 1:length(repList)){
    TErateMatrixList <- c(TErateMatrixList, list(matrix(nrow = length(boundaryLine), ncol = regionSize/binSize)))
  }
  names(TErateMatrixList) <- names(repList)
  
  for(i in 1:(regionSize/binSize)){
    print(i)
    bin <- data.frame(chr = boundaryChr, 
                      start = (boundaryLine - (regionSize/2)) + ((i - 1) * binSize), 
                      end = (boundaryLine - (regionSize/2)) + ((i) * binSize))
    bin$Known <- binSize
    sorted <- binSort(rep = repList, bins = bin,TE.names = names(repList), repType = repTypes )
    for(te in 1:length(repList)){
      TErateMatrixList[[(names(repList)[te])]][,i] <- sorted$rates[,names(repList)[te]]
    }
    
  }
  return(TErateMatrixList)
}


binScoreSort <- function(repList, bins, TE.names, repType, metadata){
  bin.gr <- GRanges(seqnames=Rle(bins$chr),
                    ranges = IRanges(start=bins$start, end = bins$end - 1))
  
  for(te in 1:length(TE.names)){
    print(TE.names[te])
    if(repType[te] == "repeats"){
      te.gr <- repList[[TE.names[te]]]
      GR <- GRanges(seqnames = Rle(te.gr$genoName),
                    ranges = IRanges(start = te.gr$genoStart, end = te.gr$genoEnd))
    }else if(repType[te] == "chromatin"){
      te.gr <- repList[[TE.names[te]]]
      GR <- GRanges(seqnames = Rle(te.gr$chrom),
                    ranges = IRanges(start = te.gr$chromStart, end = te.gr$chromEnd))
    }else if(repType[te] == "GRange"){
      GR <- repList[[te]]
    }else{
      stop("no repType")
    }
    
    # here we can do repeat coverage 
    
    OL.f <- as.matrix(findOverlaps(bin.gr,GR))
    if(length(OL.f) > 0){
      OL.agg <- aggregate( x= elementMetadata(GR[OL.f[,2]])[[metadata]] , by= list(OL.f[,1]), FUN=mean)
      bins[OL.agg[,1],TE.names[te]] <- OL.agg[,2]
      bins[is.na(bins[,TE.names[te]]),TE.names[te]] <- 0
    }else{
      bins[,TE.names[te]] <- 0
    }
    
  }
  return(bins)
}




boundaryScoreAnotation <- function(repList, binSize, regionSize, repTypes, boundaryLine, boundaryChr, metadata, TE.names){
  TErateMatrixList <- NULL
  for(i in 1:length(repList)){
    TErateMatrixList <- c(TErateMatrixList, list(matrix(nrow = length(boundaryLine), ncol = regionSize/binSize)))
  }
  names(TErateMatrixList) <- TE.names
  
  for(i in 1:(regionSize/binSize)){
    print(i)
    bin <- data.frame(chr = boundaryChr, 
                      start = (boundaryLine - (regionSize/2)) + ((i - 1) * binSize), 
                      end = (boundaryLine - (regionSize/2)) + ((i) * binSize))
    bin$Known <- binSize
    sorted <- binScoreSort(rep = repList, bins = bin,TE.names = TE.names, repType = repTypes , metadata = metadata)
    for(te in 1:length(repList)){
      TErateMatrixList[[(TE.names)[te]]][,i] <- sorted[,TE.names[te]]
    }
    
  }
  return(TErateMatrixList)
}


reuben.heatmap <- function (x, Rowv = NULL, Colv = if (symm) "Rowv" else NULL, col = heat.colors(20),
                            distfun = dist, hclustfun = hclust, 
                            reorderfun = function(d,w) reorder(d, w), add.expr, symm = FALSE, revC = identical(Colv, "Rowv"), 
                            scale = c("row", "column", "none"), na.rm = TRUE, 
                            margins = c(5, 5), ColSideColors, RowSideColors, cexRow = 0.2 + 1/log10(nr), 
                            cexCol = 0.2 + 1/log10(nc), labRow = NULL, 
                            labCol = NULL, main = NULL, xlab = NULL, ylab = NULL, keep.dendro = FALSE, 
                            verbose = getOption("verbose"), colsLabs = NULL, ...) 
{
  scale <- if (symm && missing(scale)) 
    "none"
  else match.arg(scale)
  if (length(di <- dim(x)) != 2 || !is.numeric(x)) 
    stop("'x' must be a numeric matrix")
  nr <- di[1L]
  nc <- di[2L]
  if (nr <= 1 || nc <= 1) 
    stop("'x' must have at least 2 rows and 2 columns")
  if (!is.numeric(margins) || length(margins) != 2L) 
    stop("'margins' must be a numeric vector of length 2")
  doRdend <- !identical(Rowv, NA)
  doCdend <- !identical(Colv, NA)
  if (!doRdend && identical(Colv, "Rowv")) 
    doCdend <- FALSE
  if (is.null(Rowv)) 
    Rowv <- rowMeans(x, na.rm = na.rm)
  if (is.null(Colv)) 
    Colv <- colMeans(x, na.rm = na.rm)
  if (doRdend) {
    if (inherits(Rowv, "dendrogram")) 
      ddr <- Rowv
    else {
      hcr <- hclustfun(distfun(x))
      ddr <- as.dendrogram(hcr)
      if (!is.logical(Rowv) || Rowv) 
        ddr <- reorderfun(ddr, Rowv)
    }
    if (nr != length(rowInd <- order.dendrogram(ddr))) 
      stop("row dendrogram ordering gave index of wrong length")
  }
  else rowInd <- 1L:nr
  if (doCdend) {
    if (inherits(Colv, "dendrogram")) 
      ddc <- Colv
    else if (identical(Colv, "Rowv")) {
      if (nr != nc) 
        stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
      ddc <- ddr
    }
    else {
      hcc <- hclustfun(distfun(if (symm) 
        x
        else t(x)))
      ddc <- as.dendrogram(hcc)
      if (!is.logical(Colv) || Colv) 
        ddc <- reorderfun(ddc, Colv)
    }
    if (nc != length(colInd <- order.dendrogram(ddc))) 
      stop("column dendrogram ordering gave index of wrong length")
  }
  else colInd <- 1L:nc
  x <- x[rowInd, colInd]
  
  
  
  labRow <- if (is.null(labRow)) 
    if (is.null(rownames(x))) 
      (1L:nr)[rowInd]
  else rownames(x)
  else labRow[rowInd]
  labCol <- if (is.null(labCol)) 
    if (is.null(colnames(x))) 
      (1L:nc)[colInd]
  else colnames(x)
  else labCol[colInd]
  if (scale == "row") {
    x <- sweep(x, 1L, rowMeans(x, na.rm = na.rm), check.margin = FALSE)
    sx <- apply(x, 1L, sd, na.rm = na.rm)
    x <- sweep(x, 1L, sx, "/", check.margin = FALSE)
  }
  else if (scale == "column") {
    x <- sweep(x, 2L, colMeans(x, na.rm = na.rm), check.margin = FALSE)
    sx <- apply(x, 2L, sd, na.rm = na.rm)
    x <- sweep(x, 2L, sx, "/", check.margin = FALSE)
  }
  lmat <- rbind(c(NA, 3), 2:1)
  lwid <- c(if (doRdend) 1 else 0.05, 4)
  lhei <- c((if (doCdend) 1 else 0.05) + if (!is.null(main)) 0.2 else 0, 
            4)
  if (!missing(ColSideColors)) {
    if (!is.character(ColSideColors) || length(ColSideColors) != 
        nc) 
      stop("'ColSideColors' must be a character vector of length ncol(x)")
    lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
    lhei <- c(lhei[1L], 0.2, lhei[2L])
  }
  if (!missing(RowSideColors)) {
    if (!is.character(RowSideColors) || length(RowSideColors) != 
        nr) 
      stop("'RowSideColors' must be a character vector of length nrow(x)")
    lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 
                                   1), lmat[, 2] + 1)
    lwid <- c(lwid[1L], 0.2, lwid[2L])
  }
  lmat[is.na(lmat)] <- 0
  if (verbose) {
    cat("layout: widths = ", lwid, ", heights = ", lhei, 
        "; lmat=\n")
    print(lmat)
  }
  dev.hold()
  on.exit(dev.flush())
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  layout(lmat, widths = lwid, heights = lhei, respect = TRUE)
  if (!missing(RowSideColors)) {
    par(mar = c(margins[1L], 0, 0, 0.5))
    image(rbind(if (revC) 
      nr:1L
      else 1L:nr), col = RowSideColors[rowInd], axes = FALSE)
  }
  if (!missing(ColSideColors)) {
    par(mar = c(0.5, 0, 0, margins[2L]))
    image(cbind(1L:nc), col = ColSideColors[colInd], axes = FALSE)
  }
  par(mar = c(margins[1L], 5, 0, 1))
  if (!symm || scale != "none") 
    x <- t(x)
  if (revC) {
    iy <- nr:1
    if (doRdend) 
      ddr <- rev(ddr)
    x <- x[, iy]
  }
  else iy <- 1L:nr
  
  
  x.2 <- x
  x.2[lower.tri(x.2)] <- NA
  
  if(!is.null(colsLabs)){
    colLabCol <- colsLabs[colInd]
    colLabRow <- colsLabs[rowInd]
  }else{
    colLabCol <- colLabRow <- rep("black", length(labCol))
  }
  
  
  image(1L:nc, 1L:nr, x.2, xlim = 0.5 + c(0, nc), ylim = 0.5 + 
          c(0, nr), axes = FALSE, xlab = "", ylab = "", ...,col = col)
  axis(2, 1L:nc, labels = FALSE, las = 2, line = -0.5, tick = 0, cex.axis = cexCol, col = c(1,2))
  mtext(text = labCol, side = 2,line = .5,at = 1L:nc,col = colLabCol,las = 2,adj = 1)
  if (!is.null(xlab)) 
    mtext(xlab, side = 2, line = margins[1L] - 1.25)
  axis(3, iy, labels = FALSE, las = 2, line = -0.5, tick = 0, cex.axis = cexRow)
  mtext(text = labCol, side = 3,line = .5,at = 1L:nc,col = colLabCol,las = 2,adj = 0)
  if (!is.null(ylab)) 
    mtext(ylab, side = 3, line = margins[2L] - 1.25)
  if (!missing(add.expr)) 
    eval(substitute(add.expr))
  par(mar = c(margins[1L], 0, 0, 0))
  if (doRdend) 
    plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
  else frame()
  par(mar = c(0, 0, if (!is.null(main)) 1 else 0, margins[2L]))
  if (doCdend) 
    #plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
    plot(1, type = "n", axes = F)
  else if (!is.null(main)) 
    frame()
  if (!is.null(main)) {
    par(xpd = NA)
    title(main, cex.main = 1.5 * op[["cex.main"]])
  }
  invisible(list(rowInd = rowInd, colInd = colInd, Rowv = if (keep.dendro && doRdend) ddr, Colv = if (keep.dendro && doCdend) ddc))
}


heatmap.colours <- function (x, Rowv = NULL, Colv = if (symm) "Rowv" else NULL, 
                             distfun = dist, hclustfun = hclust, reorderfun = function(d, w) reorder(d, w), 
                             add.expr, symm = FALSE, revC = identical(Colv,"Rowv"), scale = c("row", "column", "none"), na.rm = TRUE, 
                             margins = c(5, 5), ColSideColors, RowSideColors, cexRow = 0.2 + 
                               1/log10(nr), cexCol = 0.2 + 1/log10(nc), labRow = NULL, 
                             labCol = NULL, main = NULL, xlab = NULL, ylab = NULL, keep.dendro = FALSE, 
                             verbose = getOption("verbose"), ColcolsLabs = NULL, RowcolsLabs = NULL, col= col, zlim = zlim) 
{
  scale <- if (symm && missing(scale)) 
    "none"
  else match.arg(scale)
  if (length(di <- dim(x)) != 2 || !is.numeric(x)) 
    stop("'x' must be a numeric matrix")
  nr <- di[1L]
  nc <- di[2L]
  if (nr <= 1 || nc <= 1) 
    stop("'x' must have at least 2 rows and 2 columns")
  if (!is.numeric(margins) || length(margins) != 2L) 
    stop("'margins' must be a numeric vector of length 2")
  doRdend <- !identical(Rowv, NA)
  doCdend <- !identical(Colv, NA)
  if (!doRdend && identical(Colv, "Rowv")) 
    doCdend <- FALSE
  if (is.null(Rowv)) 
    Rowv <- rowMeans(x, na.rm = na.rm)
  if (is.null(Colv)) 
    Colv <- colMeans(x, na.rm = na.rm)
  if (doRdend) {
    if (inherits(Rowv, "dendrogram")) 
      ddr <- Rowv
    else {
      hcr <- hclustfun(distfun(x))
      ddr <- as.dendrogram(hcr)
      if (!is.logical(Rowv) || Rowv) 
        ddr <- reorderfun(ddr, Rowv)
    }
    if (nr != length(rowInd <- order.dendrogram(ddr))) 
      stop("row dendrogram ordering gave index of wrong length")
  }
  else rowInd <- 1L:nr
  if (doCdend) {
    if (inherits(Colv, "dendrogram")) 
      ddc <- Colv
    else if (identical(Colv, "Rowv")) {
      if (nr != nc) 
        stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
      ddc <- ddr
    }
    else {
      hcc <- hclustfun(distfun(if (symm) 
        x
        else t(x)))
      ddc <- as.dendrogram(hcc)
      if (!is.logical(Colv) || Colv) 
        ddc <- reorderfun(ddc, Colv)
    }
    if (nc != length(colInd <- order.dendrogram(ddc))) 
      stop("column dendrogram ordering gave index of wrong length")
  }
  else colInd <- 1L:nc
  x <- x[rowInd, colInd]
  labRow <- if (is.null(labRow)) 
    if (is.null(rownames(x))) 
      (1L:nr)[rowInd]
  else rownames(x)
  else labRow[rowInd]
  labCol <- if (is.null(labCol)) 
    if (is.null(colnames(x))) 
      (1L:nc)[colInd]
  else colnames(x)
  else labCol[colInd]
  if (scale == "row") {
    x <- sweep(x, 1L, rowMeans(x, na.rm = na.rm), check.margin = FALSE)
    sx <- apply(x, 1L, sd, na.rm = na.rm)
    x <- sweep(x, 1L, sx, "/", check.margin = FALSE)
  }
  else if (scale == "column") {
    x <- sweep(x, 2L, colMeans(x, na.rm = na.rm), check.margin = FALSE)
    sx <- apply(x, 2L, sd, na.rm = na.rm)
    x <- sweep(x, 2L, sx, "/", check.margin = FALSE)
  }
  lmat <- rbind(c(NA, 3), 2:1)
  lwid <- c(if (doRdend) 1 else 0.05, 4)
  lhei <- c((if (doCdend) 1 else 0.05) + if (!is.null(main)) 0.2 else 0, 
            4)
  if (!missing(ColSideColors)) {
    if (!is.character(ColSideColors) || length(ColSideColors) != 
        nc) 
      stop("'ColSideColors' must be a character vector of length ncol(x)")
    lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
    lhei <- c(lhei[1L], 0.2, lhei[2L])
  }
  if (!missing(RowSideColors)) {
    if (!is.character(RowSideColors) || length(RowSideColors) != 
        nr) 
      stop("'RowSideColors' must be a character vector of length nrow(x)")
    lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 
                                   1), lmat[, 2] + 1)
    lwid <- c(lwid[1L], 0.2, lwid[2L])
  }
  lmat[is.na(lmat)] <- 0
  if (verbose) {
    cat("layout: widths = ", lwid, ", heights = ", lhei, 
        "; lmat=\n")
    print(lmat)
  }
  dev.hold()
  on.exit(dev.flush())
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  layout(lmat, widths = lwid, heights = lhei, respect = TRUE)
  if (!missing(RowSideColors)) {
    par(mar = c(margins[1L], 0, 0, 0.5))
    image(rbind(if (revC) 
      nr:1L
      else 1L:nr), col = RowSideColors[rowInd], axes = FALSE)
  }
  if (!missing(ColSideColors)) {
    par(mar = c(0.5, 0, 0, margins[2L]))
    image(cbind(1L:nc), col = ColSideColors[colInd], axes = FALSE)
  }
  par(mar = c(margins[1L], 0, 0, margins[2L]))
  if (!symm || scale != "none") 
    x <- t(x)
  if (revC) {
    iy <- nr:1
    if (doRdend) 
      ddr <- rev(ddr)
    x <- x[, iy]
  }
  else iy <- 1L:nr
  image(1L:nc, 1L:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + 
          c(0, nr), axes = FALSE, xlab = "", ylab = "",col = col, zlim = zlim)
  axis(1, 1L:nc, labels = FALSE, las = 2, line = -0.5, tick = 0, 
       cex.axis = cexCol)
  mtext(text = labCol, side = 1,line = .5,at = 1L:nc,col = ColcolsLabs[colInd],las = 2,adj = 1)
  if (!is.null(xlab)) 
    
  if (!is.null(xlab)) 
    mtext(xlab, side = 1, line = margins[1L] - 1.25)
  axis(4, iy, labels = FALSE, las = 2, line = -0.5, tick = 0, 
       cex.axis = cexRow)
  mtext(text = labRow, side = 4,line = .5,at = 1L:nc,col = RowcolsLabs[rowInd],las = 2,adj = 0)
  if (!is.null(xlab)) 
    
  if (!is.null(ylab)) 
    mtext(ylab, side = 4, line = margins[2L] - 1.25)
  if (!missing(add.expr)) 
    eval(substitute(add.expr))
  par(mar = c(margins[1L], 0, 0, 0))
  if (doRdend) 
    plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
  else frame()
  par(mar = c(0, 0, if (!is.null(main)) 1 else 0, margins[2L]))
  if (doCdend) 
    plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
  else if (!is.null(main)) 
    frame()
  if (!is.null(main)) {
    par(xpd = NA)
    title(main, cex.main = 1.5 * op[["cex.main"]])
  }
  invisible(list(rowInd = rowInd, colInd = colInd, Rowv = if (keep.dendro &&  doRdend) ddr, 
                 Colv = if (keep.dendro && doCdend) ddc))
}





#######
# #######  Functions 
#
# binned.genome.reader
# General.G
# General.G.WG
# hostspot.G
# hotspot.G.WG
# M.j
# region.Finder
# Z.moran
# Z.moran.WG












# 
# 
# 
# A <- binned.genome.reader(bin.size=c(50000))
# 
# B <- A[[1]][10:19,]
# 
# X <- 1:nrow(B)
# 
# N <- neighbor.matrix(B,X)
# neighbors = N
# 
# 
# head(Neighbors - X.bar)*head(w)
# 
# 
# 
# Z.moran.WG(X, neighbors)
# 
# 
# General.G.WG(X,neighbors)
# 
# 
# n = length(X)
# w = matrix(data = 0, ncol = n, nrow = n)
# for(W in 1:(n-1)){
#   w[W+ 1,W] = 1
#   w[W, W +1] = 1
#   # w[W,W] = 1
# }
# 
# General.G(X,w)
# 
# 
# 
# X = 1:20
