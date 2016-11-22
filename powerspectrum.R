# This code was written by Vishwesha Guttal and then modeified by Sabiha Majumder 

powerspectrum <- function(mat) {
  
  nr <- nrow(mat)
  nc <- ncol(mat)
  
  n0x <- floor(nc/2) + 1
  n0y <- floor(nr/2) + 1
  
  # Create distance and angle matrices
  f1 <- t(replicate(nr,seq(1,nc))) - n0x
  f2 <- replicate(nc,seq(1,nr)) - n0y
  DIST <- sqrt(f1^2 + f2^2)
  ANGLE <- atan2(-f2,f1)*180/pi
  
  # Calculate DFT
  mi <- 1
  ma <- min(c(n0x,n0y))
  DISTMASK <- DIST>=mi & DIST <= ma
  
  tmp <- fft(mat)
  class(tmp) <- "matrix"
  tmpshift <- myfftshift(tmp)
  tmpshift[n0x,n0y] <- 0
  aspectr2D <- abs(tmpshift)^2 / (n0x*n0y)^4
  
  sig2 <- sum(aspectr2D[DISTMASK]) #Normalisation
  aspectr2D <- aspectr2D/sig2 #Normalisation
  
  # Now calculate r-spectrum
  STEP <- 1
  ray <- seq(mi,ma,STEP)
  rspectr <- numeric(length(ray))
  for (i in 1:length(ray))
  #for ( i in 1:10)
  {
    m <- DIST >= ray[i] - STEP/2 & DIST < ray[i] + STEP/2
    rspectr[i] <- mean(aspectr2D[m])
  }
  
  #out <- data.frame(dist  = ray,  rspec = rspectr);
  out<-rspectr
  return(out)
  
}

myfftshift <- function(X)
{
  nr=dim(X)[1]
  nc=dim(X)[2]
  shiftX = X
  if (nr != nc)
    print("Not a square matrix X")
  else
  {
    n=nc
    shift = floor(n/2)
    for (i in 1:n)
      for (j in 1:n)
      {
        a = (i+shift)%%n
        b = (j+shift)%%n
        if (a==0) a = n
        if (b==0) b = n
        shiftX[a,b] = X[i,j]
      }
  }
  return(shiftX)
}

coarse_grain <- reducedmatrix_ews <- function(mat, subsize) {
  
  N <- dim(mat)[1]
  n <- floor(N/subsize)
  
  submatrix <- matrix(nrow=subsize, ncol=subsize);
  reddata <- matrix(nrow=n, ncol=n);
  
  for (i in 1:n) {
    for (j in 1:n) {
      submatrix = mat[((i-1)*subsize+1):(i*subsize),((j-1)*subsize+1):(j*subsize)];
      reddata[i,j] = mean(as.vector(submatrix));
    }
  }
  
  return(reddata)
}
