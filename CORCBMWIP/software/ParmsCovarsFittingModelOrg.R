#
# fits single modality single reader data to ROCFIT or paired data using CORROC2, as appropriate
#
ParmsCovarsFittingModelOrg <- function(DesiredNumBins, z1, z2)
{
  I <- length(z2[,1,1])
  J <- length(z2[1,,1])
  K1 <- length(z1[1,1,])
  K2 <- length(z2[1,1,])
  # matrix order in FORTRAN code
  # COMMON/BLK7/AX,BX,AY,BY,R,RHO,T(11),U(11)
  #DIM <- 6+2*(DesiredNumBins-1)#ax,ay,bx,by,rhon,rhos,ctffsx[], ctffsy[] 
  done <- array(dim = c(I,I,J,J))
  Chi <- array(dim = c(I,I,J,J,6))
  SigmaChi <- array(dim = c(I,I,J,J,6,6))
  Ax <- array(dim = c(I,I,J,J))
  Ay <- array(dim = c(I,I,J,J))
  RhoAxAy <- array(dim = c(I,I,J,J))
  #### cycle through all combinations of modality and / or reader pairings
  for (i1 in 1:I) {
    for (i2 in 1:I) {        
      for (j1 in 1:J) {
        for (j2 in 1:J) {
          cat("########################", "i1,i2,j1,j2" = i1,i2,j1,j2,"########################\n")
          if (!is.na(done[i1,i2,j1,j2])) next
          if ((i1 ==i2) && (j1 ==j2)) {
            x1 <- z1[i1,j1,]
            x2 <- z2[i1,j1,]
            ret <- ToIntegerRatings(x1,x2,DesiredNumBins)
            K1b <- ret$fb   
            K2b <- ret$tb
            if (length(K1b) != length(K2b)) stop("should not get here")
            Kb <- array(0, dim=c(2,length(K1b)))
            Kb[1,] <- K1b
            Kb[2,] <- K2b
            ret <- RocfitR(Kb)
            if (length(ret) == 1) {
              done[i1,i2,j1,j2] <- 0 # Rocfit failed
              next 
            }
            Chi[i1,i2,j1,j2,1] <- ret$a
            Chi[i1,i2,j1,j2,2] <- ret$b 
            Chi[i1,i2,j1,j2,3] <- ret$a
            Chi[i1,i2,j1,j2,4] <- ret$b 
            Chi[i1,i2,j1,j2,5] <- 1 
            Chi[i1,i2,j1,j2,6] <- 1 
            Ax[i1,i2,j1,j2] <- ret$Az
            Ay[i1,i2,j1,j2] <- ret$Az
            RhoAxAy[i1,i2,j1,j2] <- 1
            SigmaChi[i1,i2,j1,j2,1:2,1:2] <- ret$Cov[1:2,1:2]
            done[i1,i2,j1,j2] <- 1
          } else {
            # need to apply same binning rule to both modalities
            x1 <- c(z1[i1,j1,],z1[i2,j2,]) # concatenate both modalities into one array
            x2 <- c(z2[i1,j1,],z2[i2,j2,]) # do:
            ret <- ToIntegerRatings(x1,x2,DesiredNumBins)
            z1b <- array(dim = c(2,K1));z2b <- array(dim = c(2,K2)) # to hold integer ratings
            z1b[1,] <- ret$f[1:K1]   # un-concatenate
            z1b[2,] <- ret$f[(K1+1):(K1+K1)] # do:
            z2b[1,] <- ret$t[1:K2]   # do:
            z2b[2,] <- ret$t[(K2+1):(K2+K2)] # do:
            ret <- RunCorrocOnPairedData (z1b,z2b,DesiredNumBins)
            Sys.sleep(1)
            if (ret$done == 1) {
              done[i1,i2,j1,j2] <- 1
              done[i2,i1,j2,j1] <- 1
              Chi[i1,i2,j1,j2,] <- ret$parms                    
              Chi[i2,i1,j2,j1,] <- ret$parms
              SigmaChi[i1,i2,j1,j2,,] <- ret$Cov[1:6,1:6]
              SigmaChi[i2,i1,j2,j1,,] <- ret$Cov[1:6,1:6]
              Ax[i1,i2,j1,j2] <- ret$Ax
              Ay[i1,i2,j1,j2] <- ret$Ay
              RhoAxAy[i1,i2,j1,j2] <- ret$RhoAxAy
              Ax[i2,i1,j2,j1] <- ret$Ay
              Ay[i2,i1,j2,j1] <- ret$Ax
              RhoAxAy[i2,i1,j2,j1] <- ret$RhoAxAy
            } else {
              done[i1,i2,j1,j2] <- 0
              done[i2,i1,j2,j1] <- 0
            }
          }
        }
      }
    }
  }
  
  return(list(
    done = done,
    Chi = Chi,
    SigmaChi = SigmaChi,
    Ax = Ax,
    Ay = Ay,
    RhoAxAy = RhoAxAy
  ))
  
}