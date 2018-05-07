# Three related functions:
#
# InitRandomSamples <- function( I, J, K1, K2 ) {
#
# RoeMetzZSamples <- function( I, J, K1, K2, 
#                              mu_0, delta_mu, VarStr, zSamplesInit, 
#                              isBinned, DesiredNumBins) {
#
# RoeMetzVarStr <- function( rm ) {
#

InitRandomSamples <- function( I, J, K1, K2 ) 
{
  max1 <- max( c( K1, K2 ) )
  
  R_jt <- array( rnorm( J * 2 ), c( J, 2 ) )
  C_kt <- array( rnorm( max1 * 2 ), c( max1, 2 ) )
  TR_ijt <- array( rnorm( I * J * 2 ), c( I, J, 2 ) )
  TC_ikt <- array( rnorm( I * max1 * 2 ), c( I, max1, 2 ) )
  RC_jkt <- array( rnorm( J * max1 * 2 ), c( J, max1, 2 ) )
  EP_ijkt <- array( rnorm( I * J * max1 * 2 ), c( I, J, max1, 2 ) )
  
  return( list (
    R_jt = R_jt,
    C_kt = C_kt,
    TR_ijt = TR_ijt,
    TC_ikt = TC_ikt,
    RC_jkt = RC_jkt,
    EP_ijkt = EP_ijkt
  ) )
}



RoeMetzZSamples <- function( 
  I, J, K1, K2, 
  mu_0, delta_mu, VarStr, zSamplesInit, 
  isBinned, DesiredNumBins) 
{
  
  STD_R  <- sqrt( VarStr$var_r )
  STD_TR <- sqrt( VarStr$var_tr )
  STD_C  <- sqrt( VarStr$var_c )
  STD_TC <- sqrt( VarStr$var_tc )
  STD_RC <- sqrt( VarStr$var_rc )
  STD_EP <- sqrt( VarStr$var_ep )
  
  R_jt <- zSamplesInit$R_jt
  C_kt <- zSamplesInit$C_kt
  TR_ijt <- zSamplesInit$TR_ijt
  TC_ikt <- zSamplesInit$TC_ikt
  RC_jkt <- zSamplesInit$RC_jkt
  EP_ijkt <- zSamplesInit$EP_ijkt
  
  mu <- array( 0, 2 )
  mu[ 2 ] <- mu_0
  
  tau <- array( 0, c( 2, 2 ) )
  tau[ 2, 2 ] <- delta_mu
  
  K <- c( K1, K2 )
  
  FP <- array( -Inf, c( I, J, K1 ) )
  TP <- array( -Inf, c( I, J, K2 ) )
  
  for( i in 1 : I ) {
    for( t in 1 : 2 ) {
      for( j in 1 : J ) {
        for( k in 1 : K[ t ] ) {
          temp <- 
            R_jt[ j, t ] * STD_R +
            C_kt[ k, t ] * STD_C +
            TR_ijt[ i, j, t ] * STD_TR +
            TC_ikt[ i, k, t ] * STD_TC +
            RC_jkt[ j, k, t ] * STD_RC +
            EP_ijkt[ i, j, k, t ] * STD_EP
          
          r1 <- mu[ t ] + tau[ i, t ] + temp ;				
          
          if( t == 1 ) {
            FP[ i, j, k ] <- r1
          }
          else {
            TP[ i, j, k ] <- r1 
          }
        }
      }
    }
  }
  
  #   if( isBinned ) {
  #     FPb <- FP
  #     TPb <- TP
  #     for (i in 1:I){
  #       for (j in 1:J){
  #         ret <- ToIntegerRatings.R(FP[i, j, ],TP[i, j, ],DesiredNumBins)
  #         FPb[i,j,1:K1] <- ret$f[1:K1]   
  #         TPb[i,j,1:K2] <- ret$t[1:K2]   
  #       }
  #     }
  #     FP <- FPb
  #     TP <- TPb
  #   }
  
  return( list(
    FP = FP,
    TP = TP
  ) )
}

RoeMetzVarStr <- function( rm ) {
  rmSet <- c("HL1", "HL2", "HL3", "LL1", "LL2", "LL3", "HH1", "HH2", "HH3", "LH1", "LH2", "LH3")
  if (rm %in% rmSet){
    Var <- switch(rm,
                  "HL1" = {list(
                    var_r = 0.0055,
                    var_tr = 0.0055,
                    var_c = 0.3,
                    var_tc = 0.3,
                    var_rc = 0.2,
                    var_ep = 0.2,
                    auc = 0.702,
                    varStr = rm
                  )},
                  "HL2" = {list(
                    var_r = 0.0055,
                    var_tr = 0.0055,
                    var_c = 0.3,
                    var_tc = 0.3,
                    var_rc = 0.2,
                    var_ep = 0.2,
                    auc = 0.855,
                    varStr = rm
                  )},
                  "HL3" = {list(
                    var_r = 0.0055,
                    var_tr = 0.0055,
                    var_c = 0.3,
                    var_tc = 0.3,
                    var_rc = 0.2,
                    var_ep = 0.2,
                    auc = 0.962,
                    varStr = rm
                  )},
                  "LL1" = list(
                    var_r = 0.0055,
                    var_tr = 0.0055,
                    var_c = 0.1,
                    var_tc = 0.1,
                    var_rc = 0.2,
                    var_ep = 0.6,
                    auc = 0.702,
                    varStr = rm
                  ),
                  "LL2" = {list(
                    var_r = 0.0055,
                    var_tr = 0.0055,
                    var_c = 0.1,
                    var_tc = 0.1,
                    var_rc = 0.2,
                    var_ep = 0.6,
                    auc = 0.855,
                    varStr = rm
                  )},
                  "LL3" = {list(
                    var_r = 0.0055,
                    var_tr = 0.0055,
                    var_c = 0.1,
                    var_tc = 0.1,
                    var_rc = 0.2,
                    var_ep = 0.6,
                    auc = 0.962,
                    varStr = rm
                  )},
                  "HH1" = {list(
                    var_r = 0.011,
                    var_tr = 0.011,
                    var_c = 0.3,
                    var_tc = 0.3,
                    var_rc = 0.2,
                    var_ep = 0.2,
                    auc = 0.702,
                    varStr = rm
                  )},
                  "HH2" = {list(
                    var_r = 0.03,
                    var_tr = 0.03,
                    var_c = 0.3,
                    var_tc = 0.3,
                    var_rc = 0.2,
                    var_ep = 0.2,
                    auc = 0.855,
                    varStr = rm
                  )},
                  "HH3" = {list(
                    var_r = 0.056,
                    var_tr = 0.056,
                    var_c = 0.3,
                    var_tc = 0.3,
                    var_rc = 0.2,
                    var_ep = 0.2,
                    auc = 0.962,
                    varStr = rm
                  )},
                  "LH1" = {list(
                    var_r = 0.011,
                    var_tr = 0.011,
                    var_c = 0.1,
                    var_tc = 0.1,
                    var_rc = 0.2,
                    var_ep = 0.6,
                    auc = 0.702,
                    varStr = rm
                  )},
                  "LH2" = {list(
                    var_r = 0.03,
                    var_tr = 0.03,
                    var_c = 0.1,
                    var_tc = 0.1,
                    var_rc = 0.2,
                    var_ep = 0.6,
                    auc = 0.855,
                    varStr = rm
                  )},
                  "LH3" = {list(
                    var_r = 0.056,
                    var_tr = 0.056,
                    var_c = 0.1,
                    var_tc = 0.1,
                    var_rc = 0.2,
                    var_ep = 0.6,
                    auc = 0.962,
                    varStr = rm
                  )}
    )
  }else{
    Var <- NULL
  }
  return(Var)
}

