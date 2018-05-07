SaveJafrocFile <- function( K, ROI_Ratings, isDiseased, filename) 
  
{
  N_TREATMENTS <- length( ROI_Ratings[ , 1, 1, 1, 1 ] )
  MAX_READERS <- length( ROI_Ratings[ 1, , 1, 1, 1 ] )
  K1 <- K[1]
  K2 <- K[2]
  Q <- length( ROI_Ratings[ 1, 1, 1, 1,  ] )
  
  Truth_Table <- array( 0, 0 )
  len <- 1
  for( k in 1 : K1 ) {
    elem <- array( 0, 3 )
    elem[ 1 ] <- k
    Truth_Table <- array( c( Truth_Table, elem ), c( 3, len ) )
    len <- len + 1
  }
  
  for( k in 1 : K2 ) {
    for( r in 1 : Q ) {
      if (isDiseased[k , r] == 1){
        elem <- array( 0, 3 )
        elem[ 1 ] <- k + K1 #case ID
        elem[ 2 ] <- r              #lesion ID
        elem[ 3 ] <- 0  # weight
        Truth_Table <- array( c( Truth_Table, elem ), c( 3, len ) )
        len <- len + 1
      }
    }
  }
  Truth_Table <- data.frame(cbind("CaseID" = Truth_Table[1,], "LesionID" = Truth_Table[2,], "Weights" = Truth_Table[3,]))
  write.xlsx(Truth_Table, file = filename, sheetName = "Truth", col.names = TRUE, row.names = FALSE )
  
  FP_Table <- array( 0, 0 )
  TP_Table <- array( 0, 0 )
  len <- 1
  lenp <- 1
  for (i in 1 : N_TREATMENTS) { 
    for (j in 1 : MAX_READERS) {   
      for( k in 1 : K1 ) {
        for (r in 1 : Q){ 
          elem <- array( 0, 4 )
          elem[ 1 ] <- j
          elem[ 2 ] <- i
          elem[ 3 ] <- k
          elem[ 4 ] <- ROI_Ratings[i, j, k, 1, r]
          FP_Table <- array( c( FP_Table, elem ), c( 4, len ) )
          len <- len + 1            
        }  
      }
      for( k in 1 : K2 ) {
        for (r in 1 : Q){  
          if ( isDiseased[k, r] == 0 ){
            elem <- array( 0, 4 )
            elem[ 1 ] <- j
            elem[ 2 ] <- i
            elem[ 3 ] <- k + K1
            elem[ 4 ] <- ROI_Ratings[i, j, k, 2, r]
            FP_Table <- array( c( FP_Table, elem ), c( 4, len ) )
            len <- len + 1
          }else {
            elem <- array( 0, 5 )
            elem[ 1 ] <- j
            elem[ 2 ] <- i
            elem[ 3 ] <- k + K1
            elem[ 4 ] <- r
            elem[ 5 ] <- ROI_Ratings[i, j, k, 2, r]
            TP_Table <- array( c( TP_Table, elem ), c( 5, lenp ) )
            lenp <- lenp + 1
          }
        }        
      }
    }
  } 
  FP_Table <- data.frame(cbind("ReaderID" = FP_Table[1,], "ModalityID" = FP_Table[2,], "CaseID" = FP_Table[3,]) , 
                         "FP_Rating" = signif(FP_Table[4,],6))  
  write.xlsx( FP_Table, file = filename, sheetName = "FP", col.names = TRUE, row.names = FALSE , append = TRUE)
  
  TP_Table <- data.frame(cbind("ReaderID" = TP_Table[1,], "ModalityID" = TP_Table[2,], "CaseID" = TP_Table[3,]) , 
                         "LesionID" = TP_Table[4,] , 
                         "TP_Rating" = signif(TP_Table[5,],6))  
  write.xlsx( TP_Table, file = filename, sheetName = "TP", col.names = TRUE, row.names = FALSE , append = TRUE)
}