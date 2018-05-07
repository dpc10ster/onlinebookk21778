C---------------------------------------------------------------------          
      SUBROUTINE CHECK(ICLASS,KAT,FPF,TPF,ICON,IARY)                            
C---------------------------------------------------------------------          
C                                                                               
C     PURPOSE                                                                   
C       CHECK DEGENERACY OF DATA SET.                                           
C                                                                               
C     DESCRIPTION OF PARAMETERS                                                 
C       ICLASS    -  CLASS NUMBER OF DEGENERACY                                 
C       KAT       -  NUMBER OF CATEGORIES OF DATA SET                           
C       FPF, TPF  -  OBSERVED OPERATING POINTS                                  
C       ICON,IARY -  STORAGE SPACE FOR INFORMATION TO OUPUT MESSAGE,            
C                    FOR USE IN SUBROUTINE 'DEGENE' IF NECESSARY.               
C                                                                               
      REAL  FPF(1),TPF(1)                                                       
      INTEGER  IARY(1)                                                          
      ICLASS=0                                                                  
      NN=KAT+1                                                                  
      ICON=0                                                                    
      DO 10 I=2,KAT                                                             
      T=FPF(I)*TPF(I)*(1.-FPF(I))*(1.-TPF(I))                                   
      IF (T.EQ.0.) GO TO 10                                                     
      ICON=ICON+1                                                               
      IARY(ICON)=I                                                              
10    CONTINUE                                                                  
      IF (KAT.EQ.1) GO TO 20                                                    
      IF (ICON.EQ.0) GO TO 40                                                   
      I0=IARY(1)-1                                                              
      I1=IARY(ICON)+1                                                           
      IF (KAT.EQ.2.AND.ICON.NE.0) GO TO 30                                      
      IF (ICON.GT.1) GO TO 15                                                   
      IF ((TPF(I0).NE.1.0.OR.TPF(I1).NE.0.0).AND.(FPF(I0).EQ.1.AND.             
     1    FPF(I1).EQ.0)) ICLASS=5                                               
      IF ((FPF(I0).NE.1.0.OR.FPF(I1).NE.0.0).AND.(TPF(I0).EQ.1.AND.             
     1    TPF(I1).EQ.0)) ICLASS=6                                               
      RETURN                                                                    
15    IF ((TPF(IARY(1)).EQ.TPF(IARY(ICON))).AND.                                
     1    (FPF(I0).EQ.1.AND.FPF(I1).EQ.0)) ICLASS=7                             
      IF ((TPF(I0).EQ.1.AND.TPF(I1).EQ.0).AND.                                  
     1    (FPF(IARY(1)).EQ.FPF(IARY(ICON)))) ICLASS=8                           
      RETURN                                                                    
20    ICLASS=1                                                                  
      RETURN                                                                    
30    ICLASS=2                                                                  
      RETURN                                                                    
40    ICA=0                                                                     
      ICB=0                                                                     
      DO 50 I=1,NN                                                              
      IF (FPF(I).NE.0.AND.FPF(I).NE.1) ICA=ICA+1                                
      IF (TPF(I).NE.0.AND.TPF(I).NE.1) ICB=ICB+1                                
50    CONTINUE                                                                  
      I1=0                                                                      
      I2=0                                                                      
      DO 60 I=2,KAT                                                             
      IF (FPF(I).EQ.1.AND.TPF(I).NE.1) GO TO 55                                 
      IF (FPF(I).EQ.0.AND.TPF(I).NE.0) I2=I2+1                                  
      IF (I2.EQ.1) IARY(2)=I                                                    
      GO TO 60                                                                  
55    I1=I1+1                                                                   
      IARY(1)=I                                                                 
60    CONTINUE                                                                  
      II1=0                                                                     
      II2=0                                                                     
      DO 70 I=2,KAT                                                             
      IF (FPF(I).NE.1.AND.TPF(I).EQ.1) GO TO 65                                 
      IF (FPF(I).NE.0.AND.TPF(I).EQ.0) II2=II2+1                                
      IF (II2.EQ.1) IARY(4)=I                                                   
      GO TO 70                                                                  
65    II1=II1+1                                                                 
      IARY(3)=I                                                                 
70    CONTINUE                                                                  
      IF ((ICA.NE.0.AND.ICB.NE.0).OR.((ICA.EQ.0.AND.(I1.EQ.0.OR.                
     1   I2.EQ.0)).OR.(ICB.EQ.0.AND.(II1.EQ.0.OR.II2.EQ.0))))                   
     2   GO TO 80                                                               
      IF (ICA.EQ.0.AND.I1.GE.1.AND.I2.GE.1) GO TO 75                            
      IF (ICB.EQ.0.AND.II1.GE.1.AND.II2.GE.1) ICLASS=4                          
      IARY(1)=IARY(3)                                                           
      IARY(2)=IARY(4)                                                           
      RETURN                                                                    
75    ICLASS=3                                                                  
      RETURN                                                                    
80    KK1=IARY(3)+1                                                             
      K1=IARY(1)+1                                                              
      KK2=IARY(4)-1                                                             
      K2=IARY(2)-1                                                              
      IF ((II1.NE.0.AND.FPF(KK1).EQ.0).OR.(I2.NE.0.AND.                         
     1     TPF(K2).EQ.1)) ICLASS=9                                              
      IF ((I1.NE.0.AND.TPF(K1).EQ.0).OR.(II2.NE.0.AND.                          
     1     FPF(KK2).EQ.1)) ICLASS=10                                            
      RETURN                                                                    
      END                                                                       
C---------------------------------------------------------------------          
      SUBROUTINE ESTM1(A,B,KAT,FPF,TPF,X,EMN,EMS)                               
C---------------------------------------------------------------------          
C                                                                               
C     PURPOSE                                                                   
C       OBTAIN INITIAL PARAMETER ESTIMATES ASSOCIATED WITH A MARGINAL           
C       DATA SET, FOR USE IN SUBROUTINE 'DORALF'.                               
C                                                                               
C     NOTE:                                                                     
C       THE CONTENT OF THIS SUBROUTINE IS TAKEN FROM THE PROGRAM                
C       'RSCORE II' BY DONALD D. DORFMAN (DEPT. OF PSYCHOLOGY,                  
C       THE UNIVERSITY OF IOWA) WITH ONLY MINOR CHANGES.                        
C                                                                               
C                                                                               
C     DESCRIPTION OF PARAMETERS                                                 
C       A, B      -  ESTIMATES OF 'A' AND 'B' PARAMETERS                        
C       KAT       -  NUMBER OF CATEGORIES                                       
C       FPF, TPF  -  OBSEVED OERATING POINTS                                    
C       X         -  ESTIMATES OF 'CUTOFF' PARAMETERS                           
C       EMN       -  TOTAL NUMBER OF NOISE TRIALS                               
C       EMS       -  TOTAL NUMBER OF S+N TRIALS                                 
C                                                                               
      REAL  FPF(1),TPF(1),X(1),XS(9)                                            
      Q=0.                                                                      
      DO 10 I=2,KAT                                                             
      IF (FPF(I).EQ.0.) FPF(I)=1./(2.*EMN)                                      
      IF (ABS(FPF(I)-1.).LT.1.0E-05) FPF(I)=1-(1./(2.*EMN))                     
      P=FPF(I)                                                                  
      CALL ZDEV(P,Q,D,IER)                                                      
      X(I-1)=Q                                                                  
      IF (TPF(I).EQ.0.) TPF(I)=1./(2.*EMS)                                      
      IF (ABS(TPF(I)-1.).LT.1.0E-05) TPF(I)=1.-(1./(2.*EMS))                    
      P=TPF(I)                                                                  
      CALL ZDEV(P,Q,D,IER)                                                      
      XS(I-1)=Q                                                                 
10    CONTINUE                                                                  
      IZ=KAT-2                                                                  
      DO 20 I=1,IZ                                                              
      IF (X(IZ-I+1).LE.X(IZ-I+2)) X(IZ-I+1)=X(IZ-I+2)+0.1                       
      IF (XS(IZ-I+1).LE.XS(IZ-I+2)) XS(IZ-I+1)=XS(IZ-I+2)+0.1                   
20    CONTINUE                                                                  
      SUMX=0.                                                                   
      SUMY=0.                                                                   
      SUMXY=0.                                                                  
      SUMX2=0.                                                                  
      XK=KAT-1                                                                  
      KK=KAT-1                                                                  
      DO 30 I=1,KK                                                              
      SUMX=SUMX+X(I)                                                            
      SUMY=SUMY+XS(I)                                                           
      SUMXY=SUMXY+XS(I)*X(I)                                                    
      SUMX2=SUMX2+X(I)*X(I)                                                     
30    CONTINUE                                                                  
      XMEAN=SUMX/XK                                                             
      YMEAN=SUMY/XK                                                             
      B=(XK*SUMXY-SUMX*SUMY)/(XK*SUMX2-SUMX**2)                                 
      A=YMEAN-B*XMEAN                                                           
      RETURN                                                                    
      END                                                                       
C---------------------------------------------------------------------          
      SUBROUTINE DORALF(A,B,KAT,EMN,EMS,CATN,CATS,X,VA,VB,VAB,IER)              
C---------------------------------------------------------------------          
C                                                                               
C     PURPOSE                                                                   
C       OBTAIN THE FINAL PARAMETER ESTIMATES ASSOCIATED WITH A MARGINAL         
C       DATA SET, FOR USE AS INITIAL PARAMETER ESTIMATES IN THE MAIN            
C       PROGRAM.                                                                
C                                                                               
C     NOTE:                                                                     
C       THE CONTENT OF THIS SUBROUTINE (AND THE SUBROUTINES THAT IT CALLS)      
C       ARE TAKEN FROM THE PROGRAM 'RSCORE II' BY DONALD D. DORFMAN,            
C       DEPT. OF PSYCHOLOGY, THE UNIVERSITY OF IOWA, WITH ONLY MINOR            
C       MODIFICATIONS.                                                          
C                                                                               
C     DESCRIPTION OF PARAMETERS                                                 
C       A,B       -  ESTIMATES OF 'A' AND 'B' PARAMETERS                        
C       KAT       -  NUMBER OF CATEGORIES                                       
C       EMN       -  TOTAL NUMBER OF NOISE TRIALS                               
C       EMS       -  TOTAL NUMBER OF S+N TRIALS                                 
C       CATN      -  VECTOR: NUMBER OF NOISE TRIAL RESPONSES IN EACH CATEGORY   
C       CATS      -  VECTOR: NUMBER OF S+N TRIAL RESPONSES IN EACH CATEGORY     
C       X         -  ESTIMATES OF 'CUTOFF' PARAMETERS                           
C       VA,VB,VAB -  ESTIMATES OF VARIANCES AND COVARIANCE OF 'A' AND           
C                    'B' (NOT USED IN 'CORROC')                                 
C       IER       -  RESULTING ERROR PARAMETER CODED AS FOLLOWING:              
C                  =K   'IER' CODE OBTAINED FROM MATRIX INVERSION               
C                       SUBROUTINE                                              
C                  =101  EXCEEDS 100 ITERATIONS.                                
C                                                                               
      REAL  CATN(1),CATS(1),X(1),XX(11,11),XP(11,11),                           
     1      PX(11),PXBA(11),DX(11),DXBA(11),DIFFSN(11),DIFFN(11),               
     2      PARX(11),ADJX(11),XXDUM(200)                                        
      KK=KAT-1                                                                  
      NN=KAT+1                                                                  
      SWT=0                                                                     
      LL=0                                                                      
C                                                                               
C     GET INTEGRALS AND DENSITIES                                               
C                                                                               
      DO 7045 J=1,NN                                                            
      DO 7045 K=1,NN                                                            
      XX(J,K)=0.                                                                
 7045 XP(J,K)=0.                                                                
 7050 LL=LL+1                                                                   
      PXBA(1)=0.                                                                
      PXBA(NN)=1.                                                               
      DXBA(1)=0.                                                                
      DXBA(NN)=0.                                                               
      PX(1)=0.                                                                  
      PX(NN)=1.                                                                 
      DX(1)=0.                                                                  
      DX(NN)=0.                                                                 
      IF (LL.GT.1) GO TO 7064                                                   
      DO 7055 I=2,KAT                                                           
      Y=-X(I-1)*B-A                                                             
      CALL NDTR(Y,PXBA(I),DXBA(I))                                              
 7055 CONTINUE                                                                  
      DO 7060 I=2,KAT                                                           
      Y=-X(I-1)                                                                 
      CALL NDTR(Y,PX(I),DX(I))                                                  
 7060 CONTINUE                                                                  
      DO 7062 I=1,KAT                                                           
      DIFFSN(I)=PXBA(I+1)-PXBA(I)                                               
      DIFFN(I)=PX(I+1)-PX(I)                                                    
      IF (DIFFSN(I).LE.5.0E-08) DIFFSN(I)=5.0E-08                               
      IF (DIFFN(I).LE.5.0E-08) DIFFN(I)=5.0E-08                                 
7062  CONTINUE                                                                  
      FUNLIK=0.                                                                 
      DO 7063 I=1,KAT                                                           
      FUNLIK=FUNLIK+CATN(I)*ALOG(DIFFN(I))+CATS(I)*ALOG(DIFFSN(I))              
 7063 CONTINUE                                                                  
C                                                                               
C     GET FIRST PARTIALS                                                        
C     WITH RESPECT TO A                                                         
C                                                                               
 7064 AAA=0.                                                                    
      DO 7065 I=2,KAT                                                           
 7065 AAA=AAA-DXBA(I)*(CATS(I-1)/DIFFSN(I-1)-CATS(I)/DIFFSN(I))                 
C                                                                               
C     WITH RESPECT TO B                                                         
C                                                                               
      BBB=0.                                                                    
      DO 7070 I=2,KAT                                                           
 7070 BBB=BBB-DXBA(I)*X(I-1)*(CATS(I-1)/DIFFSN(I-1)-CATS(I)/                    
     1DIFFSN(I))                                                                
C                                                                               
C     WITH RESPECT TO Z'S                                                       
C                                                                               
      DO 7075 I=2,KAT                                                           
      Q1=-DXBA(I)*B*(CATS(I-1)/DIFFSN(I-1)-CATS(I)/DIFFSN(I))                   
      Q2=-DX(I)*(CATN(I-1)/DIFFN(I-1)-CATN(I)/DIFFN(I))                         
 7075 PARX(I-1)=Q1+Q2                                                           
C                                                                               
C     GET EXPECTED SECOND AND MIXED PARTIALS                                    
C     WITH RESPECT TO A                                                         
C                                                                               
      XX(1,1)=0.                                                                
      DO 7080 I=2,KAT                                                           
 7080 XX(1,1)=XX(1,1)-DXBA(I)*((DXBA(I)-DXBA(I-1))/DIFFSN(I-1)-                 
     1(DXBA(I+1)-DXBA(I))/DIFFSN(I))                                            
      XX(1,1)=EMS*XX(1,1)                                                       
C                                                                               
C     WITH RESPECT TO B                                                         
C                                                                               
      XX(2,2)=0.                                                                
      DO 7095 I=2,KAT                                                           
      D=X(I-1)                                                                  
      IF(I.EQ.KAT)GO TO 7090                                                    
      IF(I.EQ.2)GO TO 7085                                                      
      DD=X(I-2)                                                                 
      DDD=X(I)                                                                  
      GO TO 7095                                                                
 7085 DD=0.                                                                     
      DDD=X(I)                                                                  
      GO TO 7095                                                                
 7090 DD=X(I-2)                                                                 
      DDD=0.                                                                    
 7095 XX(2,2)=XX(2,2)-DXBA(I)*X(I-1)*((DXBA(I)*D-DXBA(I-1)*DD)/                 
     1DIFFSN(I-1)-(DXBA(I+1)*DDD-DXBA(I)*D)/DIFFSN(I))                          
      XX(2,2)=EMS*XX(2,2)                                                       
C                                                                               
C     WITH RESPECT TO A AND B                                                   
C                                                                               
      XX(1,2)=0.                                                                
      DO 7100 I=2,KAT                                                           
 7100 XX(1,2)=XX(1,2)-DXBA(I)*X(I-1)*((DXBA(I)-DXBA(I-1))/                      
     1DIFFSN(I-1)-(DXBA(I+1)-DXBA(I))/DIFFSN(I))                                
      XX(1,2)=EMS*XX(1,2)                                                       
      XX(2,1)=XX(1,2)                                                           
C                                                                               
C     WITH RESPECT TO A AND Z'S                                                 
C                                                                               
      DO 7105 I=2,KAT                                                           
      XX(1,I+1)=-(EMS*B*DXBA(I)*((DXBA(I)-DXBA(I-1))/DIFFSN(I-1)                
     1-(DXBA(I+1)-DXBA(I))/DIFFSN(I)))                                          
 7105 XX(I+1,1)=XX(1,I+1)                                                       
C                                                                               
C     WITH RESPECT TO B AND Z'S                                                 
C                                                                               
      DO 7120 I=2,KAT                                                           
      XIL2=0.                                                                   
      IF(I.GT.2)XIL2=X(I-2)                                                     
      IF(I.EQ.KAT)GO TO 7110                                                    
      CHNG=X(I)                                                                 
      GO TO 7115                                                                
 7110 CHNG=0.0                                                                  
 7115 XX(2,I+1)=-(EMS*DXBA(I)*B*((DXBA(I)*X(I-1)-DXBA(I-1)*XIL2)/               
     1DIFFSN(I-1)-(DXBA(I+1)*CHNG-DXBA(I)*X(I-1))/DIFFSN(I)))                   
 7120 XX(I+1,2)=XX(2,I+1)                                                       
C                                                                               
C     WITH RESPECT TO Z'S AND MIXED WITH RESPECT TO Z'S                         
C                                                                               
      DO 7130 I=2,KAT                                                           
      IF(I.EQ.KAT)GO TO 7125                                                    
      XX(I+1,I+2)=(EMS*DXBA(I)*DXBA(I+1)*B**2)/DIFFSN(I)+(EMN*                  
     1DX(I)*DX(I+1))/DIFFN(I)                                                   
      XX(I+2,I+1)=XX(I+1,I+2)                                                   
 7125 XX(I+1,I+1)=-(EMS*(DXBA(I)*B)**2*(1./DIFFSN(I-1)+1./                      
     1DIFFSN(I))+EMN*DX(I)**2*(1./DIFFN(I-1)+1./DIFFN(I)))                      
 7130 CONTINUE                                                                  
C                                                                               
C     INVERT MATRIX                                                             
C                                                                               
      DO 7135 I=1,NN                                                            
      DO 7135 J=1,NN                                                            
 7135 XX(I,J)=-XX(I,J)                                                          
      CALL MSTR(XX,XXDUM,11,0,1)                                                
      CALL SINV(XXDUM,NN,.001,IER)                                              
      IF(IER.GE.0) GO TO 7137                                                   
      ALPHA=-0.5*ABS(ALPHA)                                                     
      GO TO 7157                                                                
 7137 CALL MSTR(XXDUM,XP,11,1,0)                                                
      IF(SWT.EQ.1)GO TO 7190                                                    
C                                                                               
C     FORM SOLUTION VECTOR                                                      
C                                                                               
      ADJA=AAA*XP(1,1)+BBB*XP(1,2)                                              
      DO 7145 I=1,KK                                                            
 7145 ADJA=ADJA+PARX(I)*XP(1,I+2)                                               
      ADJB=AAA*XP(2,1)+BBB*XP(2,2)                                              
      DO 7150 I=1,KK                                                            
 7150 ADJB=ADJB+PARX(I)*XP(2,I+2)                                               
      DO 7155 I=1,KK                                                            
      ADJX(I)=AAA*XP(I+2,1)+BBB*XP(I+2,2)                                       
      DO 7155 K=3,NN                                                            
 7155 ADJX(I)=ADJX(I)+PARX(K-2)*XP(I+2,K)                                       
C                                                                               
C     ITERATE                                                                   
C                                                                               
      ALPHA=1.0                                                                 
 7157 A=A+ALPHA*ADJA                                                            
      B=B+ALPHA*ADJB                                                            
      DO 7160 I=1,KK                                                            
 7160 X(I)=X(I)+ALPHA*ADJX(I)                                                   
      DO 7165 I=2,KAT                                                           
      Y=-X(I-1)                                                                 
      CALL NDTR(Y,PX(I),DX(I))                                                  
      Y=-X(I-1)*B-A                                                             
      CALL NDTR(Y,PXBA(I),DXBA(I))                                              
 7165 CONTINUE                                                                  
      DO 7170 I=1,KAT                                                           
      DIFFN(I)=PX(I+1)-PX(I)                                                    
      DIFFSN(I)=PXBA(I+1)-PXBA(I)                                               
      IF (DIFFN(I).LT.5.0E-8) DIFFN(I)=5.0E-8                                   
      IF (DIFFSN(I).LT.5.0E-8) DIFFSN(I)=5.0E-8                                 
 7170 CONTINUE                                                                  
      IF (ALPHA.NE.1.0) GO TO 7175                                              
C                                                                               
C     CHECK FOR MAXIMIZATION                                                    
C                                                                               
      IF (LL.EQ.1) GO TO 7175                                                   
      TNET=ABS(ADJA)+ABS(ADJB)                                                  
      DO 7173 I=1,KK                                                            
 7173 TNET=TNET+ABS(ADJX(I))                                                    
      TNET=6.0*TNET/FLOAT(1+KAT)                                                
      IF(TNET.GE..001) GO TO 7175                                               
      SWT=1                                                                     
      GO TO 7050                                                                
 7175 IF (ALPHA.EQ.1.0) FUNLI1=FUNLIK                                           
      FUNLIK=0.                                                                 
      DO 7180 I=1,KAT                                                           
      FUNLIK=FUNLIK+CATN(I)*ALOG(DIFFN(I))+CATS(I)*ALOG(DIFFSN(I))              
 7180 CONTINUE                                                                  
      IF (FUNLIK.GE.FUNLI1) GO TO 7185                                          
      FUNDIF=ABS((FUNLIK-FUNLI1)/FUNLI1)                                        
      IF (FUNDIF.LE.1.0E-5) GO TO 7185                                          
      ALPHA=-0.5*ABS(ALPHA)                                                     
      GO TO 7157                                                                
 7185 IF(LL.GT.100)GO TO 8430                                                   
      GO TO 7050                                                                
C                                                                               
C     ASSIGN THE FINAL VALUES OF A, B, AND THE CUTOFFS OBTAINED FROM            
C     THE MARGINAL DATA SET TO THE RESPECTIVE VARIABLE NAMES THAT               
C     WILL BE USED IN THE MAIN PROGRAM TO ACCOUNT FOR CORRELATIONS              
C     BETWEEN THE MARGINAL DATA SETS.  THE VALUES OF THESE NEW VARIABLES        
C     ARE FIRST USED AS INITIAL ESTIMATES IN THE MAIN PROGRAM.                  
C                                                                               
C     THE VARIANCES AND COVARIANCE OF THE 'A' AND 'B' PARAMETERS                
C     ESTIMATED FROM THE MARGINAL DATA SET ARE RETURNED ALSO, BUT               
C     ARE TO BE USED FOR COMPARATIVE PURPOSES ONLY --- WITH EACH                
C     OTHER AND WITH THE VARIANCES AND COVARIANCE OF 'A' AND 'B' FOR            
C     THE OTHER MARGINAL DATA SET, WHEN THE TWO SETS OF ROC DATA ARE            
C     TREATED AS CORRELATED IN THE MAIN PROGRAM.                                
C                                                                               
 7190 VA=XP(1,1)                                                                
      VB=XP(2,2)                                                                
      VAB=XP(1,2)                                                               
      RETURN                                                                    
C                                                                               
C     THE FOLLOWING TWO SETS OF WRITE AND FORMAT STATEMENTS WILL BE             
C     USED ONLY IF THERE ARE ANY PROBLEMS WITH ITERATION ON THE MARGINAL        
C     DATA SETS IN SUBROUTINE 'DORALF'.                                         
C                                                                               
 8430 IER=101                                                                   
      RETURN                                                                    
      END                                                                       
C---------------------------------------------------------------------          
      SUBROUTINE DEGENE(ICLASS,IARY,FPF,TPF)                               
C---------------------------------------------------------------------          
C                                                                               
C     PURPOSE                                                                   
C       OUTPUT DEGENERACY FOR EACH DEGENERACY CLASS.                            
C                                                                               
C     DESCRIPTION OF PARAMETERS                                                 
C       ICLASS    -  CLASS NUMBER OF DEGENERACY                                 
C       IARY      -  INFORMATION REGARDING DEGENERATE DATA SET,                 
C                    FROM SUBROUTINE 'CHECK'                                    
C       FPF, TPF  -  OBSERVED OPERATING POINTS                                  
C                                                                               
      INTEGER IARY(1)                                                           
      REAL  FPF(1),TPF(1)                                                       
      IF (ICLASS.EQ.1) GO TO 10                                                 
      IF (ICLASS.EQ.2) GO TO 20                                                 
      IF (ICLASS.EQ.3) GO TO 30                                                 
      IF (ICLASS.EQ.4) GO TO 40                                                 
      IF (ICLASS.EQ.5.OR.ICLASS.EQ.7) GO TO 50                                  
      IF (ICLASS.EQ.6.OR.ICLASS.EQ.8) GO TO 60                                  
      IF (ICLASS.EQ.9) GO TO 90                                                 
      write(*,5)                                                                
5     FORMAT(15X,'AND IMPLY PERFECT (BUT PERVERSE) DECISION',                   
     1       /15X,'PERFORMANCE.')                                               
      RETURN                                                                    
10    write(*,15)                                                               
15    FORMAT(15X,'AND PROVIDE NO OPERATING POINTS OFF THE (0,0)',               
     1       /15X,'AND (1,1) CORNERS.')                                         
      RETURN                                                                    
20    P=FPF(IARY(1))                                                            
      CALL ZDEV(P,QF,D,IER)                                                     
      P=TPF(IARY(1))                                                            
      CALL ZDEV(P,QT,D,IER)                                                     
      DD=QT-QF                                                                  
      write(*,25) DD                                                            
25    FORMAT(15X,'AND PROVIDE ONLY A SINGLE OPERATING POINT OFF',               
     1       /,15X,'THE (0,0) AND (1,1) CORNERS.  BINORMAL ROC CURVE',          
     2       /15X,'CANNOT BE ESTIMATED UNIQUELY.  UNIT SLOPE ROC CURVE',        
     3       /15X,'(B=1) WOULD HAVE A=',3HD'=,F8.4,'.')                         
30     write(*,35) TPF(IARY(2)),TPF(IARY(1))                                    
35    FORMAT(15X,'AND PRODUCE NO OPERATING POINTS OFF THE BORDERS',             
     1       /15X,'OF THE UNIT SQUARE.  IMPLIED EXACT-FIT BINORMAL',            
     2       /15X,'ROC CURVE IS HORIZONTAL AT UNDETERMINED HEIGHT',             
     3       /15X,'BETWEEN TPF=',F8.4,'AND TPF=',F8.4,'.')                      
      RETURN                                                                    
40    write(*,45) FPF(IARY(2)),FPF(IARY(1))                                     
45    FORMAT(15X,'AND PRODUCE NO OPERATING POINTS OFF THE BORDERS',             
     1       /15X,'OF THE UNIT SQUARE.  IMPLIED EXACT-FIT BINORMAL',            
     2       /15X,'ROC CURVE IS VERTICAL AT UNDETERMINED POSITION',             
     3       /15X,'BETWEEN FPF=',F8.4,'AND FPF=',F8.4,'.')                      
      RETURN                                                                    
50    write(*,55) TPF(IARY(1))                                                  
55    FORMAT(15X,'IMPLIED EXACT-FIT BINORMAL ROC CURVE IS HORIZONTAL',          
     1       /,15X,'AT CONSTANT TPF=',F8.4,'.')                                 
      RETURN                                                                    
60    write(*,65) FPF(IARY(1))                                                  
65    FORMAT(15X,'IMPLIED EXACT-FIT BINORMAL ROC CURVE IS VERTICAL',            
     1       /,15X,'AT CONSTANT FPF=',F8.4,'.')                                 
      RETURN                                                                    
90    write(*,95)                                                               
95    FORMAT(15X,'AND IMPLY PERFECT DECISION PERFORMANCE.')                     
      RETURN                                                                    
      END                                                                       
C------------------------------------------------------------------        
      SUBROUTINE ZDEV(P,X,D,IE)                                                 
C------------------------------------------------------------------        
C                                                                               
C     PURPOSE                                                                   
C       COMPUTES X = P**(-1) (Y), THE ARGUMENT X SUCH THAT Y=P(X)= THE          
C       PROBABILITY THAT THE RANDOM VARIABLE U, DISTRIBUTED NORMALLY            
C       N(0,1), IS LESS THAN OR EQUAL TO X.  F(X), THE ORDINATE OF THE          
C       NORMAL DENSITY, AT X, IS ALSO COMPUTED.                                 
C                                                                               
C     DESCRIPTION OF PARAMETERS                                                 
C       P     - INPUT PROBABILITY.                                              
C       X     - OUTPUT ARGUMENT SUCH THAT P=Y=THE PROBABILITY THAT U,           
C             THE RANDOM VARIABLE, IS LESS THAN OR EQUAL TO X.                  
C       D     - OUTPUT DENSITY, F(X).                                           
C       IER   - OUTPUT ERROR CODE                                               
C               = -1 IF P IS NOT IN THE INTERVAL (0,1), INCLUSIVE.              
C                     X=D=.99999E+38 IN THAT CASE.                              
C               = 0 IF THERE IS NO ERROR.  SEE REMARKS, BELOW.                  
C                                                                               
C     REMARKS                                                                   
C       MAXIMUM ERROR IS 0.00045.                                               
C       IF P=0, X IS SET TO -(10)**38.  D IS SET TO 0.                          
C       IF P=1, X IS SET TO  (10)**38.  D IS SET TO 0.                          
C     NOTE:  ORIGINAL PROGRAM SET X TO + OR -(10)**74.                          
C                                                                               
C     SUBROUTINES AND SUBPROGRAMS REQUIRED.                                     
C     NONE                                                                      
C                                                                               
C     METHOD                                                                    
C       BASED ON APPROXIMATIONS IN C. HASTINGS, APPROXIMATIONS FOR              
C       DIGITAL COMPUTERS, PRINCETON UNIV. PRESS, PRINCETON, N.J., 1955.        
C       SEE EQUATION 26.2.23, HANDBOOK OF MATHEMATICAL FUNCTIONS,               
C       ABRAMOWITZ AND STEGUN, DOVER PUBLICATIONS, INC., NEW YORK.              
C                                                                               
      IE=0                                                                      
C     X=.99999E+74                                                              
      X=.99999E+38                                                              
      D=X                                                                       
      IF(P)1,4,2                                                                
1     IE=-1                                                                     
      GO TO 12                                                                  
2     IF(P-1.0)7,5,1                                                            
C4    X=-.99999E+74                                                             
4     X=-.99999E+38                                                             
5     D=0.0                                                                     
      GO TO 12                                                                  
7     D=P                                                                       
      IF(D-0.5)9,9,8                                                            
8     D=1.0-D                                                                   
9     T2=ALOG(1.0/(D*D))                                                        
      T=SQRT(T2)                                                                
      X=T-(2.515517+0.802853*T+0.010328*T2)/(1.0+1.432788*T+0.189269            
     1*T2+0.001308*T*T2)                                                        
      IF(P-0.5)10,10,11                                                         
10    X=-X                                                                      
11    D=0.3989423*EXP(-X*X/2.0)                                                 
12    RETURN                                                                    
      END                                                                       
C------------------------------------------------------------              
      SUBROUTINE NDTR(X,P,D)                                                    
C------------------------------------------------------------              
C                                                                               
C     PURPOSE                                                                   
C       COMPUTES Y=P(X)=PROBABILITY THAT THE RANDOM VARIABLE                    
C       DISTRIBUTED NORMALLY IN (0,1), IS LESS THAN OR EQUAL                    
C       TO X.  F(X), THE ORDINATE OF THE NORMAL DENSITY AT X,                   
C       IS ALSO COMPUTED.                                                       
C                                                                               
C     DESCRIPTION OF PARAMETERS                                                 
C       X     - INPUT SCALAR FOR WHICH P(X) IS COMPUTED.                        
C       P     - OUTPUT PROBABILITY.                                             
C       D     - OUTPUT DENSITY.                                                 
C                                                                               
C     REMARKS                                                                   
C       MAXIMUM ERROR IS .0000007.                                              
C                                                                               
C     SUBROUTINE AND SUBPROGRAM REQUIRED                                        
C       NONE                                                                    
C                                                                               
C     METHOD                                                                    
C       BASED ON APPROXIMATION IN C. HASTINGS, APPROXIMATIONS                   
C       FOR DIGITAL COMPUTERS, PRINCETON UNIV. PRESS, PRINCETON,                
C       N.J., 1955.  SEE EQUATION 26.2.17, HANDBOOK OF MATHEMATICAL             
C       FUNCTIONS, ABRAMOWITZ AND STEGUN, DOVER PUBLICATION, INC.,              
C       NEW YORK.                                                               
C                                                                               
      AX=ABS(X)                                                                 
      T=1.0/(1.0+0.2316419*AX)                                                  
      IF(AX.GT.18.)GO TO 5                                                      
      D=0.3989423*EXP(-X*X/2.0)                                                 
      GO TO 6                                                                   
5     D=0.0                                                                     
6     P=1.0-D*T*((((1.330274*T-1.821256)*T+1.781478)*T-0.3565638)*T+            
     10.3193815)                                                                
      IF(X)1,2,2                                                                
1     P=1.0-P                                                                   
2     RETURN                                                                    
      END                                                                       
C------------------------------------------------------------------
      SUBROUTINE MDBNOR (X,Y,RHO,P,IER)                                         
C------------------------------------------------------------------
C
C   FUNCTION            - BIVARIATE NORMAL PROBABILITY DISTRIBUTION             
C                           FUNCTION                                            
C   USAGE               - CALL MDBNOR(X,Y,RHO,P,IER)                            
C   PARAMETERS   X      - INPUT  UPPER LIMIT OF INTEGRATION FOR THE             
C                           FIRST VARIABLE                                      
C                Y      - INPUT  UPPER LIMIT OF INTEGRATION FOR THE             
C                           SECOND VARIABLE                                     
C                RHO    - INPUT  CORRELATION COEFFICIENT                        
C                P      - OUTPUT PROBABILITY THAT THE FIRST VARIABLE            
C                           IS LESS THAN OR EQUAL TO X AND THAT THE             
C                           SECOND VARIABLE IS LESS THAN OR EQUAL TO Y          
C                IER    - ERROR PARAMETER                                       
C                         TERMINAL ERROR = 128+N                                
C                         N=1 INDICATES THE ABSOLUTE VALUE OF RHO IS            
C                             GREATER THAN OR EQUAL TO ONE                      
C   PRECISION           - SINGLE                                                
C   REQD. IMSL ROUTINES - MDTNF,UERTST                                          
C                                                                               
C                                  C0 = SQRT(.5)                                
      DATA               C0/.7071068/,C1/.999999/                               
      EPS = 0.0                                                                 
      IER = 0                                                                   
      IF (ABS(RHO) .LT. C1) GO TO 5                                             
C                                  TERMINAL - RHO OUT OF RANGE                  
      IER = 129                                                                 
      GO TO 9000                                                                
    5 F1 = 1.0/SQRT(1.0 - RHO**2)                                               
      XY = X*Y                                                                  
      IAX = 0                                                                   
      IAY = 0                                                                   
      IND = 0                                                                   
      IF (XY .EQ. 0.0) GO TO 10                                                 
      AX = F1*(Y/X - RHO)                                                       
      AY = F1*(X/Y - RHO)                                                       
      GO TO 25                                                                  
   10 IF (X .NE. 0.0) GO TO 15                                                  
      IF (Y .NE. 0.0) GO TO 20                                                  
C                                                                2 1/2          
C                                  FOR X=Y=0 AX=AY=(1-RHO)/(1-RHO )             
      AX = F1*(1.0 - RHO)                                                       
      AY = AX                                                                   
      GO TO 25                                                                  
C                                  FOR Y=0,X LESS THAN 0     TY = -1/4          
C                                  FOR Y=0,X GREATER THAN 0  TY =  1/4          
   15 TY = 0.25                                                                 
      IF (X .LT. 0.0) TY = -TY                                                  
      AX = -F1*RHO                                                              
      IND = 1                                                                   
      GO TO 25                                                                  
C                                  FOR X=0,Y LESS THAN 0     TX = -1/4          
C                                  FOR X=0,Y GREATER THAN 0  TX =  1/4          
   20 TX = 0.25                                                                 
      IF (Y .LT. 0.0) TX = -TX                                                  
      AY = -F1*RHO                                                              
      GO TO 35                                                                  
   25 IF (AX .GE. 0.0) GO TO 30                                                 
      IAX = 1                                                                   
      AX = -AX                                                                  
   30 CALL MDTNF(X,AX,EPS,TX)                                                   
      IF (IAX .NE. 0) TX = -TX                                                  
      IF (IND .NE. 0) GO TO 45                                                  
   35 IF (AY .GE. 0.0) GO TO 40                                                 
      IAY = 1                                                                   
      AY = -AY                                                                  
   40 CALL MDTNF(Y,AY,EPS,TY)                                                   
      IF (IAY .NE. 0) TY = -TY                                                  
   45 CALL NDTR(X,QX,DEN)                                                       
      CALL NDTR(Y,QY,DEN)                                                       
C                                  NOW EVALUATE P                               
      P = 0.5*(QX + QY) - TX - TY                                               
      IF (XY .LE. 0.0 .AND.(XY .NE. 0.0 .OR. X+Y .LT. 0.0)) P = P - 0.5         
      P = AMIN1(AMAX1(0.0,P),1.0)                                               
 9000 CONTINUE                                                                  
      IF (IER .NE. 0) CALL UERTST(IER,6HMDBNOR)                                 
 9005 RETURN                                                                    
      END                                                                       
C--------------------------------------------------------------------
      SUBROUTINE MDTNF  (Y,Z,EPS,T)                                             
C--------------------------------------------------------------------
C                                                                               
C   FUNCTION            - INTEGRATE T(Y,Z) FOR NON-CENTRAL T USAGE.             
C   USAGE               - CALL MDTNF(Y,Z,EPS,T)                                 
C   PARAMETERS   Y      - INPUT PARAMETER.  SEE DOCUMENTATION FOR               
C                           THE DEFINITION.                                     
C                Z      - INPUT.  INTEGRATION IS FROM 0 TO Z.                   
C                EPS    - INPUT.  ACCURACY SHOULD NOT BE LESS THAN EPS.         
C                           IF EPS=0 IS ENTERED, EPS=.000001 IS USED.           
C                T      - RESULTANT VALUE OF THE INTEGRAL.                      
C   PRECISION           - SINGLE                                                
C                                                                               
      DATA               C/.1591549/,EXPOV/174.673/                             
      EP1 = EPS                                                                 
      IF(EPS .EQ. 0.) EP1 = .000001                                             
      T = 0.0                                                                   
      B = ABS(Y)                                                                
      A = ABS(Z)                                                                
      IF(A .EQ. 0.) GO TO 35                                                    
    5 TA = ATAN(A)                                                              
      IF (A*B .LE. 4.0) GO TO 10                                                
      CALL NDTR(B,T,DEN)                                                        
      T = C*(TA+ATAN(1.0/A)) - .5*(T-.5)                                        
      GO TO 30                                                                  
C                                  APPROXIMATION FOR SMALL Y*Z                  
   10 HSQB = .5*B*B                                                             
      IF (HSQB .GT. EXPOV) GO TO 35                                             
      BEXP = EXP(-HSQB)                                                         
      ASQ = A*A                                                                 
      A4 = ASQ*ASQ                                                              
      B4 = HSQB * HSQB                                                          
      A4B4 = A4 * B4                                                            
      AHSQB = A * HSQB                                                          
      AB4 = A*B4*.5                                                             
      F = 1.0                                                                   
      SUM = 0.0                                                                 
      G = 3.0                                                                   
C                                  BEGIN SERIES EXPANSION                       
   15 G1 = G                                                                    
      BER = 0.0                                                                 
      TER = AB4                                                                 
   20 BER = BER+TER                                                             
      IF(TER .LE. BER*EP1) GO TO 25                                             
C                                  DEVELOP COEFFICIENT SERIES                   
      TER = TER*HSQB/G1                                                         
      G1 = G1+1.0                                                               
      GO TO 20                                                                  
   25 D1 = (BER+AHSQB)/F                                                        
      D2 = BER*ASQ/(F+2.0)                                                      
      D = D1-D2                                                                 
      SUM = SUM+D                                                               
      T = TA-SUM*BEXP                                                           
      AEPS = EP1*T                                                              
      AHSQB = AHSQB*A4B4/((G-1.0)*G)                                            
      AB4 = AB4*A4B4/((G +1.0)*G)                                               
      F = F+4.0                                                                 
      G = G+2.0                                                                 
C                                  SHOULD SERIES EXPANSION BE TERMINATED        
      IF (D2*BEXP .GE. AEPS) GO TO 15                                           
      T = T * C                                                                 
   30 IF(Z .LT. 0.0) T = -T                                                     
   35 RETURN                                                                    
      END                                                                       
C----------------------------------------------------------------------
      SUBROUTINE UERTST (IER,NAME)                                              
C----------------------------------------------------------------------
C                                                                               
C   FUNCTION            - ERROR MESSAGE GENERATION                              
C   USAGE               - CALL UERTST(IER,NAME)                                 
C   PARAMETERS   IER    - ERROR PARAMETER. TYPE + N  WHERE                      
C                           TYPE= 128 IMPLIES TERMINAL ERROR                    
C                                  64 IMPLIES WARNING WITH FIX                  
C                                  32 IMPLIES WARNING                           
C                           N   = ERROR CODE RELEVANT TO CALLING ROUTINE        
C                NAME   - INPUT VECTOR CONTAINING THE NAME OF THE               
C                           CALLING ROUTINE AS A SIX CHARACTER LITERAL          
C                           STRING.                                             
C                                                                               
      DIMENSION          ITYP(5,4),IBIT(4)                                      
      INTEGER*2          NAME(3)                                                
      INTEGER            WARN,WARF,TERM,PRINTR                                  
      REAL ITYP                                                                 
      EQUIVALENCE        (IBIT(1),WARN),(IBIT(2),WARF),(IBIT(3),TERM)           
      DATA     ITYP      /'WARN','ING ','    ','    ','    ',                   
     *                    'WARN','ING(','WITH',' FIX',')   ',                   
     *                    'TERM','INAL','    ','    ','    ',                   
     *                    'NON-','DEFI','NED ','    ','    '/,                  
     *         IBIT      / 32,64,128,0/                                         
      DATA     PRINTR    / 6/                                                   
      IER2=IER                                                                  
      IF (IER2 .GE. WARN) GO TO 5                                               
C                                  NON-DEFINED                                  
      IER1=4                                                                    
      GO TO 20                                                                  
   5  IF (IER2 .LT. TERM) GO TO 10                                              
C                                  TERMINAL                                     
      IER1=3                                                                    
      GO TO 20                                                                  
  10  IF (IER2 .LT. WARF) GO TO 15                                              
C                                  WARNING(WITH FIX)                            
      IER1=2                                                                    
      GO TO 20                                                                  
C                                  WARNING                                      
  15  IER1=1                                                                    
C                                  EXTRACT 'N'                                  
  20  IER2=IER2-IBIT(IER1)                                                      
C                                  PRINT ERROR MESSAGE                          
      WRITE (PRINTR,25) (ITYP(I,IER1),I=1,5),NAME,IER2,IER                      
  25  FORMAT(' *** I M S L(UERTST) ***  ',5A4,4X,3A2,4X,I2,                     
     *   ' (IER = ',I3,')')                                                     
      RETURN                                                                    
      END                                                                       
C----------------------------------------------------------------------
      SUBROUTINE MSTR(A,R,N,MSA,MSR)                                            
C----------------------------------------------------------------------
C                                                                               
C     PURPOSE                                                                   
C       CHANGE STORAGE MODE OF A MATRIX                                         
C                                                                               
C     DESCRIPTION OF PARAMETERS                                                 
C       A - NAME OF INPUT MATRIX                                                
C       R - NAME OF OUTPUT MATRIX                                               
C       N - NUMBER OF ROWS AND COLUMNS IN A AND R                               
C     MSA - ONE DIGIT NUMBER FOR STORAGE MODE OF MATRIX A                       
C             0 - GENERAL                                                       
C             1 - SYMMETRIC                                                     
C             2 - DIAGONAL                                                      
C     MSR - SAME AS MSA EXCEPT FOR MATRIX R                                     
C                                                                               
C     REMARKS                                                                   
C     MATRIX R CANNOT BE IN THE SAME LOCATION AS MATRIX A                       
C     MATRIX A MUST BE A SQUARE MATRIX                                          
C                                                                               
C     SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                             
C       LOC                                                                     
C                                                                               
C     METHOD                                                                    
C       MATRIX A IS RESTRUCTURED TO FORM MATRIX R                               
C     MSA MSR                                                                   
C      0   0  MATRIX A IS MOVED TO MATRIX R                                     
C      0   1  THE UPPER TRIANGLE OF ELEMENTS OF A GENERAL MATRIX ARE            
C             USED TO FORM A SYMMETRIC MATRIX                                   
C      0   2  THE DIAGONAL ELEMENTS OF A GENERAL MATRIX ARE USED TO FORM        
C             A DIAGONAL MATRIX                                                 
C      1   0  A SYMMETRIC MATRIX IS EXPANDED TO FORM A GENERAL MATRIX           
C      1   1  MATRIX A IS MOVED TO MATRIX R                                     
C      1   2  THE DIAGONAL ELEMENTS OF A SYMMETRIC MATRIX ARE USED TO           
C             FORM A DIAGONAL MATRIX                                            
C      2   0  A DIAGONAL MATRIX IS EXPANDED BY INSERTING MISSING ZERO           
C             ELEMENTS TO FORM A GENERAL MATRIX                                 
C      2   1  A DIAGONAL MATRIX IS EXPANDED BY INSERTING MISSING ZERO           
C             ELEMENTS TO FORM A SYMMETRIC MATRIX                               
C      2   2  MATRIX A IS MOVED TO MATRIX R                                     
C                                                                               
      DIMENSION A(1),R(1)                                                       
      DO 20 I=1,N                                                               
      DO 20 J=1,N                                                               
C                                                                               
C     IF R IS GENERAL, FORM ELEMENT                                             
C                                                                               
      IF(MSR)5,10,5                                                             
C                                                                               
C     IF IN LOWER TRIANGLE OF SYMMETRIC OR DIAGONAL R, BYPASS                   
C                                                                               
5     IF(I-J)10,10,20                                                           
10    CALL LOC(I,J,IR,N,MSR)                                                  
C                                                                               
C     IF IN UPPER AND OFF DIAGONAL OF DIAGONAL R, BYPASS                        
C                                                                               
      IF(IR)20,20,15                                                            
C                                                                               
C     OTHERWISE FORM R(I,J)                                                     
C                                                                               
15    R(IR)=0.0                                                                 
      CALL LOC(I,J,IA,N,MSA)                                                  
C                                                                               
C     IF THERE IS NO A(I,J), LEAVE R(I,J) AT 0                                  
C                                                                               
      IF(IA)20,20,18                                                            
18    R(IR)=A(IA)                                                               
20    CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
C--------------------------------------------------------------------
      SUBROUTINE LOC(I,J,IR,N,MS)                                             
C--------------------------------------------------------------------
C                                                                               
C     PURPOSE                                                                   
C       COMPUTE A VECTOR SUBSCRIPT FOR AN ELEMENT IN A MATRIX OF                
C       SPECIFIED STORAGE MODE.                                                 
C                                                                               
C     DESCRIPTION OF PARAMETERS                                                 
C       I     - ROW NUMBER OF ELEMENT                                           
C       J     - COLUMN NUMBER OF ELEMENT                                        
C       IR    - RESULTANT VECTOR SUBSCRIPT                                      
C       N     - NUMBER OF ROWS IN MATRIX                                        
C       MS    - ONE DIGIT NUMBER FOR STORAGE MODE OF MATRIX                     
C                     0 - GENERAL                                               
C                     1 - SYMMETRIC                                             
C                     2 - DIAGONAL                                              
C                                                                               
C     SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                             
C     NONE                                                                      
C                                                                               
C     METHOD                                                                    
C     MS=0    SUBSCRIPT IS COMPUTED FOR A MATRIX WITH N*M ELEMENTS IN           
C             STORAGE (GENERAL MATRIX)                                          
C     MS=1    SUBSCRIPT IS COMPUTED FOR A MATRIX WITH N*(N+1)/2 IN              
C             STORAGE (UPPER TRIANGLE OF SYMMETRIC MATRIX).  IF ELEMENT         
C             IS IN LOWER TRIANGULAR PORTION, SUBSCRIPT IS CORRESPONDING        
C             ELEMENT IN UPPER TRIANGLE.                                        
C     MS=2    SUBSCRIPT IS COMPUTED FOR A MATRIX WITH N ELEMENTS IN             
C             STORAGE (DIAGONAL ELEMENTS OF DIAGONAL MATRIX).  IF               
C             ELEMENT IS NOT ON DIAGONAL (AND THEREFORE NOT IN STORAGE),        
C             IR IS SET TO ZERO.                                                
C                                                                               
      IX=I                                                                      
      JX=J                                                                      
      IF(MS-1)10,20,30                                                          
10    IRX=N*(JX-1)+IX                                                           
      GO TO 36                                                                  
20    IF(IX-JX)22,24,24                                                         
22    IRX=IX+(JX*JX-JX)/2                                                       
      GO TO 36                                                                  
24    IRX=JX+(IX*IX-IX)/2                                                       
      GO TO 36                                                                  
30    IRX=0                                                                     
      IF(IX-JX)36,32,36                                                         
32    IRX=IX                                                                    
36    IR=IRX                                                                    
      RETURN                                                                    
      END                                                                       
C----------------------------------------------------------------------
      SUBROUTINE SINV(A,N,EPS,IER)                                              
C----------------------------------------------------------------------
C                                                                               
C     PURPOSE                                                                   
C       INVERT A GIVEN SYMMETRIC POSITIVE DEFINITE MATRIX                       
C                                                                               
C     USAGE                                                                     
C       CALL SINV(A,N,EPS,IER)                                                  
C                                                                               
C     DESCRIPTION OF PARAMETERS                                                 
C       A     - UPPER TRIANGULAR PART OF THE GIVEN SYMMETRIC POSITIVE           
C               DEFINITE N BY N COEFFICIENT MATRIX.  ON RETURN A                
C               CONTAINS THE RESULTANT UPPER TRIANGULAR MATRIX.                 
C       N     - THE NUMBER OF ROW (COLUMNS) IN GIVEN MATRIX.                    
C       EPS   - AN INPUT CONSTANT WHICH IS USED AS RELATIVE TOLERANCE           
C               FOR TEST ON LOSS OF SIGNIFICANCE.                               
C       IER   - RESULTING ERROR PARAMETER CODED AS FOLLOWS:                     
C               IER=0  - NO ERROR                                               
C               IER=-1 - NO RESULT BECAUSE OF WRONG INPUT PARAMETER N OR        
C                        BECAUSE SOME RADICAND IS NONPOSITIVE (MATRIX A         
C                        IS NOT POSITIVE DEFINITE, POSSIBLY DUE TO LOSS         
C                        OF SIGNIFICANCE)                                       
C               IER=K  - WARNING WHICH INDICATES LOSS OF SIGNIFICANCE.          
C                        THE RADICAND FORMED AT FACTORIZATION STEP K+1          
C                        WAS STILL POSITIVE BUT NO LONGER GREATER THAN          
C                        ABS(EPS*A(K+1,K+1)).                                   
C                                                                               
C     REMARKS                                                                   
C       THE UPPER TRIANGULAR PART OF GIVEN MATRIX IS ASSUMED TO BE              
C       STORED COLUMNWISE IN N*(N+1)/2 SUCCESSIVE STORAGE LOCATIONS.  IN        
C       THE SAME STORAGE LOCATIONS THE RESULTING UPPER TRIANGULAR MATRIX        
C       IS STORED COLUMNWISE TOO.                                               
C       THE PROCEDURE GIVES RESULTS IF N IS GREATER THAN 0 AND ALL              
C       CALCULATED RADICANDS ARE POSITIVE.                                      
C                                                                               
C     SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED.                            
C       MFSD                                                                    
C                                                                               
C     METHOD                                                                    
C       SOLUTION IS DONE USING THE FACTORIZATION BY SUBROUTINE MFSD.            
C                                                                               
      DIMENSION A(1)                                                            
      DOUBLE PRECISION DIN,WORK                                                 
C                                                                               
C     FACTORIZE GIVEN MATRIX BY MEANS OF SUBROUTINE MFSD                        
C                                                                               
C     A=TRANSPOSE(T)*T                                                          
C                                                                               
      CALL MFSD(A,N,EPS,IER)                                                    
      IF(IER)9,1,1                                                              
C                                                                               
C     INVERT UPPER TRIANGULAR MATRIX T                                          
C     PREPARE INVERSION-LOOP                                                    
C                                                                               
1     IPIV=N*(N+1)/2                                                            
      IND=IPIV                                                                  
C                                                                               
C     INITIALIZE INVERSION-LOOP                                                 
C                                                                               
      DO 6 I=1,N                                                                
      DIN=1.D0/DBLE(A(IPIV))                                                    
      A(IPIV)=DIN                                                               
      MIN=N                                                                     
      KEND=I-1                                                                  
      LANF=N-KEND                                                               
      IF(KEND)5,5,2                                                             
2     J=IND                                                                     
C                                                                               
C     INITIALIZE ROW-LOOP                                                       
C                                                                               
      DO 4 K=1,KEND                                                             
      WORK=0.D0                                                                 
      MIN=MIN-1                                                                 
      LHOR=IPIV                                                                 
      LVER=J                                                                    
C                                                                               
C     START INNER LOOP                                                          
C                                                                               
      DO 3 L=LANF,MIN                                                           
      LVER=LVER+1                                                               
      LHOR=LHOR+L                                                               
3     WORK=WORK+DBLE(A(LVER)*A(LHOR))                                           
C                                                                               
C     END OF INNER LOOP                                                         
C                                                                               
      A(J)=-WORK*DIN                                                            
4     J=J-MIN                                                                   
C                                                                               
C     END OF ROW-LOOP                                                           
C                                                                               
5     IPIV=IPIV-MIN                                                             
6     IND=IND-1                                                                 
C                                                                               
C     END OF INVERSION LOOP                                                     
C                                                                               
C     CALCULATE INVERSE(A) BY MEANS OF INVERSE(T)                               
C     INVERSE(A)=INVERSE(T)*TRANSPOSE(INVERSE(T))                               
C     INITIALIZE MULTIPLICATION LOOP                                            
C                                                                               
      DO 8 I=1,N                                                                
      IPIV=IPIV+I                                                               
      J=IPIV                                                                    
C                                                                               
C     INITIALIZE ROW-LOOP                                                       
C                                                                               
      DO 8 K=I,N                                                                
      WORK=0.D0                                                                 
      LHOR=J                                                                    
C                                                                               
C     START INNER LOOP                                                          
C                                                                               
      DO 7 L=K,N                                                                
      LVER=LHOR+K-I                                                             
      WORK=WORK+DBLE(A(LHOR)*A(LVER))                                           
7     LHOR=LHOR+L                                                               
C                                                                               
C     END OF INNER LOOP                                                         
C                                                                               
      A(J)=WORK                                                                 
8     J=J+K                                                                     
C                                                                               
C     END OF ROW- AND MULTIPLICATION-LOOP                                       
C                                                                               
9     RETURN                                                                    
      END                                                                       
C----------------------------------------------------------------------
      SUBROUTINE MFSD(A,N,EPS,IER)                                              
C----------------------------------------------------------------------
C                                                                               
C     PURPOSE                                                                   
C       FACTOR A GIVEN SYMMETRIC POSITIVE DEFINITE MATRIX                       
C                                                                               
C     DESCRIPTION OF PARAMETERS                                                 
C       A     - UPPER TRIANGULAR PART OF THE GIVEN SYMMETRIC POSITIVE           
C               DEFINITE N BY N COEFFICIENT MATRIX.  ON RETURN A                
C               CONTAINS THE RESULTANT UPPER TRIANGULAR MATRIX.                 
C       N     - THE NUMBER OF ROW (COLUMNS) IN GIVEN MATRIX.                    
C       EPS   - AN INPUT CONSTANT WHICH IS USED AS RELATIVE TOLERANCE           
C               FOR TEST ON LOSS OF SIGNIFICANCE.                               
C       IER   - RESULTING ERROR PARAMETER CODED AS FOLLOWS:                     
C               IER=0  - NO ERROR                                               
C               IER=-1 - NO RESULT BECAUSE OF WRONG INPUT PARAMETER N OR        
C                        BECAUSE SOME RADICAND IS NONPOSITIVE (MATRIX A         
C                        IS NOT POSITIVE DEFINITE, POSSIBLY DUE TO LOSS         
C                        OF SIGNIFICANCE)                                       
C               IER=K  - WARNING WHICH INDICATES LOSS OF SIGNIFICANCE.          
C                        THE RADICAND FORMED AT FACTORIZATION STEP K+1          
C                        WAS STILL POSITIVE BUT NO LONGER GREATER THAN          
C                        ABS(EPS*A(K+1,K+1)).                                   
C                                                                               
C     REMARKS                                                                   
C       THE UPPER TRIANGULAR PART OF GIVEN MATRIX IS ASSUMED TO BE              
C       STORED COLUMNWISE IN N*(N+1)/2 SUCCESSIVE STORAGE LOCATIONS.  IN        
C       THE SAME STORAGE LOCATIONS THE RESULTING UPPER TRIANGULAR MATRIX        
C       IS STORED COLUMNWISE TOO.                                               
C       THE PROCEDURE GIVES RESULTS IF N IS GREATER THAN 0 AND ALL              
C       CALCULATED RADICANDS ARE POSITIVE.                                      
C       THE PRODUCT OF RETURNED DIAGONAL TERMS IS EQUAL TO THE SQUARE           
C       ROOT OF THE DETERMINANT OF THE GIVEN MATRIX.                            
C                                                                               
C     SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                             
C       NONE                                                                    
C                                                                               
C     METHOD                                                                    
C       SOLUTION IS DONE USING THE SQUARE-ROOT METHOD OF CHOLESKY.              
C       THE GIVEN MATRIX IS REPRESENTED AS THE PRODUCT OF 2 TRIANGULAR          
C       MATRICES, WHERE THE LEFT HAND FACTOR IS THE TRANSPOSE OF THE            
C       RETURNED RIGHT HAND FACTOR.                                             
C                                                                               
      DIMENSION A(1)                                                            
      DOUBLE PRECISION DPIV,DSUM                                                
C                                                                               
C     TEST ON WRONG INPUT PARAMETER N                                           
C                                                                               
      IF(N-1)12,1,1                                                             
1     IER=0                                                                     
C                                                                               
C     INITIALIZE DIAGONAL-LOOP                                                  
C                                                                               
      KPIV=0                                                                    
      DO 11 K=1,N                                                               
      KPIV=KPIV+K                                                               
      IND=KPIV                                                                  
      LEND=K-1                                                                  
C                                                                               
C     CALCULATE TOLERANCE                                                       
C                                                                               
      TOL=ABS(EPS*A(KPIV))                                                      
C                                                                               
C     START FACTORIZATION-LOOP OVER K-TH ROW                                    
C                                                                               
      DO 11 I=K,N                                                               
      DSUM=0.D0                                                                 
      IF(LEND)2,4,2                                                             
C                                                                               
C     START INNER LOOP                                                          
C                                                                               
2     DO 3 L=1,LEND                                                             
      LANF=KPIV-L                                                               
      LIND=IND-L                                                                
3     DSUM=DSUM+DBLE(A(LANF)*A(LIND))                                           
C                                                                               
C     END OF INNER LOOP                                                         
C                                                                               
C     TRANSFORM ELEMENT A(IND)                                                  
C                                                                               
4     DSUM=DBLE(A(IND))-DSUM                                                    
      IF(I-K)10,5,10                                                            
C                                                                               
C     TEST FOR NEGATIVE PIVOT ELEMENT AND FOR LOSS OF SIGNIFICANCE.             
C                                                                               
5     IF(SNGL(DSUM)-TOL)6,6,9                                                   
6     IF(DSUM)12,12,7                                                           
7     IF(IER)8,8,9                                                              
8     IER=K-1                                                                   
C                                                                               
C     COMPUTE PIVOT ELEMENT                                                     
C                                                                               
9     DPIV=DSQRT(DSUM)                                                          
      A(KPIV)=DPIV                                                              
      DPIV=1.D0/DPIV                                                            
      GO TO 11                                                                  
C                                                                               
C     CALCULATE TERMS IN ROW                                                    
C                                                                               
10    A(IND)=DSUM*DPIV                                                          
11    IND=IND+I                                                                 
C                                                                               
C     END OF DIAGONAL-LOOP                                                      
C                                                                               
      RETURN                                                                    
12    IER=-1                                                                    
      RETURN                                                                    
      END                                                                       
