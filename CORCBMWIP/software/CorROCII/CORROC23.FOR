C------------------------------------------------------------------             
      SUBROUTINE SECOND                                                         
C------------------------------------------------------------------             
C                                                                               
C     COMPUTE THE EXPECTED VALUE OF EACH OF THE MINUS SECOND-ORDER              
C     PARTIAL DERIVATIVE EQUATIONS.                                             
C                                                                               
      COMMON/BLK1/VV(11,11),WW(11,11),KATI,KATJ,TCODE,FP,NNI,NNJ,KKI,           
     1  KSNX,KSNY,KKJ,KKIL,KKJL,MTRX,IDENTI(60),IDENTJ(60)                 
      COMMON/BLK2/EMN,EMS,SUMNT(11),SUMST(11),SUMNU(11),SUMSU(11),              
     1            FPFX(11),TPFX(11),FPFY(11),TPFY(11)                           
      COMMON/BLK3/TX(9),UY(9),TEXPX(11),TEXPY(11),                              
     1            TXPTUR(11,11),TXPTUP(11,11),FYXP(11,11),FXYP(11,11)
      common/blk9/FUTR(9,11),FTUR(11,9),ELR(11,11),ELP(11,11),                  
     3            CPELR(10,10),CPELP(10,10),SKP,SKR,RADP,RADR,FUNLIK            
      COMMON/BLK4/FOP(24),SOPNEG(24,24),VCOV(24,24),ESTCOR(24)
      common/blk8/CORR(24,24),XXDUM(496),ITER,CRIT,LSTAT                        
      COMMON/BLK7/AX,BX,AY,BY,R,RHO,T(11),U(11)                                 
      INTEGER TCODE                                                             
      GPI=3.14159265                                                            
C                                                                               
C     WITH RESPECT TO AX                                                        
C                                                                               
      SUM=0.                                                                    
      DO 1210 I=1,KATI                                                          
      DO 1210 J=1,KATJ                                                          
      SUM=SUM+(TEXPX(I)*(FYXP(I,J+1)-FYXP(I,J))-TEXPX(I+1)*                     
     1(FYXP(I+1,J+1)-FYXP(I+1,J)))**2/CPELP(I,J)                                
 1210 CONTINUE                                                                  
      SOPNEG(1,1)=EMS*SUM/(2.*GPI)                                              
C                                                                               
C     WITH RESPECT TO BX                                                        
C                                                                               
      SUM=0.                                                                    
      DO 1215 I=1,KATI                                                          
      DO 1215 J=1,KATJ                                                          
      SUM=SUM+(-T(I)*TEXPX(I)*(FYXP(I,J+1)-FYXP(I,J))+T(I+1)*TEXPX(I+1)*        
     1(FYXP(I+1,J+1)-FYXP(I+1,J)))**2/CPELP(I,J)                                
 1215 CONTINUE                                                                  
      SOPNEG(2,2)=EMS*SUM/(2.*GPI)                                              
C                                                                               
C     WITH RESPECT TO AY                                                        
C                                                                               
      SUM=0.                                                                    
      DO 1220 I=1,KATI                                                          
      DO 1220 J=1,KATJ                                                          
      SUM=SUM+(TEXPY(J)*(FXYP(I+1,J)-FXYP(I,J))-TEXPY(J+1)*                     
     1(FXYP(I+1,J+1)-FXYP(I,J+1)))**2/CPELP(I,J)                                
 1220 CONTINUE                                                                  
      SOPNEG(3,3)=EMS*SUM/(2.*GPI)                                              
C                                                                               
C     WITH RESPECT TO BY                                                        
C                                                                               
      SUM=0.                                                                    
      DO 1225 I=1,KATI                                                          
      DO 1225 J=1,KATJ                                                          
      SUM=SUM+(-U(J)*TEXPY(J)*(FXYP(I+1,J)-FXYP(I,J))+U(J+1)*TEXPY(J+1)*        
     1(FXYP(I+1,J+1)-FXYP(I,J+1)))**2/CPELP(I,J)                                
 1225 CONTINUE                                                                  
      SOPNEG(4,4)=EMS*SUM/(2.*GPI)                                              
C                                                                               
C     WITH RESPECT TO R=R(NOISE)                                                
C                                                                               
      SUM=0.                                                                    
      DO 1230 I=1,KATI                                                          
      DO 1230 J=1,KATJ                                                          
      SUM=SUM+(TXPTUR(I+1,J+1)+TXPTUR(I,J)-TXPTUR(I+1,J)-                       
     1TXPTUR(I,J+1))**2/CPELR(I,J)                                              
 1230 CONTINUE                                                                  
      SOPNEG(5,5)=EMN*SUM/(4.*GPI**2*SKR)                                       
C                                                                               
C     WITH RESPECT TO RHO=R(SIGNAL-PLUS-NOISE)                                  
C                                                                               
      SUM=0.                                                                    
      DO 1235 I=1,KATI                                                          
      DO 1235 J=1,KATJ                                                          
      SUM=SUM+(TXPTUP(I+1,J+1)+TXPTUP(I,J)-TXPTUP(I+1,J)-                       
     1TXPTUP(I,J+1))**2/CPELP(I,J)                                              
 1235 CONTINUE                                                                  
      SOPNEG(6,6)=EMS*SUM/(4.*GPI**2*SKP)                                       
C                                                                               
C     WITH RESPECT TO EACH OF THE T(I)                                          
C                                                                               
      DO 1245 I=1,KKI                                                           
      SUMN=0.                                                                   
      SUMS=0.                                                                   
      DO 1240 J=1,KATJ                                                          
      SUMN=SUMN+(EXP(-.5*T(I+1)**2)*(FUTR(I,J+1)-FUTR(I,J)))**2*(1./            
     1CPELR(I,J)+1./CPELR(I+1,J))                                               
      SUMS=SUMS+(TEXPX(I+1)*(FYXP(I+1,J+1)-FYXP(I+1,J)))**2*(1./                
     1CPELP(I,J)+1./CPELP(I+1,J))                                               
 1240 CONTINUE                                                                  
      SOPNEG(I+6,I+6)=(EMN*SUMN+EMS*BX**2*SUMS)/(2.*GPI)                        
 1245 CONTINUE                                                                  
C                                                                               
C     WITH RESPECT TO EACH OF THE U(J)                                          
C                                                                               
      DO 1255 J=1,KKJ                                                           
      SUMN=0.                                                                   
      SUMS=0.                                                                   
      DO 1250 I=1,KATI                                                          
      SUMN=SUMN+(EXP(-.5*U(J+1)**2)*(FTUR(I+1,J)-FTUR(I,J)))**2*(1./            
     1CPELR(I,J)+1./CPELR(I,J+1))                                               
      SUMS=SUMS+(TEXPY(J+1)*(FXYP(I+1,J+1)-FXYP(I,J+1)))**2*(1./                
     1CPELP(I,J)+1./CPELP(I,J+1))                                               
 1250 CONTINUE                                                                  
      SOPNEG(J+KKI+6,J+KKI+6)=(EMN*SUMN+EMS*BY**2*SUMS)/(2.*GPI)                
 1255 CONTINUE                                                                  
      return 
      end
c--------------------------------------------------
      subroutine secnd2
c--------------------------------------------------
      COMMON/BLK1/VV(11,11),WW(11,11),KATI,KATJ,TCODE,FP,NNI,NNJ,KKI,           
     1  KSNX,KSNY,KKJ,KKIL,KKJL,MTRX,IDENTI(60),IDENTJ(60)                 
      COMMON/BLK2/EMN,EMS,SUMNT(11),SUMST(11),SUMNU(11),SUMSU(11),              
     1            FPFX(11),TPFX(11),FPFY(11),TPFY(11)                           
      COMMON/BLK3/TX(9),UY(9),TEXPX(11),TEXPY(11),                              
     1            TXPTUR(11,11),TXPTUP(11,11),FYXP(11,11),FXYP(11,11)
      common/blk9/FUTR(9,11),FTUR(11,9),ELR(11,11),ELP(11,11),                  
     3            CPELR(10,10),CPELP(10,10),SKP,SKR,RADP,RADR,FUNLIK            
      COMMON/BLK4/FOP(24),SOPNEG(24,24),VCOV(24,24),ESTCOR(24)
      common/blk8/CORR(24,24),XXDUM(496),ITER,CRIT,LSTAT                        
      COMMON/BLK7/AX,BX,AY,BY,R,RHO,T(11),U(11)                                 
      INTEGER TCODE                                                             
      GPI=3.14159265                                                            
C                                                                               
C     MIXED WITH RESPECT TO AX AND BX                                           
C                                                                               
      SUM=0.                                                                    
      DO 1260 I=1,KATI                                                          
      DO 1260 J=1,KATJ                                                          
      SUM=SUM+(TEXPX(I)*(FYXP(I,J+1)-FYXP(I,J))-TEXPX(I+1)*                     
     1(FYXP(I+1,J+1)-FYXP(I+1,J)))*(-T(I)*TEXPX(I)*(FYXP(I,J+1)-                
     2FYXP(I,J))+T(I+1)*TEXPX(I+1)*(FYXP(I+1,J+1)-FYXP(I+1,J)))/                
     3CPELP(I,J)                                                                
 1260 CONTINUE                                                                  
      SOPNEG(1,2)=EMS*SUM/(2.*GPI)                                              
      SOPNEG(2,1)=SOPNEG(1,2)                                                   
C                                                                               
C     MIXED WITH RESPECT TO AX AND AY                                           
C                                                                               
      SUM=0.                                                                    
      DO 1265 I=1,KATI                                                          
      DO 1265 J=1,KATJ                                                          
      SUM=SUM+(TEXPX(I)*(FYXP(I,J+1)-FYXP(I,J))-TEXPX(I+1)*                     
     1(FYXP(I+1,J+1)-FYXP(I+1,J)))*(TEXPY(J)*(FXYP(I+1,J)-FXYP(I,J))-           
     2TEXPY(J+1)*(FXYP(I+1,J+1)-FXYP(I,J+1)))/CPELP(I,J)                        
 1265 CONTINUE                                                                  
      SOPNEG(1,3)=EMS*SUM/(2.*GPI)                                              
      SOPNEG(3,1)=SOPNEG(1,3)                                                   
C                                                                               
C     MIXED WITH RESPECT TO AX AND BY                                           
C                                                                               
      SUM=0.                                                                    
      DO 1270 I=1,KATI                                                          
      DO 1270 J=1,KATJ                                                          
      SUM=SUM+(TEXPX(I)*(FYXP(I,J+1)-FYXP(I,J))-TEXPX(I+1)*                     
     1(FYXP(I+1,J+1)-FYXP(I+1,J)))*(-U(J)*TEXPY(J)*(FXYP(I+1,J)-                
     2FXYP(I,J))+U(J+1)*TEXPY(J+1)*(FXYP(I+1,J+1)-FXYP(I,J+1)))/                
     3CPELP(I,J)                                                                
 1270 CONTINUE                                                                  
      SOPNEG(1,4)=EMS*SUM/(2.*GPI)                                              
      SOPNEG(4,1)=SOPNEG(1,4)                                                   
C                                                                               
C     MIXED WITH RESPECT TO AX AND RHO=R(SIGNAL-PLUS-NOISE)                     
C                                                                               
      SUM=0.                                                                    
      DO 1275 I=1,KATI                                                          
      DO 1275 J=1,KATJ                                                          
      SUM=SUM+(TEXPX(I)*(FYXP(I,J+1)-FYXP(I,J))-TEXPX(I+1)*                     
     1(FYXP(I+1,J+1)-FYXP(I+1,J)))*(TXPTUP(I+1,J+1)+TXPTUP(I,J)-                
     2TXPTUP(I+1,J)-TXPTUP(I,J+1))/CPELP(I,J)                                   
 1275 CONTINUE                                                                  
      SOPNEG(1,6)=EMS*SUM/SQRT(8.*GPI**3*SKP)                                   
      SOPNEG(6,1)=SOPNEG(1,6)                                                   
C                                                                               
C     MIXED WITH RESPECT TO AX AND EACH OF THE T(I)                             
C                                                                               
      DO 1285 I=1,KKI                                                           
      SUM=0.                                                                    
      DO 1280 J=1,KATJ                                                          
      SUM=SUM+TEXPX(I+1)*(FYXP(I+1,J+1)-FYXP(I+1,J))*((TEXPX(I)*                
     1(FYXP(I,J+1)-FYXP(I,J))-TEXPX(I+1)*(FYXP(I+1,J+1)-FYXP(I+1,J)))/          
     2CPELP(I,J)-(TEXPX(I+1)*(FYXP(I+1,J+1)-FYXP(I+1,J))-TEXPX(I+2)*            
     3(FYXP(I+2,J+1)-FYXP(I+2,J)))/CPELP(I+1,J))                                
 1280 CONTINUE                                                                  
      SOPNEG(1,I+6)=EMS*BX*SUM/(2.*GPI)                                         
      SOPNEG(I+6,1)=SOPNEG(1,I+6)                                               
 1285 CONTINUE                                                                  
C                                                                               
C     MIXED WITH RESPECT TO AX AND EACH OF THE U(J)                             
C                                                                               
      DO 1295 J=1,KKJ                                                           
      SUM=0.                                                                    
      DO 1290 I=1,KATI                                                          
      SUM=SUM+TEXPY(J+1)*(FXYP(I+1,J+1)-FXYP(I,J+1))*((TEXPX(I)*                
     1(FYXP(I,J+1)-FYXP(I,J))-TEXPX(I+1)*(FYXP(I+1,J+1)-FYXP(I+1,J)))/          
     2CPELP(I,J)-(TEXPX(I)*(FYXP(I,J+2)-FYXP(I,J+1))-TEXPX(I+1)*                
     3(FYXP(I+1,J+2)-FYXP(I+1,J+1)))/CPELP(I,J+1))                              
 1290 CONTINUE                                                                  
      SOPNEG(1,J+KKI+6)=EMS*BY*SUM/(2.*GPI)                                     
      SOPNEG(J+KKI+6,1)=SOPNEG(1,J+KKI+6)                                       
 1295 CONTINUE                                                                  
      return
      end
c----------------------------------
      subroutine secnd3
c----------------------------------
      COMMON/BLK1/VV(11,11),WW(11,11),KATI,KATJ,TCODE,FP,NNI,NNJ,KKI,           
     1  KSNX,KSNY,KKJ,KKIL,KKJL,MTRX,IDENTI(60),IDENTJ(60)                 
      COMMON/BLK2/EMN,EMS,SUMNT(11),SUMST(11),SUMNU(11),SUMSU(11),              
     1            FPFX(11),TPFX(11),FPFY(11),TPFY(11)                           
      COMMON/BLK3/TX(9),UY(9),TEXPX(11),TEXPY(11),                              
     1            TXPTUR(11,11),TXPTUP(11,11),FYXP(11,11),FXYP(11,11)
      common/blk9/FUTR(9,11),FTUR(11,9),ELR(11,11),ELP(11,11),                  
     3            CPELR(10,10),CPELP(10,10),SKP,SKR,RADP,RADR,FUNLIK            
      COMMON/BLK4/FOP(24),SOPNEG(24,24),VCOV(24,24),ESTCOR(24)
      common/blk8/CORR(24,24),XXDUM(496),ITER,CRIT,LSTAT                        
      COMMON/BLK7/AX,BX,AY,BY,R,RHO,T(11),U(11)                                 
      INTEGER TCODE                                                             
      GPI=3.14159265                                                            
C                                                                               
C     MIXED WITH RESPECT TO BX AND AY                                           
C                                                                               
      SUM=0.                                                                    
      DO 1300 I=1,KATI                                                          
      DO 1300 J=1,KATJ                                                          
      SUM=SUM+(-T(I)*TEXPX(I)*(FYXP(I,J+1)-FYXP(I,J))+T(I+1)*TEXPX(I+1)*        
     1(FYXP(I+1,J+1)-FYXP(I+1,J)))*(TEXPY(J)*(FXYP(I+1,J)-FXYP(I,J))-           
     2TEXPY(J+1)*(FXYP(I+1,J+1)-FXYP(I,J+1)))/CPELP(I,J)                        
 1300 CONTINUE                                                                  
      SOPNEG(2,3)=EMS*SUM/(2.*GPI)                                              
      SOPNEG(3,2)=SOPNEG(2,3)                                                   
C                                                                               
C     MIXED WITH RESPECT TO BX AND BY                                           
C                                                                               
      SUM=0.                                                                    
      DO 1305 I=1,KATI                                                          
      DO 1305 J=1,KATJ                                                          
      SUM=SUM+(-T(I)*TEXPX(I)*(FYXP(I,J+1)-FYXP(I,J))+T(I+1)*TEXPX(I+1)*        
     1(FYXP(I+1,J+1)-FYXP(I+1,J)))*(-U(J)*TEXPY(J)*(FXYP(I+1,J)-                
     2FXYP(I,J))+U(J+1)*TEXPY(J+1)*(FXYP(I+1,J+1)-FXYP(I,J+1)))/                
     3CPELP(I,J)                                                                
 1305 CONTINUE                                                                  
      SOPNEG(2,4)=EMS*SUM/(2.*GPI)                                              
      SOPNEG(4,2)=SOPNEG(2,4)                                                   
C                                                                               
C     MIXED WITH RESPECT TO BX AND RHO=R(SIGNAL-PLUS-NOISE)                     
C                                                                               
      SUM=0.                                                                    
      DO 1310 I=1,KATI                                                          
      DO 1310 J=1,KATJ                                                          
      SUM=SUM+(-T(I)*TEXPX(I)*(FYXP(I,J+1)-FYXP(I,J))+T(I+1)*TEXPX(I+1)*        
     1(FYXP(I+1,J+1)-FYXP(I+1,J)))*(TXPTUP(I+1,J+1)+TXPTUP(I,J)-                
     2TXPTUP(I+1,J)-TXPTUP(I,J+1))/CPELP(I,J)                                   
 1310 CONTINUE                                                                  
      SOPNEG(2,6)=EMS*SUM/SQRT(8.*GPI**3*SKP)                                   
      SOPNEG(6,2)=SOPNEG(2,6)                                                   
C                                                                               
C     MIXED WITH RESPECT TO BX AND EACH OF THE T(I)                             
C                                                                               
      DO 1320 I=1,KKI                                                           
      SUM=0.                                                                    
      DO 1315 J=1,KATJ                                                          
      SUM=SUM+TEXPX(I+1)*(FYXP(I+1,J+1)-FYXP(I+1,J))*((-T(I)*TEXPX(I)*          
     1(FYXP(I,J+1)-FYXP(I,J))+T(I+1)*TEXPX(I+1)*(FYXP(I+1,J+1)-                 
     2FYXP(I+1,J)))/CPELP(I,J)-(-T(I+1)*TEXPX(I+1)*(FYXP(I+1,J+1)-              
     3FYXP(I+1,J))+T(I+2)*TEXPX(I+2)*(FYXP(I+2,J+1)-FYXP(I+2,J)))/              
     4CPELP(I+1,J))                                                             
 1315 CONTINUE                                                                  
      SOPNEG(2,I+6)=EMS*BX*SUM/(2.*GPI)                                         
      SOPNEG(I+6,2)=SOPNEG(2,I+6)                                               
 1320 CONTINUE                                                                  
C                                                                               
C     MIXED WITH RESPECT TO BX AND EACH OF THE U(J)                             
C                                                                               
      DO 1330 J=1,KKJ                                                           
      SUM=0.                                                                    
      DO 1325 I=1,KATI                                                          
      SUM=SUM+TEXPY(J+1)*(FXYP(I+1,J+1)-FXYP(I,J+1))*((-T(I)*TEXPX(I)*          
     1(FYXP(I,J+1)-FYXP(I,J))+T(I+1)*TEXPX(I+1)*(FYXP(I+1,J+1)-                 
     2FYXP(I+1,J)))/CPELP(I,J)-(-T(I)*TEXPX(I)*(FYXP(I,J+2)-FYXP(I,J+1))        
     3+T(I+1)*TEXPX(I+1)*(FYXP(I+1,J+2)-FYXP(I+1,J+1)))/CPELP(I,J+1))           
 1325 CONTINUE                                                                  
      SOPNEG(2,J+KKI+6)=EMS*BY*SUM/(2.*GPI)                                     
      SOPNEG(J+KKI+6,2)=SOPNEG(2,J+KKI+6)                                       
 1330 CONTINUE                                                                  
      return
      end
c-----------------------------------
      subroutine secnd4
c-----------------------------------
      COMMON/BLK1/VV(11,11),WW(11,11),KATI,KATJ,TCODE,FP,NNI,NNJ,KKI,           
     1  KSNX,KSNY,KKJ,KKIL,KKJL,MTRX,IDENTI(60),IDENTJ(60)                 
      COMMON/BLK2/EMN,EMS,SUMNT(11),SUMST(11),SUMNU(11),SUMSU(11),              
     1            FPFX(11),TPFX(11),FPFY(11),TPFY(11)                           
      COMMON/BLK3/TX(9),UY(9),TEXPX(11),TEXPY(11),                              
     1            TXPTUR(11,11),TXPTUP(11,11),FYXP(11,11),FXYP(11,11)
      common/blk9/FUTR(9,11),FTUR(11,9),ELR(11,11),ELP(11,11),                  
     3            CPELR(10,10),CPELP(10,10),SKP,SKR,RADP,RADR,FUNLIK            
      COMMON/BLK4/FOP(24),SOPNEG(24,24),VCOV(24,24),ESTCOR(24)
      common/blk8/CORR(24,24),XXDUM(496),ITER,CRIT,LSTAT                        
      COMMON/BLK7/AX,BX,AY,BY,R,RHO,T(11),U(11)                                 
      INTEGER TCODE                                                             
      GPI=3.14159265                                                            
C                                                                               
C     MIXED WITH RESPECT TO AY AND BY                                           
C                                                                               
      SUM=0.                                                                    
      DO 1335 I=1,KATI                                                          
      DO 1335 J=1,KATJ                                                          
      SUM=SUM+(TEXPY(J)*(FXYP(I+1,J)-FXYP(I,J))-TEXPY(J+1)*                     
     1(FXYP(I+1,J+1)-FXYP(I,J+1)))*(-U(J)*TEXPY(J)*(FXYP(I+1,J)-                
     2FXYP(I,J))+U(J+1)*TEXPY(J+1)*(FXYP(I+1,J+1)-FXYP(I,J+1)))/                
     3CPELP(I,J)                                                                
 1335 CONTINUE                                                                  
      SOPNEG(3,4)=EMS*SUM/(2.*GPI)                                              
      SOPNEG(4,3)=SOPNEG(3,4)                                                   
C                                                                               
C     MIXED WITH RESPECT TO AY AND RHO=R(SIGNAL-PLUS-NOISE)                     
C                                                                               
      SUM=0.                                                                    
      DO 1340 I=1,KATI                                                          
      DO 1340 J=1,KATJ                                                          
      SUM=SUM+(TEXPY(J)*(FXYP(I+1,J)-FXYP(I,J))-TEXPY(J+1)*                     
     1(FXYP(I+1,J+1)-FXYP(I,J+1)))*(TXPTUP(I+1,J+1)+TXPTUP(I,J)-                
     2TXPTUP(I+1,J)-TXPTUP(I,J+1))/CPELP(I,J)                                   
 1340 CONTINUE                                                                  
      SOPNEG(3,6)=EMS*SUM/SQRT(8.*GPI**3*SKP)                                   
      SOPNEG(6,3)=SOPNEG(3,6)                                                   
C                                                                               
C     MIXED WITH RESPECT TO AY AND EACH OF THE T(I)                             
C                                                                               
      DO 1350 I=1,KKI                                                           
      SUM=0.                                                                    
      DO 1345 J=1,KATJ                                                          
      SUM=SUM+TEXPX(I+1)*(FYXP(I+1,J+1)-FYXP(I+1,J))*((TEXPY(J)*                
     1(FXYP(I+1,J)-FXYP(I,J))-TEXPY(J+1)*(FXYP(I+1,J+1)-FXYP(I,J+1)))/          
     2CPELP(I,J)-(TEXPY(J)*(FXYP(I+2,J)-FXYP(I+1,J))-TEXPY(J+1)*                
     3(FXYP(I+2,J+1)-FXYP(I+1,J+1)))/CPELP(I+1,J))                              
 1345 CONTINUE                                                                  
      SOPNEG(3,I+6)=EMS*BX*SUM/(2.*GPI)                                         
      SOPNEG(I+6,3)=SOPNEG(3,I+6)                                               
 1350 CONTINUE                                                                  
C                                                                               
C     MIXED WITH RESPECT TO AY AND EACH OF THE U(J)                             
C                                                                               
      DO 1360 J=1,KKJ                                                           
      SUM=0.                                                                    
      DO 1355 I=1,KATI                                                          
      SUM=SUM+TEXPY(J+1)*(FXYP(I+1,J+1)-FXYP(I,J+1))*((TEXPY(J)*                
     1(FXYP(I+1,J)-FXYP(I,J))-TEXPY(J+1)*(FXYP(I+1,J+1)-FXYP(I,J+1)))/          
     2CPELP(I,J)-(TEXPY(J+1)*(FXYP(I+1,J+1)-FXYP(I,J+1))-TEXPY(J+2)*            
     3(FXYP(I+1,J+2)-FXYP(I,J+2)))/CPELP(I,J+1))                                
 1355 CONTINUE                                                                  
      SOPNEG(3,J+KKI+6)=EMS*BY*SUM/(2.*GPI)                                     
      SOPNEG(J+KKI+6,3)=SOPNEG(3,J+KKI+6)                                       
 1360 CONTINUE                                                                  
      return
      end
c----------------------------------------
      subroutine secnd5
c----------------------------------------
      COMMON/BLK1/VV(11,11),WW(11,11),KATI,KATJ,TCODE,FP,NNI,NNJ,KKI,           
     1  KSNX,KSNY,KKJ,KKIL,KKJL,MTRX,IDENTI(60),IDENTJ(60)                 
      COMMON/BLK2/EMN,EMS,SUMNT(11),SUMST(11),SUMNU(11),SUMSU(11),              
     1            FPFX(11),TPFX(11),FPFY(11),TPFY(11)                           
      COMMON/BLK3/TX(9),UY(9),TEXPX(11),TEXPY(11),                              
     1            TXPTUR(11,11),TXPTUP(11,11),FYXP(11,11),FXYP(11,11)
      common/blk9/FUTR(9,11),FTUR(11,9),ELR(11,11),ELP(11,11),                  
     3            CPELR(10,10),CPELP(10,10),SKP,SKR,RADP,RADR,FUNLIK            
      COMMON/BLK4/FOP(24),SOPNEG(24,24),VCOV(24,24),ESTCOR(24)
      common/blk8/CORR(24,24),XXDUM(496),ITER,CRIT,LSTAT                        
      COMMON/BLK7/AX,BX,AY,BY,R,RHO,T(11),U(11)                                 
      INTEGER TCODE                                                             
      GPI=3.14159265                                                            
C                                                                               
C     MIXED WITH RESPECT TO BY AND RHO=R(SIGNAL-PLUS-NOISE)                     
C                                                                               
      SUM=0.                                                                    
      DO 1365 I=1,KATI                                                          
      DO 1365 J=1,KATJ                                                          
      SUM=SUM+(-U(J)*TEXPY(J)*(FXYP(I+1,J)-FXYP(I,J))+U(J+1)*TEXPY(J+1)*        
     1(FXYP(I+1,J+1)-FXYP(I,J+1)))*(TXPTUP(I+1,J+1)+TXPTUP(I,J)-                
     2TXPTUP(I+1,J)-TXPTUP(I,J+1))/CPELP(I,J)                                   
 1365 CONTINUE                                                                  
      SOPNEG(4,6)=EMS*SUM/SQRT(8.*GPI**3*SKP)                                   
      SOPNEG(6,4)=SOPNEG(4,6)                                                   
C                                                                               
C     MIXED WITH RESPECT TO BY AND EACH OF THE T(I)                             
C                                                                               
      DO 1375 I=1,KKI                                                           
      SUM=0.                                                                    
      DO 1370 J=1,KATJ                                                          
      SUM=SUM+TEXPX(I+1)*(FYXP(I+1,J+1)-FYXP(I+1,J))*((-U(J)*TEXPY(J)*          
     1(FXYP(I+1,J)-FXYP(I,J))+U(J+1)*TEXPY(J+1)*(FXYP(I+1,J+1)-                 
     2FXYP(I,J+1)))/CPELP(I,J)-(-U(J)*TEXPY(J)*(FXYP(I+2,J)-FXYP(I+1,J))        
     3+U(J+1)*TEXPY(J+1)*(FXYP(I+2,J+1)-FXYP(I+1,J+1)))/CPELP(I+1,J))           
 1370 CONTINUE                                                                  
      SOPNEG(4,I+6)=EMS*BX*SUM/(2.*GPI)                                         
      SOPNEG(I+6,4)=SOPNEG(4,I+6)                                               
 1375 CONTINUE                                                                  
C                                                                               
C     MIXED WITH RESPECT TO BY AND EACH OF THE U(J)                             
C                                                                               
      DO 1385 J=1,KKJ                                                           
      SUM=0.                                                                    
      DO 1380 I=1,KATI                                                          
      SUM=SUM+TEXPY(J+1)*(FXYP(I+1,J+1)-FXYP(I,J+1))*((-U(J)*TEXPY(J)*          
     1(FXYP(I+1,J)-FXYP(I,J))+U(J+1)*TEXPY(J+1)*(FXYP(I+1,J+1)-                 
     2FXYP(I,J+1)))/CPELP(I,J)-(-U(J+1)*TEXPY(J+1)*(FXYP(I+1,J+1)-              
     3FXYP(I,J+1))+U(J+2)*TEXPY(J+2)*(FXYP(I+1,J+2)-FXYP(I,J+2)))/              
     4CPELP(I,J+1))                                                             
 1380 CONTINUE                                                                  
      SOPNEG(4,J+KKI+6)=EMS*BY*SUM/(2.*GPI)                                     
      SOPNEG(J+KKI+6,4)=SOPNEG(4,J+KKI+6)                                       
 1385 CONTINUE                                                                  
C                                                                               
C     MIXED WITH RESPECT TO R=R(NOISE) AND EACH OF THE T(I)                     
C                                                                               
      DO 1395 I=1,KKI                                                           
      SUM=0.                                                                    
      DO 1390 J=1,KATJ                                                          
      SUM=SUM+EXP(-.5*T(I+1)**2)*(FUTR(I,J+1)-FUTR(I,J))*                       
     1((TXPTUR(I+1,J+1)+TXPTUR(I,J)-TXPTUR(I+1,J)-TXPTUR(I,J+1))/               
     2CPELR(I,J)-(TXPTUR(I+2,J+1)+TXPTUR(I+1,J)-TXPTUR(I+2,J)-                  
     3TXPTUR(I+1,J+1))/CPELR(I+1,J))                                            
 1390 CONTINUE                                                                  
      SOPNEG(5,I+6)=EMN*SUM/SQRT(8.*GPI**3*SKR)                                 
      SOPNEG(I+6,5)=SOPNEG(5,I+6)                                               
 1395 CONTINUE                                                                  
C                                                                               
C     MIXED WITH RESPECT TO R=R(NOISE) AND EACH OF THE U(J)                     
C                                                                               
      DO 1405 J=1,KKJ                                                           
      SUM=0.                                                                    
      DO 1400 I=1,KATI                                                          
      SUM=SUM+EXP(-.5*U(J+1)**2)*(FTUR(I+1,J)-FTUR(I,J))*                       
     1((TXPTUR(I+1,J+1)+TXPTUR(I,J)-TXPTUR(I+1,J)-TXPTUR(I,J+1))/               
     2CPELR(I,J)-(TXPTUR(I+1,J+2)+TXPTUR(I,J+1)-TXPTUR(I+1,J+1)-                
     3TXPTUR(I,J+2))/CPELR(I,J+1))                                              
 1400 CONTINUE                                                                  
      SOPNEG(5,J+KKI+6)=EMN*SUM/SQRT(8.*GPI**3*SKR)                             
      SOPNEG(J+KKI+6,5)=SOPNEG(5,J+KKI+6)                                       
 1405 CONTINUE                                                                  
      return
      end
c-------------------------------
      subroutine secnd6
c-------------------------------
      COMMON/BLK1/VV(11,11),WW(11,11),KATI,KATJ,TCODE,FP,NNI,NNJ,KKI,           
     1  KSNX,KSNY,KKJ,KKIL,KKJL,MTRX,IDENTI(60),IDENTJ(60)                 
      COMMON/BLK2/EMN,EMS,SUMNT(11),SUMST(11),SUMNU(11),SUMSU(11),              
     1            FPFX(11),TPFX(11),FPFY(11),TPFY(11)                           
      COMMON/BLK3/TX(9),UY(9),TEXPX(11),TEXPY(11),                              
     1            TXPTUR(11,11),TXPTUP(11,11),FYXP(11,11),FXYP(11,11)
      common/blk9/FUTR(9,11),FTUR(11,9),ELR(11,11),ELP(11,11),                  
     3            CPELR(10,10),CPELP(10,10),SKP,SKR,RADP,RADR,FUNLIK            
      COMMON/BLK4/FOP(24),SOPNEG(24,24),VCOV(24,24),ESTCOR(24)
      common/blk8/CORR(24,24),XXDUM(496),ITER,CRIT,LSTAT                        
      COMMON/BLK7/AX,BX,AY,BY,R,RHO,T(11),U(11)                                 
      INTEGER TCODE
      GPI=3.14159265
C
      DO 1415 I=1,KKI                                                           
      SUM=0.                                                                    
      DO 1410 J=1,KATJ                                                          
      SUM=SUM+TEXPX(I+1)*(FYXP(I+1,J+1)-FYXP(I+1,J))*((TXPTUP(I+1,J+1)+         
     1TXPTUP(I,J)-TXPTUP(I+1,J)-TXPTUP(I,J+1))/CPELP(I,J)-                      
     2(TXPTUP(I+2,J+1)+TXPTUP(I+1,J)-TXPTUP(I+2,J)-TXPTUP(I+1,J+1))/            
     3CPELP(I+1,J))                                                             
 1410 CONTINUE                                                                  
      SOPNEG(6,I+6)=EMS*BX*SUM/SQRT(8.*GPI**3*SKP)                              
      SOPNEG(I+6,6)=SOPNEG(6,I+6)                                               
 1415 CONTINUE                                                                  
C                                                                               
C     MIXED WITH RESPECT TO RHO=R(SIGNAL-PLUS-NOISE) AND EACH OF THE U(J)       
C                                                                               
      DO 1425 J=1,KKJ                                                           
      SUM=0.                                                                    
      DO 1420 I=1,KATI                                                          
      SUM=SUM+TEXPY(J+1)*(FXYP(I+1,J+1)-FXYP(I,J+1))*((TXPTUP(I+1,J+1)+         
     1TXPTUP(I,J)-TXPTUP(I+1,J)-TXPTUP(I,J+1))/CPELP(I,J)-                      
     2(TXPTUP(I+1,J+2)+TXPTUP(I,J+1)-TXPTUP(I+1,J+1)-TXPTUP(I,J+2))/            
     3CPELP(I,J+1))                                                             
 1420 CONTINUE                                                                  
      SOPNEG(6,J+KKI+6)=EMS*BY*SUM/SQRT(8.*GPI**3*SKP)                          
      SOPNEG(J+KKI+6,6)=SOPNEG(6,J+KKI+6)                                       
 1425 CONTINUE                                                                  
C                                                                               
C     MIXED WITH RESPECT TO T(I) AND T(I+1) FOR I=1,KKI-1                       
C                                                                               
      DO 1435 I=1,KKIL                                                          
      SUMN=0.                                                                   
      SUMS=0.                                                                   
      DO 1430 J=1,KATJ                                                          
      SUMN=SUMN+(-EXP(-.5*T(I+1)**2)*(FUTR(I,J+1)-FUTR(I,J)))*(EXP(-.5*         
     1T(I+2)**2)*(FUTR(I+1,J+1)-FUTR(I+1,J)))/CPELR(I+1,J)                      
      SUMS=SUMS+(-TEXPX(I+1)*(FYXP(I+1,J+1)-FYXP(I+1,J)))*(TEXPX(I+2)*          
     1(FYXP(I+2,J+1)-FYXP(I+2,J)))/CPELP(I+1,J)                                 
 1430 CONTINUE                                                                  
      SOPNEG(I+6,I+7)=(EMN*SUMN+EMS*BX**2*SUMS)/(2.*GPI)                        
      SOPNEG(I+7,I+6)=SOPNEG(I+6,I+7)                                           
 1435 CONTINUE                                                                  
C                                                                               
C     MIXED WITH RESPECT TO EACH OF THE T(I) AND EACH OF THE U(J)               
C                                                                               
      DO 1440 I=1,KKI                                                           
      DO 1440 J=1,KKJ                                                           
      SUMN=EXP(-.5*U(J+1)**2)*(FTUR(I+1,J)-FTUR(I,J))*(EXP(-.5*                 
     1T(I+1)**2)*(FUTR(I,J+1)-FUTR(I,J))/CPELR(I,J)-EXP(-.5*T(I+1)**2)*         
     2(FUTR(I,J+2)-FUTR(I,J+1))/CPELR(I,J+1))-EXP(-.5*U(J+1)**2)*               
     3(FTUR(I+2,J)-FTUR(I+1,J))*(EXP(-.5*T(I+1)**2)*(FUTR(I,J+1)-               
     4FUTR(I,J))/CPELR(I+1,J)-EXP(-.5*T(I+1)**2)*(FUTR(I,J+2)-                  
     5FUTR(I,J+1))/CPELR(I+1,J+1))                                              
      SUMS=BY*TEXPY(J+1)*(FXYP(I+1,J+1)-FXYP(I,J+1))*(BX*TEXPX(I+1)*            
     1(FYXP(I+1,J+1)-FYXP(I+1,J))/CPELP(I,J)-BX*TEXPX(I+1)*                     
     2(FYXP(I+1,J+2)-FYXP(I+1,J+1))/CPELP(I,J+1))-BY*TEXPY(J+1)*                
     3(FXYP(I+2,J+1)-FXYP(I+1,J+1))*(BX*TEXPX(I+1)*(FYXP(I+1,J+1)-              
     4FYXP(I+1,J))/CPELP(I+1,J)-BX*TEXPX(I+1)*(FYXP(I+1,J+2)-                   
     5FYXP(I+1,J+1))/CPELP(I+1,J+1))                                            
      SOPNEG(I+6,J+KKI+6)=(EMN*SUMN+EMS*SUMS)/(2.*GPI)                          
      SOPNEG(J+KKI+6,I+6)=SOPNEG(I+6,J+KKI+6)                                   
 1440 CONTINUE                                                                  
C                                                                               
C     MIXED WITH RESPECT TO U(J) AND U(J+1)FOR J=1,KKJ-1                        
C                                                                               
      DO 1450 J=1,KKJL                                                          
      SUMN=0.                                                                   
      SUMS=0.                                                                   
      DO 1445 I=1,KATI                                                          
      SUMN=SUMN+(-EXP(-.5*U(J+1)**2)*(FTUR(I+1,J)-FTUR(I,J)))*(EXP(-.5*         
     1U(J+2)**2)*(FTUR(I+1,J+1)-FTUR(I,J+1)))/CPELR(I,J+1)                      
      SUMS=SUMS+(-TEXPY(J+1)*(FXYP(I+1,J+1)-FXYP(I,J+1)))*(TEXPY(J+2)*          
     1(FXYP(I+1,J+2)-FXYP(I,J+2)))/CPELP(I,J+1)                                 
 1445 CONTINUE                                                                  
      SOPNEG(J+KKI+6,J+KKI+7)=(EMN*SUMN+EMS*BY**2*SUMS)/(2.*GPI)                
      SOPNEG(J+KKI+7,J+KKI+6)=SOPNEG(J+KKI+6,J+KKI+7)                           
 1450 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
C                                                                               
C     THE EXPECTED VALUES OF ALL OTHER MINUS SECOND-ORDER PARTIAL               
C     DERIVATIVES ARE ZERO.  THESE ARE THE MIXED PARTIALS WITH                  
C     RESPECT TO R AND EACH OF AX,BX,AY,BY AND RHO; THE MIXED PARTIALS          
C     WITH RESPECT TO T(I) AND T(I') WHERE ABS(I-I')>1; THE MIXED PARTIALS      
C     WITH RESPECT TO U(J) AND U(J') WHERE ABS(J-J')>1; AND THE MIXED           
C     PARTIALS WITH RESPECT TO R AND RHO.  THESE HAVE ALREADY BEEN SET          
C     TO ZERO WHEN ALL MATRIX ELEMENTS WERE INITIALIZED TO ZERO IN THE          
C     1125 DO-LOOP.                                                             
C                                                                               
C---------------------------------------------------------------------          
      SUBROUTINE ITERAT                                                         
C---------------------------------------------------------------------          
C                                                                               
C     NOW THAT THE VALUES OF ALL MATRIX ELEMENTS HAVE BEEN COMPUTED,            
C     THE NEXT STEP IS TO INVERT THE MATRIX.  WHEN THE ITERATIVE                
C     PROCEDURE HAS CONVERGED, THE INVERTED MATRIX OF MINUS SECOND-ORDER        
C     PARTIAL DERIVATIVES WILL BE EQUIVALENT TO THE VARIANCE-COVARIANCE         
C     MATRIX.                                                                   
C                                                                               
      COMMON/BLK1/VV(11,11),WW(11,11),KATI,KATJ,TCODE,FP,NNI,NNJ,KKI,           
     1  KSNX,KSNY,KKJ,KKIL,KKJL,MTRX,IDENTI(60),IDENTJ(60)                 
      COMMON/BLK2/EMN,EMS,SUMNT(11),SUMST(11),SUMNU(11),SUMSU(11),              
     1            FPFX(11),TPFX(11),FPFY(11),TPFY(11)                           
      COMMON/BLK3/TX(9),UY(9),TEXPX(11),TEXPY(11),                              
     1            TXPTUR(11,11),TXPTUP(11,11),FYXP(11,11),FXYP(11,11)
      common/blk9/FUTR(9,11),FTUR(11,9),ELR(11,11),ELP(11,11),                  
     3            CPELR(10,10),CPELP(10,10),SKP,SKR,RADP,RADR,FUNLIK            
      COMMON/BLK4/FOP(24),SOPNEG(24,24),VCOV(24,24),ESTCOR(24)
      common/blk8/CORR(24,24),XXDUM(496),ITER,CRIT,LSTAT                        
      COMMON/BLK7/AX,BX,AY,BY,R,RHO,T(11),U(11)                                 
      INTEGER TCODE                                                             
      IFLAG=0                                                                   
      CALL MSTR(SOPNEG,XXDUM,24,0,1)                                            
      CALL SINV(XXDUM,MTRX,.001,IER)                                            
      IF(IER.GE.0) GO TO 1454                                                   
      ALPHA=-0.5*ABS(ALPHA)                                                     
      IFLAG=1                                                                   
      GO TO 1463                                                                
 1454 IF(IER.GT.0)write(*,1455)IER,ITER                                         
 1455 FORMAT(1X,'LOSS OF SIGNIFICANCE. STEP ',I5,'+1   ---   ON THE ',          
     1I3,'TH ITERATION'///)                                                     
      CALL MSTR(XXDUM,VCOV,24,1,0)                                              
      IF(LSTAT.EQ.1)LSTAT=2                                                     
      IF(LSTAT.EQ.2)RETURN                                                      
C                                                                               
C     FORM SOLUTION VECTOR                                                      
C                                                                               
      DO 1460 K=1,MTRX                                                          
      ESTCOR(K)=0.                                                              
      DO 1460 L=1,MTRX                                                          
      ESTCOR(K)=ESTCOR(K)+VCOV(K,L)*FOP(L)                                      
 1460 CONTINUE                                                                  
C                                                                               
C     ITERATE                                                                   
C                                                                               
 1462 ALPHA=1.0                                                                 
 1463 AX=AX+ALPHA*ESTCOR(1)                                                     
      BX=BX+ALPHA*ESTCOR(2)                                                     
      AY=AY+ALPHA*ESTCOR(3)                                                     
      BY=BY+ALPHA*ESTCOR(4)                                                     
      IF (IFLAG.EQ.1) GO TO 1464                                                
      R=R+ALPHA*ESTCOR(5)                                                       
      RHO=RHO+ALPHA*ESTCOR(6)                                                   
      GO TO 1466                                                                
 1464 R=R2+ALPHA*ESTCOR(5)                                                      
      RHO=RHO2+ALPHA*ESTCOR(6)                                                  
 1466 R2=R                                                                      
      RHO2=RHO                                                                  
      DO 1465 I=1,KKI                                                           
      T(I+1)=T(I+1)+ALPHA*ESTCOR(I+6)                                           
 1465 CONTINUE                                                                  
      DO 1470 J=1,KKJ                                                           
      U(J+1)=U(J+1)+ALPHA*ESTCOR(J+KKI+6)                                       
 1470 CONTINUE                                                                  
      DO 1471 I=2,KATI                                                          
      TX(I-1)=T(I)*BX-AX                                                        
      CALL NDTR(TX(I-1),ELP(I,NNJ),D)                                           
 1471 CALL NDTR(T(I),ELR(I,NNJ),D)                                              
      DO 1472 J=2,KATJ                                                          
      UY(J-1)=U(J)*BY-AY                                                        
      CALL NDTR(UY(J-1),ELP(NNI,J),D)                                           
 1472 CALL NDTR(U(J),ELR(NNI,J),D)                                              
C                                                                               
C     COMPUTE THE TERMS USED IN CALCULATION OF LIKELIHOOD FUNCTION              
C                                                                               
      RHO1=RHO                                                                  
      R1=R                                                                      
      IF (RHO1.GT.0.9999) RHO1=0.9999                                           
      IF (RHO1.LT.-0.9999) RHO1=-0.9999                                         
      IF (R1.GT.0.9999) R1=0.9999                                               
      IF (R1.LT.-0.9999) R1=-0.9999                                             
      DO 1473 I=2,KATI                                                          
      DO 1473 J=2,KATJ                                                          
      CALL MDBNOR(TX(I-1),UY(J-1),RHO1,ELP(I,J),IER)                            
 1473 CALL MDBNOR(T(I),U(J),R1,ELR(I,J),IER)                                    
      DO 1474 I=1,KATI                                                          
      DO 1474 J=1,KATJ                                                          
      CPELR(I,J)=ELR(I+1,J+1)+ELR(I,J)-ELR(I+1,J)-ELR(I,J+1)                    
      CPELP(I,J)=ELP(I+1,J+1)+ELP(I,J)-ELP(I+1,J)-ELP(I,J+1)                    
      IF(CPELR(I,J).LT.5.0E-08) CPELR(I,J)=5.0E-08                              
      IF(CPELP(I,J).LT.5.0E-08) CPELP(I,J)=5.0E-08                              
 1474 CONTINUE                                                                  
      IF (ALPHA.NE.1.0) GO TO 1476                                              
C                                                                               
C     CHECK FOR MAXIMIZATION                                                    
C                                                                               
      IF (ITER.EQ.1) GO TO 1476                                                 
      CRIT=0.                                                                   
      DO 1475 I=1,MTRX                                                          
      CRIT=CRIT+ABS(ESTCOR(I))                                                  
 1475 CONTINUE                                                                  
      CRIT=14.0*CRIT/FLOAT(MTRX)                                                
      IF(CRIT.GE..001) GO TO 1476                                               
      LSTAT=1                                                                   
      RETURN                                                                    
C                                                                               
C     CHECK IF LIKELIHOOD FUNCTION VALUES ARE INCREASING, OR THE                
C     DIFFERENCE BETWEEN LIKELIHOOD FUNCTION VALUES ARE QUITE                   
C     SMALL; IF NOT, HALVE THE CORRECTION VECTOR.                               
C                                                                               
 1476 IF (ALPHA.EQ.1.0) FUNLI1=FUNLIK                                           
      FUNLIK=VLOGL(KATI,KATJ,VV,WW,CPELR,CPELP)                                 
      IF (FUNLIK.GE.FUNLI1) GO TO 1478                                          
      FUNDIF=ABS((FUNLIK-FUNLI1)/FUNLI1)                                        
      IF (FUNDIF.LE.1.0E-5) GO TO 1478                                          
      ALPHA=-0.5*ABS(ALPHA)                                                     
      GO TO 1463                                                                
 1478 RETURN                                                                    
      END                                                                       
C----------------------------------------------------------------               
      FUNCTION VLOGL(KATI,KATJ,VV,WW,CPELR,CPELP)                               
C----------------------------------------------------------------               
C
C     COMPUTE THE VALUE OF LOG L                                                
C                                                                               
      DIMENSION VV(11,11),WW(11,11),CPELR(10,10),CPELP(10,10)                   
 1490 SUMN=0.                                                                   
      SUMS=0.                                                                   
      DO 1500 I=1,KATI                                                          
      DO 1500 J=1,KATJ                                                          
      SUMN=SUMN+VV(J,I)*ALOG(CPELR(I,J))                                        
      SUMS=SUMS+WW(J,I)*ALOG(CPELP(I,J))                                        
 1500 CONTINUE                                                                  
      VLOGL=SUMN+SUMS                                                           
      RETURN                                                                    
      END                                                                       
