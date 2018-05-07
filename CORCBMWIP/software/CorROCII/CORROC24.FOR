C-------------------------------------------------------------------            
      SUBROUTINE OUTRSL                                                         
C-------------------------------------------------------------------            
C                                                                               
C     COMPLETE THE PRINTED OUTPUT PAGE                                          
C                                                                               
      COMMON/BLK1/VV(11,11),WW(11,11),KATI,KATJ,TCODE,FP,NNI,NNJ,KKI,           
     1  KSNX,KSNY,KKJ,KKIL,KKJL,MTRX,IDENTI(60),IDENTJ(60)                 
      COMMON/BLK3/TX(9),UY(9),TEXPX(11),TEXPY(11),                              
     1            TXPTUR(11,11),TXPTUP(11,11),FYXP(11,11),FXYP(11,11)
      common/blk9/FUTR(9,11),FTUR(11,9),ELR(11,11),ELP(11,11),                  
     3            CPELR(10,10),CPELP(10,10),SKP,SKR,RADP,RADR,FUNLIK            
      COMMON/BLK4/FOP(24),SOPNEG(24,24),VCOV(24,24),ESTCOR(24)
      common/blk8/CORR(24,24),XXDUM(496),ITER,CRIT,LSTAT                        
      COMMON/BLK7/AX,BX,AY,BY,R,RHO,T(11),U(11)                                 
      common/io/nin,nout,nnew
      INTEGER TCODE                                                             
      write(*,1500)ITER
      write(*,1505)AX,BX,AY,BY,R,RHO                                       
 1500 FORMAT(//9X,'PROCEDURE CONVERGES AFTER ',I3,' ITERATIONS'/)
 1505 FORMAT(16X,'FINAL ESTIMATES OF THE PARAMETERS:'/
     2    6X,'AX=',F7.4,5X,'BX=',F7.4,5X,'AY=',F7.4,5X,'BY=',F7.4/
     3    6X,'R(NEGATIVE CASES)=',F7.4,5X,'R(POSITIVE CASES)=',F7.4)          
      write(*,1506) (T(I),I=2,KATI)                                             
      write(*,1507) (U(I),I=2,KATJ)                                             
 1506 FORMAT(6X,'T(I)',9(F6.3,1x))
 1507 FORMAT(6X,'U(J)',9(F6.3,1x))
      FUNLIK=VLOGL(KATI,KATJ,VV,WW,CPELR,CPELP)                                 
      write(*,1515)FUNLIK                                                       
 1515 FORMAT(/6X,'LOGL =',F10.3)                                               
      write(*,1520)                                                             
 1520 FORMAT(//3X,'VARIANCE-COVARIANCE MATRIX')
c     +            'CORRELATION PARAMETERS ONLY:'//)
	write (*, 1528) MTRX
	do 1527 I = 1,MTRX
	do 1527 J = 1,MTRX
      write(*,1521) VCOV(1,J)                                        
 1521 FORMAT(E14.6)
 1527	continue

 1528	format(I10)

c      write(*,1522) (VCOV(2,J),J=1,6)
c 1522 FORMAT(6X,'BX  ',1X,6(2X,F7.4))
c      write(*,1523) (VCOV(3,J),J=1,6)
c 1523 FORMAT(6X,'AY  ',1X,6(2X,F7.4))
c      write(*,1524) (VCOV(4,J),J=1,6)
c 1524 FORMAT(6X,'BY  ',1X,6(2X,F7.4))
c      write(*,1525) (VCOV(5,J),J=1,6)
c 1525 FORMAT(6X,'R(-)',1X,6(2X,F7.4))
c      write(*,1526) (VCOV(6,J),J=1,6)
c 1526 FORMAT(6X,'R(+)',1X,6(2X,F7.4))
	write(*,*)
	write(*,*)
	write(*,*)
	write(*,*)
c
      if(nout.ne.6)then
	write(nout,*)
	write(nout,*)
	write(nout,*)
        write(nout,1505)AX,BX,AY,BY,R,RHO
      endif
      RETURN                                                                    
      END                                                                       
C-------------------------------------------------------------------            
      SUBROUTINE TPFVAL                                                         
C-------------------------------------------------------------------            
C                                                                               
C     CALCULATE TPF VALUES FOR CONDITION X AND CONDITION Y ON THE FITTED        
C     ROC CURVES, AT SELECTED FPF VALUES.                                       
C                                                                               
      COMMON/BLK4/FOP(24),SOPNEG(24,24),VCOV(24,24),ESTCOR(24)
      common/blk8/CORR(24,24),XXDUM(496),ITER,CRIT,LSTAT                        
      COMMON/BLK5/FPVAL(26),TPVALX(26),TPVALY(26)                               
      COMMON/BLK7/AX,BX,AY,BY,R,RHO,T(11),U(11)                                 
      DO 1605 I=1,26                                                            
      CALL ZDEV(FPVAL(I),Q,D,IER)                                               
      DEVX=AX+BX*Q                                                              
      DEVY=AY+BY*Q                                                              
      CALL NDTR(DEVX,TPVALX(I),D)                                               
      CALL NDTR(DEVY,TPVALY(I),D)                                               
 1605 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
C-----------------------------------------------------------------------        
      SUBROUTINE TEST
C-----------------------------------------------------------------------        
C                                                                               
C     COMPUTE AREA('A SUB Z') FOR CONDITION X AND CONDITION Y RESPECTIVELY;     
C     ALSO CALCULATE THE VARIANCES AND COVARIANCES OF AREAS.                    
C                                                                               
      COMMON/BLK1/VV(11,11),WW(11,11),KATI,KATJ,TCODE,FP,NNI,NNJ,KKI,           
     1  KSNX,KSNY,KKJ,KKIL,KKJL,MTRX,IDENTI(60),IDENTJ(60)                 
      COMMON/BLK4/FOP(24),SOPNEG(24,24),VCOV(24,24),ESTCOR(24)
      common/blk8/CORR(24,24),XXDUM(496),ITER,CRIT,LSTAT                        
      COMMON/BLK7/AX,BX,AY,BY,R,RHO,T(11),U(11)                                 
      COMMON/BLK18/AZX,AZY,VAZX,VAZY,COVAZ,RAREA,XSTAT,PLEVEL,
     1            ZSTAT,ZLEVEL,TPSTAT,TPLEVL                                    
      INTEGER TCODE                                                             
      RDATA=((VCOV(1,2)+VCOV(3,4))-(VCOV(1,4)+VCOV(2,3)))/                      
     1SQRT((VCOV(1,1)+VCOV(3,3)-2.*VCOV(1,3))*(VCOV(2,2)+VCOV(4,4)-2.*          
     2VCOV(2,4)))                                                               
      WORKX=BX*BX+1                                                             
      WORKY=BY*BY+1                                                             
      ADEVX=AX/SQRT(WORKX)                                                      
      ADEVY=AY/SQRT(WORKY)                                                      
      CALL NDTR(ADEVX,AZX,DENX)                                                 
      CALL NDTR(ADEVY,AZY,DENY)                                                 
      VAZX=(DENX*DENX)*(VCOV(1,1)/WORKX+((AX*BX)**2)*VCOV(2,2)/                 
     1     (WORKX**3)-2.*AX*BX*VCOV(1,2)/(WORKX*WORKX))                         
      VAZY=(DENY*DENY)*(VCOV(3,3)/WORKY+((AY*BY)**2)*VCOV(4,4)/                 
     1     (WORKY**3)-2.*AY*BY*VCOV(3,4)/(WORKY*WORKY))                         
      COVAZ=(DENX*DENY)*(VCOV(1,3)/SQRT(WORKX*WORKY)-AY*BY*VCOV(1,4)/           
     1     SQRT(WORKX*(WORKY**3))-AX*BX*VCOV(2,3)/SQRT(WORKY*(WORKX**3))        
     2      +AX*BX*AY*BY*VCOV(2,4)/SQRT((WORKX*WORKY)**3))                      
      RAREA=COVAZ/SQRT(VAZX*VAZY)                                               
C                                                                               
C     OBTAIN THE RESULTS OF REQUESTED STATISTICAL TEST                          
C                                                                               
      IF (TCODE.EQ.2) GO TO 1549                                                
      IF (TCODE.EQ.3) GO TO 1550                                                
C                                                                               
C     OBTAIN THE BIVARIATE CHI-SQUARE TEST STATISTIC VALUE FOR THE              
C     SIGNIFICANCE OF THE DIFFERENCE BETWEEN THE TWO CORRELATED ROC CURVES      
C                                                                               
      XSTAT=1./(1.-RDATA**2)*((AX-AY)**2/(VCOV(1,1)+VCOV(3,3)-2.*               
     1VCOV(1,3))+(BX-BY)**2/(VCOV(2,2)+VCOV(4,4)-2.*VCOV(2,4))-2.*RDATA*        
     2(AX-AY)*(BX-BY)/(SQRT(VCOV(1,1)+VCOV(3,3)-2.*VCOV(1,3))*                  
     3SQRT(VCOV(2,2)+VCOV(4,4)-2.*VCOV(2,4))))                                  
      PLEVEL=EXP(-XSTAT/2.)                                                     
      RETURN                                                                    
C                                                                               
C     CALCULATE THE AREA Z-SCORE TEST STATISTIC VALUE                           
C                                                                               
 1549 ZSTAT=(AZX-AZY)/SQRT(VAZX+VAZY-2.*COVAZ)                                  
      CALL NDTR(ZSTAT,P,D)                                                      
      IF (P.LT.0.5) ZLEVEL=2.*P                                                 
      IF (P.GE.0.5) ZLEVEL=2.*(1.-P)                                            
      RETURN                                                                    
C                                                                               
C     COMPUTE THE TRUE-POSITIVE-FRACTION TEST STATISTIC VALUE AT                
C     FALSE-POSITIVE FRACTION SPECIFIED BY USER.                                
C                                                                               
 1550 P=1-FP                                                                    
      CALL ZDEV(P,Z,D,IER)                                                      
      XTPFX=Z*BX-AX                                                             
      YTPFY=Z*BY-AY                                                             
      VTPX=VCOV(1,1)+Z*Z*VCOV(2,2)-                                             
     1     2.*Z*VCOV(1,2)                                                       
      VTPY=VCOV(3,3)+Z*Z*VCOV(4,4)-                                             
     1     2.*Z*VCOV(3,4)                                                       
      COVTP=VCOV(1,3)-Z*VCOV(1,4)-Z*VCOV(2,3)+                                  
     1      Z*Z*VCOV(2,4)                                                       
      TPSTAT=(XTPFX-YTPFY)/SQRT(VTPX+VTPY-2.*COVTP)                             
      CALL NDTR(TPSTAT,P,D)                                                     
      IF (P.LT.0.5) TPLEVL=2.*P                                                 
      IF (P.GT.0.5) TPLEVL=2.*(1.-P)                                            
      RETURN                                                                    
      END                                                                       
C--------------------------------------------------------------------           
      SUBROUTINE INFORM
C--------------------------------------------------------------------           
C                                                                               
C     OUTPUT THE TPF VALUES AT SELECTED FPF VALUES, AND ALSO                    
C     AREAS, VARIANCES, COVARIANCE, AND CORRELATION.                            
C                                                                               
      COMMON/BLK1/VV(11,11),WW(11,11),KATI,KATJ,TCODE,FP,NNI,NNJ,KKI,           
     1  KSNX,KSNY,KKJ,KKIL,KKJL,MTRX,IDENTI(60),IDENTJ(60)                 
      COMMON/BLK5/FPVAL(26),TPVALX(26),TPVALY(26)                               
      COMMON/BLK7/AX,BX,AY,BY,R,RHO,T(11),U(11)                                 
      COMMON/BLK18/AZX,AZY,VAZX,VAZY,COVAZ,RAREA,XSTAT,PLEVEL,
     1            ZSTAT,ZLEVEL,TPSTAT,TPLEVL                                    
      common/io/nin,nout,nnew
      INTEGER TCODE                                                             
      write(*,1605)
 1605 format(1h1)
      write(*,1610)                                                             
 1610 FORMAT(7X,1H','PLOTS',1H',' OF THE FITTED BINORMAL',
     +      ' ROC CURVES:',
     1      //3X,'**************************************************',
     2      ///12X,'FPF',7X,'TPF FOR',8X,'TPF FOR'/21X,'CONDITION X',
     3      4X,'CONDITION Y'/11X,'-----',4X,'-------------',2X,               
     4      '-------------'/)                                                  
      DO 1630 I=1,26                                                            
      write(*,1625) FPVAL(I),TPVALX(I),TPVALY(I)                                
 1625 FORMAT(11X,F5.3,6X,F6.4,9X,F6.4)                                          
 1630 CONTINUE
	stdazx=sqrt(vazx)
	stdazy=sqrt(vazy)
      write(*,1635) AZX,AZY,stdazx,stdazy,RAREA
 1635 FORMAT(////12X,'AREA(X)= ',F7.4,5X,'AREA(Y)= ',F7.4/12X,               
     1 'STD DEV OF AREA(X)= ',F7.4,5X,'STD DEV OF AREA(Y)= ',
     2  F7.4/
     3 12X,'CORRELATION OF AREA(X) AND AREA(Y) =',F7.4)
C                                                                               
C     OUTPUT THE RESULTS OF THE REQUESTED STATISTICAL TEST                      
C                                                                               
 1639 write(*,1640)                                                             
 1640 FORMAT(/////13X,'STATISTICAL SIGNIFICANCE OF THE DIFFERENCE',/           
     1      6X,'BETWEEN THE TWO CORRELATED ROC CURVES ACCORDING TO',           
     2      ' THE SELECTED TEST:'//
     3  2X,'******************************************************',
     4      '***********************')
      IF (TCODE.EQ.2) GO TO 1642                                                
      IF (TCODE.EQ.3) GO TO 1648                                                
      write(*,1641) XSTAT,PLEVEL                                                
 1641 FORMAT(//4X'THE COMPUTED CORRELATED ',1H','BIVARIATE CHI-',             
     1       'SQUARE TEST',1H',' STATISTIC VALUE IS',F8.4/4X,                 
     2       'WITH A CORRESPONDING P-LEVEL OF',F7.4,'.')                        
      go to 2600                                                                
 1642 ONEPL=0.5*ZLEVEL                                                          
      write(*,1645) ZSTAT,ZLEVEL,ONEPL                                          
 1645 FORMAT(//4X'THE COMPUTED CORRELATED ',1H','AREA TEST',1H',              
     1' STATISTIC VALUE IS',2X,F8.4/4X,'WITH A CORRESPONDING TWO-',           
     2'TAILED P-LEVEL OF',F7.4,/4x,'AND ONE-TAILED P-LEVEL OF',F7.4,'.')
      go to 2600                                                                
 1648 ONEPL=0.5*TPLEVL                                                          
      write(*,1650) FP,TPSTAT,TPLEVL,ONEPL                                      
 1650 FORMAT(//4X,'AT FALSE-POSITIVE-FRACTION EQUAL TO ',F4.2,', '
     1    /4X,'THE COMPUTED CORRELATED ',1H', 'TRUE-POSITIVE-FRACTION',        
     2      ' TEST',1H'/4x,'STATISTIC VALUE IS',F8.4,' WITH A CORRES',        
     3     'PONDING TWO-TAILED P-LEVEL OF',F7.4,/4x,'AND ONE-TAILED P-',
     4      'LEVEL OF',F7.4,'.')                                                
c
 2600 if(nout.eq.3)then
        write(nout,2610) AZX,stdazx,AZY,stdazy
 2610   format(/6X,'AREA INDEX (X) =',F7.4,3X,'WITH STD. DEV. =',F7.4/
     +         6X,'AREA INDEX (Y) =',F7.4,3X,'WITH STD. DEV. =',F7.4///)
	write(nout,1610)
	DO 2630 I=1,26
	  write(nout,1625) FPVAL(I),TPVALX(I),TPVALY(I)
 2630	  CONTINUE
      endif
      RETURN                                                                    
      END                                                                       
