C------------------------------------------------------------------             
      SUBROUTINE HEADLN
C------------------------------------------------------------------             
C                                                                               
C     CALCULATE THE MARGINAL DATA FOR THE TWO EXPERIMENTAL CONDITIONS,          
C     RESPECTIVELY, BY SUMMING THE 'NOISE' AND 'SINGLE-PLUS-NOISE' MATRICES     
C     IN ROWS (FOR MARGINAL DATA OF CONDITION X) AND IN COLUMNS (FOR            
C     MARGINAL DATA OF CONDITION Y).                                            
C                                                                               
      REAL RVV(11,11),RWW(11,11),RSUMNT(11),RSUMST(11),RSUMNU(11),              
     1     RSUMSU(11)                                                           
      COMMON/BLK1/VV(11,11),WW(11,11),KATI,KATJ,TCODE,FP,NNI,NNJ,KKI,           
     1  KSNX,KSNY,KKJ,KKIL,KKJL,MTRX,IDENTI(60),IDENTJ(60)
      COMMON/BLK2/EMN,EMS,SUMNT(11),SUMST(11),SUMNU(11),SUMSU(11),              
     1            FPFX(11),TPFX(11),FPFY(11),TPFY(11)                           
      common/io/nin,nout,nnew
      INTEGER TCODE                                                             
      DO 40 I=1,KATI                                                          
      SUMNT(I)=0.                                                               
      SUMST(I)=0.                                                               
 40   CONTINUE                                                                  
      DO 45 J=1,KATJ                                                          
      SUMNU(J)=0.                                                               
      SUMSU(J)=0.                                                               
 45   CONTINUE                                                                  
      DO 50 I=1,KATI                                                          
      DO 50 J=1,KATJ                                                          
      SUMNT(I)=SUMNT(I)+VV(J,I)                                                 
      SUMST(I)=SUMST(I)+WW(J,I)                                                 
 50   CONTINUE                                                                  
      DO 55 J=1,KATJ                                                          
      DO 55 I=1,KATI                                                          
      SUMNU(J)=SUMNU(J)+VV(J,I)                                                 
      SUMSU(J)=SUMSU(J)+WW(J,I)                                                 
 55   CONTINUE                                                                  
      EMN=0.                                                                    
      EMS=0.                                                                    
      DO 60 I=1,KATI                                                          
      EMN=EMN+SUMNT(I)                                                          
      EMS=EMS+SUMST(I)                                                          
 60   CONTINUE                                                                  
C                                                                               
C     PRINT PAGE HEADING, IDENTIFICATION NAMES, AND RATING-DATA MATRICES        
C     FOR NOISE CASES AND SIGNAL-PLUS-NOISE CASES.                             
C                                                                               
      write(*,70)                                                             
 70   FORMAT(1H1,26X,'C O R R O C 2 :'//,
     + 5X,'M A X I M U M   L I K E L I H O O D   E S T I '           
     1,'M A T I O N'/19X,'O F   T H E   P A R A M E T E R S'/5X,               
     2'O F   T H E   B I V A R I A T E   B I N O R M A L   ',                   
     3'M O D E L'/9X,'F O R   C O R R E L A T E D   R A T I N G   ',           
     4'D A T A')                                                                
      write(*,71)                                                             
 71   FORMAT(//6X,'STATISTICAL TEST TO BE EMPLOYED:')                          
      IF(TCODE.EQ.1)then
	write(*,72)
      else IF(TCODE.EQ.2)then
	write(*,73)                                               
      else 
	write(*,74)FP                                             
      endif
 72   FORMAT(10X,'CHI-SQUARE TEST FOR THE DIFFERENCE BETWEEN ',
     1  '(AX,BX) AND (AY,BY)')                                                  
 73   FORMAT(10X,'AREA (A SUB Z) TEST')                                         
 74   FORMAT(10X,'TPF TEST AT FPF = ',F7.4)                                     
      write(*,76)(IDENTI(I),I=1,60)
 76   FORMAT(//6X,'CONDITION X: ',60A1)
      write(*,77)KATI,KSNX
 77   format(6X,'DATA COLLECTED IN ',I2,' CATEGORIES'/
     2 6X,'WITH CATEGORY ',I2,' REPRESENTING STRONGEST EVIDENCE ',              
     3'OF POSITIVITY'/6x,'(E.G., THAT ABNORMALITY IS PRESENT).')
      write(*,78)(IDENTJ(J),J=1,60)
 78   format(//6X,'CONDITION Y: ',60A1)
      write(*,77)KATJ,KSNY
      write(*,80)(I,I=1,KATI)
 80   format(///6X,'RATING-DATA MATRIX FOR ACTUALLY NEGATIVE CASES:'//
     2          20X,'CONDITION X RATING'/3X,'CONDITION',2X,10(I2,4X))
      write(*,81)                                                             
 81   FORMAT(3X,'Y RATING')                                                    
      DO 86 J=1,KATJ                                                          
      JOPP=NNJ-J                                                                
      write(*,85)JOPP,(VV(JOPP,I),I=1,KATI),SUMNU(JOPP)                       
 85   FORMAT(7X,I2,4X,11(F4.0,2X),/)                                           
 86   CONTINUE                                                                  
      write(*,87) (SUMNT(I),I=1,KATI),EMN                                     
 87   FORMAT(/6X,'SUMX',3X,11(F4.0,2X)/)                                     
      write(*,88) (I,I=1,KATI)                                                
 88   FORMAT(//6X,'RATING-DATA MATRIX FOR ACTUALLY POSITIVE CASES:',
     1       //20X,'CONDITION X RATING'/3X,'CONDITION',2X,10(I2,4X))
      write(*,81)                                                             
      DO 89 J=1,KATJ                                                          
      JOPP=NNJ-J                                                                
      write(*,85)JOPP,(WW(JOPP,I),I=1,KATI),SUMSU(JOPP)                       
 89   CONTINUE                                                                  
      write(*,87) (SUMST(I),I=1,KATI),EMS                                     
      IF(KSNX.EQ.KATI)GO TO 200                                                
      DO 110 I=1,KATI                                                          
      DO 110 J=1,KATJ                                                          
      RVV(J,I)=VV(J,KATI+1-I)                                                   
 110  RWW(J,I)=WW(J,KATI+1-I)                                                   
      DO 120 I=1,KATI                                                          
      DO 120 J=1,KATJ                                                          
      VV(J,I)=RVV(J,I)                                                          
 120  WW(J,I)=RWW(J,I)                                                          
      DO 130 I=1,KATI                                                          
      RSUMNT(KATI-I+1)=SUMNT(I)                                                 
 130  RSUMST(KATI-I+1)=SUMST(I)                                                 
      DO 140 I=1,KATI                                                          
      SUMNT(I)=RSUMNT(I)                                                        
 140  SUMST(I)=RSUMST(I)                                                        
C                                                                               
 200  IF(KSNY.EQ.KATJ)GO TO 300                                                
      DO 210 I=1,KATI                                                          
      DO 210 J=1,KATJ                                                          
      RVV(J,I)=VV(KATJ+1-J,I)                                                   
 210  RWW(J,I)=WW(KATJ+1-J,I)                                                   
      DO 220 I=1,KATI                                                          
      DO 220 J=1,KATJ                                                          
      VV(J,I)=RVV(J,I)                                                          
 220  WW(J,I)=RWW(J,I)                                                          
      DO 230 J=1,KATJ                                                          
      RSUMNU(KATJ-J+1)=SUMNU(J)                                                 
 230  RSUMSU(KATJ-J+1)=SUMSU(J)                                                 
      DO 240 J=1,KATJ                                                          
      SUMNU(J)=RSUMNU(J)                                                        
 240  SUMSU(J)=RSUMSU(J)                                                        
c
 300  if(nout.ne.6)then
	write(nout,70)
	write(nout,76)(identi(i),i=1,60)
	write(nout,78)(identj(i),i=1,60)
      endif
      RETURN                                                                    
      END                                                                       
C------------------------------------------------------------------             
      SUBROUTINE CLAPSE                                                         
C-----------------------------------------------------------------              
C                                                                               
C       COLLAPSE DATA                                                           
C                                                                               
      COMMON/BLK1/VV(11,11),WW(11,11),KATI,KATJ,TCODE,FP,NNI,NNJ,KKI,           
     1  KSNX,KSNY,KKJ,KKIL,KKJL,MTRX,IDENTI(60),IDENTJ(60)                 
      COMMON/BLK2/EMN,EMS,SUMNT(11),SUMST(11),SUMNU(11),SUMSU(11),              
     1            FPFX(11),TPFX(11),FPFY(11),TPFY(11)                           
      INTEGER TCODE                                                             
      KATI1=KATI                                                                
      DO 1005 I=1,KATI                                                          
      IT=KATI-I+1                                                               
      IF (SUMNT(IT).NE.0.OR.SUMST(IT).NE.0) GO TO 1005                          
      KATI1=KATI1-1                                                             
      IF (IT.GT.KATI1) GO TO 1005                                               
      DO 1004 J=IT,KATI1                                                        
      SUMNT(J)=SUMNT(J+1)                                                       
      SUMST(J)=SUMST(J+1)                                                       
      DO 1003 IX=1,KATJ                                                         
      VV(IX,J)=VV(IX,J+1)                                                       
      WW(IX,J)=WW(IX,J+1)                                                       
 1003 CONTINUE                                                                  
 1004 CONTINUE                                                                  
 1005 CONTINUE                                                                  
      KATI=KATI1                                                                
      KATJ1=KATJ                                                                
      DO 1015 I=1,KATJ                                                          
      IT=KATJ-I+1                                                               
      IF (SUMNU(IT).NE.0.OR.SUMSU(IT).NE.0) GO TO 1015                          
      KATJ1=KATJ1-1                                                             
      IF (IT.GT.KATJ1) GO TO 1015                                               
      DO 1010 J=IT,KATJ1                                                        
      SUMNU(J)=SUMNU(J+1)                                                       
      SUMSU(J)=SUMSU(J+1)                                                       
      DO 1007 IY=1,KATI                                                         
      VV(J,IY)=VV(J+1,IY)                                                       
      WW(J,IY)=WW(J+1,IY)                                                       
 1007 CONTINUE                                                                  
 1010 CONTINUE                                                                  
 1015 CONTINUE                                                                  
      KATJ=KATJ1                                                                
      KKI=KATI-1                                                                
      KKJ=KATJ-1                                                                
      NNI=KKI+2                                                                 
      NNJ=KKJ+2                                                                 
      KKIL=KKI-1                                                                
      KKJL=KKJ-1                                                                
      RETURN                                                                    
      END                                                                       
C---------------------------------------------------------------                
      SUBROUTINE TPFFPF                                                         
C---------------------------------------------------------------                
C                                                                               
C     CALCULATE OBSERVED TRUE POSITIVE FRACTIONS AND                            
C     FALSE POSITIVE FRACTIONS AND OUTPUT THE OBSERVED                          
C     OPERATING POINTS                                                          
C                                                                               
      COMMON/BLK1/VV(11,11),WW(11,11),KATI,KATJ,TCODE,FP,NNI,NNJ,KKI,           
     1  KSNX,KSNY,KKJ,KKIL,KKJL,MTRX,IDENTI(60),IDENTJ(60)                 
      COMMON/BLK2/EMN,EMS,SUMNT(11),SUMST(11),SUMNU(11),SUMSU(11),              
     1            FPFX(11),TPFX(11),FPFY(11),TPFY(11)                           
      common/io/nin,nout,nnew
      INTEGER TCODE                                                             
      FPFX(NNI)=0.                                                              
      TPFX(NNI)=0.                                                              
      FPFX(KATI)=SUMNT(KATI)                                                    
      TPFX(KATI)=SUMST(KATI)                                                    
      DO 1020 I=2,KATI                                                          
      FPFX(KATI-I+1)=FPFX(KATI-I+2)+SUMNT(KATI-I+1)                             
      TPFX(KATI-I+1)=TPFX(KATI-I+2)+SUMST(KATI-I+1)                             
 1020 CONTINUE                                                                  
      FPFY(NNJ)=0.                                                              
      TPFY(NNJ)=0.                                                              
      FPFY(KATJ)=SUMNU(KATJ)                                                    
      TPFY(KATJ)=SUMSU(KATJ)                                                    
      DO 1025 I=2,KATJ                                                          
      FPFY(KATJ-I+1)=FPFY(KATJ-I+2)+SUMNU(KATJ-I+1)                             
      TPFY(KATJ-I+1)=TPFY(KATJ-I+2)+SUMSU(KATJ-I+1)                             
 1025 CONTINUE                                                                  
      DO 1030 I=1,KATI                                                          
      FPFX(I)=FPFX(I)/EMN                                                       
      TPFX(I)=TPFX(I)/EMS                                                       
 1030 CONTINUE                                                                  
      DO 1035 I=1,KATJ                                                          
      FPFY(I)=FPFY(I)/EMN                                                       
      TPFY(I)=TPFY(I)/EMS                                                       
 1035 CONTINUE                                                                  
      write(*,1090)                                                             
 1090 FORMAT(/12X,'OBSERVED OPERATING POINTS FOR CONDITION X:')
      write(*,1091) (FPFX(NNI-I+1),I=1,NNI)                                     
 1091 FORMAT(6X,' FPF: ',11(F5.3,1X))                                          
      write(*,1092) (TPFX(NNI-I+1),I=1,NNI)                                     
 1092 FORMAT(6X,' TPF: ',11(F5.3,1X))                                          
      write(*,1093)                                                             
 1093 FORMAT(/12X,'OBSERVED OPERATING POINTS FOR CONDITION Y:')                
      write(*,1091) (FPFY(NNJ-I+1),I=1,NNJ)                                     
      write(*,1092) (TPFY(NNJ-I+1),I=1,NNJ)                                     
c
      if(nout.ne.6)then
	write(nout,1090)
	write(nout,1091) (FPFX(NNI-I+1),I=1,NNI)
	write(nout,1092) (TPFX(NNI-I+1),I=1,NNI)
	write(nout,1093)
	write(nout,1091) (FPFY(NNJ-I+1),I=1,NNJ)
	write(nout,1092) (TPFY(NNJ-I+1),I=1,NNJ)
      endif
      RETURN                                                                    
      END                                                                       
C------------------------------------------------------------------             
      SUBROUTINE DEGMSG(KEYDEG)                                                 
C------------------------------------------------------------------             
C                                                                               
C     WRITE OUT DEGENERACY MESSAGE FOR CONDITION X OR CONDITION Y               
C                                                                               
      COMMON/BLK2/EMN,EMS,SUMNT(11),SUMST(11),SUMNU(11),SUMSU(11),              
     1            FPFX(11),TPFX(11),FPFY(11),TPFY(11)                           
      COMMON/BLK6/ICLASS,JCLASS,IERX,IERY,ICON,JCON,IARY(11),JARY(11)           
C                                                                               
      IF (ICLASS.EQ.0) GO TO 1120                                               
      IF (ICLASS.GE.5.AND.ICLASS.LE.8) write(*,1105)                            
      IF (ICLASS.LT.5.OR.ICLASS.GT.8) write(*,1110)                             
 1105 FORMAT(///,15X,'MARGINAL DATA FOR CONDITION X ARE DEGENERATE.')           
 1110 FORMAT(///,15X,'MARGINAL DATA FOR CONDITION X ARE DEGENERATE')            
      CALL DEGENE(ICLASS,IARY,FPFX,TPFX)                                   
 1120 IF (JCLASS.EQ.0) GO TO 1135                                               
      IF (JCLASS.GE.5.AND.JCLASS.LE.8) write(*,1125)                            
      IF (JCLASS.LT.5.OR.JCLASS.GT.8) write(*,1130)                             
 1125 FORMAT(///,15X,'MARGINAL DATA FOR CONDITION Y ARE DEGENERATE.')           
 1130 FORMAT(///,15X,'MARGINAL DATA FOR CONDITION Y ARE DEGENERATE')            
      CALL DEGENE(JCLASS,JARY,FPFY,TPFY)                                   
 1135 write(*,1140)                                                             
 1140 FORMAT(//,15X,'EXECUTION OF THIS BIVARIATE DATA SET ABORTED.')            
      KEYDEG=1                                                                  
      RETURN                                                                    
      END                                                                       
C--------------------------------------------------------------------           
      SUBROUTINE ITRMSG                                                         
C--------------------------------------------------------------------           
C                                                                               
C     WRITE OUT MESSAGE OF MATRIX INVERSION PROBLEM OR ITERATIONS               
C     EXCEEDING MAXIMUM RESTRICTION FOR MARGINAL DATA SETS                      
C                                                                               
      COMMON/BLK6/ICLASS,JCLASS,IERX,IERY,ICON,JCON,IARY(11),JARY(11)           
      IF (IERX.EQ.0) GO TO 1205                                                 
      IF (IERX.GT.100) write(*,1203)                                            
      IF (IERX.LE.100) write(*,1204) IERX                                       
 1203 FORMAT(6X,'PROCEDURE DOES NOT CONVERGE -- 100 ITERATIONS --'/
     1       6x,'THIS WAS WITH MARGINAL DATA SET FOR CONDITION X')
 1204 FORMAT(6X,'PROBLEM WITH MATRIX INVERSION ON ',I3,'TH ITERATION'/
     2       6x,' WITH MARGINAL DATA SET FOR CONDITION X')
      IF (IERY.EQ.0) RETURN                                                     
 1205 IF (IERY.GT.100) write(*,1206)                                            
      IF (IERY.LE.100) write(*,1207) IERY                                       
 1206 FORMAT(6X,'PROCEDURE DOES NOT CONVERGE -- 100 ITERATIONS --'/
     1       6x,'THIS WAS WITH MARGINAL DATA SET FOR CONDITION Y')
 1207 FORMAT(6X,'PROBLEM WITH MATRIX INVERSION ON ',I3,'TH ITERATION'/
     2       6x,' WITH MARGINAL DATA SET FOR CONDITION Y')
      RETURN                                                                    
      END                                                                       
C-----------------------------------------------------------------              
      SUBROUTINE CORCOE                                                         
C-----------------------------------------------------------------              
C                                                                               
C     CALCULATE INITIAL ESTIMATES FOR R=R(NOISE) AND RHO=R(SIGNAL-PLUS-NOISE)   
C     BY COMPUTING PEARSON PRODUCT MOMENT CORRELATION COEFFICIENTS FOR THE      
C     NOISE RATING-DATA MATRIX AND THE SIGNAL-PLUS-NOISE RATING-DATA MATRIX,    
C     RESPECTIVELY.                                                             
C                                                                               
      COMMON/BLK1/VV(11,11),WW(11,11),KATI,KATJ,TCODE,FP,NNI,NNJ,KKI,           
     1  KSNX,KSNY,KKJ,KKIL,KKJL,MTRX,IDENTI(60),IDENTJ(60)                 
      COMMON/BLK2/EMN,EMS,SUMNT(11),SUMST(11),SUMNU(11),SUMSU(11),              
     1            FPFX(11),TPFX(11),FPFY(11),TPFY(11)                           
      COMMON/BLK7/AX,BX,AY,BY,R,RHO,T(11),U(11)                                 
      INTEGER TCODE                                                             
      SUMRX=0.                                                                  
      SUMPX=0.                                                                  
      SUMRY=0.                                                                  
      SUMPY=0.                                                                  
      SUMRX2=0.                                                                 
      SUMPX2=0.                                                                 
      SUMRY2=0.                                                                 
      SUMPY2=0.                                                                 
      SUMRXY=0.                                                                 
      SUMPXY=0.                                                                 
      DO 1065 I=1,KATI                                                          
      SUMRX=SUMRX+SUMNT(I)*I                                                    
      SUMPX=SUMPX+SUMST(I)*I                                                    
      SUMRX2=SUMRX2+SUMNT(I)*I**2                                               
      SUMPX2=SUMPX2+SUMST(I)*I**2                                               
 1065 CONTINUE                                                                  
      DO 1070 J=1,KATJ                                                          
      SUMRY=SUMRY+SUMNU(J)*J                                                    
      SUMPY=SUMPY+SUMSU(J)*J                                                    
      SUMRY2=SUMRY2+SUMNU(J)*J**2                                               
      SUMPY2=SUMPY2+SUMSU(J)*J**2                                               
 1070 CONTINUE                                                                  
      DO 1075 I=1,KATI                                                          
      DO 1075 J=1,KATJ                                                          
      SUMRXY=SUMRXY+VV(J,I)*J*I                                                 
      SUMPXY=SUMPXY+WW(J,I)*J*I                                                 
 1075 CONTINUE                                                                  
      R=(EMN*SUMRXY-SUMRX*SUMRY)/SQRT((EMN*SUMRX2-SUMRX**2)*(EMN*SUMRY2         
     1-SUMRY**2))                                                               
      RHO=(EMS*SUMPXY-SUMPX*SUMPY)/SQRT((EMS*SUMPX2-SUMPX**2)*(EMS*             
     1SUMPY2-SUMPY**2))                                                         
      RETURN                                                                    
      END                                                                       
C-------------------------------------------------------------------            
      SUBROUTINE OUTINI(VAX,VBX,VABX,VAY,VBY,VABY)                              
C-------------------------------------------------------------------            
C                                                                               
C     PRINT INITIAL ESTIMATES OF ALL OF THE MODEL PARAMETERS TO BE              
C     USED NEXT IN CALCULATING PARAMETERS OF THE BIVARIATE MODEL BY             
C     THE 'METHOD OF SCORING'.                                                  
C                                                                               
      COMMON/BLK1/VV(11,11),WW(11,11),KATI,KATJ,TCODE,FP,NNI,NNJ,KKI,           
     1  KSNX,KSNY,KKJ,KKIL,KKJL,MTRX,IDENTI(60),IDENTJ(60)                 
      COMMON/BLK7/AX,BX,AY,BY,R,RHO,T(11),U(11)                                 
      INTEGER TCODE                                                             
C                                                                               
      write(*,1095) AX,BX,AY,BY,R,RHO                                           
 1095 FORMAT(//16X,'INITIAL ESTIMATES OF THE PARAMETERS:'/
     16X,'AX=',F7.4,5X,'BX=',F7.4,5X,'AY=',F7.4,5X,'BY=',F7.4/
     26X,'R(NEGATIVE CASES)=',F7.4,5X,'R(POSITIVE CASES)=',F7.4)
      write(*,1096) (T(I),I=2,KATI)                                             
 1096 FORMAT(6X,'T(I)',9(F6.3,1x))                                             
      write(*,1097) (U(I),I=2,KATJ)                                             
 1097 FORMAT(6X,'U(J)',9(F6.3,1x))                                             
      CORX=VABX/SQRT(VAX*VBX)                                                   
      CORY=VABY/SQRT(VAY*VBY)                                                   
c     write(*,1098) VAX,VBX,VABX,CORX,VAY,VBY,VABY,CORY                         
c1098 FORMAT(6X,'(VAX=',F6.4,2X,'VBX=',F6.4,2X,'COVABX=',F7.4,2X,               
c    1'CORABX=',F7.4/7X,'VAY=',F6.4,2X,'VBY=',F6.4,2X,'COVABY=',F7.4,2X,        
c    2'CORABY=',F7.4,')')                                                       
      RETURN                                                                    
      END                                                                       
C---------------------------------------------------------------------          
      SUBROUTINE INIVAR                                                         
C---------------------------------------------------------------------          
C                                                                               
C     INITIALIZE VARIABLES                                                      
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
      INTEGER TCODE                                                             
C                                                                               
      TEXPX(1)=0.                                                               
      TEXPX(NNI)=0.                                                             
      TEXPY(1)=0.                                                               
      TEXPY(NNI)=0.                                                             
      DO 1100 J=1,NNJ                                                           
      FYXP(1,J)=0.                                                              
      FYXP(NNI,J)=0.                                                            
      ELR(1,J)=0.                                                               
      ELP(1,J)=0.                                                               
      TXPTUR(1,J)=0.                                                            
      TXPTUR(NNI,J)=0.                                                          
      TXPTUP(1,J)=0.                                                            
      TXPTUP(NNI,J)=0.                                                          
 1100 CONTINUE                                                                  
      DO 1105 I=2,KATI                                                          
      FYXP(I,1)=0.                                                              
      FYXP(I,NNJ)=1.                                                            
      TXPTUR(I,1)=0.                                                            
      TXPTUR(I,NNJ)=0.                                                          
      TXPTUP(I,1)=0.                                                            
      TXPTUP(I,NNJ)=0.                                                          
      FUTR(I-1,NNJ)=1.                                                          
      FUTR(I-1,1)=0.                                                            
 1105 CONTINUE                                                                  
      DO 1110 I=1,NNI                                                           
      FXYP(I,1)=0.                                                              
      FXYP(I,NNJ)=0.                                                            
 1110 CONTINUE                                                                  
      DO 1115 J=2,KATJ                                                          
      FXYP(1,J)=0.                                                              
      FXYP(NNI,J)=1.                                                            
      FTUR(NNI,J-1)=1.                                                          
      FTUR(1,J-1)=0.                                                            
 1115 CONTINUE                                                                  
      DO 1120 I=2,NNI                                                           
      ELP(I,1)=0.                                                               
      ELR(I,1)=0.                                                               
 1120 CONTINUE                                                                  
      ELP(NNI,NNJ)=1.                                                           
      ELR(NNI,NNJ)=1.                                                           
      MTRX=NNI+NNJ+2                                                            
      DO 1125 I=1,MTRX                                                          
      DO 1125 J=1,MTRX                                                          
      SOPNEG(I,J)=0.                                                            
      VCOV(I,J)=0.                                                              
      CORR(I,J)=0.                                                              
 1125 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
C---------------------------------------------------------------------          
      SUBROUTINE TERMS                                                          
C---------------------------------------------------------------------          
C                                                                               
C     COMPUTE THE MAJOR TERMS THAT OCCUR IN THE DERIVATIVE EXPRESSIONS.         
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
      IF (RHO.GT.0.9999) RHO=0.9999                                             
      IF (R.GT.0.9999) R=0.9999                                                 
      IF (RHO.LT.-0.9999) RHO=-0.9999                                           
      IF (R.LT.-0.9999) R=-0.9999                                               
      SKP=1.-RHO**2                                                             
      SKR=1.-R**2                                                               
      RADP=SQRT(SKP)                                                            
      RADR=SQRT(SKR)                                                            
      DO 1135 I=2,KATI                                                          
      TX(I-1)=T(I)*BX-AX                                                        
      TEXPX(I)=EXP(-.5*TX(I-1)**2)                                              
 1135 CONTINUE                                                                  
      DO 1140 J=2,KATJ                                                          
      UY(J-1)=U(J)*BY-AY                                                        
      TEXPY(J)=EXP(-.5*UY(J-1)**2)                                              
 1140 CONTINUE                                                                  
      DO 1145 I=1,KKI                                                           
      DO 1145 J=2,KATJ                                                          
      ARG=(U(J)-R*T(I+1))/RADR                                                  
      CALL NDTR(ARG,FUTR(I,J),D)                                                
 1145 CONTINUE                                                                  
      DO 1150 J=1,KKJ                                                           
      DO 1150 I=2,KATI                                                          
      ARG=(T(I)-R*U(J+1))/RADR                                                  
      CALL NDTR(ARG,FTUR(I,J),D)                                                
 1150 CONTINUE                                                                  
      DO 1155 I=2,KATI                                                          
      DO 1155 J=2,KATJ                                                          
      TXPTUR(I,J)=EXP(-1./(2.*SKR)*(T(I)**2-2.*R*T(I)*U(J)+U(J)**2))            
      TXPTUP(I,J)=EXP(-1./(2.*SKP)*(TX(I-1)**2-2.*RHO*TX(I-1)*UY(J-1)+          
     1UY(J-1)**2))                                                              
      ARG1=(UY(J-1)-RHO*TX(I-1))/RADP                                           
      CALL NDTR(ARG1,FYXP(I,J),D)                                               
      ARG2=(TX(I-1)-RHO*UY(J-1))/RADP                                           
      CALL NDTR(ARG2,FXYP(I,J),D)                                               
 1155 CONTINUE                                                                  
      IF (ITER.GT.1) RETURN                                                     
      DO 1156 I=2,KATI                                                          
      CALL NDTR(TX(I-1),ELP(I,NNJ),D)                                           
 1156 CALL NDTR (T(I),ELR(I,NNJ),D)                                             
      DO 1157 J=2,KATJ                                                          
      CALL NDTR(UY(J-1),ELP(NNI,J),D)                                           
 1157 CALL NDTR(U(J),ELR(NNI,J),D)                                              
      DO 1158 I=2,KATI                                                          
      DO 1158 J=2,KATJ                                                          
      CALL MDBNOR(TX(I-1),UY(J-1),RHO,ELP(I,J),IER)                             
 1158 CALL MDBNOR(T(I),U(J),R,ELR(I,J),IER)                                     
      DO 1159 I=1,KATI                                                          
      DO 1159 J=1,KATJ                                                          
      CPELR(I,J)=ELR(I+1,J+1)+ELR(I,J)-ELR(I+1,J)-ELR(I,J+1)                    
      CPELP(I,J)=ELP(I+1,J+1)+ELP(I,J)-ELP(I+1,J)-ELP(I,J+1)                    
      IF(CPELR(I,J).LT.5.0E-08) CPELR(I,J)=5.0E-08                              
      IF(CPELP(I,J).LT.5.0E-08) CPELP(I,J)=5.0E-08                              
 1159 CONTINUE                                                                  
      FUNLIK=VLOGL(KATI,KATJ,VV,WW,CPELR,CPELP)                                 
      write(*,1160)FUNLIK                                                       
 1160 FORMAT(/6X,'LOGL =',F10.3)                                               
      RETURN                                                                    
      END                                                                       
C---------------------------------------------------------------------          
      SUBROUTINE FIRST                                                          
C---------------------------------------------------------------------          
C                                                                               
C     COMPUTE THE VALUE OF EACH OF THE FIRST-ORDER PARTIAL DERIVATIVE           
C     EXPRESSIONS                                                               
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
 1162 SUM=0.                                                                    
      DO 1163 I=1,KATI                                                          
      DO 1163 J=1,KATJ                                                          
      SUM=SUM+WW(J,I)*(TEXPX(I)*(FYXP(I,J+1)-FYXP(I,J))-TEXPX(I+1)*             
     1(FYXP(I+1,J+1)-FYXP(I+1,J)))/CPELP(I,J)                                   
 1163 CONTINUE                                                                  
      FOP(1)=SUM/SQRT(2.*GPI)                                                   
C                                                                               
C     WITH RESPECT TO BX                                                        
C                                                                               
      SUM=0.                                                                    
      DO 1165 I=1,KATI                                                          
      DO 1165 J=1,KATJ                                                          
      SUM=SUM+WW(J,I)*(-T(I)*TEXPX(I)*(FYXP(I,J+1)-FYXP(I,J))+T(I+1)*           
     1TEXPX(I+1)*(FYXP(I+1,J+1)-FYXP(I+1,J)))/CPELP(I,J)                        
 1165 CONTINUE                                                                  
      FOP(2)=SUM/SQRT(2.*GPI)                                                   
C                                                                               
C     WITH RESPECT TO AY                                                        
C                                                                               
      SUM=0.                                                                    
      DO 1170 I=1,KATI                                                          
      DO 1170 J=1,KATJ                                                          
      SUM=SUM+WW(J,I)*(TEXPY(J)*(FXYP(I+1,J)-FXYP(I,J))-TEXPY(J+1)*             
     1(FXYP(I+1,J+1)-FXYP(I,J+1)))/CPELP(I,J)                                   
 1170 CONTINUE                                                                  
      FOP(3)=SUM/SQRT(2.*GPI)                                                   
C                                                                               
C     WITH RESPECT TO BY                                                        
C                                                                               
      SUM=0.                                                                    
      DO 1175 I=1,KATI                                                          
      DO 1175 J=1,KATJ                                                          
      SUM=SUM+WW(J,I)*(-U(J)*TEXPY(J)*(FXYP(I+1,J)-FXYP(I,J))+U(J+1)*           
     1TEXPY(J+1)*(FXYP(I+1,J+1)-FXYP(I,J+1)))/CPELP(I,J)                        
 1175 CONTINUE                                                                  
      FOP(4)=SUM/SQRT(2.*GPI)                                                   
C                                                                               
C     WITH RESPECT TO R=R(NOISE)                                                
C                                                                               
      SUM=0.                                                                    
      DO 1180 I=1,KATI                                                          
      DO 1180 J=1,KATJ                                                          
      SUM=SUM+VV(J,I)*(TXPTUR(I+1,J+1)+TXPTUR(I,J)-TXPTUR(I+1,J)-               
     1TXPTUR(I,J+1))/CPELR(I,J)                                                 
 1180 CONTINUE                                                                  
      FOP(5)=SUM/(2.*GPI*RADR)                                                  
C                                                                               
C     WITH RESPECT TO RHO=R(SIGNAL-PLUS-NOISE)                                  
C                                                                               
      SUM=0.                                                                    
      DO 1185 I=1,KATI                                                          
      DO 1185 J=1,KATJ                                                          
      SUM=SUM+WW(J,I)*(TXPTUP(I+1,J+1)+TXPTUP(I,J)-TXPTUP(I+1,J)-               
     1TXPTUP(I,J+1))/CPELP(I,J)                                                 
 1185 CONTINUE                                                                  
      FOP(6)=SUM/(2.*GPI*RADP)                                                  
C                                                                               
C     WITH RESPECT TO EACH OF THE T(I)                                          
C                                                                               
      DO 1195 I=1,KKI                                                           
      SUMN=0.                                                                   
      SUMS=0.                                                                   
      DO 1190 J=1,KATJ                                                          
      SUMN=SUMN+EXP(-.5*T(I+1)**2)*(FUTR(I,J+1)-FUTR(I,J))*(VV(J,I)/            
     1CPELR(I,J)-VV(J,I+1)/CPELR(I+1,J))                                        
      SUMS=SUMS+TEXPX(I+1)*(FYXP(I+1,J+1)-FYXP(I+1,J))*(WW(J,I)/                
     1CPELP(I,J)-WW(J,I+1)/CPELP(I+1,J))                                        
 1190 CONTINUE                                                                  
      FOP(I+6)=(SUMN+BX*SUMS)/SQRT(2.*GPI)                                      
 1195 CONTINUE                                                                  
C                                                                               
C     WITH RESPECT TO EACH OF THE U(J)                                          
C                                                                               
      DO 1205 J=1,KKJ                                                           
      SUMN=0.                                                                   
      SUMS=0.                                                                   
      DO 1200 I=1,KATI                                                          
      SUMN=SUMN+EXP(-.5*U(J+1)**2)*(FTUR(I+1,J)-FTUR(I,J))*(VV(J,I)/            
     1CPELR(I,J)-VV(J+1,I)/CPELR(I,J+1))                                        
      SUMS=SUMS+TEXPY(J+1)*(FXYP(I+1,J+1)-FXYP(I,J+1))*(WW(J,I)/                
     1CPELP(I,J)-WW(J+1,I)/CPELP(I,J+1))                                        
 1200 CONTINUE                                                                  
      FOP(J+KKI+6)=(SUMN+BY*SUMS)/SQRT(2.*GPI)                                  
 1205 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
