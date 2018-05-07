C--------------------------------------------------------------------
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
      COMMON/BLK5/FPVAL(26),TPVALX(26),TPVALY(26)                               
      COMMON/BLK6/ICLASS,JCLASS,IERX,IERY,ICON,JCON,IARY(11),JARY(11)           
      COMMON/BLK7/AX,BX,AY,BY,R,RHO,T(11),U(11)                                 
      common/io/nin,nout,nnew
      INTEGER TCODE                                                             
      REAL XI(11),XJ(11)                                                        
      iend=0
 5    IERX=0                                                                    
      IERY=0                                                                    
      ICLASS=0                                                                  
      JCLASS=0                                                                  
      call iofile(iend)
      IF (iend.NE.0) GO TO 100
      CALL READIN
      CALL HEADLN
      CALL CLAPSE                                                               
      CALL TPFFPF                                                               
C                                                                               
C     CLASSIFY DATA SET DEGENERACY, IF ANY.  IF DATA ARE NOT DEGENERATE,        
C     GET INITIAL ESTIMATES OF THE AX,BX, T(I), AY, BY, AND U(J)                
C     PARAMETERS FOR THE BIVARIATE CALCULATION BY RUNNING SUBROUTINE            
C     'DORALF' ON THE MARGINAL DATA SETS.                                       
C                                                                               
      CALL CHECK(ICLASS,KATI,FPFX,TPFX,ICON,IARY)                               
      CALL CHECK(JCLASS,KATJ,FPFY,TPFY,JCON,JARY)                               
      KEYDEG=0                                                                  
      IF(ICLASS.NE.0.OR.JCLASS.NE.0)CALL DEGMSG(KEYDEG)                         
      IF(KEYDEG.EQ.1)then
	if(nout.ne.6)close(nout,status='delete')
	go to 5
      endif
      CALL ESTM1(AX,BX,KATI,FPFX,TPFX,XI,EMN,EMS)                               
      CALL DORALF(AX,BX,KATI,EMN,EMS,SUMNT,SUMST,XI,VAX,VBX,VABX,IERX)          
      IF (IERX.ne.0)then
	CALL ITRMSG
	if(nout.ne.6)close(nout,status='delete')
	GO TO 5
      endif
      CALL ESTM1(AY,BY,KATJ,FPFY,TPFY,XJ,EMN,EMS)                               
      CALL DORALF(AY,BY,KATJ,EMN,EMS,SUMNU,SUMSU,XJ,VAY,VBY,VABY,IERY)          
      IF (IERY.ne.0)then
	CALL ITRMSG
	if(nout.ne.6)close(nout,status='delete')
	GO TO 5
      endif
C                                                                               
C     THIS PART CALCULATES MAXIMUM LIKELIHOOD ESTIMATES OF THE PARAMETERS       
C     OF THE BIVARIATE MODEL FROM CORRELATED PAIRS OF RATINGS.                  
C     INITIAL ESTIMATES OF AX, BX, AY, BY AND THE T'S AND U'S WERE              
C     OBTAINED ABOVE BY USING SUBROUTINE 'DORALF' ON THE MARGINAL DATA          
C     SETS.  INITIAL ESTIMATES OF THE CORRELATION COEFFICIENTS OF THE           
C     UNDERLYING BIVARIATE NOISE AND SIGNAL-PLUS-NOISE DISTRIBUTIONS            
C     WILL BE OBTAINED BY CALCULATING A PRODUCT MOMENT CORRELATION              
C     COEFFICIENT FROM THE CORRESPONDING RATING-DATA MATRIX.                    
C                                                                               
20    ITER=0                                                                    
      LSTAT=0                                                                   
C                                                                               
C     THE LSTAT VARIABLE IS A DUMMY VARIABLE FOR THE STATUS OF THE              
C     ITERATIVE PROCEDURE.  0-INCOMPLETE , 1-COMPLETE                           
C                                                                               
      T(1)=0.                                                                   
      U(1)=0.                                                                   
      T(NNI)=0.                                                                 
      U(NNJ)=0.                                                                 
      DO 30 I=1,KKI                                                             
      T(I+1)=-XI(I)                                                             
30    CONTINUE                                                                  
      DO 40 I=1,KKJ                                                             
      U(I+1)=-XJ(I)                                                             
40    CONTINUE                                                                  
      CALL CORCOE                                                               
      CALL OUTINI(VAX,VBX,VABX,VAY,VBY,VABY)                                    
      CALL INIVAR                                                               
50    ITER=ITER+1                                                               
      CALL TERMS                                                                
      if(iter.eq.1)write(*,*)'--- program is running; please wait ---'
      CALL FIRST                                                                
      CALL SECOND                                                               
      call secnd2
      call secnd3
      call secnd4
      call secnd5
      call secnd6
      CALL ITERAT                                                               
      IF (ITER.GT.100) then
	write(*,95) CRIT
95	FORMAT(1X,'PROCEDURE DOES NOT CONVERGE AFTER 100 ITERATIONS'/
     1            'ON THE LAST ITERATION CRIT=',F7.5///)
	if(nout.ne.6)close(nout,status='delete')
	GO TO 5
      endif
      IF (LSTAT.NE.2) GO TO 50                                                  
C                                                                               
C     GET CORRELATION MATRIX ON FINAL ITERATION                                 
C                                                                               
      DO 85 I=1,MTRX                                                            
      DO 85 J=1,MTRX                                                            
      CORR(I,J)=VCOV(I,J)/SQRT(VCOV(I,I)*VCOV(J,J))                             
85    CONTINUE                                                                  
      CALL TPFVAL                                                               
      CALL TEST
      CALL OUTRSL                                                               
      CALL INFORM
      if(nout.ne.6)close(nout,status='keep')
      GO TO 5                                                                   
100   stop
      end
c---------------------------------------
      subroutine iofile(iend)
c---------------------------------------
c
c     assign input/output files
c
      COMMON/PASS/LSTRING(80),LINE(80),LENGTH
      character*12 infname,outfname,newfname
      logical*1 ex
      common/io/nin,nout,nnew
      nnew=0
      nin=0
      nout=0
c
c     information about input file
c
      write(*,50)
   50 format(1h1)
  100 write(*,*)'Do you want to use data from a previously created ', 
     +          'input file for the next run?'
      write(*,*)'(Y/N, or Q to quit)'
      read(5,101)(line(i),i=1,80)
  101 FORMAT(80A1)
      call getwrd(ierror)
      if(ierror.ne.0)then
	write(*,*)'   --- no data entered ---' 
	go to 100
      endif
c
c     input file exists
c
      if(lstring(1).eq.'Y' .or. lstring(1).eq.'y')then
  110	write(*,*)'enter INPUT file name:'
	read(5,120)infname
  120	format(a)
        ex=.false.
	inquire(file=infname,exist=ex)
	if(ex .eqv. .false.)then
	  write(*,125)infname
  125     format(1x,'The input file (',a,') CANNOT be found in your ',
     +              'current subdirectory!')
	  go to 110
	endif
        nin=2
	open(nin,file=infname,status='old')
c
c     enter from keyboard
c
      else if(lstring(1).eq.'N' .or. lstring(1).eq.'n')then
c
c     ask for saving input file
c
  130   write(*,*)'Do you want to create a file to store the ', 
     +            'data you will input? (Y/N)'
        read(*,101)(line(i),i=1,80)
	call getwrd(ierror)
	if(ierror.ne.0)then
	  write(*,*)'   --- no data entered ---' 
	  go to 130
	endif
        if(lstring(1).eq.'Y' .or. lstring(1).eq.'y')then
c
c     save the input file
c
  140  	  write(*,*)'Assign a name for this input data file:'
	  write(*,150)
  150	  format(2x,'(Filename must take the form FILENAME.',
     +              'EXTENSION, where FILENAME is 1 to 8'/
     +              '  characters, EXTENSION is 1 to 3 ',
     +              'characters, and all characters must be'/
     +              '  alphanumeric. No blanks are allowed ',
     +              'in the filename.)')
	  read(5,120)newfname
	  ex=.false.
	  inquire(file=newfname,exist=ex)
	  if(ex .eqv. .true.)then
  160	    write(*,230)newfname
            read(5,101)(line(i),i=1,80)
            call getwrd(ierror)
	    if(ierror.ne.0)then
	      write(*,*)'   --- no data entered ---'
	      go to 160
	    endif
            if(lstring(1).eq.'Y' .or. lstring(1).eq.'y')then
c
c      scratch the old input file
c
	      open(11,file=newfname,status='old')
	      close(11,status='delete')
	    else if(lstring(1).eq.'N' .or. lstring(1).eq.'n')then
              go to 140
	    else
	      write(*,*)'   --- invalid input ---'
	      go to 160
	    endif
          endif 
c
c      create a file for saving input data
c
	  nnew=12
	  open(nnew,file=newfname,status='new')
	else if(lstring(1).ne.'N' .and. lstring(1).ne.'n')then
	  write(*,*)'   --- invalid input ---'
	  go to 130
	endif
	nin=5
      else if(lstring(1).eq.'Q' .or. lstring(1).eq.'q')then
c
c     quit the program
c
        iend=1
	return
      else
	write(*,*)'   --- invalid input ---'
	go to 100
      endif
c
c     information for output file
c
  200 write(*,*)'Do you want to save an output file for plotting ',
     +          'the ROC curve? (Y/N)'
      read(5,101)(line(i),i=1,80)
      call getwrd(ierror)
      if(ierror.ne.0)then
	write(*,*)'   --- no data entered ---'
	go to 200
      endif
c
c     save ploting output file
c
      if(lstring(1).eq.'Y' .or. lstring(1).eq.'y')then
  210	write(*,*)'enter OUTPUT file name:'
	write(*,150)
	read(5,120)outfname
	ex=.false.
	inquire(file=outfname,exist=ex)
c
c	if output file name exists 
c
	if(ex .eqv. .true.)then
  220      write(*,230)outfname
  230      format(1x,'The file name ( ',a,') already exists;',
     +     ' do you want to OVERWRITE it?'/1x,
     +     '(Y/N).  If you enter ''Y'', the old file will be ',
     +     'overwritten and lost!!!')
           read(5,101)(line(i),i=1,80)
           call getwrd(ierror)
           if(ierror.ne.0)then
             write(*,*)'   --- no data entered ---' 
             go to 220
           endif
c
c	scratch old file
c
           if(lstring(1).eq.'Y' .or. lstring(1).eq.'y')then
             open(11,file=outfname,status='old')
             close(11,status='delete')
c
c	allow to reenter a new output file name
c
           else if(lstring(1).eq.'N' .or. lstring(1).eq.'n')then
	     go to 210
           else
             write(*,*)'   --- invalid input ---'
	     go to 220
	   endif
	endif
c
c	output file name does not exist, creat a file with the
c	assigned name
c
	nout=3
	open(nout,file=outfname,status='new')
c
c     no ploting output file requested
c
      else if(lstring(1).eq.'N' .or. lstring(1).eq.'n')then
	nout=6
      else
	write(*,*)'   --- invalid input ---'
	go to 200
      endif
      return
      end
C-----------------------------------------------
      BLOCK DATA                                                                
C-----------------------------------------------
      COMMON/BLK5/FPVAL(26),TPVALX(26),TPVALY(26)                               
      DATA   FPVAL/0.005,0.01,0.02,0.03,0.04,0.05,                              
     1             0.06,0.07,0.08,0.09,0.10,0.11,                               
     2             0.12,0.13,0.14,0.15,0.20,0.25,                               
     3             0.30,0.40,0.50,0.60,0.70,0.80,                               
     4             0.90,0.95/                                                   
      END                                                                       
C-----------------------------------------------
      SUBROUTINE READIN
C-----------------------------------------------
C
C     READ INPUT DATA
C
      COMMON/BLK1/VV(11,11),WW(11,11),KATI,KATJ,TCODE,FP,NNI,NNJ,KKI,           
     1  KSNX,KSNY,KKJ,KKIL,KKJL,MTRX,IDENTI(60),IDENTJ(60)
      COMMON/BLK2/EMN,EMS,SUMNT(11),SUMST(11),SUMNU(11),SUMSU(11),              
     1            FPFX(11),TPFX(11),FPFY(11),TPFY(11)                           
      common/io/nin,nout,nnew
      INTEGER TCODE
C
C	READ DATA DESCRIPTION, TOTAL NUMBER OF CATEGORY, AND 
C       THE CATEGORY REPRESENTING THE STRONGEST EVIDENCE OF 
C       SIGNAL
C
        if(nin.eq.5)then
	  write(*,70)
 70	  FORMAT(1H1,26X,'C O R R O C 2 :'//,
     +    5X,'M A X I M U M   L I K E L I H O O D   E S T I '           
     1    ,'M A T I O N'/19X,'O F   T H E   P A R A M E T E R S'/5X,
     2    'O F   T H E   B I V A R I A T E   B I N O R M A L   ',
     3    'M O D E L'/9X,'F O R   C O R R E L A T E D   R A T I N G   ',
     4    'D A T A'/)
          write(*,*)
	  write(*,*)'For the first condition (''condition X'') --'
        endif
	CALL READHL(IDENTI,KATI,KSNX)
	if(nin.eq.5)then
          write(*,*)
	  write(*,*)'For the second condition (''condition Y''): --'
	endif
	CALL READHL(IDENTJ,KATJ,KSNY)
C
C       READ RATING-DATA PAIRS AND GENERATE RATING-DATA MATRICES
C       FOR NOISE AND SIGNAL-PLUS-NOISE TRIALS
C
        if(nin.eq.5)then
          write(*,*)
	  write(*,*)'For each actually NEGATIVE trial, enter on a ',
     +              'single line: the rating'
          write(*,75)
  75      format(1x,'from condition X, one or more blanks, and the ',
     +              'rating from condition Y.')
          write(*,*)'After all actually NEGATIVE trials have been ',
     +              'input, enter an asterisk (*)'
          write(*,76)
  76      format(1x,'in response to the next prompt.')
        endif
	CALL READQ(1,MN,KATI,KATJ,VV)
        if(nin.eq.5)then
	  if(nnew.eq.12)write(nnew,77)
  77      format(2x,'*')
	  write(*,*)'*** end of actually negative trials ***'
          write(*,*)
	  write(*,*)'For each actually POSITIVE trial, enter on a ', 
     +              'single line: the rating'
          write(*,75)
	  write(*,*)'After all actually POSITIVE trials have been ',
     +              'input, enter an asterisk (*)'
          write(*,76)
        endif
	CALL READQ(2,MS,KATI,KATJ,WW)
	if(nin.eq.5)write(*,*)'*** end of actually positive trials ***'
	if(nnew.eq.12)write(nnew,77)
80	NNI=KATI+1
	NNJ=KATJ+1
C
C	READ "TCODE"
C
	CALL READT(TCODE,FP)
199     RETURN
	END
C--------------------------------------------
      SUBROUTINE READHL(ITITLE,KAT,KSN)
C--------------------------------------------
C
C     THIS SUBROUTINE READS IN A FREE-TEXT DESCRIPTION OF THE DATA,
C     TOTAL NUMBER OF CATEGORY AND A CODE WHICH INDICATES WHETHER 
C     CATEGORY 1 OR CATEGORY KAT IS ASSOCIATED WITH ACTUALLY 
C     POSITIVE TRIALS.
C                                                                               
      COMMON/PASS/LSTRING(80),LINE(80),LENGTH
      common/io/nin,nout,nnew
      DIMENSION ITITLE(60)
      DATA IBLANK/1H /
C
      KAT=0
      KSN=0
      if(nin.eq.5)write(*,*)'Enter data description (up to 60 ',
     +                      'characters):'
      read(nin,100)(LINE(I),I=1,80)
 100  FORMAT(80A1)
      DO 120 I=1,80
	IF(LINE(I).EQ.IBLANK)GO TO 120
	IPNT=I
	GO TO 130
 120	CONTINUE
      IPNT=1
 130  DO 150 J=IPNT,60
        ITITLE(J-IPNT+1)=LINE(J)
 150	CONTINUE
      if(nnew.eq.12)write(nnew,155)(ititle(i),i=1,60)
 155  format(60a1)
C
 200  if(nin.eq.5)write(*,*)'Enter total number of categories ',
     +                      '(3 to 10):'
      read(nin,100)(LINE(I),I=1,80)
      CALL GETWRD(IERROR)
      IF(IERROR.ne.0)then
	write(*,*)'   --- no data entered ---'
	go to 200
      endif
      CALL TONUM(LENGTH,LSTRING,RVALUE,IERROR)
      IF(IERROR.ne.0)then
	write(*,*)'   --- invalid form of input data ---'
	go to 200
      endif
      DO 202 I=3,10
         IF(RVALUE.EQ.I)GO TO 203
 202     CONTINUE
      write(*,*)'   --- invalid form of input data ---'
      go to 200
 203  call getwrd(notget)
      if(notget.eq.0)then
        write(*,*)'   --- invalid total number of input data ---'
	go to 200
      endif
      KAT=RVALUE
      if(nnew.eq.12)write(nnew,210)kat
 210  format(2x,i6)
C
 300  if(nin.eq.5)then
        write(*,*)'Input the category number that ',
     +            'represents the strongest evidence of positivity'
        write(*,*)'(e.g., that abnormality is present):'
      endif
      read(nin,100)(LINE(I),I=1,80)
      CALL GETWRD(IERROR)
      IF(IERROR.ne.0)then
	write(*,*)'   --- no data entered ---'
	go to 300
      endif
      CALL TONUM(LENGTH,LSTRING,RVALUE,IERROR)
      IF(IERROR.ne.0 .or. (RVALUE.NE.1.AND.RVALUE.NE.KAT))then
	write(*,*)'   --- invalid form of input data ---'
	go to 300
      endif
      call getwrd(notget)
      if(notget.eq.0)then
        write(*,*)'   --- invalid total number of input data ---'
	go to 300
      endif
      KSN=RVALUE
      if(nnew.eq.12)write(nnew,210)ksn
      RETURN
      END
C--------------------------------------------
      SUBROUTINE READQ(kkk,NUM,MAXX,MAXY,QQ)
C--------------------------------------------
C                                                                               
C     THIS SUBROUTINE READS IN A SEQUENCE OF INPUT DATA PAIR IN FREE 
C     FORMAT. THE ONLY FORMAT REQUIREMENTS ARE THAT  (1)
C     ANY TWO INPUT VALUES MUST BE SEPARATED BY AT LEAST ONE                    
C     IBLANK COLUMN AND (2) THE INPUT DATA SEQUENCE MUST BE TERMINATED 
C     BY AN NON NUMERICAL CHARACTER
C                                                                               
      COMMON/PASS/LSTRING(80),LINE(80),LENGTH
      common/io/nin,nout,nnew
      DIMENSION QQ(11,11)
      DATA IBLANK/1H /,ISTAR/1H*/
C                                                                               
      NUM=1
      DO 50 I=1,11
      DO 50 J=1,11
         QQ(I,J)=0
  50     CONTINUE
C
 100  if(nin.eq.5)then
      if(kkk.eq.1)then
	write(*,105)num
 105	format(2x,'actually negative trial #',i3,': ')
      else
	write(*,106)num
 106	format(2x,'actually positive trial #',i3,': ')
      endif
      endif
      read(nin,110)(LINE(I),I=1,80)
 110  FORMAT(80A1)
C
C     READ RATING-DATA FOR X-CONDITION
C
      CALL GETWRD(IERROR)
      IF(IERROR.ne.0)then
	write(*,*)'   --- no data entered ---'
	go to 100
      endif
      IF(LSTRING(1).EQ.ISTAR)then
	num=num-1
	return
      endif
      CALL TONUM(LENGTH,LSTRING,RVALUE,IERROR)
      IF(IERROR.ne.0)then
	write(*,*)'   --- invalid form of input data ---'
	go to 100
      endif
      DO 121 I=1,MAXX
         IF(RVALUE.EQ.I)GO TO 122
 121     CONTINUE
      write(*,*)'   --- invalid form of input data ---'
      go to 100
 122  IX=RVALUE
C
C     READ RATING-DATA FOR Y-CONDITION
C
      CALL GETWRD(IERROR)
      IF(IERROR.ne.0)then
	write(*,*)'   --- invalid total number of input data ---'
	go to 100
      endif
      CALL TONUM(LENGTH,LSTRING,RVALUE,IERROR)
      IF(IERROR.ne.0)then
	write(*,*)'   --- invalid form of input data ---'
	go to 100
      endif
      DO 123 I=1,MAXY
         IF(RVALUE.EQ.I)GO TO 124
 123     CONTINUE
      write(*,*)'   --- invalid form of input data ---'
      go to 100
 124  call getwrd(notget)
      if(notget.eq.0)then
        write(*,*)'   --- invalid total number of input data ---'
	go to 100
      endif
      IY=RVALUE
      QQ(IY,IX)=QQ(IY,IX)+1
      if(nnew.eq.12)write(nnew,130)ix,iy
 130  format(2x,i2,2x,i2)
      num=num+1
      GO TO 100
      END                                                                       
C--------------------------------------------
	SUBROUTINE READT(ICODE,FP)
C--------------------------------------------
C
C     READ STATISTICAL TEST CODE
C
      COMMON/PASS/LSTRING(80),LINE(80),LENGTH
      common/io/nin,nout,nnew
      DATA IB1/1HB/,IB2/1Hb/,IA1/1HA/,IA2/1Ha/,IT1/1HT/,IT2/1Ht/
C
  50  if(nin.eq.5)write(*,*)'Enter statistical test desired ',
     +                      '(bivariate/area/tpf):'
      read(nin,100)(LINE(I),I=1,80)
 100  FORMAT(80A1)
      CALL GETWRD(IERROR)
      IF(IERROR.ne.0)then
	write(*,*)'   --- no data entered ---'
	go to 50
      endif
      IF(LSTRING(1).EQ.IB1 .OR. LSTRING(1).EQ.IB2)then
	ICODE=1
      else IF(LSTRING(1).EQ.IA1 .OR. LSTRING(1).EQ.IA2)then
	ICODE=2
      else IF(LSTRING(1).EQ.IT1 .OR. LSTRING(1).EQ.IT2)then
	ICODE=3
      else
	write(*,*)'   --- invalid form of input data ---'
	go to 50
      endif
      if(nnew.eq.12)write(nnew,110)(lstring(i),i=1,20)
 110  format(20a1)
      IF(ICODE.lt.3)then
	if(nin.ne.5)close(nin,status='keep')
	if(nnew.eq.12)close(nnew,status='keep')
	RETURN
      endif
C
C     READ IN FPF VALUE IF ICODE=3
C
 200  if(nin.eq.5)then
        write(*,*)'TPF test is selected; please input the ',
     +          'FPF value (>0 and <1) at which'
        write(*,*)'the TPFs are to be compared:'
      endif
      read(nin,100)(LINE(I),I=1,80)
      CALL GETWRD(IERROR)
      IF(IERROR.ne.0)then
	write(*,*)'   --- no data entered ---'
	go to 200
      endif
      CALL TONUM(LENGTH,LSTRING,RVALUE,IERROR)
      IF(IERROR.ne.0)then
	write(*,*)'   --- invalid form of input data ---'
	go to 200
      endif
      FP=RVALUE
      IF(FP.GE.1.OR.FP.LE.0)then
	write(*,*)'   --- invalid FPF value ---'
	go to 200
      endif
      if(nnew.eq.12)write(nnew,210)fp
 210  format(2x,f7.4)
      if(nnew.eq.12)close(nnew,status='keep')
      if(nin.ne.5)close(nin,status='keep')
      RETURN
      END
C--------------------------------------------
      SUBROUTINE GETWRD(NOTGET)
C--------------------------------------------
C
C     GET A STRING FROM THE INPUT LINE
C
      COMMON /PASS/LSTRING(80),LINE(80),LENGTH
      DATA IBLANK/1H /
C
      NOTGET=0
      LENGTH=0
      I=1
      DO 10 J=1,80
         LSTRING(J)=IBLANK
10       CONTINUE
C
20    IF(I .GT. 80) GO TO 90
      IF(LINE(I) .NE. IBLANK) GO TO 40
      I= I+1
      GO TO 20
C
40    II= 1
50    CONTINUE
      LSTRING(II)= LINE(I)
      LENGTH=II
      I= I+1
      IF(I .GT. 80) GO TO 60
      IF(LINE(I).EQ.IBLANK) GO TO 70
      II= II+1
      GO TO 50
C
60    DO 65 J=1,80
         LINE(J)=IBLANK
65    CONTINUE
      RETURN
C
70    DO 80 J=I,80
         LINE(J-I+1)=LINE(J)
80       CONTINUE
      INIT=80-I+2
      DO 85 J=INIT,80
         LINE(J)=IBLANK
85       CONTINUE
      RETURN
C
90    NOTGET=1
      RETURN
      END
C-----------------------------------------------------------
      SUBROUTINE TONUM(LEN,LINPUT,RNUM,IERR)
C-----------------------------------------------------------
C
C     CONVERT CHARACTER STRING TO REAL NUMBER
C
      DIMENSION LINPUT(80)
      DATA IZERO/1H0/,IDEC/1H./
C
      IERR= 0
      IF(LEN .EQ. 0) GO TO 120
C
C     CHECK INVALID CHARACTER (I.E., ANY CHARACTER EXCEPT 0-9, AND PERIOD)
C
      IP= 0
      DO 20 I= 1,LEN
      IK= IORD(LINPUT(I)) - IORD(IZERO)
      IF(IK .GE. 0 .AND. IK .LE. 9) GO TO 20
      IF(LINPUT(I) .NE. IDEC) GO TO 120
      IP= IP+1
      IF(IP .GT. 1 .OR. (IP.EQ.1.AND.LEN.EQ.1)) GO TO 120
20    CONTINUE
C
C     CONVERT TO NUMERICAL STRING
C
      IU= 0
      TOTAL= 0
      K= 1
40    IF(LINPUT(K) .EQ. IDEC) GO TO 70
      IK=IORD(LINPUT(K))-IORD(IZERO)
      IF(IK.LT.0 .OR. IK.GT.9)GO TO 120
      TOTAL= TOTAL*10 + IK
      K= K+1
      IF(K.GT.LEN) GO TO 100
      GO TO 40
C
70    K= K+1
      IF(K.GT.LEN) GO TO 100      
      IK=IORD(LINPUT(K))-IORD(IZERO)
      IF(IK.LT.0 .OR. IK.GT.9)GO TO 120
      TOTAL= TOTAL*10 + IK
      IU= IU+1
      GO TO 70
C
100   RNUM = TOTAL/(10.0**IU)
      RETURN
C
120   IERR= 1
200   RETURN
      END
C----------------------------------------------
      INTEGER FUNCTION IORD(LCH)
C----------------------------------------------
C
C      CHANGE FROM CHARACTER TO ASCII CODE
C
      DIMENSION LSTD(95)
      DATA LSTD/1H ,1H!,1H",1H#,1H$,1H%,1H&,1H',1H(,1H),1H*,1H+,1H-,
     +	        1H,,1H.,1H/,1H0,1H1,1H2,1H3,1H4,1H5,1H6,1H7,1H8,1H9,
     +          1H:,1H;,1H<,1H=,1H>,1H?,1H@,1HA,1HB,1HC,1HD,1HE,1HF,
     +          1HG,1HH,1HI,1HJ,1HK,1HL,1HM,1HN,1HO,1HP,1HQ,1HR,1HS,
     +          1HT,1HU,1HV,1HW,1HX,1HY,1HZ,1H[,1H\,1H],1H^,1H_,1H`,
     +          1Ha,1Hb,1Hc,1Hd,1He,1Hf,1Hg,1Hh,1Hi,1Hj,1Hk,1Hl,1Hm,
     +          1Hn,1Ho,1Hp,1Hq,1Hr,1Hs,1Ht,1Hu,1Hv,1Hw,1Hx,1Hy,1Hz,
     +          1H{,1H|,1H},1H~/
      DO 20 I= 1,95
      IF(LSTD(I).NE.LCH) GO TO 20
	IORD=I
	RETURN
20    CONTINUE
50    write(*,60)
60    FORMAT(2X,'--- CHARACTER NOT FOUND ---')
      RETURN
      END
