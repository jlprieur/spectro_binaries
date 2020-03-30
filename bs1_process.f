      SUBROUTINEBS1_PROCESS(OBJECT,PERIOD,SIGMA)
      IMPLICITNONE
      INTEGERIDIM
      PARAMETER(IDIM=120)
      REALPOI(IDIM),TOBS(IDIM),VIT(IDIM),RES(IDIM)
      REALERR(6),EQM(6),CORR(6),W1,W2
      INTEGERLU_IN,LU_OUT,LU_OUT_TXT,NP,X1,Y1
      INTEGERNORB,NN,NITER,I,IOPT,STATUS
      REALEE,OMEGA,KA,VZERO,PI,PERIOD,SIGMA
      REALTZERO,ASINI,FDEM,ERR_ASINI,ERR_FDEM
      DOUBLEPRECISIONPP,T
      CHARACTEROBJECT*40,OUT_RESULTS*60,OUT_LATEX*60
      CHARACTERBUFFER*80
      LOGICALPUBLI,INTERACTIVE
      PI=3.1415926535
   10 FORMAT(A)
      INTERACTIVE=(OBJECT(1:16).NE.'MODE-AUTOMATIQUE').AND.(OBJECT(1:16)
     +.NE.'mode-automatique')
      LU_IN=7
      OPEN(LU_IN,FILE='VR_BS1.DAT',STATUS='OLD')
      NORB=1
 1000 FORMAT(F10.6,F10.3,2F10.5,F10.2,F6.2,1X,2I1)
      READ(LU_IN,1000)PP,TZERO,OMEGA,EE,KA,VZERO,X1,Y1
      IF(INTERACTIVE)THEN
      WRITE(6,*)' P=',PP,' T0=',TZERO,' Omega=',OMEGA,' E=',EE,' KA=',KA
      WRITE(6,*)' V0=',VZERO
      ENDIF
      T=TZERO
      IF(X1.NE.1)X1=0
      IF(Y1.NE.1.AND.Y1.NE.2)Y1=0
      IF(INTERACTIVE)WRITE(6,*)'Option: X1=',X1,' Y1=',Y1
      IF(Y1.EQ.2)THEN
      EE=0.
      OMEGA=0.
      ENDIF
      IF(.NOT.INTERACTIVE.AND.PERIOD.GT.0.)THEN
      PP=PERIOD
      WRITE(6,*)' OK: Mode de recherche automatique: Période=',PP
      ENDIF
      READ(LU_IN,*)NN,NITER
      IF(INTERACTIVE)WRITE(6,*)'Nombre de mesures:',NN,' Nombre d''itera
     +tions',NITER
      DO80000I=1,NN
 1200 FORMAT(F9.3,F7.1,F6.2)
      READ(LU_IN,10)BUFFER
      POI(I)=-1.0
      READ(BUFFER,*)TOBS(I),VIT(I)
      READ(BUFFER,*,ERR=181,END=181)W1,W2,POI(I)
  181 IF(POI(I).EQ.-1.0)POI(I)=1.0
80000 CONTINUE
      CLOSE(LU_IN)
      LU_OUT_TXT=8
      OUT_RESULTS=OBJECT(1:INDEX(OBJECT,' ')-1)//'.txt'
      IF(INTERACTIVE)OPEN(LU_OUT_TXT,FILE=OUT_RESULTS,STATUS='UNKNOWN')
      NP=6-X1-Y1
      IF(INTERACTIVE)WRITE(LU_OUT_TXT,1500)PP,T,OMEGA,EE,KA,VZERO
 1500 FORMAT(13X,'ELEMENTS PROVISOIRES',/,16X,'P =',F15.7,/,16X,'T =',F1
     +1.3,/,12X,'OMEGA =',F13.5,/,16X,'E =',F13.5,/,16X,'K =',F10.2,/,15
     +X,'VO =',F10.2,/,/)
      CALLBS1_FIT(POI,TOBS,VIT,RES,ERR,EQM,CORR,PP,T,TZERO,OMEGA,EE,KA,V
     +ZERO,ASINI,FDEM,ERR_ASINI,ERR_FDEM,LU_OUT_TXT,NP,X1,Y1,NN,NITER,IN
     +TERACTIVE,STATUS)
      IF(STATUS.NE.0)THEN
      SIGMA=1000.
      IF(INTERACTIVE)THEN
      WRITE(6,*)"BS1_PROCESS/Errorwhileinvertingmatrix"
      WRITE(LU_OUT_TXT,*)"BS1_PROCESS/Errorwhileinvertingmatrix"
      ENDIF
      CLOSE(LU_OUT_TXT)
      RETURN
      ENDIF
      CALLRMS_RESIDUS(RES,POI,NN,SIGMA)
      IF(INTERACTIVE)THEN
      WRITE(6,3000)SIGMA
      WRITE(LU_OUT_TXT,3000)SIGMA
      ENDIF
 3000 FORMAT('Ecart-type moyen des residus (Sigma O-C): ',F5.2,' km s$^{
     +-1}$',/)
      CLOSE(LU_OUT_TXT)
      IF(INTERACTIVE)THEN
      WRITE(6,*)' Courbe: 0=aucune  1=orbite  2=residus '
      WRITE(6,*)'         3=les deux  4=mode pour publi. ?'
      READ(5,*)IOPT
      PUBLI=(IOPT.EQ.4)
      IF(IOPT.EQ.1)THEN
      CALLJLP_PLOT_BS1(TOBS,VIT,POI,T,PP,NN,EE,OMEGA,KA,VZERO,OBJECT,PUB
     +LI)
      ELSEIF(IOPT.EQ.2)THEN
      CALLJLP_PLOT_RES1(TOBS,RES,POI,NN,OBJECT,PUBLI)
      ELSEIF(IOPT.EQ.3.OR.IOPT.EQ.4)THEN
      CALLJLP_PLOT_BS1(TOBS,VIT,POI,T,PP,NN,EE,OMEGA,KA,VZERO,OBJECT,PUB
     +LI)
      CALLJLP_PLOT_RES1(TOBS,RES,POI,NN,OBJECT,PUBLI)
      ENDIF
      OUT_LATEX=OBJECT(1:INDEX(OBJECT,' ')-1)//'.tex'
      CALLSB1_TO_LATEX(PP,T,OMEGA,EE,KA,VZERO,ASINI,FDEM,SIGMA,ERR,ERR_A
     +SINI,ERR_FDEM,TOBS,VIT,RES,POI,NN,OUT_LATEX,OBJECT)
      WRITE(6,1300)CORR(1),PP,ERR(1),CORR(2),T,ERR(2),CORR(3)*180/PI,OME
     +GA*180/PI,ERR(3)*180/PI,CORR(3),OMEGA,ERR(3),CORR(4),EE,ERR(4),COR
     +R(5),KA,ERR(5),CORR(6),VZERO,ERR(6),ASINI,ERR_ASINI,FDEM,ERR_FDEM
 1300 FORMAT(1X,'CORRECTIONS',15X,'NOUVEAUX ELEMENTS',15X,'ERREURS STD.'
     +,/,3X,E9.3,17X,'P =',F15.7,13X,E9.3,/,3X,E9.3,17X,'T =',F11.3,17X,
     +E9.3,/,3X,E9.3,13X,'OMEGA =',F13.5,15X,E9.3,4X,'(deg.)',/,3X,E9.3,
     +13X,'OMEGA =',F13.5,15X,E9.3,4X,'(rad.)',/,3X,E9.3,17X,'E =',F13.5
     +,15X,E9.3,/,3X,E9.3,17X,'K =',F10.2,18X,E9.3,/,3X,E9.3,16X,'V0 =',
     +F10.2,18X,E9.3,/,24X,'A.SINI =',6X,E10.4,12X,E9.3,/,26X,'F(M) =',6
     +X,E10.4,12X,E9.3,/,/)
      LU_OUT=8
      OPEN(LU_OUT,FILE='RESID.DAT',STATUS='UNKNOWN')
      WRITE(LU_OUT,*)NN
      WRITE(LU_OUT,1200)(TOBS(I),RES(I),POI(I),I=1,NN)
      CLOSE(LU_OUT)
      WRITE(6,*)' Fichiers de sortie: ',OUT_RESULTS(1:INDEX(OUT_RESULTS,
     +' ')-1),', ',OUT_LATEX(1:INDEX(OUT_LATEX,' ')-1),', RESID.DAT,'
      IF(IOPT.EQ.3)THEN
      WRITE(6,*)' et graphes dans: ',OBJECT(1:INDEX(OBJECT,' ')-1),'.ps 
     +et ',OBJECT(1:INDEX(OBJECT,' ')-1),'_resi.ps'
      ELSEIF(IOPT.EQ.1)THEN
      WRITE(6,*)' et graphe dans: ',OBJECT(1:INDEX(OBJECT,' ')-1),'.ps'
      ELSEIF(IOPT.EQ.2)THEN
      WRITE(6,*)' et graphe dans: ',OBJECT(1:INDEX(OBJECT,' ')-1),'_resi
     +.ps'
      ENDIF
      ENDIF
      PERIOD=PP
      RETURN
      END
      SUBROUTINEJLP_PLOT_BS1(TOBS,VIT,POI,T,PP,NN,EE,OMEGA,KA,VZERO,OBJE
     +CT,PUBLI)
      IMPLICITNONE
      INTEGERIDIM,KCUR
      PARAMETER(IDIM=1024,KCUR=3)
      REALTOBS(*),VIT(*),POI(*),PHI_MINI
      DOUBLEPRECISIONT,PP,PHI
      REALPI,EE,OMEGA,KA,VZERO,VIC,VV
      INTEGERITER,I,I1,I2,NN,KK,NPTS(KCUR)
      REALXPLOT(IDIM,KCUR),YPLOT(IDIM,KCUR)
      CHARACTERCHAR1*30,CHAR2*30,CHAR3*40,OBJECT*40
      CHARACTERPLOTDEV*40,IN_FILE*40,IN_COMMENTS*80
      CHARACTER*4NCHAR(KCUR)
      LOGICALPUBLI
      PI=3.1415926535
      PHI_MINI=0.1
      I1=0
      I2=0
      DO80001ITER=1,3
      DO80002I=1,NN
      PHI=DMOD((TOBS(I)-T)/PP,1.D0)
      IF(PHI.LT.0)PHI=1.+PHI
      IF(ITER.EQ.1)THEN
      IF(POI(I).EQ.1)THEN
      I1=I1+1
      XPLOT(I1,1)=PHI
      YPLOT(I1,1)=VIT(I)
      ELSE
      I2=I2+1
      XPLOT(I2,2)=PHI
      YPLOT(I2,2)=VIT(I)
      ENDIF
      ELSEIF(ITER.EQ.2)THEN
      IF(PHI.LT.PHI_MINI)THEN
      PHI=PHI+1.
      IF(POI(I).EQ.1)THEN
      I1=I1+1
      XPLOT(I1,1)=PHI
      YPLOT(I1,1)=VIT(I)
      ELSE
      I2=I2+1
      XPLOT(I2,2)=PHI
      YPLOT(I2,2)=VIT(I)
      ENDIF
      ENDIF
      ELSEIF(ITER.EQ.3)THEN
      IF(PHI.GT.(1.-PHI_MINI))THEN
      PHI=PHI-1.
      IF(POI(I).EQ.1)THEN
      I1=I1+1
      XPLOT(I1,1)=PHI
      YPLOT(I1,1)=VIT(I)
      ELSE
      I2=I2+1
      XPLOT(I2,2)=PHI
      YPLOT(I2,2)=VIT(I)
      ENDIF
      ENDIF
      ENDIF
80002 CONTINUE
80001 CONTINUE
      NPTS(1)=I1
      NPTS(2)=I2
      NCHAR(1)='913'
      NCHAR(2)='213'
      DO80003I=1,1000
      PHI=-PHI_MINI+(1+2.*PHI_MINI)*REAL(I)/REAL(1000)
      CALLVITESSE_CALCULEE(PHI,VIC,VV,EE,OMEGA,KA,VZERO)
      XPLOT(I,3)=PHI
      YPLOT(I,3)=VIC
80003 CONTINUE
      NCHAR(3)='L0'
      NPTS(3)=1000
      KK=3
      IF(PUBLI)THEN
      WRITE(PLOTDEV,102)OBJECT(1:INDEX(OBJECT,' ')-1),'.ps'
  102 FORMAT('&LANDSCAPE/',A,A)
      CHAR3=' '
      ELSE
      WRITE(PLOTDEV,103)OBJECT(1:INDEX(OBJECT,' ')-1),'.ps'
  103 FORMAT('&landscape/',A,A)
      CHAR3=OBJECT
      ENDIF
      WRITE(6,*)' Output to: ',PLOTDEV
      CHAR1='Phase'
      CHAR2='Radial velocity (km/s)'
      IN_FILE=' '
      IN_COMMENTS=' '
      CALLNEWPLOT(XPLOT,YPLOT,NPTS,IDIM,KK,CHAR1,CHAR2,CHAR3,NCHAR,PLOTD
     +EV,IN_FILE,IN_COMMENTS)
      RETURN
      END
      SUBROUTINEVITESSE_CALCULEE(PHI,VIC,VV,EE,OMEGA,KA,VZERO)
      IMPLICITNONE
      DOUBLEPRECISIONPHI
      REALPI,M,EMME,ANOM,ANOM_OLD,VV,VIC
      REALEE,OMEGA,KA,VZERO
      INTEGERI
      PI=3.1415926535
      M=2.*PI*PHI
      EMME=AMOD(M,2.*PI)
      ANOM_OLD=EMME
      DO80004I=1,1000
      ANOM=EMME+EE*SIN(ANOM_OLD)
      IF(ABS(ANOM-ANOM_OLD).LT.1.E-6)GOTO12
      ANOM_OLD=ANOM
80004 CONTINUE
   12 VV=2.*ATAN(TAN(ANOM/2.)*SQRT((1+EE)/(1-EE)))
      VIC=VZERO+KA*(EE*COS(OMEGA)+COS(VV+OMEGA))
      RETURN
      END
      SUBROUTINERMS_RESIDUS(RES,POI,NN,SIGMA)
      IMPLICITNONE
      INTEGERI,NN
      REALSUM,SUMSQ
      REALRES(*),POI(*),SIGMA
      SUMSQ=0.
      SUM=0.
      DO80005I=1,NN
      SUMSQ=SUMSQ+POI(I)*RES(I)*RES(I)
      SUM=SUM+POI(I)
80005 CONTINUE
      SIGMA=SQRT(SUMSQ/SUM)
      RETURN
      END
      SUBROUTINEMAT_INV(A,D,NP,INTERACTIVE,STATUS)
      IMPLICITNONE
      INTEGERNP,STATUS
      REALA(6,6),C(6,6),D(6,6),G,H
      INTEGERI1,J1,I,J,K,KK,L,L1
      LOGICALINTERACTIVE
      DO10I=1,NP
      DO10J=1,NP
      D(I,J)=0
      D(I,I)=1.
   10 CONTINUE
      DO116K=1,NP
      DO11I1=1,NP
      DO11J1=1,NP
      C(I1,J1)=A(I1,J1)
   11 CONTINUE
      IF(A(K,K).EQ.0)THEN
      if(INTERACTIVE)write(6,*)' Pivot null pour K=',K
      if(INTERACTIVE)write(6,*)' (JLP: j''ai des doutes sur la suite...)
     +'
      KK=K+1
      L1=-1
      DO12L=KK,NP
      IF(A(L,K).NE.0)THEN
      L1=L
      GOTO13
      ENDIF
   12 CONTINUE
      IF(L1.EQ.-1)THEN
      STATUS=-1
      RETURN
      ENDIF
   13 DO14J1=1,NP
      C(K,J1)=A(L1,J1)
      C(L1,J1)=A(K,J1)
      A(K,J1)=C(K,J1)
      A(L1,J1)=C(L1,J1)
      G=D(L1,J1)
      H=D(K,J1)
      D(K,J1)=G
      D(L1,J1)=H
   14 CONTINUE
      ENDIF
      DO61I=1,NP
      IF(A(I,K).NE.0.AND.I.NE.K)THEN
      DO16J=1,NP
      A(I,J)=A(I,J)/C(I,K)-A(K,J)/C(K,K)
      D(I,J)=D(I,J)/C(I,K)-D(K,J)/C(K,K)
   16 CONTINUE
      ENDIF
   61 CONTINUE
  116 CONTINUE
      DO17I=1,NP
      DO17J=1,NP
      D(I,J)=D(I,J)/A(I,I)
   17 CONTINUE
      STATUS=0
      RETURN
      END
      SUBROUTINEJLP_PLOT_RES1(TOBS,RES1,POI1,NN,OBJECT,PUBLI)
      IMPLICITNONE
      INTEGERDIM,KCUR
      PARAMETER(DIM=1024,KCUR=1)
      REALTOBS(120),RES1(120),POI1(120)
      INTEGERKK,NPTS(KCUR),I,NN,I1
      REALXPLOT(DIM,KCUR),YPLOT(DIM,KCUR)
      CHARACTERCHAR1*30,CHAR2*30,CHAR3*40,OBJECT*40
      CHARACTERPLOTDEV*40,IN_FILE*40,IN_COMMENTS*80
      CHARACTER*4NCHAR(KCUR)
      LOGICALPUBLI
      I1=1
      DO80006I=1,NN
      IF(POI1(I).NE.0.)THEN
      XPLOT(I1,1)=TOBS(I)
      YPLOT(I1,1)=RES1(I)
      I1=I1+1
      ENDIF
80006 CONTINUE
      NCHAR(1)='913'
      NPTS(1)=I1-1
      KK=1
      IF(PUBLI)THEN
      WRITE(PLOTDEV,102)OBJECT(1:INDEX(OBJECT,' ')-1),'_resi.ps'
  102 FORMAT('&LANDSCAPE/',A,A)
      CHAR3=' '
      ELSE
      WRITE(PLOTDEV,103)OBJECT(1:INDEX(OBJECT,' ')-1),'_resi.ps'
  103 FORMAT('&landscape/',A,A)
      CHAR3=OBJECT
      ENDIF
      WRITE(6,*)' Output to: ',PLOTDEV
      CHAR1='Date'
      CHAR2='Residuals (km/s)'
      IN_FILE=' '
      IN_COMMENTS=' '
      CALLNEWPLOT(XPLOT,YPLOT,NPTS,DIM,KK,CHAR1,CHAR2,CHAR3,NCHAR,PLOTDE
     +V,IN_FILE,IN_COMMENTS)
      RETURN
      END
      SUBROUTINESB1_TO_LATEX(PP,T,OMEGA,EE,KA,VZERO,ASINI,FDEM,SIGMA,ERR
     +,ERR_ASINI,ERR_FDEM,TOBS,VIT,RES,POI,NN,OUT_LATEX,OBJECT)
      IMPLICITNONE
      DOUBLEPRECISIONPP,T,PHI
      INTEGERNN
      REALPOI(NN),TOBS(NN),VIT(NN),RES(NN)
      REALOMEGA,EE,KA,VZERO,ASINI,FDEM,SIGMA,ERR(6)
      REALERR_ASINI,ERR_FDEM
      CHARACTEROUT_LATEX*40,OBJECT*40
      INTEGERI,LU_OUT
      LU_OUT=8
      OPEN(LU_OUT,FILE=OUT_LATEX,STATUS='UNKNOWN')
      WRITE(LU_OUT,3000)OBJECT(1:INDEX(OBJECT,' ')-1)
 3000 FORMAT('\\documentclass{article}',/,'\\usepackage{graphicx}',/,'\\
     +voffset=-1cm',/,'\\hoffset=-3cm',/,'\\textwidth=15.6cm',/,'\\texth
     +eight=26cm',/,'\\newcommand{\\nodata}{\\ldots}',/,/,'\\begin{docum
     +ent}',/,/,'\\centerline{\\large \\bf ',A,'}',/,/,'\\vskip 1cm',/,'
     +\\tabcolsep=1mm',/,'\\begin{tabular}{lcccccccccccc}',/,'\\hline',/
     +,'Name & $P$ & $T_0$ (JD)& $\\omega$ & $e$ &',' $K_1$ & $K_2$ & $V
     +_0$ ','& $a_1 \\sin i$ & $a_2 \\sin i$',/,'& $f(m_1)$ & $f(m_2)$ &
     + $\\sigma_{(O-C)}$ \\\\',/,'& days & 2400000+ & deg. & & km s$^{-1
     +}$ ','& km s$^{-1}$ & km s$^{-1}$','& Gm & Gm & M$_\\odot$ & M$_\\
     +odot$ & km s$^{-1}$ \\\\',/,'\\hline')
      WRITE(LU_OUT,3001)OBJECT(1:INDEX(OBJECT,' ')-1),PP,T,OMEGA*180/3.1
     +4159,EE,KA,VZERO,ASINI/1.e6,FDEM,SIGMA
 3001 FORMAT(A,' & ',F13.6,' & ',F13.3,' & ',F5.1,' & ',F8.3,' & ',F8.2,
     +' & \\nodata ',' & ',F8.2,' & ',F8.2,' & \\nodata ',' & ',F10.6,' 
     +& \\nodata ',' & ',F5.2,' \\\\')
      WRITE(LU_OUT,3002)ERR(1),ERR(2),ERR(3)*180/3.14159,ERR(4),ERR(5),E
     +RR(6),ERR_ASINI/1.e6,ERR_FDEM
 3002 FORMAT(' & $\\pm',F10.6,'$ & $\\pm',F8.3,'$ & $\\pm',F5.1,'$ & $\\
     +pm',F8.3,'$ & $\\pm',F8.2,'$ & \\nodata & $\\pm',F8.2,'$ & $\\pm',
     +F8.2,'$ & \\nodata & $\\pm',F10.6,'$ & \\nodata & \\\\')
      WRITE(LU_OUT,3003)
 3003 FORMAT('\\end{tabular}',/,'\\vskip 1cm')
      WRITE(LU_OUT,3004)OBJECT(1:INDEX(OBJECT,' ')-1),OBJECT(1:INDEX(OBJ
     +ECT,' ')-1)
 3004 FORMAT('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%',/,'\\
     +begin{figure}[h]',/,'\\centering\\includegraphics*[width=12cm]{',A
     +,'.ps}',/,'\\caption{',A,': radial velocity curve.}',/,'\\end{figu
     +re}',/)
      WRITE(LU_OUT,3005)OBJECT(1:INDEX(OBJECT,' ')-1),OBJECT(1:INDEX(OBJ
     +ECT,' ')-1)
 3005 FORMAT('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%',/,'\\
     +begin{figure}[h]',/,'\\centering\\includegraphics*[width=12cm]{',A
     +,'_resi.ps}',/,'\\caption{',A,': radial velocity residuals.}',/,'\
     +\end{figure}',/)
      WRITE(LU_OUT,3009)OBJECT(1:INDEX(OBJECT,' ')-1)
 3009 FORMAT('\\twocolumn',/,/,'\\footnotesize',/,'\\begin{tabular}{rrrr
     +}',/,'\\hline',/,'\\multicolumn{4}{c}{',A,'}\\\\',/,'\\hline',/,'D
     +ate (JD) & Cycle\\hfil & $RV$\\hfil & $(O-C)$\\hfil \\\\',/,' 2400
     +000+ & & km s$^{-1}$ & km s$^{-1}$ \\\\',/,'\\hline')
      DO80007I=1,NN
      IF(POI(I).NE.0.)THEN
      PHI=(TOBS(I)-T)/PP
      WRITE(LU_OUT,3010)TOBS(I),PHI,VIT(I),RES(I)
 3010 FORMAT(F12.2,' & ',F12.2,' & ',F8.1,' & ',F8.1,' \\\\')
      ENDIF
80007 CONTINUE
      WRITE(LU_OUT,3011)
 3011 FORMAT('\\hline',/,'\\end{tabular}',/)
      WRITE(LU_OUT,3012)OBJECT(1:INDEX(OBJECT,' ')-1),OBJECT(1:INDEX(OBJ
     +ECT,' ')-1)
 3012 FORMAT('\\bigskip',/,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     +%%%%%%%',/,'\\begin{figure}[h]',/,'\\centering\\includegraphics*[w
     +idth=12cm]{',A,'_per.ps}',/,'\\caption{',A,': periodogram.}',/,'\\
     +end{figure}',/)
      WRITE(LU_OUT,3013)
 3013 FORMAT(/,'\\end{document}')
      CLOSE(LU_OUT)
      RETURN
      END
      SUBROUTINEBS1_FIT(POI,TOBS,VIT,RES,ERR,EQM,CORR,PP,T,TZERO,OMEGA,E
     +E,KA,VZERO,ASINI,FDEM,ERR_ASINI,ERR_FDEM,LU_OUT_TXT,NP,X1,Y1,NN,NI
     +TER,INTERACTIVE,STATUS)
      IMPLICITNONE
      INTEGERIDIM
      PARAMETER(IDIM=120)
      REALPOI(IDIM),TOBS(IDIM),VIT(IDIM),RES(IDIM)
      REALAA(6,IDIM),B(6),DELT(6),ERR(6),EQM(6),CORR(6)
      REALA(6,6),D(6,6),SOM_POI
      INTEGERLU_OUT_TXT,NP,X1,Y1,F,DELTA(6)
      INTEGERNN,NITER,I,J,K,ITER,STATUS
      REALEE,OMEGA,KA,VZERO,PI
      REALTZERO,SOME,ASINI,DIST
      REALVV,VIC,S,Q,SOMB,SOM,SOMD,FDEM,ERR_ASINI,ERR_FDEM
      DOUBLEPRECISIONPP,T,PHI
      LOGICALINTERACTIVE
      ITER=0
      PI=3.1415926535
      DIST=1.
      SOM_POI=0.
      DO80008I=1,NN
      SOM_POI=SOM_POI+POI(I)
80008 CONTINUE
   89 DO80009J=1,NP
      B(J)=0
80009 CONTINUE
      DO9J=1,NP
      DO9K=1,NP
      A(J,K)=0
    9 CONTINUE
      IF(ITER.EQ.0.OR.ITER.EQ.NITER)THEN
      IF(INTERACTIVE)WRITE(LU_OUT_TXT,1600)
 1600 FORMAT('OBS.',2X,'DATE (JJ)',1X,'PHASE',3X,'VR OBS.',2X,'VR CALC.'
     +,3X,'O-C',4X,'POIDS')
      ENDIF
      DO53I=1,NN
      PHI=DMOD((TOBS(I)-T)/PP,1.D0)
      IF(PHI.LT.0)PHI=1.+PHI
      CALLVITESSE_CALCULEE(PHI,VIC,VV,EE,OMEGA,KA,VZERO)
      RES(I)=VIT(I)-VIC
      IF(ITER.EQ.0.OR.ITER.EQ.NITER)THEN
      IF(INTERACTIVE)WRITE(LU_OUT_TXT,1400)I,TOBS(I),PHI,VIT(I),VIC,RES(
     +I),POI(I)
 1400 FORMAT(I3,3X,F9.3,2X,F4.3,3(2X,F7.2),3X,F5.2)
      ENDIF
      IF(ITER.EQ.NITER.OR.DIST.LT.1.E-7)GOTO53
      S=SIN(VV+OMEGA)
      Q=S*(1+EE*COS(VV))**2/(1-EE**2)**1.5
      IF(X1.EQ.0)THEN
      AA(1,I)=(2.*PI/PP**2)*(TOBS(I)-T)*KA*Q
      ENDIF
      AA(2-X1,I)=KA*Q*(2.*PI/PP)
      IF(Y1.NE.2)THEN
      AA(3-X1,I)=-KA*(EE*SIN(OMEGA)+S)
      IF(Y1.EQ.0)THEN
      AA(4-X1,I)=KA*(COS(OMEGA)-S*SIN(VV)*(2.+EE*COS(VV))/(1-EE**2))
      ENDIF
      ENDIF
      AA(5-X1-Y1,I)=EE*COS(OMEGA)+COS(VV+OMEGA)
      AA(6-X1-Y1,I)=1
      DO80010J=1,NP
      B(J)=B(J)+POI(I)*AA(J,I)*RES(I)
80010 CONTINUE
      DO3J=1,NP
      DO3K=J,NP
      A(J,K)=A(J,K)+POI(I)*AA(J,I)*AA(K,I)
    3 CONTINUE
   53 CONTINUE
      IF(ITER.EQ.NITER.OR.DIST.LT.1.E-7)GOTO81
      DO4J=2,NP
      DO4K=1,J-1
      A(J,K)=A(K,J)
    4 CONTINUE
      CALLMAT_INV(A,D,NP,INTERACTIVE,STATUS)
      IF(STATUS.NE.0)THEN
      RETURN
      ENDIF
      DO80011I=1,NP
      DELT(I)=0
80011 CONTINUE
      DO6I=1,NP
      DO6J=1,NP
      DELT(I)=DELT(I)+D(I,J)*B(J)
    6 CONTINUE
      SOMB=0
      DO80012I=1,NN
      SOMB=SOMB+POI(I)*RES(I)**2
80012 CONTINUE
      SOMB=SOMB*NN/SOM_POI
      SOMD=0
      SOME=SOMB-SOMD
      SOM=SQRT(SOME/(NN-NP))
      DO80013I=1,NP
      EQM(I)=SOM*SQRT(D(I,I))
80013 CONTINUE
      DIST=0
      F=0
      DO80014I=1,6
      DELTA(I)=0
80014 CONTINUE
      DELTA(1)=X1
      DELTA(4)=Y1
      IF(Y1.EQ.2)THEN
      DELTA(3)=1
      DELTA(4)=1
      ENDIF
      DO801I=1,6
      F=F+DELTA(I)
      J=I-F+DELTA(I)
      ERR(I)=(1-DELTA(I))*EQM(J)
      CORR(I)=(1-DELTA(I))*DELT(J)
      DIST=DIST+ABS(CORR(I))
  801 CONTINUE
      PP=PP+CORR(1)
      T=TZERO+DMOD(T-TZERO+CORR(2),PP)
      OMEGA=AMOD(OMEGA+CORR(3),2.*PI)
      IF(OMEGA.LT.0)OMEGA=OMEGA+2.*PI
      IF(EE+CORR(4).LT.0.)THEN
      EE=EE/2.
      ELSE
      EE=EE+CORR(4)
      ENDIF
      IF(EE.GE.1.)EE=0.9
      KA=KA+CORR(5)
      IF(KA.LT.0.)KA=-KA
      VZERO=VZERO+CORR(6)
      ASINI=0
      FDEM=0
      IF(ABS(EE).LT.1.)THEN
      ASINI=43200.*PP*KA*SQRT(1.-EE**2)/PI
      FDEM=0.00000010385*(1.-EE**2)**1.5*PP*KA**3
      ENDIF
      ERR_ASINI=ASINI*(ERR(1)/PP+ERR(5)/KA+ERR(4)*EE/(1.-EE**2))
      ERR_FDEM=FDEM*(3.*EE*ERR(4)/(1.-EE**2)+ERR(1)/PP+3.*ERR(5)/KA)
      ITER=ITER+1
      IF(INTERACTIVE)THEN
      WRITE(6,2200)ITER,DIST
      WRITE(LU_OUT_TXT,2200)ITER,DIST
 2200 FORMAT(3X,I2,'eme iteration',' dist=',E10.3)
      IF(X1.EQ.1)WRITE(LU_OUT_TXT,1700)
 1700 FORMAT(3X,'PERIODE FIXEE')
      IF(Y1.EQ.1)WRITE(LU_OUT_TXT,1800)
 1800 FORMAT(3X,'EXCENTRICITE FIXEE')
      IF(Y1.EQ.2)WRITE(LU_OUT_TXT,1801)
 1801 FORMAT(3X,'ORBITE CIRCULAIRE')
      WRITE(LU_OUT_TXT,1300)CORR(1),PP,ERR(1),CORR(2),T,ERR(2),CORR(3)*1
     +80./PI,OMEGA*180./PI,ERR(3)*180./PI,CORR(3),OMEGA,ERR(3),CORR(4),E
     +E,ERR(4),CORR(5),KA,ERR(5),CORR(6),VZERO,ERR(6),ASINI,ERR_ASINI,FD
     +EM,ERR_FDEM
 1300 FORMAT(1X,'CORRECTIONS',15X,'NOUVEAUX ELEMENTS',15X,'ERREURS STD.'
     +,/,3X,E9.3,17X,'P =',F15.7,13X,E9.3,/,3X,E9.3,17X,'T =',F11.3,17X,
     +E9.3,/,3X,E9.3,13X,'OMEGA =',F13.5,15X,E9.3,4X,'(deg.)',/,3X,E9.3,
     +13X,'OMEGA =',F13.5,15X,E9.3,4X,'(rad.)',/,3X,E9.3,17X,'E =',F13.5
     +,15X,E9.3,/,3X,E9.3,17X,'K =',F10.2,18X,E9.3,/,3X,E9.3,16X,'V0 =',
     +F10.2,18X,E9.3,/,24X,'A.SINI =',6X,E10.4,12X,E9.3,/,26X,'F(M) =',6
     +X,E10.4,12X,E9.3,/,/)
      ENDIF
      IF(ABS(EE).GE.1.)THEN
      WRITE(6,2000)
 2000 FORMAT(5X,'BS1_FIT/FATAL ERROR/ ABS(EE).GE.1,',' POURSUITE DU CALC
     +UL IMPOSSIBLE')
      IF(INTERACTIVE)THEN
      WRITE(LU_OUT_TXT,2000)
      CLOSE(LU_OUT_TXT)
      ENDIF
      STATUS=-2
      RETURN
      ELSE
      GOTO89
      ENDIF
   81 STATUS=0
      RETURN
      END
