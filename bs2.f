      PROGRAMBS2
      IMPLICITNONE
      INTEGERIDIM
      PARAMETER(IDIM=120)
      REALTOBS(IDIM),VIT1(IDIM),VIT2(IDIM),POI1(IDIM),POI2(IDIM)
      REALRES1(IDIM),RES2(IDIM)
      REALAA(7,200),B(7),DELT(7),ERR(7),EQM(7),CORR(7)
      REALA(7,7),D(7,7)
      DOUBLEPRECISIONPHI,PP,T,TZERO
      INTEGERNP,X1,Y1,F,DELTA(7),LU_IN,LU_OUT,IT,NN,ITER,NITER
      REALOMEGA,EE,PI,VV,VIC1,VIC2,OME
      REALKA,K1,K2,SIGMA1,SIGMA2,SIGMA3,ENNE,NU,Q,S
      REALFDEM1,FDEM2,A1SINI,A2SINI,SOM,SOMB,SOMD,SOME
      REALDIST,VZERO,ERR_A1SINI,ERR_A2SINI
      REALERR_FDEM1,ERR_FDEM2,POI_0,SOM_POI
      INTEGERI,J,K,L,J1,J2,JJ,N2,IOPT
      CHARACTEROBJECT*20,OUT_RESULTS*60,OUT_LATEX*60
      LOGICALPUBLI
      PI=3.1415926535
      WRITE(6,*)' Programme bs2: version du 18/10/2005'
      WRITE(6,*)' Titre (sans blancs, car utilise comme prefixe):'
      READ(5,10)OBJECT
   10 FORMAT(A)
      LU_IN=7
      OPEN(LU_IN,FILE='VR_BS2.DAT',STATUS='OLD')
      LU_OUT=8
      OUT_RESULTS=OBJECT(1:INDEX(OBJECT,' ')-1)//'.txt'
      OPEN(LU_OUT,FILE=OUT_RESULTS,STATUS='UNKNOWN')
   88 READ(LU_IN,1000)PP,TZERO,OMEGA,EE,K1,K2,VZERO,X1,Y1
 1000 FORMAT(F10.7,F10.3,2F10.4,2F10.2,F6.2,2X,2I1)
      WRITE(6,1001)PP,TZERO,OMEGA*180./PI,EE,K1,K2,VZERO
 1001 FORMAT(' P=',F12.7,' T0=',F10.3,' Omega=',F9.3,'(deg) E=',F9.3,' K
     +1=',F9.2,' K2=',F9.2,'V0=',F9.2)
      T=TZERO
      WRITE(6,*)'Option: X1=',X1,' Y1=',Y1
      T=TZERO
      IF(Y1.EQ.2)THEN
      EE=0.
      OMEGA=0.
      ENDIF
  881 READ(LU_IN,1100)NN,NITER
 1100 FORMAT(I2,6X,I2)
      WRITE(6,*)'Nombre de mesures:',NN,' Nombre d''iterations',NITER
      READ(LU_IN,1200)(TOBS(I),VIT1(I),POI1(I),VIT2(I),POI2(I),I=1,NN)
 1200 FORMAT((F9.3,2(F7.1,F6.2,7X)))
      CLOSE(LU_IN)
      SOM_POI=0.
      DO80000I=1,NN
      SOM_POI=SOM_POI+POI1(I)+POI2(I)
80000 CONTINUE
      NP=7-X1-Y1
      WRITE(LU_OUT,1500)PP,T,OMEGA,EE,K1,K2,VZERO
 1500 FORMAT(13X,'ELEMENTS PROVISOIRES',/,16X,'P =',F15.7,/,16X,'T =',F1
     +1.3,/,12X,'OMEGA =',F13.5,/,16X,'E =',F13.5,/,15X,'K1 =',F10.2,/,1
     +5X,'K2 =',F10.2,/,15X,'V0 =',F10.2,/,/)
      DIST=1.
      N2=NN*2
      DO80001I=1,NP
      DO80002J=1,N2
      AA(I,J)=0
80002 CONTINUE
80001 CONTINUE
      ITER=0
      DO80003IT=0,NITER
      DO80004J=1,NP
      B(J)=0
80004 CONTINUE
      DO80005J=1,NP
      DO80006K=1,NP
      A(J,K)=0
80006 CONTINUE
80005 CONTINUE
      NU=0
      IF(ITER.EQ.0.OR.ITER.EQ.NITER)THEN
      WRITE(LU_OUT,1600)
 1600 FORMAT('OBS.',1X,'DATE (JJ)',2X,'PHASE',1X,'V_OBS1',2X,'V_CALC1',1
     +X,'(O-C)1',1X,'POIDS1',1X,'V_OBS2',1X,'V_CALC2',1X,'(O-C)2',1X,'PO
     +IDS2',/)
      ENDIF
      ENNE=2.*PI/PP
      DO53I=1,NN
      PHI=DMOD((TOBS(I)-T)/PP,1.D0)
      IF(PHI.LT.0)PHI=1.+PHI
      CALLVITESSE_CALCULEE_BS2(PHI,VIC1,VIC2,VV,EE,OMEGA,K1,K2,VZERO)
      RES1(I)=VIT1(I)-VIC1
      RES2(I)=VIT2(I)-VIC2
      IF(ITER.EQ.0.OR.ITER.EQ.NITER)THEN
      IF(POI1(I).EQ.0)THEN
      WRITE(LU_OUT,2000)I,TOBS(I),PHI,VIC1,VIT2(I),VIC2,RES2(I),POI2(I)
 2000 FORMAT(I3,2X,F9.3,2X,F4.3,9X,F7.2,14X,3(1X,F7.2),1X,F5.2)
      ELSE
      IF(POI2(I).NE.0)THEN
      WRITE(LU_OUT,1400)I,TOBS(I),PHI,VIT1(I),VIC1,RES1(I),POI1(I),VIT2(
     +I),VIC2,RES2(I),POI2(I)
 1400 FORMAT(I3,2X,F9.3,2X,F4.3,1X,3(F7.2,1X),F5.2,1X,3(F7.2,1X),F5.2)
      ELSE
      WRITE(LU_OUT,2100)I,TOBS(I),PHI,VIT1(I),VIC1,RES1(I),POI1(I),VIC2
 2100 FORMAT(I3,2X,F9.3,2X,F4.3,1X,3(F7.2,1X),F5.2,9X,F7.2)
      ENDIF
      ENDIF
      ENDIF
      IF(ITER.EQ.NITER.OR.DIST.LT.1.E-7)GOTO53
      DO28K=1,2
      L=I+(K-1)*NN
      IF(K.EQ.1)THEN
      IF(POI1(I).EQ.0)GOTO28
      NU=NU+1
      KA=K1
      OME=OMEGA
      J1=5
      J2=6
      ELSE
      IF(POI2(I).EQ.0)GOTO28
      NU=NU+1
      KA=K2
      OME=OMEGA+PI
      J1=6
      J2=5
      ENDIF
      S=SIN(VV+OME)
      Q=S*(1+EE*COS(VV))**2/(1-EE**2)**1.5
      IF(X1.EQ.0)THEN
      AA(1,L)=ENNE*(TOBS(I)-T)*KA*Q/PP
      ENDIF
      AA(2-X1,L)=KA*ENNE*Q
      IF(Y1.NE.2)THEN
      AA(3-X1,L)=-KA*(EE*SIN(OME)+S)
      IF(Y1.EQ.0)THEN
      AA(4-X1,L)=KA*(COS(OME)-S*SIN(VV)*(2+EE*COS(VV))/(1-EE**2))
      ENDIF
      ENDIF
      AA(J1-X1-Y1,L)=EE*COS(OME)+COS(VV+OME)
      AA(J2-X1-Y1,L)=0
      AA(NP,L)=1.
   28 CONTINUE
      DO80007J=1,NP
      B(J)=B(J)+POI1(I)*AA(J,I)*RES1(I)+POI2(I)*AA(J,I+NN)*RES2(I)
80007 CONTINUE
      DO80008J=1,NP
      DO80009K=J,NP
      A(J,K)=A(J,K)+POI1(I)*AA(J,I)*AA(K,I)+POI2(I)*AA(J,I+NN)*AA(K,I+NN
     +)
80009 CONTINUE
80008 CONTINUE
   53 CONTINUE
      IF(ITER.EQ.NITER.OR.DIST.LT.1.E-7)GOTO81
      DO80010J=2,NP
      JJ=J-1
      DO80011K=1,JJ
      A(J,K)=A(K,J)
80011 CONTINUE
80010 CONTINUE
      CALLMAT_INV7(A,D,NP)
      DO80012I=1,NP
      DELT(I)=0
80012 CONTINUE
      DO6I=1,NP
      DO6J=1,NP
      DELT(I)=DELT(I)+D(I,J)*B(J)
    6 CONTINUE
      SOMB=0
      DO7I=1,NN
      SOMB=SOMB+POI1(I)*RES1(I)**2+POI2(I)*RES2(I)**2
    7 CONTINUE
      SOMB=SOMB*NU/SOM_POI
      SOMD=0
      SOME=SOMB-SOMD
      SOM=SQRT(SOME/(NU-NP))
      DO80I=1,NP
      EQM(I)=SOM*SQRT(D(I,I))
   80 CONTINUE
      DIST=0
      F=0
      DO80013I=1,7
      DELTA(I)=0
80013 CONTINUE
      DELTA(1)=X1
      DELTA(4)=Y1
      IF(Y1.EQ.2)THEN
      DELTA(3)=1
      DELTA(4)=1
      ENDIF
      DO80014I=1,7
      F=F+DELTA(I)
      J=I-F+DELTA(I)
      ERR(I)=(1-DELTA(I))*EQM(J)
      CORR(I)=(1-DELTA(I))*DELT(J)
      DIST=DIST+ABS(CORR(I))
80014 CONTINUE
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
      K1=K1+CORR(5)
      IF(K1.LT.0.)K1=-K1
      K2=K2+CORR(6)
      IF(K2.LT.0.)K2=-K2
      VZERO=VZERO+CORR(7)
      A1SINI=0
      A2SINI=0
      FDEM1=0
      FDEM2=0
      IF(ABS(EE).LT.1.)THEN
      A1SINI=43200.*PP*K1*SQRT(1-EE**2)/PI
      A2SINI=43200.*PP*K2*SQRT(1-EE**2)/PI
      FDEM1=0.00000010385*(K1+K2)**2*K2*PP*(1-EE**2)**1.5
      FDEM2=0.00000010385*(K1+K2)**2*K1*PP*(1-EE**2)**1.5
      ENDIF
      ERR_A1SINI=A1SINI*(ERR(1)/PP+ERR(5)/K1+ERR(4)*EE/(1.-EE**2))
      ERR_A2SINI=A2SINI*(ERR(1)/PP+ERR(6)/K2+ERR(4)*EE/(1.-EE**2))
      ERR_FDEM1=FDEM1*(3.*EE*ERR(4)/(1.-EE**2)+ERR(1)/PP+ERR(5)/K1+2.*(E
     +RR(5)+ERR(6))/(K1+K2))
      ERR_FDEM2=FDEM2*(3.*EE*ERR(4)/(1.-EE**2)+ERR(1)/PP+ERR(6)/K2+2.*(E
     +RR(5)+ERR(6))/(K1+K2))
      ITER=ITER+1
      WRITE(LU_OUT,2200)ITER
 2200 FORMAT(/,5X,I2,'E ITERATION')
      IF(X1.EQ.1)WRITE(LU_OUT,1700)
 1700 FORMAT(5X,'PERIODE FIXEE')
      IF(Y1.EQ.1)WRITE(LU_OUT,1800)
 1800 FORMAT(5X,'EXCENTRICITE FIXEE')
      IF(Y1.EQ.2)WRITE(LU_OUT,1801)
 1801 FORMAT(5X,'ORBITE CIRCULAIRE ')
      WRITE(LU_OUT,1300)CORR(1),PP,ERR(1),CORR(2),T,ERR(2),CORR(3)*180./
     +PI,OMEGA*180./PI,ERR(3)*180./PI,CORR(3),OMEGA,ERR(3),CORR(4),EE,ER
     +R(4),CORR(5),K1,ERR(5),CORR(6),K2,ERR(6),CORR(7),VZERO,ERR(7),A1SI
     +NI,ERR_A1SINI,A2SINI,ERR_A2SINI,FDEM1,ERR_FDEM1,FDEM2,ERR_FDEM2
 1300 FORMAT(5X,'CORRECTIONS',15X,'NOUVEAUX ELEMENTS',15X,'ERREURS STD.'
     +,/,3X,E9.3,17X,'P =',F15.7,13X,E9.3,/,3X,E9.3,17X,'T =',F11.3,17X,
     +E9.3,/,3X,E9.3,13X,'OMEGA =',F13.5,15X,E9.3,4X,'(deg.)',/,3X,E9.3,
     +13X,'OMEGA =',F13.5,15X,E9.3,4X,'(rad.)',/,3X,E9.3,17X,'E =',F13.5
     +,15X,E9.3,/,3X,E9.3,16X,'K1 =',F10.2,18X,E9.3,/,3X,E9.3,16X,'K2 ='
     +,F10.2,18X,E9.3,/,3X,E9.3,16X,'V0 =',F10.2,18X,E9.3,/,23X,'A1.SINI
     + =',6X,E10.4,12X,E9.3,/,23X,'A2.SINI =',6X,E10.4,12X,E9.3,/,22X,'M
     +1.SIN3I =',6X,E10.4,12X,E9.3,/,22X,'M2.SIN3I =',6X,E10.4,12X,E9.3,
     +/)
      IF(ABS(EE).GE.1.)THEN
      WRITE(LU_OUT,2300)
 2300 FORMAT(5X,'ABS(E).GE.1, Poursuite du calcul impossible!')
      CLOSE(LU_OUT)
      STOP
      ENDIF
80003 CONTINUE
   81 CONTINUE
      WRITE(6,*)' Courbe:  0=aucune  1=orbite  2=residus'
      WRITE(6,*)'          3=les deux  4=mode pour publi. ?'
      READ(5,*)IOPT
      PUBLI=(IOPT.EQ.4)
      IF(IOPT.EQ.1)THEN
      CALLJLP_PLOT_BS2(TOBS,VIT1,VIT2,POI1,POI2,T,PP,NN,EE,OMEGA,K1,K2,V
     +ZERO,OBJECT,PUBLI)
      ELSEIF(IOPT.EQ.2)THEN
      CALLJLP_PLOT_RES2(TOBS,RES1,RES2,POI1,POI2,NN,OBJECT,PUBLI)
      ELSEIF(IOPT.EQ.3.OR.IOPT.EQ.4)THEN
      CALLJLP_PLOT_BS2(TOBS,VIT1,VIT2,POI1,POI2,T,PP,NN,EE,OMEGA,K1,K2,V
     +ZERO,OBJECT,PUBLI)
      CALLJLP_PLOT_RES2(TOBS,RES1,RES2,POI1,POI2,NN,OBJECT,PUBLI)
      ENDIF
      POI_0=1
      DO80015I=1,NN
      IF(POI1(I).NE.1)POI_0=POI1(I)
80015 CONTINUE
      IF(POI_0.NE.1)THEN
      DO80016I=1,2
      CALLRMS_RESIDUS_PARTIELS(RES1,POI1,NN,SIGMA1,POI_0)
      CALLRMS_RESIDUS_PARTIELS(RES2,POI2,NN,SIGMA2,POI_0)
      WRITE(LU_OUT,3001)POI_0,SIGMA1,SIGMA2
 3001 FORMAT(/,'Résidus partiels pour poids=',F6.3,/,'Sigma O-C: Etoile 
     +primaire (1)   ',F5.2,' km s$^{-1}$',/,'Sigma O-C: Etoile secondai
     +re (2) ',F5.2,' km s$^{-1}$')
      POI_0=1.0
80016 CONTINUE
      ENDIF
      CALLRMS_RESIDUS_NP(RES1,POI1,NN,SIGMA1)
      CALLRMS_RESIDUS_NP(RES2,POI2,NN,SIGMA2)
      CALLRMS_RESIDUS_NP_2(RES1,RES2,POI1,POI2,NN,SIGMA3)
      WRITE(6,3000)SIGMA1,SIGMA2,SIGMA3
      WRITE(LU_OUT,3002)SIGMA1,SIGMA2,SIGMA3
 3002 FORMAT(/,'Sigma O-C (non pondéré): Etoile primaire (1)   ',F5.2,' 
     +km s$^{-1}$',/,'Sigma O-C (non pondéré): Etoile secondaire (2) ',F
     +5.2,' km s$^{-1}$',/,'Sigma O-C (non pondéré): global: ',F5.2,' km
     + s$^{-1}$',/)
      CALLRMS_RESIDUS(RES1,POI1,NN,SIGMA1)
      CALLRMS_RESIDUS(RES2,POI2,NN,SIGMA2)
      CALLRMS_RESIDUS_2(RES1,RES2,POI1,POI2,NN,SIGMA3)
      WRITE(6,3000)SIGMA1,SIGMA2,SIGMA3
      WRITE(LU_OUT,3000)SIGMA1,SIGMA2,SIGMA3
 3000 FORMAT(/,'Sigma O-C (pondéré): Etoile primaire (1)   ',F5.2,' km s
     +$^{-1}$',/,'Sigma O-C (pondéré): Etoile secondaire (2) ',F5.2,' km
     + s$^{-1}$',/,'Sigma O-C (pondéré): global: ',F5.2,' km s$^{-1}$',/
     +)
      CLOSE(LU_OUT)
      OUT_LATEX=OBJECT(1:INDEX(OBJECT,' ')-1)//'.tex'
      CALLSB2_TO_LATEX(PP,T,OMEGA,EE,K1,K2,VZERO,A1SINI,A2SINI,FDEM1,FDE
     +M2,SIGMA1,SIGMA2,SIGMA3,ERR,ERR_A1SINI,ERR_A2SINI,ERR_FDEM1,ERR_FD
     +EM2,TOBS,VIT1,VIT2,RES1,RES2,POI1,POI2,NN,OUT_LATEX,OBJECT)
      WRITE(6,1300)CORR(1),PP,ERR(1),CORR(2),T,ERR(2),CORR(3)*180./PI,OM
     +EGA*180./PI,ERR(3)*180./PI,CORR(3),OMEGA,ERR(3),CORR(4),EE,ERR(4),
     +CORR(5),K1,ERR(5),CORR(6),K2,ERR(6),CORR(7),VZERO,ERR(7),A1SINI,ER
     +R_A1SINI,A2SINI,ERR_A2SINI,FDEM1,ERR_FDEM1,FDEM2,ERR_FDEM2
      OPEN(LU_OUT,FILE='BS2_GRAV.DAT',STATUS='UNKNOWN')
      DO80017I=1,NN
      IF(POI1(I)*POI2(I).NE.0.)THEN
      VIC1=(VIT1(I)-VZERO)/K1+(VIT2(I)-VZERO)/K2
      VIC1=VIC1*(K1+K2)/2.
      WRITE(LU_OUT,1201)TOBS(I),VIC1,1.0
 1201 FORMAT(F9.3,F7.1,F6.2)
      ENDIF
80017 CONTINUE
      CLOSE(LU_OUT)
      WRITE(6,*)' Fichiers de sortie: ',OUT_RESULTS(1:INDEX(OUT_RESULTS,
     +' ')-1),', ',OUT_LATEX(1:INDEX(OUT_LATEX,' ')-1),', '
      WRITE(6,*)' mouvement du centre de gravite dans "BS2_GRAV.DAT",'
      IF(IOPT.EQ.3)THEN
      WRITE(6,*)' et graphes dans: ',OBJECT(1:INDEX(OBJECT,' ')-1),'.ps 
     +et ',OBJECT(1:INDEX(OBJECT,' ')-1),'_resi.ps'
      ELSEIF(IOPT.EQ.1)THEN
      WRITE(6,*)' et graphe dans: ',OBJECT(1:INDEX(OBJECT,' ')-1),'.ps'
      ELSEIF(IOPT.EQ.2)THEN
      WRITE(6,*)' et graphe dans: ',OBJECT(1:INDEX(OBJECT,' ')-1),'_resi
     +.ps'
      ENDIF
      STOP
      END
      SUBROUTINEJLP_PLOT_BS2(TOBS,VIT1,VIT2,POI1,POI2,T,PP,NN,EE,OMEGA,K
     +1,K2,VZERO,OBJECT,PUBLI)
      PARAMETER(IDIM=1024,KCUR=4)
      REALTOBS(120),VIT1(120),VIT2(120),POI1(120),POI2(120)
      DOUBLEPRECISIONT,PP,PHI
      REALPI,EE,OMEGA,K1,K2,VZERO,VV
      INTEGERKK,NPTS(KCUR),I1,I2
      REALXPLOT(IDIM,KCUR),YPLOT(IDIM,KCUR)
      CHARACTERCHAR1*30,CHAR2*30,CHAR3*40,OBJECT*20
      CHARACTERPLOTDEV*40,IN_FILE*40,IN_COMMENTS*80
      CHARACTER*4NCHAR(KCUR)
      LOGICALPUBLI
      I1=1
      I2=1
      DO80018I=1,NN
      PHI=DMOD((TOBS(I)-T)/PP,1.D0)
      IF(PHI.LT.0)PHI=1.+PHI
      IF(POI1(I).NE.0.)THEN
      XPLOT(I1,1)=PHI
      YPLOT(I1,1)=VIT1(I)
      I1=I1+1
      ENDIF
      IF(POI2(I).NE.0.)THEN
      XPLOT(I2,2)=PHI
      YPLOT(I2,2)=VIT2(I)
      I2=I2+1
      ENDIF
80018 CONTINUE
      NCHAR(1)='913'
      NCHAR(2)='813'
      NPTS(1)=I1-1
      NPTS(2)=I2-1
      PI=3.1415926535
      DO80019I=1,1000
      PHI=REAL(I)/REAL(1000)
      CALLVITESSE_CALCULEE_BS2(PHI,VIC1,VIC2,VV,EE,OMEGA,K1,K2,VZERO)
      XPLOT(I,3)=PHI
      YPLOT(I,3)=VIC1
      XPLOT(I,4)=PHI
      YPLOT(I,4)=VIC2
80019 CONTINUE
      NCHAR(3)='L0'
      NCHAR(4)='L1'
      NPTS(3)=1000
      NPTS(4)=1000
      KK=4
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
      SUBROUTINEVITESSE_CALCULEE_BS2(PHI,VIC1,VIC2,VV,EE,OMEGA,K1,K2,VZE
     +RO)
      IMPLICITNONE
      DOUBLEPRECISIONPHI
      REALM,PI,EMME,ANOM,ANOM_OLD,VV
      REALEE,OMEGA,K1,K2,VZERO,VIC1,VIC2
      INTEGERI
      PI=3.1415926535
      M=2.*PI*PHI
      EMME=AMOD(M,2.*PI)
      ANOM_OLD=EMME
      DO80020I=1,1000
      ANOM=EMME+EE*SIN(ANOM_OLD)
      IF(ABS(ANOM-ANOM_OLD).LT.1.E-6)GOTO12
      ANOM_OLD=ANOM
80020 CONTINUE
   12 VV=2.*ATAN(TAN(ANOM/2.)*SQRT((1+EE)/(1-EE)))
      VIC1=VZERO+K1*(EE*COS(OMEGA)+COS(VV+OMEGA))
      VIC2=VZERO-K2*(EE*COS(OMEGA)+COS(VV+OMEGA))
      RETURN
      END
      SUBROUTINERMS_RESIDUS_NP(RES,POI,NN,SIGMA)
      INTEGERI
      REALSUM,SUMSQ
      REALRES(120),POI(120)
      SUMSQ=0.
      SUM=0.
      DO80021I=1,NN
      IF(POI(I).GT.0)THEN
      SUMSQ=SUMSQ+RES(I)*RES(I)
      SUM=SUM+1.0
      ENDIF
80021 CONTINUE
      IF(SUM.NE.0.)THEN
      SIGMA=SQRT(SUMSQ/SUM)
      ELSE
      SIGMA=0.
      ENDIF
      RETURN
      END
      SUBROUTINERMS_RESIDUS(RES,POI,NN,SIGMA)
      INTEGERI
      REALSUM,SUMSQ
      REALRES(120),POI(120)
      SUMSQ=0.
      SUM=0.
      DO80022I=1,NN
      SUMSQ=SUMSQ+POI(I)*RES(I)*RES(I)
      SUM=SUM+POI(I)
80022 CONTINUE
      IF(SUM.NE.0.)THEN
      SIGMA=SQRT(SUMSQ/SUM)
      ELSE
      SIGMA=0.
      ENDIF
      RETURN
      END
      SUBROUTINERMS_RESIDUS_PARTIELS(RES,POI,NN,SIGMA,POI_0)
      INTEGERI
      REALSUM,SUMSQ,POI_0
      REALRES(120),POI(120)
      SUMSQ=0.
      SUM=0.
      DO80023I=1,NN
      IF(POI_0.EQ.POI(I))THEN
      SUMSQ=SUMSQ+POI(I)*RES(I)*RES(I)
      SUM=SUM+POI(I)
      ENDIF
80023 CONTINUE
      IF(SUM.NE.0.)THEN
      SIGMA=SQRT(SUMSQ/SUM)
      ELSE
      SIGMA=0.
      ENDIF
      RETURN
      END
      SUBROUTINERMS_RESIDUS_2(RES1,RES2,POI1,POI2,NN,SIGMA)
      INTEGERI
      REALSUM,SUMSQ
      REALRES1(120),RES2(120),POI1(120),POI2(120)
      SUMSQ=0.
      SUM=0.
      DO80024I=1,NN
      SUMSQ=SUMSQ+POI1(I)*RES1(I)*RES1(I)+POI2(I)*RES2(I)*RES2(I)
      SUM=SUM+POI1(I)+POI2(I)
80024 CONTINUE
      IF(SUM.NE.0.)THEN
      SIGMA=SQRT(SUMSQ/SUM)
      ELSE
      SIGMA=0.
      ENDIF
      RETURN
      END
      SUBROUTINERMS_RESIDUS_NP_2(RES1,RES2,POI1,POI2,NN,SIGMA)
      INTEGERI
      REALSUM,SUMSQ
      REALRES1(120),RES2(120),POI1(120),POI2(120)
      SUMSQ=0.
      SUM=0.
      DO80025I=1,NN
      IF(POI1(I).GT.0)THEN
      SUMSQ=SUMSQ+RES1(I)*RES1(I)
      SUM=SUM+1.0
      ENDIF
80025 CONTINUE
      DO80026I=1,NN
      IF(POI2(I).GT.0)THEN
      SUMSQ=SUMSQ+RES2(I)*RES2(I)
      SUM=SUM+1.0
      ENDIF
80026 CONTINUE
      IF(SUM.NE.0.)THEN
      SIGMA=SQRT(SUMSQ/SUM)
      ELSE
      SIGMA=0.
      ENDIF
      RETURN
      END
      SUBROUTINEMAT_INV7(A,D,NP)
      IMPLICITNONE
      INTEGERNP
      REALA(7,7),C(7,7),D(7,7),G,H
      INTEGERI1,J1,I,J,K,KK,L,L1
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
      write(6,*)' Pivot null pour K=',K
      write(6,*)' (JLP: j''ai des doutes sur la suite...)'
      KK=K+1
      DO12L=KK,NP
      IF(A(L,K).NE.0)THEN
      L1=L
      GOTO13
      ENDIF
   12 CONTINUE
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
      RETURN
      END
      SUBROUTINEJLP_PLOT_RES2(TOBS,RES1,RES2,POI1,POI2,NN,OBJECT,PUBLI)
      IMPLICITNONE
      INTEGERDIM,KCUR
      PARAMETER(DIM=1024,KCUR=2)
      REALTOBS(120),RES1(120),RES2(120),POI1(120),POI2(120)
      INTEGERKK,NPTS(KCUR),I,NN,I1,I2
      REALXPLOT(DIM,KCUR),YPLOT(DIM,KCUR)
      CHARACTERCHAR1*30,CHAR2*30,CHAR3*40,OBJECT*20
      CHARACTERPLOTDEV*40,IN_FILE*40,IN_COMMENTS*80
      CHARACTER*4NCHAR(KCUR)
      LOGICALPUBLI
      I1=1
      I2=1
      DO80027I=1,NN
      IF(POI1(I).NE.0.)THEN
      XPLOT(I1,1)=TOBS(I)
      YPLOT(I1,1)=RES1(I)
      I1=I1+1
      ENDIF
      IF(POI2(I).NE.0.)THEN
      XPLOT(I2,2)=TOBS(I)
      YPLOT(I2,2)=RES2(I)
      I2=I2+1
      ENDIF
80027 CONTINUE
      NCHAR(1)='910'
      NCHAR(2)='810'
      NPTS(1)=I1-1
      NPTS(2)=I2-1
      KK=2
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
      SUBROUTINESB2_TO_LATEX(PP,T,OMEGA,EE,K1,K2,VZERO,A1SINI,A2SINI,FDE
     +M1,FDEM2,SIGMA1,SIGMA2,SIGMA3,ERR,ERR_A1SINI,ERR_A2SINI,ERR_FDEM1,
     +ERR_FDEM2,TOBS,VIT1,VIT2,RES1,RES2,POI1,POI2,NN,OUT_LATEX,OBJECT)
      IMPLICITNONE
      DOUBLEPRECISIONPP,T,PHI
      INTEGERNN
      REALTOBS(NN),POI1(NN),POI2(NN),VIT1(NN),VIT2(NN)
      REALRES1(NN),RES2(NN),ERR(7)
      REALOMEGA,EE,K1,K2,VZERO,A1SINI,A2SINI,FDEM1,FDEM2
      REALSIGMA1,SIGMA2,SIGMA3
      REALERR_A1SINI,ERR_A2SINI,ERR_FDEM1,ERR_FDEM2
      CHARACTEROUT_LATEX*40,OBJECT*20
      INTEGERI,LU_OUT
      LU_OUT=8
      OPEN(LU_OUT,FILE=OUT_LATEX,STATUS='UNKNOWN')
      WRITE(LU_OUT,3000)OBJECT(1:INDEX(OBJECT,' ')-1)
 3000 FORMAT('\\documentclass{article}',/,'\\usepackage{graphicx}',/,'\\
     +voffset=-1cm',/,'\\hoffset=-4cm',/,'\\textwidth=15.6cm',/,'\\texth
     +eight=26cm',/,'\\newcommand{\\nodata}{\\ldots}',/,/,'\\begin{docum
     +ent}',/,/,'\\centerline{\\large \\bf ',A,'}',/,/,'\\vskip 1cm',/,'
     +\\tabcolsep=1mm',/,'\\begin{tabular}{lcccccccccccccc}',/,'\\hline'
     +,/,'Name & $P$ & $T_0$ (JD)& $\\omega$ & $e$ &',' $K_1$ & $K_2$ & 
     +$V_0$ ','& $a_1 \\sin i$ & $a_2 \\sin i$',/,'& $f(m_1)$ & $f(m_2)$
     + & $\\sigma_{1(O-C)}$ ',/,'& $\\sigma_{2(O-C)}$ & $\\sigma_{(O-C)}
     +$ \\\\',/,'& days & 2400000+ & deg. & & km s$^{-1}$ &',' km s$^{-1
     +}$ & km s$^{-1}$ ','& Gm & Gm & M$_\\odot$ & M$_\\odot$ &',' km s$
     +^{-1}$ & km s$^{-1}$ ','& km s$^{-1}$ \\\\',/,'\\hline')
      WRITE(LU_OUT,3001)OBJECT(1:INDEX(OBJECT,' ')-1),PP,T,OMEGA*180/3.1
     +4159,EE,K1,K2,VZERO,A1SINI/1.e6,A2SINI/1.e6,FDEM1,FDEM2,SIGMA1,SIG
     +MA2,SIGMA3
 3001 FORMAT(A,' & ',F12.5,' & ',F12.2,' & ',F5.1,' & ',F8.3,' & ',F8.2,
     +' & ',F8.2,' & ',F8.2,/,' & ',F8.2,' & ',F8.2,' & ',F8.3,' & ',F8.
     +3,' & ',F5.2,' & ',F5.2,' & ',F5.2,' \\\\')
      WRITE(LU_OUT,3002)ERR(1),ERR(2),ERR(3)*180/3.14159,ERR(4),ERR(5),E
     +RR(6),ERR(7),ERR_A1SINI/1.e6,ERR_A2SINI/1.e6,ERR_FDEM1,ERR_FDEM2
 3002 FORMAT(' & $\\pm',F10.5,'$ & $\\pm',F8.2,'$ & $\\pm',F5.1,'$ & $\\
     +pm',F8.3,'$ & $\\pm',F8.2,'$ ',/,'& $\\pm',F8.2,'$ & $\\pm',F8.2,'
     +$ & $\\pm',F8.2,'$ & $\\pm',F8.2,'$ & $\\pm',F8.3,'$ & $\\pm',F8.2
     +,'$ & & & \\\\')
      WRITE(LU_OUT,3003)
 3003 FORMAT('\\end{tabular}',/,'\\vskip 1cm')
      WRITE(LU_OUT,3004)OBJECT(1:INDEX(OBJECT,' ')-1),OBJECT(1:INDEX(OBJ
     +ECT,' ')-1)
 3004 FORMAT('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%',/,'\\
     +begin{figure}[h]',/,'\\centering\\includegraphics*[width=12cm]{',A
     +,'.ps}',/,'\\caption{',A,': radial velocity curve.',/,'Filled circ
     +les: primary component, open circles: secondary.}',/,'\\end{figure
     +}',/)
      WRITE(LU_OUT,3005)OBJECT(1:INDEX(OBJECT,' ')-1),OBJECT(1:INDEX(OBJ
     +ECT,' ')-1)
 3005 FORMAT('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%',/,'\\
     +begin{figure}[h]',/,'\\centering\\includegraphics*[width=12cm]{',A
     +,'_resi.ps}',/,'\\caption{',A,': radial velocity residuals.',/,'Fi
     +lled circles: primary component, open circles: secondary.}',/,'\\e
     +nd{figure}',/)
      WRITE(LU_OUT,3009)OBJECT(1:INDEX(OBJECT,' ')-1)
 3009 FORMAT('\\twocolumn',/,/,'\\footnotesize',/,'\\begin{tabular}{rrcr
     +cr}',/,'\\hline',/,'\\multicolumn{5}{c}{',A,'}\\\\',/,'\\hline',/,
     +'Date (JD) & Cycle & $RV_1$ & $(O\\!-\\!C)_1$ ','& $RV_2$ & $(O\\!
     +-\\!C)_2$ \\\\',/,' 2400000+ & & km s$^{-1}$ & km s$^{-1}$ &',' km
     + s$^{-1}$ & km s$^{-1}$\\\\',/,'\\hline')
      DO80028I=1,NN
      PHI=(TOBS(I)-T)/PP
      IF(POI1(I).NE.0..AND.POI2(I).NE.0.)THEN
      WRITE(LU_OUT,3020)TOBS(I),PHI,VIT1(I),RES1(I),VIT2(I),RES2(I)
 3020 FORMAT(F12.2,' & ',F8.2,' & ',F8.1,' & ',F8.1,' & ',F8.1,' & ',F8.
     +1,' \\\\')
      ELSEIF(POI1(I).EQ.0.)THEN
      WRITE(LU_OUT,3021)TOBS(I),PHI,VIT2(I),RES2(I)
 3021 FORMAT(F12.2,' & ',F8.2,' & \\nodata & \\nodata & ',F8.1,' & ',F8.
     +1,' \\\\')
      ELSE
      WRITE(LU_OUT,3022)TOBS(I),PHI,VIT1(I),RES1(I)
 3022 FORMAT(F12.2,' & ',F8.2,' & ',F8.1,' & ',F8.1,' & \\nodata & ','\\
     +nodata \\\\')
      ENDIF
80028 CONTINUE
      WRITE(LU_OUT,3011)
 3011 FORMAT('\\hline',/,'\\end{tabular}',/)
      WRITE(LU_OUT,3012)OBJECT(1:INDEX(OBJECT,' ')-1),OBJECT(1:INDEX(OBJ
     +ECT,' ')-1)
 3012 FORMAT('\\bigskip',/,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     +%%%%%%%',/,'\\begin{figure}[h]',/,'\\centering\\includegraphics*[w
     +idth=10cm]{',A,'_per.ps}',/,'\\caption{',A,': periodogram.}',/,'\\
     +end{figure}',/)
      WRITE(LU_OUT,3013)
 3013 FORMAT(/,'\\end{document}')
      CLOSE(LU_OUT)
      RETURN
      END
