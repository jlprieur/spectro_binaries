C************************************************************************
C     AMELIORATION DES ELEMENTS ORBITAUX                      *******
C     POUR LES BINAIRES A UN SPECTRE ET 3 CORPS               * BS4 *
C                                                             *******
C     Comme BS1 mais avec des parametres supplementaires:
C        PP3: periode du 3eme corps
C        TT3: T0 passage au periastre du 3eme corps
C        K3
C     
C Historique:
C  Cree a partir de BS3 (version 15/05/2002) 
C
C JLP
C Version 16/05/2002 
C************************************************************************
      PROGRAM BS4
      IMPLICIT NONE
      INTEGER IDIM,NPAR
      PARAMETER (IDIM=120,NPAR=12)
      REAL TOBS(IDIM),VIT1(IDIM),POI1(IDIM),RES1(IDIM)
      REAL AA(NPAR,200),B(NPAR),DELT(NPAR),ERR(NPAR),EQM(NPAR),CORR(NPAR)
      REAL A(NPAR,NPAR),D(NPAR,NPAR)
      DOUBLE PRECISION PHI,PHI3,PP,T,PP3,T3,TZERO,TZERO3
      INTEGER NP,X1,Y1,X3,Y3,F,DELTA(NPAR),LU_IN,LU_OUT,IT,NN,ITER,NITER
      REAL OMEGA,OMEGA3,EE,EE3,PI,VV,VV3,VIC1,WW1
      REAL K1,K3,SIGMA1,SUMSQ1,Q,S,Q3,S3,SOM_POI
      REAL FDEM1,FDEM3,A1SINI,A3SINI,SOM,SOMB,SOMD,SOME
      REAL DIST,VZERO,ERR_A1SINI,ERR_A3SINI,ERR_FDEM1,ERR_FDEM3
      INTEGER I,J,K,JJ,IOPT
      LOGICAL ORBIT1_IS_FIXED,ORBIT3_IS_FIXED
      CHARACTER OBJECT*20,OUT_RESULTS*60,OUT_LATEX*60
      PI=3.1415926535

      WRITE(6,*)' Titre (sans blancs, car utilise comme prefixe):'
      READ(5,10) OBJECT
10    FORMAT(A)

C
C Ouverture du fichier en entree nomme VR_BS4.DAT;
      LU_IN=7
      OPEN(LU_IN,FILE='VR_BS4.DAT',STATUS='OLD')
      WRITE(6,42)
 42   FORMAT(' Menu:',/,' 0 = fit both orbits',/,
     1 ' 1 = fit orbit1 only',/,' 3 = fit orbit3 only')
      READ(5,*) I
      IF(I.EQ.1)THEN
       ORBIT1_IS_FIXED=.FALSE.
       ORBIT3_IS_FIXED=.TRUE.
      ELSE IF(I.EQ.3)THEN
       ORBIT1_IS_FIXED=.TRUE.
       ORBIT3_IS_FIXED=.FALSE.
      ELSE
       ORBIT1_IS_FIXED=.FALSE.
       ORBIT3_IS_FIXED=.FALSE.
      ENDIF
C*******************************************************************************
C Format of VR_BS4.DAT
C Periode (jours)
C Passage au periastre ou noeud ascendant si circulaire
C angle omega (noeuds/absides ou 0 si circulaire)
C eccentricite
C k1 (semi amplitude de la courbe de vitesses radiales de la primaire en km/s)
C V0 vitesse du centre de gravite du systeme en km/s
C Parametres optionnels: X1=1 si periode fixeee Y1=2 si orbite circulaire
C 14.2081   48682.860 0.        0.        42.5      2.8      1    2
C
C Nouvelle ligne, pour le 3eme corps:
C Periode (jours)
C Passage au periastre ou noeud ascendant si circulaire
C angle omega 
C eccentricite 
C k3 (semi amplitude de la courbe de vitesses radiales de la secondaire en km/s)
C 0 
C 0 
C Parametres optionnels: X1=1 si periode fixeee Y1=2 si orbite circulaire
C 1700.     50701.    0.        0.        2.7        0.        0.      1    2 
C
C Nombre de mesures de vitesses radiales (dans ce fichier)
C Nombre d'iterations (pour converger)
C 42       5
C Date en jour Julien (-2000000), vitesse radiale de la primaire, poids de la mesure
C vitesse radiale de la secondaire, poids de la mesure
C*******************************************************************************
C Creation du fichier en sortie OBJECT.txt
      LU_OUT=8
      OUT_RESULTS=OBJECT(1:INDEX(OBJECT,' ')-1)//'.txt'
      OPEN(LU_OUT,FILE=OUT_RESULTS,STATUS='UNKNOWN')
C Format des elements provisoires de l'orbite:
 1000 FORMAT(F10.6,F10.3,2F10.5,F10.2,F6.2,1X,2I1)
   88 READ(LU_IN,1000)PP,TZERO,OMEGA,EE,K1,VZERO,X1,Y1
      WRITE(6,1001)PP,TZERO,OMEGA*180./PI,EE,K1,VZERO
1001  FORMAT(' P=',F12.7,' T0=',F10.3,' Omega=',F9.3,'(deg) E=',F9.3,
     1       ' K1=',F9.2,'V0=',F9.2)
      T=TZERO
      IF(X1.NE.1)X1=0
      IF(Y1.NE.1.AND.Y1.NE.2)Y1=0
      WRITE(6,*)'Option: X1=',X1,' Y1=',Y1
C X1=1: periode fixe
C Y1=1: ellipticite fixe
C Y1=2: orbite circulaire
      IF(Y1.EQ.2)THEN
       EE=0.
       OMEGA=0.
      ENDIF
      IF(EE.GE.1.OR.EE3.GE.1)THEN
        WRITE(6,*)'Erreur fatale: EE=',EE,' et EE3=',EE3
      ENDIF

      READ(LU_IN,1000)PP3,TZERO3,OMEGA3,EE3,K3,WW1,X3,Y3
      WRITE(6,1002)PP3,TZERO3,OMEGA3*180./PI,EE3,K3
1002  FORMAT(' P3=',F12.5,' T3=',F10.3,' Omega3=',F9.3,'(deg) E3=',F9.3,
     1       ' K3=',F9.2)
      T3=TZERO3
      IF(X3.NE.1)X3=0
      IF(Y3.NE.1.AND.Y3.NE.2)Y3=0
      WRITE(6,*)'Option: X3=',X3,' Y3=',Y3
C X3=1: periode fixe
C Y3=1: ellipticite fixe
C Y3=2: orbite circulaire
      IF(Y3.EQ.2)THEN
       EE3=0.
       OMEGA3=0.
      ENDIF

C Nombre de mesures, nombre d'iterations:
881   READ(LU_IN,*)NN,NITER
      WRITE(6,*) 'Nombre de mesures:',NN,' Nombre d''iterations',NITER

      DO I=1,NN
C      READ(LU_IN,1200)TOBS(I),VIT1(I),POI1(I)
C 1200 FORMAT(F9.3,F7.1,F6.2)
CCCCCCCCCCCCCCCCC
C Sans poids:
       READ(LU_IN,*)TOBS(I),VIT1(I)
       POI1(I)=1.0
CCCCCCCCCCCCCCCCC
C Avec poids:
C       READ(LU_IN,*)TOBS(I),VIT1(I),POI1(I)
       IF(POI1(I).LE.0.)THEN
        WRITE(6,*)' Erreur fatale: poids aberrant: ',POI1(I)
        STOP
       ENDIF
      END DO
      CLOSE(LU_IN)

C Calcul de la somme des poids:
      SOM_POI=0.
      DO I=1,NN
        SOM_POI=SOM_POI+POI1(I)
      END DO

C Nombre de parametres
C 1=PP 2=T 3=OMEGA 4=EE 5=K1 6=VZERO 7=PP3 8=T3 9=OMEGA3 10=EE3 11=K3
C Before:
C 1=PP 2=T 3=OMEGA 4=EE 5=K1 6=K2 7=VZERO 8=PP3 9=T3 10=OMEGA3 11=EE3 12=K3
C X=1: periode fixe
C Y=1: ellipticite fixe
C Y=2: orbite circulaire

      IF(ORBIT3_IS_FIXED)THEN
        NP=6-X1-Y1
      ELSE
        NP=11-X1-Y1-X3-Y3
      ENDIF 

      WRITE(LU_OUT,1500)PP,T,OMEGA,EE,K1,PP3,T3,K3,VZERO
1500  FORMAT(13X,'ELEMENTS PROVISOIRES',/,16X,'P =',F15.7,/,16X,'T =',
     1 F11.3,/,12X,'OMEGA =',F13.5,/,16X,'E =',F13.5,/,16X,'K1 =',
     1 F10.2,/,16X,'P3 =',F15.7,/,16X,'T3 =',F11.3,/,
     1 16X,'K3 =',F10.2,/,16X,'V0 =',F10.2,/,/)
C
      DIST=1.
      DO I=1,NP
        DO J=1,NN
         AA(I,J)=0
        END DO
      END DO

C Boucle principale:
      ITER=0
      DO IT=0,NITER
C Remise a zero:
        DO J=1,NP
         B(J)=0
        END DO
        DO J=1,NP
          DO K=1,NP
            A(J,K)=0
          END DO
        END DO
C Sortie des RES(I) = O-C
C seulement pour les elements provisoires et pour la derniere iteration:
      IF(ITER.EQ.0.OR.ITER.EQ.NITER)THEN
        WRITE(LU_OUT,1600)
1600    FORMAT('OBS.',1X,'DATE (JJ)',2X,'PHASE',1X,'V_OBS1',2X,'V_CALC1',
     1  1X,'(O-C)1',1X,'POIDS1')
      ENDIF
      DO 53 I=1,NN
      PHI=DMOD((TOBS(I)-T)/PP,1.D0)
      IF(PHI.LT.0)PHI=1.+PHI
      PHI3=DMOD((TOBS(I)-T3)/PP3,1.D0)
      IF(PHI3.LT.0)PHI3=1.+PHI3
C
C Calcul des coefficients des equations normales
C
C Calcul des residus:
      CALL VITESSE_CALCULEE_BS4(PHI,PHI3,VIC1,VV,VV3,
     1                          EE,EE3,OMEGA,OMEGA3,K1,K3,VZERO)
      RES1(I)=VIT1(I)-VIC1
C Sortie des RES(I) = O-C
C seulement pour les elements provisoires et pour la derniere iteration:
      IF(ITER.EQ.0.OR.ITER.EQ.NITER)THEN
          WRITE(LU_OUT,1400)I,TOBS(I),PHI,VIT1(I),VIC1,
     1                   RES1(I),POI1(I)
1400      FORMAT(I3,2X,F9.3,2X,F4.3,3(2X,F7.2),3X,F5.2)
      ENDIF
C Aiguillage vers la sortie de la boucle secondaire:
      IF(ITER.EQ.NITER.OR.DIST.LT.1.E-7)GOTO 53

C Calcul des elements de la matrice AA(NP=11-X1-Y1-X3-Y3,NN)

C 1=PP 2=T 3=OMEGA 4=EE 5=K1 6=VZERO 7=PP3 8=T3 9=OMEGA3 10=EE3 11=K3
C Before:
C 1=PP 2=T 3=OMEGA 4=EE 5=K1 6=K2 7=VZERO 8=PP3 9=T3 10=OMEGA3 11=EE3 12=K3
C Calcul avec parametres precedents:
        S=SIN(VV+OMEGA)
        Q=S*(1+EE*COS(VV))**2/(1-EE**2)**1.5
        S3=SIN(VV3+OMEGA3)
        Q3=S3*(1+EE3*COS(VV3))**2/(1-EE3**2)**1.5
C Correction de (-2.*PI/PP**2) par rapport a p83 (coeff "D"),
C puisque qu'on utilise la periode PP au lieu de n (moyen mouvement)
C X=1: periode fixe
        IF(X1.EQ.0)THEN
C Coeff "D" corrige: 
          AA(1,I)=(2.*PI/PP**2)*(TOBS(I)-T)*K1*Q
        ENDIF
C Coeff "F" p83
        AA(2-X1,I)=K1*Q*(2.*PI/PP)
C Y=0: ellipticite libre
C Y=1: ellipticite fixe
C Y=2: orbite circulaire
        IF(Y1.NE.2)THEN
C Coeff "B" p 83 (omega)
          AA(3-X1,I)=-K1*(EE*SIN(OMEGA)+S)
          IF(Y1.EQ.0)THEN
C Coeff "C" p 83 (e)
            AA(4-X1,I)=K1*(COS(OMEGA)-S*SIN(VV)*(2+EE*COS(VV))/(1-EE**2))
          ENDIF
        ENDIF
C 1=PP 2=T 3=OMEGA 4=EE 5=K1 6=VZERO 7=PP3 8=T3 9=OMEGA3 10=EE3 11=K3
C Before:
C 1=PP 2=T 3=OMEGA 4=EE 5=K1 6=K2 7=VZERO 8=PP3 9=T3 10=OMEGA3 11=EE3 12=K3
C Coeff "A" p 83 (k)
        AA(5-X1-Y1,I)=EE*COS(OMEGA)+COS(VV+OMEGA)
        AA(6-X1-Y1,I)=1.
C X3=1 si periode P3 fixe
        IF(X3.EQ.0)THEN
C Coeff "D" corrige: (P3) 
           AA(7-X1-Y1,I)=(2.*PI/PP3**2)*(TOBS(I)-T3)*K3*Q3
        ENDIF
C Coeff "F" p83 (T3)
        AA(8-X1-Y1-X3,I)=K3*Q3*(2.*PI/PP3)
        IF(Y3.NE.2)THEN
C Coeff "B" p 83 (omega)
          AA(9-X1-Y1-X3,I)=-K3*(EE3*SIN(OMEGA3)+S3)
          IF(Y3.EQ.0)THEN
C Coeff "C" p 83 (e)
            AA(10-X1-Y1-X3,I)=K3*(COS(OMEGA3)-S3*SIN(VV3)
     1                     *(2.+EE3*COS(VV3))/(1.-EE3**2))
          ENDIF
        ENDIF
C Coeff "A" p 83 (K3)
        AA(11-X1-Y1-X3-Y3,I)=EE3*COS(OMEGA3)+COS(VV3+OMEGA3)

C Calcul du 2eme membre de l'equation normale
      DO J=1,NP
        B(J)=B(J)+POI1(I)*AA(J,I)*RES1(I)
      END DO
C Calcul du 1er membre de l'equation normale
      DO J=1,NP
        DO K=J,NP
         A(J,K)=A(J,K)+POI1(I)*AA(J,I)*AA(K,I)
        END DO
      END DO
C Fin de la boucle secondaire
   53 CONTINUE
C Aiguillage vers la sortie de la boucle principale: 
      IF(ITER.EQ.NITER.OR.DIST.LT.1.E-7)GOTO 81

C Symmetrise la matrice A
C (puisque par construction elle est symmetrique, mais
C on n'avait pas encore affecte l'autre partie de cette matrice)
      DO J=2,NP
        JJ=J-1
        DO K=1,JJ
          A(J,K)=A(K,J)
        END DO
      END DO

C Inversion de la matrice:
C en sortie: D=matrice inverse
      CALL MAT_INV12(A,D,NP)

C
C Nouveaux elements - erreurs

C Calcul des corrections DELT
      DO I=1,NP
        DELT(I)=0
      END DO
      DO 6 I=1,NP
      DO 6 J=1,NP
      DELT(I)=DELT(I)+D(I,J)*B(J)
    6 CONTINUE

C Calcul des erreurs EQM
      SOMB=0
      DO 7 I=1,NN
        SOMB=SOMB+POI1(I)*RES1(I)**2
    7 CONTINUE
C BOF JLP2004: normalisation des poids a 1:
      SOMB=SOMB*NN/SOM_POI
C Ajout d'un terme tenant compte des corrections actuelles:
C (Il s'annule lorsqu'on converge, et peut conduire
C a reduire les erreurs s'il est positif,
C alors je le neutralise):
      SOMD=0
C      DO I=1,NP
C        SOMD=SOMD+B(I)*DELT(I)
C      END DO
C EOF JLP2004
      SOME=SOMB-SOMD
      SOM=SQRT(SOME/(NN-NP))
      DO 80 I=1,NP
        EQM(I)=SOM*SQRT(D(I,I))
   80 CONTINUE

C Nouveaux elements:
      DIST=0
      F=0
      DO I=1,12
        DELTA(I)=0
      END DO
C X1=1 periode fixe
      DELTA(1)=X1
C Y1=1 ellipticite fixe 
      DELTA(4)=Y1
C Y1=2 orbite circulaire 
      IF(Y1.EQ.2)THEN
        DELTA(3)=1
        DELTA(4)=1
      ENDIF
C X3=1 periode fixe
      DELTA(8)=X3
C Y3=1 ellipticite fixe 
      DELTA(10)=Y3
C Y3=2 orbite circulaire 
      IF(Y3.EQ.2)THEN
        DELTA(9)=1
        DELTA(10)=1
      ENDIF

C 1=PP 2=T 3=OMEGA 4=EE 5=K1 6=VZERO 7=PP3 8=T3 9=OMEGA3 10=EE3 11=K3
C Before:
C 1=PP 2=T 3=OMEGA 4=EE 5=K1 6=K2 7=VZERO 8=PP3 9=T3 10=OMEGA3 11=EE3 12=K3
C Remise en forme pour les 11 parametres:
      DO I=1,11
        F=F+DELTA(I)
        J=I-F+DELTA(I)
        ERR(I)=(1-DELTA(I))*EQM(J)
        CORR(I)=(1-DELTA(I))*DELT(J)
        DIST=DIST+ABS(CORR(I))
      END DO

      IF(.NOT.ORBIT1_IS_FIXED)THEN
        PP=PP+CORR(1)
        T=TZERO+DMOD(T-TZERO+CORR(2),PP)
        OMEGA=AMOD(OMEGA+CORR(3),2.*PI)
        IF(OMEGA.LT.0)OMEGA=OMEGA+2.*PI
C JLP2002: j'ajoute une contrainte sur la positivite de EE:
        IF(EE+CORR(4).LT.0.)THEN
         EE=EE/2.
C Surtout ne pas rajouter OMEGA=OMEGA+PI, comme on pourrait le penser!
C (diverge avec HD100054B_fixed)
        ELSE 
         EE=EE+CORR(4)
        ENDIF
        IF(EE.GE.1.)EE=0.9
        K1=K1+CORR(5)
        IF(K1.LT.0)K1=-K1
        VZERO=VZERO+CORR(6)
      ENDIF

      IF(.NOT.ORBIT3_IS_FIXED)THEN
        PP3=PP3+CORR(7)
        T3=TZERO3+DMOD(T3-TZERO3+CORR(8),PP3)
        OMEGA3=AMOD(OMEGA3+CORR(9),2.*PI)
        IF(OMEGA3.LT.0)OMEGA3=OMEGA3+2.*PI
C JLP2002: j'ajoute une contrainte sur la positivite de EE3:
        IF(EE3+CORR(10).LT.0.)THEN
         EE3=EE3/2.
C Surtout ne pas rajouter OMEGA3=OMEGA3+PI, comme on pourrait le penser!
C (diverge avec HD100054B_fixed)
        ELSE 
         EE3=EE3+CORR(10)
        ENDIF
        IF(EE3.GE.1.)EE3=0.9
        K3=K3+CORR(11)
        IF(K3.LT.0)K3=-K3
      ENDIF

      A1SINI=0
      FDEM1=0
      IF(ABS(EE).LT.1.)THEN
        A1SINI=43200.*PP*K1*SQRT(1-EE**2)/PI
        FDEM1=0.00000010385*K1**3*PP*(1-EE**2)**1.5
      ENDIF
      ERR_A1SINI=A1SINI*(ERR(1)/PP + ERR(5)/K1 + ERR(4)*EE/(1.-EE**2))
      ERR_FDEM1=FDEM1*(3.*EE*ERR(4)/(1.-EE**2) + ERR(1)/PP
     1         + 3.*ERR(5)/K1)
      A3SINI=0
      FDEM3=0
      IF(ABS(EE3).LT.1.)THEN
        A3SINI=43200.*PP3*K3*SQRT(1.-EE3**2)/PI
        FDEM3=0.00000010385*PP3*(K3**3)*((1.-EE3**2)**1.5)
      ENDIF
C 1=PP 2=T 3=OMEGA 4=EE 5=K1 6=VZERO 7=PP3 8=T3 9=OMEGA3 10=EE3 11=K3
      ERR_A3SINI=A3SINI*(ERR(7)/PP3 + ERR(11)/K3 
     1           + ERR(10)*EE3/(1.-EE3**2))
      ERR_FDEM3=FDEM3*(3.*EE3*ERR(10)/(1.-EE3**2) + ERR(7)/PP3
     1         + 3.*ERR(11)/K3)
      ITER=ITER+1
      WRITE(LU_OUT,2200)ITER
2200  FORMAT('*******',I2,'eme iteration ********')
      IF(X1.EQ.1)WRITE(LU_OUT,1700)
1700  FORMAT('Periode fixee pour la primaire')
      IF(Y1.EQ.1)WRITE(LU_OUT,1800)
1800  FORMAT('Excentricite fixee pour la primaire')
      IF(Y1.EQ.2)WRITE(LU_OUT,1801)
1801  FORMAT('Orbite circulaire pour la primaire')
      IF(X3.EQ.1)WRITE(LU_OUT,1703)
1703  FORMAT('Periode fixee pour le 3eme corps')
      IF(Y3.EQ.1)WRITE(LU_OUT,1803)
1803  FORMAT('Excentricite fixee pour le 3eme corps')
      IF(Y3.EQ.2)WRITE(LU_OUT,1804)
1804  FORMAT('Orbite circulaire pour le 3eme corps')
      WRITE(LU_OUT,1300)CORR(1),PP,ERR(1),CORR(2),T,ERR(2),
     1 CORR(3)*180./PI,OMEGA*180./PI,ERR(3)*180./PI,
     1 CORR(3),OMEGA,ERR(3),
     1 CORR(4),EE,ERR(4),CORR(5),K1,ERR(5),
     1 CORR(6),VZERO,ERR(6),A1SINI,ERR_A1SINI,FDEM1,ERR_FDEM1
1300  FORMAT(5X,'CORRECTIONS',15X,'NOUVEAUX ELEMENTS',
     1 15X,'ERREURS STD.',/,
     1 3X,E9.3,17X,'P =',F15.7,13X,E9.3,/,
     1 3X,E9.3,17X,'T =',F11.3,17X,E9.3,/,
     1 3X,E9.3,13X,'OMEGA =',F13.5,15X,E9.3,4X,'(deg.)',/,
     1 3X,E9.3,13X,'OMEGA =',F13.5,15X,E9.3,4X,'(rad.)',/,
     1 3X,E9.3,17X,'E =',F13.5,15X,E9.3,/,
     1 3X,E9.3,17X,'K =',F10.2,18X,E9.3,/,
     1 3X,E9.3,16X,'V0 =',F10.2,18X,E9.3,/,
     1 24X,'A.SINI =',6X,E10.4,12X,E9.3,/,
     1 26X,'F(M) =',6X,E10.4,12X,E9.3,/)
C 1=PP 2=T 3=OMEGA 4=EE 5=K1 6=VZERO 7=PP3 8=T3 9=OMEGA3 10=EE3 11=K3
C Before:
C 1=PP 2=T 3=OMEGA 4=EE 5=K1 6=K2 7=VZERO 8=PP3 9=T3 10=OMEGA3 11=EE3 12=K3
      WRITE(LU_OUT,2001)CORR(7),PP3,ERR(7),CORR(8),T3,ERR(8),
     1 CORR(9)*180./PI,OMEGA3*180./PI,ERR(9)*180./PI,
     1 CORR(9),OMEGA3,ERR(9),
     1 CORR(10),EE3,ERR(10),CORR(11),K3,ERR(11),
     1 A3SINI,ERR_A3SINI,FDEM3,ERR_FDEM3
2001  FORMAT(
     1 3X,E9.3,16X,'P3 =',F15.7,13X,E9.3,/,
     1 3X,E9.3,16X,'T3 =',F11.3,17X,E9.3,/,
     1 3X,E9.3,12X,'OMEGA3 =',F13.5,15X,E9.3,4X,'(deg.)',/,
     1 3X,E9.3,12X,'OMEGA3 =',F13.5,15X,E9.3,4X,'(rad.)',/,
     1 3X,E9.3,16X,'E3 =',F13.5,15X,E9.3,/,
     1 3X,E9.3,16X,'K3 =',F13.5,15X,E9.3,/,
     1 23X,'A3.SINI =',6X,E10.4,12X,E9.3,/,
     1 26X,'F(M) =',6X,E10.4,12X,E9.3,/,/)
C Erreur grave
      IF(ABS(EE).GE.1.)THEN
          WRITE(LU_OUT,2300)
          WRITE(6,2300)
2300      FORMAT(5X,'ABS(EE).GE.1, Poursuite du calcul impossible!')
          CLOSE(LU_OUT)
          STOP
        ENDIF
C Fin de la boucle principale:
      END DO
C Sortie normale:
  81  CONTINUE 
C Visualisation des courbes:
      WRITE(6,*) ' Courbe: 0=aucune  1=courte periode'
      WRITE(6,*) ' 2=longue periode  3=residus  4=toutes ?'
      READ(5,*) IOPT
      IF(IOPT.EQ.1)THEN
C Courte period:
        CALL JLP_PLOT_BS4(TOBS,VIT1,POI1,
     1                  T,T3,PP,PP3,NN,EE,EE3,OMEGA,OMEGA3,
     1                  K1,K3,VZERO,OBJECT)
      ELSE IF(IOPT.EQ.2)THEN
C Longue periode:
        CALL JLP_PLOT_BS4_3(TOBS,VIT1,POI1,
     1                    T,T3,PP,PP3,NN,EE,EE3,OMEGA,OMEGA3,
     1                    K1,K3,VZERO,OBJECT) 
      ELSE IF(IOPT.EQ.3)THEN 
C Residus:
        CALL JLP_PLOT_RES(TOBS,RES1,POI1,NN,OBJECT)
      ELSE IF(IOPT.EQ.4.OR.IOPT.EQ.5)THEN 
        CALL JLP_PLOT_BS4(TOBS,VIT1,POI1,
     1                  T,T3,PP,PP3,NN,EE,EE3,OMEGA,OMEGA3,
     1                  K1,K3,VZERO,OBJECT)
        CALL JLP_PLOT_BS4_3(TOBS,VIT1,POI1,
     1                    T,T3,PP,PP3,NN,EE,EE3,OMEGA,OMEGA3,
     1                    K1,K3,VZERO,OBJECT) 
        CALL JLP_PLOT_RES(TOBS,RES1,POI1,NN,OBJECT)
      ENDIF
C RMS des residus:
      CALL RMS_RESIDUS(RES1,POI1,NN,SIGMA1,SUMSQ1)
      WRITE(6,3000)SIGMA1,SUMSQ1
      WRITE(LU_OUT,3000)SIGMA1,SUMSQ1
3000  FORMAT(/,'Sigma O-C: Etoile primaire (1)   ',F5.2,' km/s',
     1       /,' Sum of squares: ',E12.5,/)
C 1=PP 2=T 3=OMEGA 4=EE 5=K1 6=VZERO 7=PP3 8=T3 9=OMEGA3 10=EE3 11=K3
C Latex format:
      WRITE(LU_OUT,3001)PP,T,OMEGA*180/3.14159,EE,K1,VZERO,
     1                  A1SINI/1.e6,FDEM1,SIGMA1
      WRITE(LU_OUT,3002)ERR(1),ERR(2),ERR(3)*180/3.14159,ERR(4),
     1                  ERR(5),ERR(6),ERR_A1SINI/1.e6,ERR_FDEM1
3001  FORMAT(' & ',F10.3,' & ',F12.2,' & ',F5.1,' & ',F8.3,
     1       ' & ',F8.2,' & ',F8.2,' & ',F8.2,
     1       ' & ',F8.3,' & ',F5.2,' \\\\')
3002  FORMAT(' & $\\pm',F8.3,'$ & $\\pm',F8.2,'$ & $\\pm',
     1       F5.1,'$ & $\\pm',F8.3,'$ & $\\pm',F8.2,'$ & $\\pm',
     1       F8.2,'$ & $\\pm',F8.2,'$ & $\\pm',F8.3,'$ & \\\\')
      WRITE(LU_OUT,3001)PP3,T3,OMEGA3*180/3.14159,EE3,K3,VZERO,
     1                  A3SINI/1.e6,FDEM3,SIGMA1
      WRITE(LU_OUT,3002)ERR(7),ERR(8),ERR(9)*180/3.14159,ERR(10),
     1                  ERR(11),ERR(6),ERR_A3SINI/1.e6,ERR_FDEM3
      CLOSE(LU_OUT)
C Latex format:
      OUT_LATEX=OBJECT(1:INDEX(OBJECT,' ')-1)//'.tex'
      CALL SB4_TO_LATEX(PP,PP3,T,T3,OMEGA,OMEGA3,EE,EE3,
     1   K1,K3,VZERO,A1SINI,A3SINI,
     1   FDEM1,FDEM3,SIGMA1,
     1   ERR,ERR_A1SINI,ERR_A3SINI,
     1   ERR_FDEM1,ERR_FDEM3,TOBS,VIT1,
     1   RES1,POI1,NN,
     1   OUT_LATEX,OBJECT)
C Information sur l'ecran:
      WRITE(6,1300)CORR(1),PP,ERR(1),CORR(2),T,ERR(2),
     1 CORR(3)*180./PI,OMEGA*180./PI,ERR(3)*180./PI,
     1 CORR(3),OMEGA,ERR(3),
     1 CORR(4),EE,ERR(4),CORR(5),K1,ERR(5),CORR(6),
     1 VZERO,ERR(6),A1SINI,ERR_A1SINI,FDEM1,ERR_FDEM1
      WRITE(6,2001)CORR(7),PP3,ERR(7),CORR(8),T3,ERR(8),
     1 CORR(9)*180./PI,OMEGA3*180./PI,ERR(9)*180./PI,
     1 CORR(9),OMEGA3,ERR(9),
     1 CORR(10),EE3,ERR(10),CORR(11),
     1 K3,ERR(11),A3SINI,ERR_A3SINI,FDEM3,ERR_FDEM3
      WRITE(6,*)' Fichiers de sortie: ',
     1         OUT_RESULTS(1:INDEX(OUT_RESULTS,' ')-1),', ',
     1         OUT_LATEX(1:INDEX(OUT_LATEX,' ')-1),', '
      IF(IOPT.EQ.4.OR.IOPT.EQ.5)THEN
        WRITE(6,*)' et graphes dans: ',
     1    OBJECT(1:INDEX(OBJECT,' ')-1),'_short.ps,',
     1    OBJECT(1:INDEX(OBJECT,' ')-1),'_long.ps'
        WRITE(6,*)'    et ',OBJECT(1:INDEX(OBJECT,' ')-1),'_resi.ps,'
      ELSE IF(IOPT.EQ.1)THEN
        WRITE(6,*)' et graphe dans: ',OBJECT(1:INDEX(OBJECT,' ')-1),
     1            '_short.ps'
      ELSE IF(IOPT.EQ.2)THEN
        WRITE(6,*)' et graphe dans: ',OBJECT(1:INDEX(OBJECT,' ')-1),
     1            '_long.ps'
      ELSE IF(IOPT.EQ.3)THEN
        WRITE(6,*)' et graphe dans: ',OBJECT(1:INDEX(OBJECT,' ')-1),
     1            '_resi.ps'
      ENDIF
      STOP
      END
C*********************************************************************
C Pour visualiser les points et la courbe (courte periode)
C Mode simple
C*********************************************************************
      SUBROUTINE JLP_PLOT_BS44(TOBS,VIT1,POI1,
     1                        T,T3,PP,PP3,NN,EE,EE3,OMEGA,OMEGA3,
     1                        K1,K3,VZERO,OBJECT)
      IMPLICIT NONE
      INTEGER DIM,KCUR
      PARAMETER (DIM=1024,KCUR=4)
      REAL TOBS(120),VIT1(120),POI1(120)
      DOUBLE PRECISION T,T3,PP,PP3,PHI,PHI3
      REAL PI,EE,EE3,OMEGA,OMEGA3,K1,K3,VZERO,VV,VV3
      REAL VIC1
      INTEGER KK,NPTS(KCUR),I,NN,I1,I2
      REAL XPLOT(DIM,KCUR),YPLOT(DIM,KCUR)
      CHARACTER CHAR1*30,CHAR2*30,CHAR3*40,OBJECT*20
      CHARACTER PLOTDEV*40,IN_FILE*40,IN_COMMENTS*80
      CHARACTER*4 NCHAR(KCUR)

C On charge les mesures de poids 1 en K=1
      I1 = 1
      I2 = 1
      DO I=1,NN
        PHI=DMOD((TOBS(I)-T)/PP,1.D0)
        IF(PHI.LT.0)PHI=1.+PHI
        PHI3=DMOD((TOBS(I)-T3)/PP3,1.D0)
        IF(PHI3.LT.0)PHI3=1.+PHI3
C Calcule la correction de vitesse: K3*COS(OMEGA3+VV3)
C Pour cela on appelle VITESSE_CALCULEE_BS4 avec K1=VZERO=0.
        CALL VITESSE_CALCULEE_BS4(PHI,PHI3,VIC1,VV,VV3,
     1                            EE,EE3,OMEGA,OMEGA3,0.,K3,0.)
        IF(POI1(I).EQ.1)THEN
           XPLOT(I1,1)=PHI
           YPLOT(I1,1)=VIT1(I)-VIC1
           I1=I1+1
        ELSE
           XPLOT(I2,2)=PHI
           YPLOT(I2,2)=VIT1(I)-VIC1
           I2=I2+1
        ENDIF
      END DO
C Cercles pleins (primaire) 
      NCHAR(1)='913'
      NPTS(1)=I1-1

C Croix 
      NCHAR(2)='513'
      NPTS(2)=I2-1

C On charge la courbe en K=3
      PI=3.1415926535
      DO I=1,1000
        PHI=FLOAT(I)/1000.
        PHI3=0.
        CALL VITESSE_CALCULEE_BS4(PHI,PHI3,VIC1,VV,VV3,
     1                            EE,EE3,OMEGA,OMEGA3,K1,0.,VZERO)
        XPLOT(I,3)=PHI
        YPLOT(I,3)=VIC1
      END DO
C
      NCHAR(3)='L0'
      NPTS(3)=1000
      KK=3
      WRITE(PLOTDEV,103)OBJECT(1:INDEX(OBJECT,' ')-1),'_short.ps'
103   FORMAT('&landscape/',A,A)
      WRITE(CHAR1,11) 'Phase',CHAR(0)
      WRITE(CHAR2,11) 'Radial velocity (km/s)',CHAR(0)
      WRITE(CHAR3,11) OBJECT,CHAR(0)
11    FORMAT(A,A1)
      IN_FILE=' '
      IN_COMMENTS=' '
      CALL NEWPLOT(XPLOT,YPLOT,NPTS,DIM,KK,CHAR1,CHAR2,
     1  CHAR3,NCHAR,PLOTDEV,IN_FILE,IN_COMMENTS)

      RETURN
      END
C*********************************************************************
C Pour visualiser les points et la courbe (courte periode)
C Mode avec cadre plus large (cycle de -0.1 à 1.1) 
C*********************************************************************
      SUBROUTINE JLP_PLOT_BS4(TOBS,VIT1,POI1,
     1                        T,T3,PP,PP3,NN,EE,EE3,OMEGA,OMEGA3,
     1                        K1,K3,VZERO,OBJECT)
      IMPLICIT NONE
      INTEGER DIM,KCUR
      PARAMETER (DIM=1024,KCUR=4)
      REAL TOBS(120),VIT1(120),POI1(120)
      DOUBLE PRECISION T,T3,PP,PP3,PHI,PHI3,PHI_MINI
      REAL PI,EE,EE3,OMEGA,OMEGA3,K1,K3,VZERO,VV,VV3
      REAL VIC1
      INTEGER KK,NPTS(KCUR),I,NN,I1,I2,ITER
      REAL XPLOT(DIM,KCUR),YPLOT(DIM,KCUR)
      CHARACTER CHAR1*30,CHAR2*30,CHAR3*40,OBJECT*20
      CHARACTER PLOTDEV*40,IN_FILE*40,IN_COMMENTS*80
      CHARACTER*4 NCHAR(KCUR)

C JLP2006/ Modification suggested by Imbert/Griffin,
C -0.1,1.1
      PHI_MINI=0.1

C On charge les mesures de poids 1 en K=1
      I1 = 1
      I2 = 1
      DO ITER=1,3
      DO I=1,NN
        PHI=DMOD((TOBS(I)-T)/PP,1.D0)
        IF(PHI.LT.0)PHI=1.+PHI
        PHI3=DMOD((TOBS(I)-T3)/PP3,1.D0)
        IF(PHI3.LT.0)PHI3=1.+PHI3
C Calcule la correction de vitesse: K3*COS(OMEGA3+VV3)
C Pour cela on appelle VITESSE_CALCULEE_BS4 avec K1=VZERO=0.
        CALL VITESSE_CALCULEE_BS4(PHI,PHI3,VIC1,VV,VV3,
     1                            EE,EE3,OMEGA,OMEGA3,0.,K3,0.)
       IF(ITER.EQ.1)THEN
          IF(POI1(I).EQ.1)THEN
             XPLOT(I1,1)=PHI
             YPLOT(I1,1)=VIT1(I)-VIC1
             I1=I1+1
          ELSE
             XPLOT(I2,2)=PHI
             YPLOT(I2,2)=VIT1(I)-VIC1
             I2=I2+1
          ENDIF
C JLP2006/ Modification suggested by Imbert/Griffin,
C to display the peaks of highly eccentric orbits in one piece:
       ELSEIF(ITER.EQ.2)THEN
         IF(PHI .LT. PHI_MINI) THEN
          PHI = PHI + 1.
          IF(POI1(I).EQ.1)THEN
             XPLOT(I1,1)=PHI
             YPLOT(I1,1)=VIT1(I)-VIC1
             I1=I1+1
          ELSE
             XPLOT(I2,2)=PHI
             YPLOT(I2,2)=VIT1(I)-VIC1
             I2=I2+1
          ENDIF
         ENDIF
C JLP2006/ Modification suggested by Griffin,
C to display the peaks of highly eccentric orbits in one piece:
        ELSEIF(ITER.EQ.3)THEN
         IF(PHI .GT. (1.-PHI_MINI)) THEN
          PHI = PHI - 1.
          IF(POI1(I).EQ.1)THEN
             XPLOT(I1,1)=PHI
             YPLOT(I1,1)=VIT1(I)-VIC1
             I1=I1+1
          ELSE
             XPLOT(I2,2)=PHI
             YPLOT(I2,2)=VIT1(I)-VIC1
             I2=I2+1
          ENDIF
         ENDIF
        ENDIF
C EOF loop on I:
      END DO
C EOF loop on ITER:
      END DO
C Cercles pleins (primaire) 
      NCHAR(1)='913'
      NPTS(1)=I1-1

C Croix 
      NCHAR(2)='513'
      NPTS(2)=I2-1

C On charge la courbe en K=3
      PI=3.1415926535
      DO I=1,1000
C Phase PHI entre -PHI_MINI et 1+PHI_MINI:
        PHI = -PHI_MINI + (1+2.*PHI_MINI)*REAL(I)/1000.
        PHI3=0.
        CALL VITESSE_CALCULEE_BS4(PHI,PHI3,VIC1,VV,VV3,
     1                            EE,EE3,OMEGA,OMEGA3,K1,0.,VZERO)
        XPLOT(I,3)=PHI
        YPLOT(I,3)=VIC1
      END DO
C
      NCHAR(3)='L0'
      NPTS(3)=1000
      KK=3
      WRITE(PLOTDEV,103)OBJECT(1:INDEX(OBJECT,' ')-1),'_short.ps'
103   FORMAT('&landscape/',A,A)
      WRITE(CHAR1,11) 'Phase',CHAR(0)
      WRITE(CHAR2,11) 'Radial velocity (km/s)',CHAR(0)
      WRITE(CHAR3,11) OBJECT,CHAR(0)
11    FORMAT(A,A1)
      IN_FILE=' '
      IN_COMMENTS=' '
      CALL NEWPLOT(XPLOT,YPLOT,NPTS,DIM,KK,CHAR1,CHAR2,
     1  CHAR3,NCHAR,PLOTDEV,IN_FILE,IN_COMMENTS)

      RETURN
      END
C*********************************************************************
C Pour visualiser les residus 
C*********************************************************************
      SUBROUTINE JLP_PLOT_RES(TOBS,RES1,POI1,NN,OBJECT)
      IMPLICIT NONE
      INTEGER DIM,KCUR
      PARAMETER (DIM=1024,KCUR=4)
      REAL TOBS(120),RES1(120),POI1(120)
      INTEGER KK,NPTS(KCUR),I,NN,I1,I2
      REAL XPLOT(DIM,KCUR),YPLOT(DIM,KCUR)
      CHARACTER CHAR1*30,CHAR2*30,CHAR3*40,OBJECT*20
      CHARACTER PLOTDEV*40,IN_FILE*40,IN_COMMENTS*80
      CHARACTER*4 NCHAR(KCUR)

C On charge les residus en K=1 et K=2:
C On charge les mesures de poids 1 en K=1
      I1 = 1
      I2 = 1
      DO I=1,NN
        IF(POI1(I).EQ.1.)THEN
           XPLOT(I1,1)=TOBS(I)
           YPLOT(I1,1)=RES1(I)
           I1=I1+1
        ELSE 
           XPLOT(I2,2)=TOBS(I)
           YPLOT(I2,2)=RES1(I)
           I2=I2+1
        ENDIF
      END DO
C Cercles pleins (primaire) 
      NCHAR(1)='913'
      NPTS(1)=I1-1
C Croix  
      NCHAR(2)='513'
      NPTS(2)=I2-1

      KK=2
C For papers: '&LANDSCAPE' (no rotation)
C For work  : '&landscape' (rotation)
      WRITE(PLOTDEV,103)OBJECT(1:INDEX(OBJECT,' ')-1),'_resi.ps'
103   FORMAT('&landscape/',A,A)
      WRITE(CHAR1,11) 'Date',CHAR(0)
      WRITE(CHAR2,11) 'Residuals (km/s)',CHAR(0)
      WRITE(CHAR3,11) OBJECT,CHAR(0)
11    FORMAT(A,A1)
      IN_FILE=' '
      IN_COMMENTS=' '
      CALL NEWPLOT(XPLOT,YPLOT,NPTS,DIM,KK,CHAR1,CHAR2,
     1  CHAR3,NCHAR,PLOTDEV,IN_FILE,IN_COMMENTS)

      RETURN
      END
C*********************************************************************
C Pour visualiser la courbe du 3eme corps
C Mode simple
C*********************************************************************
      SUBROUTINE JLP_PLOT_BS4_33(TOBS,VIT1,POI1,
     1                          T,T3,PP,PP3,NN,EE,EE3,OMEGA,OMEGA3,
     1                          K1,K3,VZERO,OBJECT)
      IMPLICIT NONE
      INTEGER DIM,KCUR
      PARAMETER (DIM=1024,KCUR=4)
      REAL TOBS(120),VIT1(120),POI1(120)
      DOUBLE PRECISION T,T3,PP,PP3,PHI,PHI3
      REAL PI,EE,EE3,OMEGA,OMEGA3,K1,K3,VZERO,VV,VV3,VIC1
      INTEGER KK,NPTS(KCUR),I,NN,I1,I2
      REAL XPLOT(DIM,KCUR),YPLOT(DIM,KCUR)
      CHARACTER CHAR1*30,CHAR2*30,CHAR3*40,OBJECT*20
      CHARACTER PLOTDEV*40,IN_FILE*40,IN_COMMENTS*80
      CHARACTER*4 NCHAR(KCUR)

C On charge les mesures en K=1 et K=2
      I1 = 1
      I2 = 1
      DO I=1,NN
        PHI=DMOD((TOBS(I)-T)/PP,1.D0)
        IF(PHI.LT.0)PHI=1.+PHI
        PHI3=DMOD((TOBS(I)-T3)/PP3,1.D0)
        IF(PHI3.LT.0)PHI3=1.+PHI3
C On calcule la vitesse radiale du centre de gravite
C Pour cela on appelle VITESSE_CALCULEE_BS4 avec K3=0.
        CALL VITESSE_CALCULEE_BS4(PHI,0.D0,VIC1,VV,VV3,
     1                            EE,EE3,OMEGA,OMEGA3,K1,0.,VZERO)
        IF(POI1(I).EQ.1.)THEN
           XPLOT(I1,1)=PHI3
           YPLOT(I1,1)=VIT1(I)-VIC1+VZERO
           I1=I1+1
        ELSE 
           XPLOT(I2,2)=PHI3
           YPLOT(I2,2)=VIT1(I)-VIC1+VZERO
           I2=I2+1
        ENDIF
      END DO
C Cercles pleins 
      NCHAR(1)='913'
      NPTS(1)=I1-1
C Croix 
      NCHAR(2)='513'
      NPTS(2)=I2-1

      write(6,*)'npts(1)=',npts(1),'npts(2)',npts(2)

C On charge la courbe en K=3 
      PI=3.1415926535
      DO I=1,1000
        PHI3=FLOAT(I)/1000.
        CALL VITESSE_CALCULEE_BS4(0.D0,PHI3,VIC1,VV,VV3,
     1                            EE,EE3,OMEGA,OMEGA3,0.,K3,VZERO)
        XPLOT(I,3)=PHI3
        YPLOT(I,3)=VIC1
      END DO
C
      NCHAR(3)='L0'
      NPTS(3)=1000
      KK=3
C For papers: '&LANDSCAPE' (no rotation)
C For work  : '&landscape' (rotation)
      WRITE(PLOTDEV,103)OBJECT(1:INDEX(OBJECT,' ')-1),'_long.ps'
103   FORMAT('&landscape/',A,A)
      WRITE(CHAR1,11) 'Phase',CHAR(0)
      WRITE(CHAR2,11) 'Radial velocity (km/s)',CHAR(0)
      WRITE(CHAR3,11) OBJECT,CHAR(0)
11    FORMAT(A,A1)
      IN_FILE=' '
      IN_COMMENTS=' '
      CALL NEWPLOT(XPLOT,YPLOT,NPTS,DIM,KK,CHAR1,CHAR2,
     1  CHAR3,NCHAR,PLOTDEV,IN_FILE,IN_COMMENTS)

      RETURN
      END
C*********************************************************************
C Pour visualiser la courbe du 3eme corps
C Mode avec cadre plus grand (cycle -0.1 1.1) 
C*********************************************************************
      SUBROUTINE JLP_PLOT_BS4_3(TOBS,VIT1,POI1,
     1                          T,T3,PP,PP3,NN,EE,EE3,OMEGA,OMEGA3,
     1                          K1,K3,VZERO,OBJECT)
      IMPLICIT NONE
      INTEGER DIM,KCUR
      PARAMETER (DIM=1024,KCUR=4)
      REAL TOBS(120),VIT1(120),POI1(120),PHI_MINI
      DOUBLE PRECISION T,T3,PP,PP3,PHI,PHI3
      REAL PI,EE,EE3,OMEGA,OMEGA3,K1,K3,VZERO,VV,VV3,VIC1
      INTEGER ITER,KK,NPTS(KCUR),I,NN,I1,I2
      REAL XPLOT(DIM,KCUR),YPLOT(DIM,KCUR)
      CHARACTER CHAR1*30,CHAR2*30,CHAR3*40,OBJECT*20
      CHARACTER PLOTDEV*40,IN_FILE*40,IN_COMMENTS*80
      CHARACTER*4 NCHAR(KCUR)

C JLP2006/ Modification suggested by Imbert/Griffin,
C -0.1,1.1
      PHI_MINI=0.1

C On charge les mesures en K=1 et K=2
      I1 = 1
      I2 = 1
      DO ITER=1,3
      DO I=1,NN
        PHI=DMOD((TOBS(I)-T)/PP,1.D0)
        IF(PHI.LT.0)PHI=1.+PHI
        PHI3=DMOD((TOBS(I)-T3)/PP3,1.D0)
        IF(PHI3.LT.0)PHI3=1.+PHI3
C On calcule la vitesse radiale du centre de gravite
C Pour cela on appelle VITESSE_CALCULEE_BS4 avec K3=0.
        CALL VITESSE_CALCULEE_BS4(PHI,0.D0,VIC1,VV,VV3,
     1                            EE,EE3,OMEGA,OMEGA3,K1,0.,VZERO)
       IF(ITER.EQ.1)THEN
        IF(POI1(I).EQ.1.)THEN
           XPLOT(I1,1)=PHI3
           YPLOT(I1,1)=VIT1(I)-VIC1+VZERO
           I1=I1+1
        ELSE 
           XPLOT(I2,2)=PHI3
           YPLOT(I2,2)=VIT1(I)-VIC1+VZERO
           I2=I2+1
        ENDIF
C JLP2006/ Modification suggested by Imbert/Griffin,
C to display the peaks of highly eccentric orbits in one piece:
       ELSEIF(ITER.EQ.2)THEN
         IF(PHI3 .LT. PHI_MINI) THEN
          PHI3 = PHI3 + 1.
          IF(POI1(I).EQ.1.)THEN
              XPLOT(I1,1)=PHI3
              YPLOT(I1,1)=VIT1(I)-VIC1+VZERO
              I1=I1+1
            ELSE
              XPLOT(I2,2)=PHI3
              YPLOT(I2,2)=VIT1(I)-VIC1+VZERO
              I2=I2+1
            ENDIF
          ENDIF
C JLP2006/ Modification suggested by Griffin,
C to display the peaks of highly eccentric orbits in one piece:
        ELSEIF(ITER.EQ.3)THEN
          IF(PHI3 .GT. (1.-PHI_MINI)) THEN
            PHI3 = PHI3 - 1.
            IF(POI1(I).EQ.1)THEN
              XPLOT(I1,1)=PHI3
              YPLOT(I1,1)=VIT1(I)-VIC1+VZERO
              I1=I1+1
            ELSE
              XPLOT(I2,2)=PHI3
              YPLOT(I2,2)=VIT1(I)-VIC1+VZERO
              I2=I2+1
            ENDIF
          ENDIF
       ENDIF
      END DO
      END DO
C Cercles pleins 
      NCHAR(1)='913'
      NPTS(1)=I1-1
C Croix 
      NCHAR(2)='513'
      NPTS(2)=I2-1

      write(6,*)'npts(1)=',npts(1),'npts(2)',npts(2)

C On charge la courbe en K=3 
      PI=3.1415926535
      DO I=1,1000
C Phase PHI3 entre -PHI_MINI et 1+PHI_MINI:
        PHI3 = -PHI_MINI + (1+2.*PHI_MINI)*REAL(I)/1000.
        CALL VITESSE_CALCULEE_BS4(0.D0,PHI3,VIC1,VV,VV3,
     1                            EE,EE3,OMEGA,OMEGA3,0.,K3,VZERO)
        XPLOT(I,3)=PHI3
        YPLOT(I,3)=VIC1
      END DO
C
      NCHAR(3)='L0'
      NPTS(3)=1000
      KK=3
C For papers: '&LANDSCAPE' (no rotation)
C For work  : '&landscape' (rotation)
      WRITE(PLOTDEV,103)OBJECT(1:INDEX(OBJECT,' ')-1),'_long.ps'
103   FORMAT('&landscape/',A,A)
      WRITE(CHAR1,11) 'Phase',CHAR(0)
      WRITE(CHAR2,11) 'Radial velocity (km/s)',CHAR(0)
      WRITE(CHAR3,11) OBJECT,CHAR(0)
11    FORMAT(A,A1)
      IN_FILE=' '
      IN_COMMENTS=' '
      CALL NEWPLOT(XPLOT,YPLOT,NPTS,DIM,KK,CHAR1,CHAR2,
     1  CHAR3,NCHAR,PLOTDEV,IN_FILE,IN_COMMENTS)

      RETURN
      END
C*************************************************************
C Resolution de l'equation de Kepler
C Calcul des vitesses radiales
C
C Equation de Kepler: u - e sin u = n(t-t0)
C u est l'anomalie excentrique (relative au centre de l'ellipse)
C v est l'anomalie vraie (relative au foyer de l'ellipse)
C n mouvement moyen: n=2*PI/PP
C
C
C INPUT:
C PHI: phase entre 0 et 1:
C OMEGA: angle (GN,GP): centre sur le foyer G de l'ellipse
C        entre le noeud ascendant N et le periastre P,
C        mesure dans le plan de l'orbite
C EE: ellipticite
C K1: parametre K1 (demi-amplitude en VR)
C VZERO: vitesse radiale du centre de gravite
C
C OUTPUT:
C VV: anomalie vraie
C VIC1: vitesse radiale de la primaire
C*************************************************************
      SUBROUTINE VITESSE_CALCULEE_BS4(PHI,PHI3,VIC1,VV,VV3,
     1                                EE,EE3,OMEGA,OMEGA3,K1,K3,VZERO)
      IMPLICIT NONE
      DOUBLE PRECISION PHI,PHI3
      REAL M,PI,EMME,ANOM,ANOM_OLD,VV,VV3
      REAL EE,EE3,OMEGA,OMEGA3,K1,K3,VZERO,VIC1
      INTEGER I
      
      PI=3.1415926535
      M=2.*PI*PHI
      EMME=AMOD(M,2.*PI)
C M ou EMME: anomalie moyenne 
C ANOM_OLD : anomalie excentrique a l'iteration n-1 
C ANOM : anomalie excentrique a l'iteration n 
      ANOM_OLD=EMME
      DO I=1,1000
        ANOM=EMME+EE*SIN(ANOM_OLD)
        IF(ABS(ANOM-ANOM_OLD).LT.1.E-6)GOTO 12
        ANOM_OLD=ANOM
      ENDDO
C Passage de anomalie excentrique ANOM a l'anomalie vraie VV
12    VV=2.*ATAN(TAN(ANOM/2.)*SQRT((1+EE)/(1-EE)))

C Same for VV3:
      M=2.*PI*PHI3
      EMME=AMOD(M,2.*PI)
      ANOM_OLD=EMME
      DO I=1,1000
        ANOM=EMME+EE3*SIN(ANOM_OLD)
        IF(ABS(ANOM-ANOM_OLD).LT.1.E-6)GOTO 14
        ANOM_OLD=ANOM
      ENDDO
14    VV3=2.*ATAN(TAN(ANOM/2.)*SQRT((1+EE3)/(1-EE3)))

C Calcul de la vitesse radiale (formule fondamentale, p81):
      VIC1=VZERO+K1*(EE*COS(OMEGA)+COS(VV+OMEGA))
     1     +K3*(EE3*COS(OMEGA3)+COS(VV3+OMEGA3))

      RETURN
      END
C*************************************************************
C Calcul de l'ecart-type des residus pour une composante
C (en fait: ecart quadradtique moyen = root mean square)
C et non pas un ecart-type qui tient compte du nombre de
C degres de liberte.
C*************************************************************
      SUBROUTINE RMS_RESIDUS(RES,POI,NN,SIGMA,SUMSQ)
      INTEGER I
      REAL SUM,SUMSQ
      REAL RES(120),POI(120)
      SUMSQ=0.
      SUM=0.
      DO I=1,NN
       SUMSQ=SUMSQ+POI(I)*RES(I)*RES(I)
       SUM=SUM+POI(I)
      END DO
      SIGMA=SQRT(SUMSQ/SUM)
      RETURN
      END
C************************************************************************
C Inversion de la matrice A de rang R
C Same as for BS2 or BS1 but with 12*12 elements 
C Methode du pivot (non optimisee)
C
C INPUT:
C A(NP,NP): matrice a inverser
C
C OUTPUT:
C D(NP,NP): matrice inverse
C************************************************************************
      SUBROUTINE MAT_INV12(A,D,NP)
      IMPLICIT NONE
      INTEGER NP,NPAR
      PARAMETER(NPAR=12)
      REAL A(NPAR,NPAR),C(NPAR,NPAR),D(NPAR,NPAR),G,H
      INTEGER I1,J1,I,J,K,KK,L,L1

C D: Initialise D par la matrice identite
      DO 10 I=1,NP
        DO 10 J=1,NP
          D(I,J)=0
          D(I,I)=1.
   10 CONTINUE
      DO 116 K=1,NP

C Sauve matrice A dans C (utile pour permutation)
      DO 11 I1=1,NP
        DO 11 J1=1,NP
         C(I1,J1)=A(I1,J1)
   11 CONTINUE

C Cas ou le pivot est null:
C JLP: d'apres moi, c'est faux, car on ne peut pas faire de
C permutation dans A, sans le faire aussi dans B
      IF(A(K,K).EQ.0)THEN
        write(6,*)' Pivot null pour K=',K
        write(6,*)' (JLP: j''ai des doutes sur la suite...)'
        KK=K+1
        DO 12 L=KK,NP
          IF(A(L,K).NE.0)THEN
            L1=L
            GO TO 13
          ENDIF
   12   CONTINUE
   13   DO 14 J1=1,NP
C Permutation L1,K des pivots dans A et C:
          C(K,J1)=A(L1,J1)
          C(L1,J1)=A(K,J1)
          A(K,J1)=C(K,J1)
          A(L1,J1)=C(L1,J1)
C Permutation L1,K dans la matrice inverse:
          G=D(L1,J1)
          H=D(K,J1)
          D(K,J1)=G
          D(L1,J1)=H
   14   CONTINUE
      ENDIF
      DO 61 I=1,NP
      IF(A(I,K).NE.0.AND.I.NE.K)THEN
        DO 16 J=1,NP
          A(I,J)=A(I,J)/C(I,K)-A(K,J)/C(K,K)
          D(I,J)=D(I,J)/C(I,K)-D(K,J)/C(K,K)
   16   CONTINUE
      ENDIF
   61 CONTINUE
  116 CONTINUE

C  D=Matrice inverse:
      DO 17 I=1,NP
       DO 17 J=1,NP
        D(I,J)=D(I,J)/A(I,I)
   17 CONTINUE
      RETURN
      END
C*****************************************************************
C Latex format:
C*****************************************************************
      SUBROUTINE SB4_TO_LATEX(PP,PP3,T,T3,OMEGA,OMEGA3,EE,EE3,
     1   K1,K3,VZERO,A1SINI,A3SINI,
     1   FDEM1,FDEM3,SIGMA1,
     1   ERR,ERR_A1SINI,ERR_A3SINI,
     1   ERR_FDEM1,ERR_FDEM3,TOBS,VIT1,
     1   RES1,POI1,NN,
     1   OUT_LATEX,OBJECT)
      IMPLICIT NONE
      DOUBLE PRECISION PP,PP3,T,T3,PHI,PHI3
      INTEGER NN
      REAL TOBS(NN),POI1(NN),VIT1(NN)
      REAL RES1(NN),ERR(11)
      REAL OMEGA,EE,K1,VZERO,A1SINI,FDEM1
      REAL SIGMA1
      REAL ERR_A1SINI,ERR_FDEM1
      REAL OMEGA3,EE3,K3,A3SINI,FDEM3,ERR_FDEM3,ERR_A3SINI
      CHARACTER OUT_LATEX*40,OBJECT*20
      INTEGER I,LU_OUT
      LU_OUT=8
      OPEN(LU_OUT,FILE=OUT_LATEX,STATUS='UNKNOWN')
      WRITE(LU_OUT,3000)OBJECT(1:INDEX(OBJECT,' ')-1)
3000  FORMAT('\\documentclass{article}',/,
     1 '\\usepackage{graphicx}',/,
     1 '\\voffset=-1cm',/,
     1 '\\hoffset=-4cm',/,
     1 '\\textwidth=15.6cm',/,
     1 '\\textheight=26cm',/,
     1 '\\newcommand{\\nodata}{\\ldots}',/,/,
     1 '\\begin{document}',/,/,
     1 '\\centerline{\\large \\bf ',A,'}',/,/,
     1 '\\vskip 1cm',/,
     1 '\\tabcolsep=1mm',/,
     1 '\\begin{tabular}{lccccccccc}',/,
     1 '\\hline',/,
     1 'Name & $P$ & $T_0$ (JD)& $\\omega$ & $e$ &',
     1 ' $K_1$ & $V_0$ ',
     1 '& $a_1 \\sin i$ ',/,
     1 '& $f(m_1)$ & $\\sigma_{(O-C)}$ \\\\',/,
     1 '& days & 2400000+ & deg. & & km/s & km/s ',
     1 '& Gm & M$_\\odot$ & km/s \\\\',/,'\\hline')
      WRITE(LU_OUT,3001)OBJECT(1:INDEX(OBJECT,' ')-1),
     1                  PP,T,OMEGA*180/3.14159,EE,K1,VZERO,
     1                  A1SINI/1.e6,FDEM1,SIGMA1
3001  FORMAT(A,' & ',F12.5,' & ',F12.2,' & ',F5.1,' & ',F8.3,
C K1:
     1       ' & ',F8.2,
C V0
     1       ' & ',F8.2,
C A1 SINI:
     1       ' & ',F8.2,
C F(M1):
     1       ' & ',F8.3,
C SIGMA1:
     1       ' & ',F5.2,' \\\\')
      WRITE(LU_OUT,3002)ERR(1),ERR(2),ERR(3)*180/3.14159,ERR(4),
     1                  ERR(5),ERR(6),ERR_A1SINI/1.e6,ERR_FDEM1
3002  FORMAT(' & $\\pm',F10.5,'$ & $\\pm',F8.2,'$ & $\\pm',
     1       F5.1,'$ & $\\pm',F8.3,'$ & $\\pm',F8.2,
     1       '$ & $\\pm',F8.2,'$ & $\\pm',F8.2,
     1       '$ & $\\pm',F8.3,'$ & \\\\')
      WRITE(LU_OUT,3006)PP3,T3,OMEGA3*180/3.14159,EE3,K3,
     1                  A3SINI/1.e6,FDEM3
3006  FORMAT('outer orbit & ',F12.5,' & ',F12.2,' & ',F5.1,' & ',F8.3,
C K1:
     1       ' & ',F8.2,
C V0
     1       ' & \\nodata',
C A1 SINI:
     1       ' & ',F8.2,
C F(M1):
     1       ' & ',F8.3,
C SIGMA1:
     1       ' & \\nodata \\\\')
C 1=PP 2=T 3=OMEGA 4=EE 5=K1 6=VZERO 7=PP3 8=T3 9=OMEGA3 10=EE3 11=K3
      WRITE(LU_OUT,3007)ERR(7),ERR(8),ERR(9)*180/3.14159,ERR(10),
     1                  ERR(11),ERR_A3SINI/1.e6,ERR_FDEM3
3007  FORMAT(' & $\\pm',F10.5,'$ & $\\pm',F8.2,'$ & $\\pm',
     1       F5.1,'$ & $\\pm',F8.3,'$ & $\\pm',F8.2,
     1       '$ & \\nodata & $\\pm',F8.2,
     1       '$ & $\\pm',F8.3,'$ & \\\\')
      WRITE(LU_OUT,3003)
3003  FORMAT('\\hline',/,'\\end{tabular}',/,'\\vskip 1cm')
C*****************************************
C Figures: 
C*****************************************
      WRITE(LU_OUT,3004)OBJECT(1:INDEX(OBJECT,' ')-1),
     1 OBJECT(1:INDEX(OBJECT,' ')-1)
3004  FORMAT('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%',/,
     1 '\\begin{figure}[h]',/,
     1 '\\centering\\includegraphics*[width=12cm]{',A,'_short.ps}',/,
     1 '\\caption{',A,
     1 ': radial velocity curve of the short period orbit.}',/,
     1 '\\end{figure}',/)

      WRITE(LU_OUT,30041)OBJECT(1:INDEX(OBJECT,' ')-1),
     1 OBJECT(1:INDEX(OBJECT,' ')-1)
30041  FORMAT('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%',/,
     1 '\\begin{figure}[h]',/,
     1 '\\centering\\includegraphics*[width=12cm]{',A,'_long.ps}',/,
     1 '\\caption{',A,
     1 ': radial velocity curve of the long period orbit.}',/,
     1 '\\end{figure}',/)

      WRITE(LU_OUT,3005)OBJECT(1:INDEX(OBJECT,' ')-1),
     1 OBJECT(1:INDEX(OBJECT,' ')-1)
3005  FORMAT('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%',/,
     1 '\\begin{figure}[h]',/,
     1 '\\centering\\includegraphics*[width=12cm]{',A,'_resi.ps}',/,
     1 '\\caption{',A,': radial velocity residuals.',/,
     1 'Filled circles: primary component, open circles: secondary.}',/,
     1 '\\end{figure}',/)

C*****************************************
C O-C
C*****************************************
      WRITE(LU_OUT,3009) OBJECT(1:INDEX(OBJECT,' ')-1)
3009  FORMAT('\\twocolumn',/,/,
     1 '\\footnotesize',/,
     1 '\\begin{tabular}{rrrrr}',/,
     1 '\\hline',/,
     1 '\\multicolumn{5}{c}{',A,'}\\\\',/,
     1 '\\hline',/,
     1 'Date (JD) & Phase 1\\hfil & Phase 3\\hfil & $RV_O$ ',
     1 '& $RV_{(O\\!-\\!C)}$ \\\\',/,
     1 ' 2400000+ & & & km/s & km/s \\\\',/,
     1 '\\hline')
      DO I=1,NN
        IF(POI1(I).NE.0.)THEN
C Cycle count instead of phase (cf. Griffin)
          PHI=(TOBS(I)-T)/PP
          PHI3=(TOBS(I)-T3)/PP3
          WRITE(LU_OUT,3020)TOBS(I),PHI,PHI3,VIT1(I),RES1(I)
3020      FORMAT(F12.3,' & ',F12.2,' & ',F12.2,' & ',
     1           F8.1,' & ',F8.1,' \\\\')
        ENDIF
      ENDDO
      WRITE(LU_OUT,3011)
3011  FORMAT('\\end{tabular}',/)

C      WRITE(LU_OUT,3012)OBJECT(1:INDEX(OBJECT,' ')-1),
C     1 OBJECT(1:INDEX(OBJECT,' ')-1)
3012  FORMAT('\\bigskip',/,
     1 '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%',/,
     1 '\\begin{figure}[h]',/,
     1 '\\centering\\includegraphics*[width=10cm]{',A,'_per.ps}',/,
     1 '\\caption{',A,': periodogram.}',/,
     1 '\\end{figure}',/)

      WRITE(LU_OUT,3013)
3013  FORMAT(/,'\\end{document}')
      CLOSE(LU_OUT)
      RETURN
      END
