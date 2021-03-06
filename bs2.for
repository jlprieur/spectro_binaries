C************************************************************************
C     AMELIORATION DES ELEMENTS ORBITAUX                         *******
C     POUR LES BINAIRES A DEUX SPECTRES                          * BS2 *
C                                                                *******
C     CONSIGNES POUR UTILISATION DU PROGRAMME : NE PAS OUBLIER DE METTRE
C     UN POIDS POUR CHAQUE COMPOSANTE MESUREE. L'ABSENCE DE POIDS SIGNI-
C     FIE QUE LA COMPOSANTE N'A PAS ETE MESUREE. X=1 : PERIODE FIXEE,
C     Y=1 : EXCENTRICITE FIXEE, Y=2 : ORBITE CIRCULAIRE.
C Historique:
C  Cree par Robert Nadal en juillet 1978          
C  Revise par Jean-Louis Prieur en octobre 2000 (sortie graphique)
C  Revu en profondeur en mai 2002
C
C JLP
C Version 18/10/2005 
C************************************************************************
      PROGRAM BS2
      IMPLICIT NONE
      INTEGER IDIM
      PARAMETER (IDIM=120)
      REAL TOBS(IDIM),VIT1(IDIM),VIT2(IDIM),POI1(IDIM),POI2(IDIM)
      REAL RES1(IDIM),RES2(IDIM)
      REAL AA(7,200),B(7),DELT(7),ERR(7),EQM(7),CORR(7)
      REAL A(7,7),D(7,7)
      DOUBLE PRECISION PHI,PP,T,TZERO
      INTEGER NP,X1,Y1,F,DELTA(7),LU_IN,LU_OUT,IT,NN,ITER,NITER
      REAL OMEGA,EE,PI,VV,VIC1,VIC2,OME
      REAL KA,K1,K2,SIGMA1,SIGMA2,SIGMA3,ENNE,NU,Q,S
      REAL FDEM1,FDEM2,A1SINI,A2SINI,SOM,SOMB,SOMD,SOME
      REAL DIST,VZERO,ERR_A1SINI,ERR_A2SINI
      REAL ERR_FDEM1,ERR_FDEM2,POI_0,SOM_POI
      INTEGER I,J,K,L,J1,J2,JJ,N2,IOPT
      CHARACTER OBJECT*20,OUT_RESULTS*60,OUT_LATEX*60
      LOGICAL PUBLI
      PI=3.1415926535

      WRITE(6,*) ' Programme bs2: version du 18/10/2005'
      WRITE(6,*)' Titre (sans blancs, car utilise comme prefixe):'
      READ(5,10) OBJECT
10    FORMAT(A)

C Ouverture du fichier en entree nomme VR.DAT;
      LU_IN=7
      OPEN(LU_IN,FILE='VR_BS2.DAT',STATUS='OLD')
C*******************************************************************************
C Format of VR.DAT
C Periode (jours)
C Passage au periastre ou noeud ascendant si circulaire
C angle omega (noeuds/absides ou 0 si circulaire)
C eccentricite
C k1 (semi amplitude de la courbe de vitesses radiales de la primaire en km/s)
C k2 (semi amplitude de la courbe de vitesses radiales de la secondaire en km/s)
C V0 vitesse du centre de gravite du systeme en km/s
C Parametres optionnels: X1=1 si periode fixeee Y1=2 si orbite circulaire
C 14.2081   48682.860 0.        0.        42.5      49.1      2.8      1  2
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
88    READ(LU_IN,1000)PP,TZERO,OMEGA,EE,K1,K2,VZERO,X1,Y1
1000  FORMAT(F10.7,F10.3,2F10.4,2F10.2,F6.2,2X,2I1)
      WRITE(6,1001)PP,TZERO,OMEGA*180./PI,EE,K1,K2,VZERO
1001  FORMAT(' P=',F12.7,' T0=',F10.3,' Omega=',F9.3,'(deg) E=',F9.3,
     1       ' K1=',F9.2,' K2=',F9.2,'V0=',F9.2)
      T=TZERO
      WRITE(6,*)'Option: X1=',X1,' Y1=',Y1
      T=TZERO
C X1=0: periode libre
C X1=1: periode fixe
C Y1=0: periode libre
C Y1=1: ellipticite fixe
C Y1=2: orbite circulaire
      IF(Y1.EQ.2)THEN
       EE=0.
       OMEGA=0.
      ENDIF


C Nombre de mesures, nombre d'iterations:
881   READ(LU_IN,1100)NN,NITER
1100  FORMAT(I2,6X,I2)
      WRITE(6,*) 'Nombre de mesures:',NN,' Nombre d''iterations',NITER

      READ(LU_IN,1200)(TOBS(I),VIT1(I),POI1(I),VIT2(I),POI2(I),I=1,NN)
1200  FORMAT((F9.3,2(F7.1,F6.2,7X)))
      CLOSE(LU_IN)

C Calcul de la somme des poids:
C (Corrig� en 2005)
      SOM_POI=0.
      DO I=1,NN
        SOM_POI=SOM_POI+POI1(I)+POI2(I)
      END DO

C Nombre de parametres
C 1=PP 2=T 3=OMEGA 4=EE 5=K1 6=K2 7=VZERO
C X1=1: periode fixe
C Y1=1: ellipticite fixe
C Y1=2: orbite circulaire
      NP=7-X1-Y1
      WRITE(LU_OUT,1500)PP,T,OMEGA,EE,K1,K2,VZERO
1500  FORMAT(13X,'ELEMENTS PROVISOIRES',/,16X,'P =',F15.7,/,16X,'T =',
     1F11.3,/,12X,'OMEGA =',F13.5,/,16X,'E =',F13.5,/,15X,'K1 =',F10.2,/,15X,
     2'K2 =',F10.2,/,15X,'V0 =',F10.2,/,/)
C
      DIST=1.
      N2=NN*2
      DO I=1,NP
        DO J=1,N2
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
      NU=0
C Sortie des RES(I) = O-C
C seulement pour les elements provisoires et pour la derniere iteration:
      IF(ITER.EQ.0.OR.ITER.EQ.NITER)THEN
        WRITE(LU_OUT,1600)
1600    FORMAT('OBS.',1X,'DATE (JJ)',2X,'PHASE',1X,'V_OBS1',2X,'V_CALC1',
     1  1X,'(O-C)1',1X,'POIDS1',1X,'V_OBS2',1X,'V_CALC2',1X,'(O-C)2',
     2  1X,'POIDS2',/)
      ENDIF
      ENNE=2.*PI/PP
      DO 53 I=1,NN
      PHI=DMOD((TOBS(I)-T)/PP,1.D0)
      IF(PHI.LT.0)PHI=1.+PHI
C
C Calcul des coefficients des equations normales
C
C Calcul des residus:
      CALL VITESSE_CALCULEE_BS2(PHI,VIC1,VIC2,VV,
     1                          EE,OMEGA,K1,K2,VZERO)
      RES1(I)=VIT1(I)-VIC1
      RES2(I)=VIT2(I)-VIC2
C Sortie des RES(I) = O-C
C seulement pour les elements provisoires et pour la derniere iteration:
      IF(ITER.EQ.0.OR.ITER.EQ.NITER)THEN
C Affichage en verifiant qu'un des poids n'est pas nul:
      IF(POI1(I).EQ.0)THEN
        WRITE(LU_OUT,2000)I,TOBS(I),PHI,VIC1,VIT2(I),VIC2,RES2(I),POI2(I)
2000    FORMAT(I3,2X,F9.3,2X,F4.3,9X,F7.2,14X,3(1X,F7.2),1X,F5.2)
      ELSE
        IF(POI2(I).NE.0)THEN
          WRITE(LU_OUT,1400)I,TOBS(I),PHI,VIT1(I),VIC1,
     1                   RES1(I),POI1(I),VIT2(I),
     1                   VIC2,RES2(I),POI2(I)
1400      FORMAT(I3,2X,F9.3,2X,F4.3,1X,3(F7.2,1X),F5.2,1X,3(F7.2,1X),F5.2)
        ELSE
          WRITE(LU_OUT,2100)I,TOBS(I),PHI,VIT1(I),VIC1,RES1(I),POI1(I),VIC2
2100      FORMAT(I3,2X,F9.3,2X,F4.3,1X,3(F7.2,1X),F5.2,9X,F7.2)
        ENDIF
      ENDIF
      ENDIF
C Aiguillage vers la sortie de la boucle secondaire:
      IF(ITER.EQ.NITER.OR.DIST.LT.1.E-7)GOTO 53

C Calcul des elements de la matrice AA(7,NN)
C Deux cas: V1 et V2, avec K1 et K2, et omega et omega+PI:
      DO 28 K=1,2
C Ligne L de la matrice AA:
C AA contient d'abord les lignes relatives a V1 puis a V2:
        L=I+(K-1)*NN
C Cas K=1:
        IF(K.EQ.1)THEN
          IF(POI1(I).EQ.0)GOTO 28
          NU=NU+1
          KA=K1
          OME=OMEGA
          J1=5
          J2=6
C Cas K=2:
        ELSE
          IF(POI2(I).EQ.0)GOTO 28
          NU=NU+1
          KA=K2
          OME=OMEGA+PI
          J1=6
          J2=5
        ENDIF

C 1=PP 2=T 3=OMEGA 4=EE 5=K1 6=K2 7=VZERO
C Calcul avec parametres precedents:
        S=SIN(VV+OME)
        Q=S*(1+EE*COS(VV))**2/(1-EE**2)**1.5
C Correction de (-2.*PI/PP**2) par rapport a p83 (coeff "D"),
C puisque qu'on utilise la periode PP au lieu de n (moyen mouvement)
C X=1: periode fixe
        IF(X1.EQ.0)THEN
          AA(1,L)=ENNE*(TOBS(I)-T)*KA*Q/PP
        ENDIF
C Coeff "F" p83
        AA(2-X1,L)=KA*ENNE*Q
C Y1=0: ellipticite libre
C Y1=1: ellipticite fixe
C Y1=2: orbite circulaire
        IF(Y1.NE.2)THEN
C Coeff "B" p 83
          AA(3-X1,L)=-KA*(EE*SIN(OME)+S)
          IF(Y1.EQ.0)THEN
C Coeff "C" p 83
            AA(4-X1,L)=KA*(COS(OME)-S*SIN(VV)*(2+EE*COS(VV))/(1-EE**2))
          ENDIF
        ENDIF
C Coeff "A" p 83
        AA(J1-X1-Y1,L)=EE*COS(OME)+COS(VV+OME)
        AA(J2-X1-Y1,L)=0
        AA(NP,L)=1.
   28 CONTINUE

C Calcul du 2eme membre de l'equation normale
      DO J=1,NP
        B(J)=B(J)+POI1(I)*AA(J,I)*RES1(I)+POI2(I)*AA(J,I+NN)*RES2(I)
      END DO
C Calcul du 1er membre de l'equation normale
      DO J=1,NP
        DO K=J,NP
         A(J,K)=A(J,K)+POI1(I)*AA(J,I)*AA(K,I)
     1          +POI2(I)*AA(J,I+NN)*AA(K,I+NN)
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
      CALL MAT_INV7(A,D,NP)

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
        SOMB=SOMB+POI1(I)*RES1(I)**2+POI2(I)*RES2(I)**2
    7 CONTINUE
C BOF JLP2004: normalisation des poids a 1:
      SOMB=SOMB*NU/SOM_POI
C Ajout d'un terme tenant compte des corrections actuelles:
C (Il s'annule lorsqu'on converge, et peut conduire
C a reduire les erreurs s'il est positif,
C alors je le neutralise): 
C (Nota: j'ai verifie en 2005 que cela ne changeait rien)
      SOMD=0
C      DO I=1,NP
C        SOMD=SOMD+B(I)*DELT(I)
C      END DO
C EOF JLP2004
      SOME=SOMB-SOMD
      SOM=SQRT(SOME/(NU-NP))
      DO 80 I=1,NP
        EQM(I)=SOM*SQRT(D(I,I))
   80 CONTINUE

C Nouveaux elements:
      DIST=0
      F=0
      DO I=1,7
        DELTA(I)=0
      END DO
C X1=0 periode libre
C X1=1 periode fixe
      DELTA(1)=X1
C Y1=0 ellipticite libre
C Y1=1 ellipticite fixe
      DELTA(4)=Y1
C Y1=2 orbite circulaire
      IF(Y1.EQ.2)THEN
        DELTA(3)=1
        DELTA(4)=1
      ENDIF

      DO I=1,7
        F=F+DELTA(I)
        J=I-F+DELTA(I)
        ERR(I)=(1-DELTA(I))*EQM(J)
        CORR(I)=(1-DELTA(I))*DELT(J)
        DIST=DIST+ABS(CORR(I))
      END DO

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
      ERR_A1SINI=A1SINI*(ERR(1)/PP + ERR(5)/K1 + ERR(4)*EE/(1.-EE**2))
      ERR_A2SINI=A2SINI*(ERR(1)/PP + ERR(6)/K2 + ERR(4)*EE/(1.-EE**2))
      ERR_FDEM1=FDEM1*(3.*EE*ERR(4)/(1.-EE**2) + ERR(1)/PP
     1         + ERR(5)/K1 + 2.*(ERR(5) + ERR(6))/(K1 + K2))
      ERR_FDEM2=FDEM2*(3.*EE*ERR(4)/(1.-EE**2) + ERR(1)/PP
     1         + ERR(6)/K2 + 2.*(ERR(5) + ERR(6))/(K1 + K2))
      ITER=ITER+1
      WRITE(LU_OUT,2200)ITER
2200  FORMAT(/,5X,I2,'E ITERATION')
      IF(X1.EQ.1)WRITE(LU_OUT,1700)
1700  FORMAT(5X,'PERIODE FIXEE')
      IF(Y1.EQ.1)WRITE(LU_OUT,1800)
1800  FORMAT(5X,'EXCENTRICITE FIXEE')
      IF(Y1.EQ.2)WRITE(LU_OUT,1801)
1801  FORMAT(5X,'ORBITE CIRCULAIRE ')
      WRITE(LU_OUT,1300)CORR(1),PP,ERR(1),CORR(2),T,ERR(2),
     1 CORR(3)*180./PI,OMEGA*180./PI,ERR(3)*180./PI,
     1 CORR(3),OMEGA,ERR(3),
     1 CORR(4),EE,ERR(4),CORR(5),K1,ERR(5),CORR(6),K2,ERR(6),
     1 CORR(7),VZERO,ERR(7),A1SINI,ERR_A1SINI,A2SINI,ERR_A2SINI,
     1 FDEM1,ERR_FDEM1,FDEM2,ERR_FDEM2
1300  FORMAT(5X,'CORRECTIONS',15X,'NOUVEAUX ELEMENTS',
     1 15X,'ERREURS STD.',/,
     1 3X,E9.3,17X,'P =',F15.7,13X,E9.3,/,
     1 3X,E9.3,17X,'T =',F11.3,17X,E9.3,/,
     1 3X,E9.3,13X,'OMEGA =',F13.5,15X,E9.3,4X,'(deg.)',/,
     1 3X,E9.3,13X,'OMEGA =',F13.5,15X,E9.3,4X,'(rad.)',/,
     1 3X,E9.3,17X,'E =',F13.5,15X,E9.3,/,
     1 3X,E9.3,16X,'K1 =',F10.2,18X,E9.3,/,
     1 3X,E9.3,16X,'K2 =',F10.2,18X,E9.3,/,
     1 3X,E9.3,16X,'V0 =',F10.2,18X,E9.3,/,
     1 23X,'A1.SINI =',6X,E10.4,12X,E9.3,/,
     1 23X,'A2.SINI =',6X,E10.4,12X,E9.3,/,
     1 22X,'M1.SIN3I =',6X,E10.4,12X,E9.3,/,
     1 22X,'M2.SIN3I =',6X,E10.4,12X,E9.3,/)
C Erreur grave
      IF(ABS(EE).GE.1.)THEN
          WRITE(LU_OUT,2300)
2300      FORMAT(5X,'ABS(E).GE.1, Poursuite du calcul impossible!')
          CLOSE(LU_OUT)
          STOP
        ENDIF
C Fin de la boucle principale:
      END DO
C Sortie normale:
  81  CONTINUE 
C Visualisation des courbes:
      WRITE(6,*)' Courbe:  0=aucune  1=orbite  2=residus'
      WRITE(6,*)'          3=les deux  4=mode pour publi. ?'
      READ(5,*) IOPT
      PUBLI=(IOPT.EQ.4)
      IF(IOPT.EQ.1)THEN
C Vitesse/Phase
        CALL JLP_PLOT_BS2(TOBS,VIT1,VIT2,POI1,POI2,T,PP,NN,EE,
     1                    OMEGA,K1,K2,VZERO,OBJECT,PUBLI)
      ELSE IF(IOPT.EQ.2)THEN
C Residus:
        CALL JLP_PLOT_RES2(TOBS,RES1,RES2,POI1,POI2,NN,OBJECT,
     1                     PUBLI)
      ELSE IF(IOPT.EQ.3.OR.IOPT.EQ.4)THEN
        CALL JLP_PLOT_BS2(TOBS,VIT1,VIT2,POI1,POI2,T,PP,NN,EE,
     1                    OMEGA,K1,K2,VZERO,OBJECT,PUBLI)
        CALL JLP_PLOT_RES2(TOBS,RES1,RES2,POI1,POI2,NN,OBJECT,
     1                     PUBLI)
      ENDIF

C Calcul de l'ecart-type des residus partiels
      POI_0=1
      DO I=1,NN
        IF(POI1(I).NE.1)POI_0=POI1(I)
      ENDDO
      IF(POI_0.NE.1)THEN
        DO I=1,2
        CALL RMS_RESIDUS_PARTIELS(RES1,POI1,NN,SIGMA1,POI_0)
        CALL RMS_RESIDUS_PARTIELS(RES2,POI2,NN,SIGMA2,POI_0)
        WRITE(LU_OUT,3001)POI_0,SIGMA1,SIGMA2
3001  FORMAT(/,'R�sidus partiels pour poids=',F6.3,/,
     1       'Sigma O-C: Etoile primaire (1)   ',F5.2,' km s$^{-1}$',/,
     1       'Sigma O-C: Etoile secondaire (2) ',F5.2,' km s$^{-1}$')
C Deuxi�me iteration avec poids=1.0
        POI_0=1.0
        ENDDO
      ENDIF
C Calcul de l'ecart-type des residus (non pond�r�s)
      CALL RMS_RESIDUS_NP(RES1,POI1,NN,SIGMA1)
      CALL RMS_RESIDUS_NP(RES2,POI2,NN,SIGMA2)
      CALL RMS_RESIDUS_NP_2(RES1,RES2,POI1,POI2,NN,SIGMA3)
      WRITE(6,3000)SIGMA1,SIGMA2,SIGMA3
      WRITE(LU_OUT,3002)SIGMA1,SIGMA2,SIGMA3
3002  FORMAT(/,'Sigma O-C (non pond�r�): Etoile primaire (1)   ',
     1       F5.2,' km s$^{-1}$',/,
     1         'Sigma O-C (non pond�r�): Etoile secondaire (2) ',
     1       F5.2,' km s$^{-1}$',/,
     1       'Sigma O-C (non pond�r�): global: ',F5.2,' km s$^{-1}$',/)
C Calcul de l'ecart-type des residus (pond�r�s)
      CALL RMS_RESIDUS(RES1,POI1,NN,SIGMA1)
      CALL RMS_RESIDUS(RES2,POI2,NN,SIGMA2)
      CALL RMS_RESIDUS_2(RES1,RES2,POI1,POI2,NN,SIGMA3)
      WRITE(6,3000)SIGMA1,SIGMA2,SIGMA3
      WRITE(LU_OUT,3000)SIGMA1,SIGMA2,SIGMA3
3000  FORMAT(/,'Sigma O-C (pond�r�): Etoile primaire (1)   ',
     1       F5.2,' km s$^{-1}$',/,
     1         'Sigma O-C (pond�r�): Etoile secondaire (2) ',
     1       F5.2,' km s$^{-1}$',/,
     1       'Sigma O-C (pond�r�): global: ',F5.2,' km s$^{-1}$',/)
      CLOSE(LU_OUT)
C Latex format:
      OUT_LATEX=OBJECT(1:INDEX(OBJECT,' ')-1)//'.tex'
      CALL SB2_TO_LATEX(PP,T,OMEGA,EE,K1,K2,VZERO,A1SINI,A2SINI,
     1   FDEM1,FDEM2,SIGMA1,SIGMA2,SIGMA3,ERR,ERR_A1SINI,ERR_A2SINI,
     1   ERR_FDEM1,ERR_FDEM2,TOBS,VIT1,VIT2,RES1,RES2,POI1,POI2,NN,
     1   OUT_LATEX,OBJECT)
C Information sur l'ecran:
      WRITE(6,1300)CORR(1),PP,ERR(1),CORR(2),T,ERR(2),
     1 CORR(3)*180./PI,OMEGA*180./PI,ERR(3)*180./PI,
     1 CORR(3),OMEGA,ERR(3),
     1 CORR(4),EE,ERR(4),CORR(5),K1,ERR(5),CORR(6),K2,ERR(6),
     1 CORR(7),VZERO,ERR(7),A1SINI,ERR_A1SINI,A2SINI,ERR_A2SINI,
     1 FDEM1,ERR_FDEM1,FDEM2,ERR_FDEM2
C Sauvegarde du mouvement du centre de gravite:
      OPEN(LU_OUT,FILE='BS2_GRAV.DAT',STATUS='UNKNOWN')
      DO I=1,NN
        IF(POI1(I)*POI2(I).NE.0.)THEN
          VIC1=(VIT1(I)-VZERO)/K1 + (VIT2(I)-VZERO)/K2
          VIC1=VIC1*(K1+K2)/2.
          WRITE(LU_OUT,1201)TOBS(I),VIC1,1.0
1201      FORMAT(F9.3,F7.1,F6.2)
        ENDIF
      END DO
      CLOSE(LU_OUT)
      WRITE(6,*)' Fichiers de sortie: ',
     1         OUT_RESULTS(1:INDEX(OUT_RESULTS,' ')-1),', ',
     1         OUT_LATEX(1:INDEX(OUT_LATEX,' ')-1),', '
      WRITE(6,*)' mouvement du centre de gravite dans "BS2_GRAV.DAT",'
      IF(IOPT.EQ.3)THEN
        WRITE(6,*)' et graphes dans: ',OBJECT(1:INDEX(OBJECT,' ')-1)
     1            ,'.ps et ',OBJECT(1:INDEX(OBJECT,' ')-1),'_resi.ps'
      ELSE IF(IOPT.EQ.1)THEN
        WRITE(6,*)' et graphe dans: ',OBJECT(1:INDEX(OBJECT,' ')-1),
     1            '.ps'
      ELSE IF(IOPT.EQ.2)THEN
        WRITE(6,*)' et graphe dans: ',OBJECT(1:INDEX(OBJECT,' ')-1),
     1            '_resi.ps'
      ENDIF
      STOP
      END
C*********************************************************************
C Pour visualiser les points et la courbe
C*********************************************************************
      SUBROUTINE JLP_PLOT_BS2(TOBS,VIT1,VIT2,POI1,POI2,
     1      T,PP,NN,EE,OMEGA,K1,K2,VZERO,OBJECT,PUBLI)
      PARAMETER (IDIM=1024,KCUR=4)
      REAL TOBS(120),VIT1(120),VIT2(120),POI1(120),POI2(120)
      DOUBLE PRECISION T,PP,PHI
      REAL PI,EE,OMEGA,K1,K2,VZERO,VV
      INTEGER KK,NPTS(KCUR),I1,I2
      REAL XPLOT(IDIM,KCUR),YPLOT(IDIM,KCUR)
      CHARACTER CHAR1*30,CHAR2*30,CHAR3*40,OBJECT*20
      CHARACTER PLOTDEV*40,IN_FILE*40,IN_COMMENTS*80
      CHARACTER*4 NCHAR(KCUR)
      LOGICAL PUBLI

C On charge les mesures en K=1
      I1 = 1
      I2 = 1
      DO I=1,NN
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
      END DO
C Cercles pleins (primaire) 
      NCHAR(1)='913'
C Cercles vides (secondaire)
      NCHAR(2)='813'
      NPTS(1)=I1-1
      NPTS(2)=I2-1

C On charge la courbe en K=2 
      PI=3.1415926535
      DO I=1,1000
        PHI = REAL(I)/REAL(1000)
        CALL VITESSE_CALCULEE_BS2(PHI,VIC1,VIC2,VV,
     1                            EE,OMEGA,K1,K2,VZERO)
        XPLOT(I,3)=PHI
        YPLOT(I,3)=VIC1
        XPLOT(I,4)=PHI
        YPLOT(I,4)=VIC2
      END DO
      NCHAR(3)='L0'
      NCHAR(4)='L1'
      NPTS(3)=1000
      NPTS(4)=1000
      KK=4
C For papers: '&LANDSCAPE'
C For work (with title): '&landscape'
      IF(PUBLI)THEN
        WRITE(PLOTDEV,102)OBJECT(1:INDEX(OBJECT,' ')-1),'.ps'
102     FORMAT('&LANDSCAPE/',A,A)
        CHAR3=' '
      ELSE
        WRITE(PLOTDEV,103)OBJECT(1:INDEX(OBJECT,' ')-1),'.ps'
103     FORMAT('&landscape/',A,A)
        CHAR3=OBJECT
      ENDIF
      WRITE(6,*) ' Output to: ',PLOTDEV
      CHAR1='Phase'
      CHAR2='Radial velocity (km/s)'
      IN_FILE=' '
      IN_COMMENTS=' '
      CALL NEWPLOT(XPLOT,YPLOT,NPTS,IDIM,KK,CHAR1,CHAR2,
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
C K1, K2: parametres K1 et K2
C VZERO: vitesse radiale du centre de gravite
C
C OUTPUT:
C VV: anomalie vraie
C VIC1, VIC2: vitesses radiales
C*************************************************************
      SUBROUTINE VITESSE_CALCULEE_BS2(PHI,VIC1,VIC2,VV,
     1                                EE,OMEGA,K1,K2,VZERO)
      IMPLICIT NONE
      DOUBLE PRECISION PHI
      REAL M,PI,EMME,ANOM,ANOM_OLD,VV
      REAL EE,OMEGA,K1,K2,VZERO,VIC1,VIC2
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

C Calcul de la vitesse radiale (formule fondamentale, p81):
      VIC1=VZERO+K1*(EE*COS(OMEGA)+COS(VV+OMEGA))
      VIC2=VZERO-K2*(EE*COS(OMEGA)+COS(VV+OMEGA))

      RETURN
      END
C*************************************************************
C Calcul de l'ecart-type des residus (non pond�r�) pour une composante
C*************************************************************
      SUBROUTINE RMS_RESIDUS_NP(RES,POI,NN,SIGMA)
      INTEGER I
      REAL SUM,SUMSQ
      REAL RES(120),POI(120)
      SUMSQ=0.
      SUM=0.
      DO I=1,NN
       IF(POI(I).GT.0)THEN
         SUMSQ=SUMSQ+RES(I)*RES(I)
         SUM=SUM+1.0
       ENDIF
      END DO
      IF(SUM.NE.0.)THEN
         SIGMA=SQRT(SUMSQ/SUM)
      ELSE
         SIGMA=0.
      ENDIF
      RETURN
      END
C*************************************************************
C Calcul de l'ecart-type des residus (pond�r�) pour une composante
C*************************************************************
      SUBROUTINE RMS_RESIDUS(RES,POI,NN,SIGMA)
      INTEGER I
      REAL SUM,SUMSQ
      REAL RES(120),POI(120)
      SUMSQ=0.
      SUM=0.
      DO I=1,NN
       SUMSQ=SUMSQ+POI(I)*RES(I)*RES(I)
       SUM=SUM+POI(I)
      END DO
      IF(SUM.NE.0.)THEN
         SIGMA=SQRT(SUMSQ/SUM)
      ELSE
         SIGMA=0.
      ENDIF
      RETURN
      END
C*************************************************************
C Calcul de l'ecart-type des residus pour une composante
C et un jeu de mesures donn\'e (rep\'er\'e par les poids)
C*************************************************************
      SUBROUTINE RMS_RESIDUS_PARTIELS(RES,POI,NN,SIGMA,POI_0)
      INTEGER I
      REAL SUM,SUMSQ,POI_0
      REAL RES(120),POI(120)
      SUMSQ=0.
      SUM=0.
      DO I=1,NN
       IF(POI_0.EQ.POI(I))THEN
          SUMSQ=SUMSQ+POI(I)*RES(I)*RES(I)
          SUM=SUM+POI(I)
       ENDIF
      END DO
      IF(SUM.NE.0.)THEN
         SIGMA=SQRT(SUMSQ/SUM)
      ELSE
         SIGMA=0.
      ENDIF
      RETURN
      END
C*************************************************************
C Calcul de l'ecart-type global des residus (pond�r�)
C*************************************************************
      SUBROUTINE RMS_RESIDUS_2(RES1,RES2,POI1,POI2,NN,SIGMA)
      INTEGER I
      REAL SUM,SUMSQ
      REAL RES1(120),RES2(120),POI1(120),POI2(120)
      SUMSQ=0.
      SUM=0.
      DO I=1,NN
       SUMSQ=SUMSQ+POI1(I)*RES1(I)*RES1(I)+POI2(I)*RES2(I)*RES2(I)
       SUM=SUM+POI1(I)+POI2(I)
      END DO
      IF(SUM.NE.0.)THEN
         SIGMA=SQRT(SUMSQ/SUM)
      ELSE
         SIGMA=0.
      ENDIF
      RETURN
      END
C*************************************************************
C Calcul de l'ecart-type global des residus (non pond�r�)
C*************************************************************
      SUBROUTINE RMS_RESIDUS_NP_2(RES1,RES2,POI1,POI2,NN,SIGMA)
      INTEGER I
      REAL SUM,SUMSQ
      REAL RES1(120),RES2(120),POI1(120),POI2(120)
      SUMSQ=0.
      SUM=0.
      DO I=1,NN
       IF(POI1(I).GT.0)THEN
         SUMSQ=SUMSQ+RES1(I)*RES1(I)
         SUM=SUM+1.0
       ENDIF
      END DO
      DO I=1,NN
       IF(POI2(I).GT.0)THEN
         SUMSQ=SUMSQ+RES2(I)*RES2(I)
         SUM=SUM+1.0
       ENDIF
      END DO
      IF(SUM.NE.0.)THEN
         SIGMA=SQRT(SUMSQ/SUM)
      ELSE
         SIGMA=0.
      ENDIF
      RETURN
      END
C************************************************************************
C Inversion de la matrice A de rang R
C Same as for BS1 but with 7,7 instead of 6,6
C Methode du pivot (non optimisee)
C
C INPUT:
C A(NP,NP): matrice a inverser
C
C OUTPUT:
C D(NP,NP): matrice inverse
C************************************************************************
      SUBROUTINE MAT_INV7(A,D,NP)
      IMPLICIT NONE
      INTEGER NP
      REAL A(7,7),C(7,7),D(7,7),G,H
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

C*********************************************************************
C Pour visualiser les residus des 2 courbes
C*********************************************************************
      SUBROUTINE JLP_PLOT_RES2(TOBS,RES1,RES2,POI1,POI2,NN,OBJECT,
     1                          PUBLI)
      IMPLICIT NONE
      INTEGER DIM,KCUR
      PARAMETER (DIM=1024,KCUR=2)
      REAL TOBS(120),RES1(120),RES2(120),POI1(120),POI2(120)
      INTEGER KK,NPTS(KCUR),I,NN,I1,I2
      REAL XPLOT(DIM,KCUR),YPLOT(DIM,KCUR)
      CHARACTER CHAR1*30,CHAR2*30,CHAR3*40,OBJECT*20
      CHARACTER PLOTDEV*40,IN_FILE*40,IN_COMMENTS*80
      CHARACTER*4 NCHAR(KCUR)
      LOGICAL PUBLI

C On charge les residus en K=1 et K=2:
      I1 = 1
      I2 = 1
      DO I=1,NN
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
      END DO
C Cercles pleins (primaire)
      NCHAR(1)='910'
C Cercles vides (secondaire)
      NCHAR(2)='810'
      NPTS(1)=I1-1
      NPTS(2)=I2-1

      KK=2
C For papers: '&LANDSCAPE'
C For work (with title): '&landscape'
      IF(PUBLI)THEN
        WRITE(PLOTDEV,102)OBJECT(1:INDEX(OBJECT,' ')-1),'_resi.ps'
102     FORMAT('&LANDSCAPE/',A,A)
        CHAR3=' '
      ELSE
        WRITE(PLOTDEV,103)OBJECT(1:INDEX(OBJECT,' ')-1),'_resi.ps'
103     FORMAT('&landscape/',A,A)
        CHAR3=OBJECT
      ENDIF
      WRITE(6,*) ' Output to: ',PLOTDEV
      CHAR1='Date'
      CHAR2='Residuals (km/s)'
      IN_FILE=' '
      IN_COMMENTS=' '
      CALL NEWPLOT(XPLOT,YPLOT,NPTS,DIM,KK,CHAR1,CHAR2,
     1  CHAR3,NCHAR,PLOTDEV,IN_FILE,IN_COMMENTS)

      RETURN
      END
C*****************************************************************
C Latex format:
C*****************************************************************
      SUBROUTINE SB2_TO_LATEX(PP,T,OMEGA,EE,K1,K2,VZERO,A1SINI,A2SINI,
     1   FDEM1,FDEM2,SIGMA1,SIGMA2,SIGMA3,ERR,ERR_A1SINI,ERR_A2SINI,
     1   ERR_FDEM1,ERR_FDEM2,TOBS,VIT1,VIT2,RES1,RES2,POI1,POI2,NN,
     1   OUT_LATEX,OBJECT)
      IMPLICIT NONE
      DOUBLE PRECISION PP,T,PHI
      INTEGER NN
      REAL TOBS(NN),POI1(NN),POI2(NN),VIT1(NN),VIT2(NN)
      REAL RES1(NN),RES2(NN),ERR(7)
      REAL OMEGA,EE,K1,K2,VZERO,A1SINI,A2SINI,FDEM1,FDEM2
      REAL SIGMA1,SIGMA2,SIGMA3
      REAL ERR_A1SINI,ERR_A2SINI,ERR_FDEM1,ERR_FDEM2
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
     1 '\\begin{tabular}{lcccccccccccccc}',/,
     1 '\\hline',/,
     1 'Name & $P$ & $T_0$ (JD)& $\\omega$ & $e$ &',
     1 ' $K_1$ & $K_2$ & $V_0$ ',
     1 '& $a_1 \\sin i$ & $a_2 \\sin i$',/,
     1 '& $f(m_1)$ & $f(m_2)$ & $\\sigma_{1(O-C)}$ ',/,
     1 '& $\\sigma_{2(O-C)}$ & $\\sigma_{(O-C)}$ \\\\',/,
     1 '& days & 2400000+ & deg. & & km s$^{-1}$ &',
     1 ' km s$^{-1}$ & km s$^{-1}$ ',
     1 '& Gm & Gm & M$_\\odot$ & M$_\\odot$ &',
     1 ' km s$^{-1}$ & km s$^{-1}$ ',
     1 '& km s$^{-1}$ \\\\',/,'\\hline')
      WRITE(LU_OUT,3001)OBJECT(1:INDEX(OBJECT,' ')-1),
     1  PP,T,OMEGA*180/3.14159,EE,K1,K2,VZERO,A1SINI/1.e6,A2SINI/1.e6,
     1  FDEM1,FDEM2,SIGMA1,SIGMA2,SIGMA3
3001  FORMAT(A,' & ',F12.5,' & ',F12.2,' & ',F5.1,' & ',F8.3,
C K1, K2:
     1       ' & ',F8.2,' & ',F8.2, 
C V0
     1       ' & ',F8.2,/,
C A1 SINI, A2 SINI:
     1       ' & ',F8.2,' & ',F8.2,
C F(M1), F(M2):
     1       ' & ',F8.3,' & ',F8.3,
     1       ' & ',F5.2,' & ',F5.2,' & ',F5.2,' \\\\')
      WRITE(LU_OUT,3002)ERR(1),ERR(2),ERR(3)*180/3.14159,ERR(4),
     1  ERR(5),ERR(6),ERR(7),
     1  ERR_A1SINI/1.e6,ERR_A2SINI/1.e6,ERR_FDEM1,ERR_FDEM2
3002  FORMAT(' & $\\pm',F10.5,'$ & $\\pm',F8.2,'$ & $\\pm',
     1       F5.1,'$ & $\\pm',F8.3,'$ & $\\pm',F8.2,
     1       '$ ',/,'& $\\pm',F8.2,'$ & $\\pm',F8.2,'$ & $\\pm',F8.2,
     1       '$ & $\\pm',F8.2,'$ & $\\pm',F8.3,'$ & $\\pm',F8.2,
     1       '$ & & & \\\\')
      WRITE(LU_OUT,3003)
3003  FORMAT('\\end{tabular}',/,'\\vskip 1cm')
C*****************************************
C Figures: 
C*****************************************
      WRITE(LU_OUT,3004)OBJECT(1:INDEX(OBJECT,' ')-1),
     1 OBJECT(1:INDEX(OBJECT,' ')-1)
3004  FORMAT('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%',/,
     1 '\\begin{figure}[h]',/,
     1 '\\centering\\includegraphics*[width=12cm]{',A,'.ps}',/,
     1 '\\caption{',A,': radial velocity curve.',/,
     1 'Filled circles: primary component, open circles: secondary.}',/,
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
     1 '\\begin{tabular}{rrcrcr}',/,
     1 '\\hline',/,
     1 '\\multicolumn{5}{c}{',A,'}\\\\',/,
     1 '\\hline',/,
     1 'Date (JD) & Cycle & $RV_1$ & $(O\\!-\\!C)_1$ ',
     1 '& $RV_2$ & $(O\\!-\\!C)_2$ \\\\',/,
     1 ' 2400000+ & & km s$^{-1}$ & km s$^{-1}$ &',
     1 ' km s$^{-1}$ & km s$^{-1}$\\\\',/,
     1 '\\hline')
      DO I=1,NN
C          PHI=DMOD((TOBS(I)-T)/PP,1.D0)
C          IF(PHI.LT.0)PHI=1.+PHI
C Cycle
          PHI=(TOBS(I)-T)/PP
        IF(POI1(I).NE.0..AND.POI2(I).NE.0.)THEN
          WRITE(LU_OUT,3020)TOBS(I),PHI,VIT1(I),RES1(I),
     1                      VIT2(I),RES2(I)
3020      FORMAT(F12.2,' & ',F8.2,
     1           ' & ',F8.1,' & ',F8.1,' & ',F8.1,' & ',
     1           F8.1,' \\\\')
        ELSE IF(POI1(I).EQ.0.)THEN
          WRITE(LU_OUT,3021)TOBS(I),PHI,VIT2(I),RES2(I)
3021      FORMAT(F12.2,' & ',F8.2,
     1           ' & \\nodata & \\nodata & ',F8.1,' & ',
     1           F8.1,' \\\\')
        ELSE
          WRITE(LU_OUT,3022)TOBS(I),PHI,VIT1(I),RES1(I)
3022      FORMAT(F12.2,' & ',F8.2,
     1           ' & ',F8.1,' & ',F8.1,' & \\nodata & ',
     1           '\\nodata \\\\')
        ENDIF
      ENDDO
      WRITE(LU_OUT,3011)
3011  FORMAT('\\hline',/,'\\end{tabular}',/)

      WRITE(LU_OUT,3012)OBJECT(1:INDEX(OBJECT,' ')-1),
     1 OBJECT(1:INDEX(OBJECT,' ')-1)
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
