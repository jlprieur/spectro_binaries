C************************************************************************
C     AMELIORATION DES ELEMENTS ORBITAUX                         0000000
C     POUR LES BINAIRES A UN SPECTRE                             0 bs1 0
C                                                                0000000
C     Options:P fixe (X=1), e fixe (Y=1), orbite circulaire (Y=2)
C     (X at position 58, Y at position 59)
C     Pour l'option orb. circul., fixer omega=0 et entrer pour T
C     l'epoque de passage au noeud ascendant.
C
C Historique:
C  Cree par Robert Nadal en juillet 1978 
C  Revise par Jean-Louis Prieur en octobre 2000 (sortie graphique)
C  Revu en profondeur en mai 2002
C
C JLP
C Version 30/05/2005
C************************************************************************
      PROGRAM BS1 
      IMPLICIT NONE
      INTEGER IDIM
      PARAMETER (IDIM=120)
      REAL POI(IDIM),TOBS(IDIM),VIT(IDIM),RES(IDIM)
      REAL ERR(6),EQM(6),CORR(6),W1,W2
      INTEGER LU_IN,LU_OUT,LU_OUT_TXT,NP,X1,Y1
      INTEGER NORB,NN,NITER,I,IOPT
      REAL EE,OMEGA,KA,VZERO,PI,SIGMA
      REAL TZERO,ASINI,FDEM,ERR_ASINI,ERR_FDEM
      DOUBLE PRECISION PP,T
      CHARACTER OBJECT*20,OUT_RESULTS*60,OUT_LATEX*60
      CHARACTER BUFFER*80
      LOGICAL PUBLI,INTERACTIVE
      PI=3.1415926535

      WRITE(6,*)' BS1            JLP/Version 30/05/2005'
      WRITE(6,*)' Titre (sans blancs, car utilisé comme préfixe):'
      READ(5,10) OBJECT 
10    FORMAT(A)
      INTERACTIVE=(OBJECT(1:15).NE.'NON-INTERACTIVE')
     1            .AND.(OBJECT(1:15).NE.'non-interactive')

C Ouverture du fichier en entree nomme VR_BS1.DAT;
      LU_IN=7
      OPEN(LU_IN,FILE='VR_BS1.DAT',STATUS='OLD')
C**************************************************************************
C Format of VR.DAT
C Periode (jours)
C Passage au periastre ou noeud ascendant si circulaire
C angle omega (noeuds/absides ou 0 si circulaire)
C eccentricite
C k1 (semi amplitude de la courbe de vitesses radiales de la primaire en km/s)
C V0 vitesse du centre de gravite du systeme en km/s
C Parametres optionnels: X1=1 si periode fixeee Y1=2 si orbite circulaire
C 14.2081   48682.860 0.        0.        42.5      2.8      1     2
C Nombre de mesures de vitesses radiales (dans ce fichier)
C Nombre d'iterations (pour converger)
C 42       5
C Date en jour Julien (-2000000), vitesse radiale de la primaire, poids de la mesure
C vitesse radiale de la secondaire, poids de la mesure
C**************************************************************************
      NORB=1
C Format des elements provisoires de l'orbite: 
 1000 FORMAT(F10.6,F10.3,2F10.5,F10.2,F6.2,1X,2I1)
      READ(LU_IN,1000)PP,TZERO,OMEGA,EE,KA,VZERO,X1,Y1
      WRITE(6,*)' P=',PP,' T0=',TZERO,' Omega=',OMEGA,' E=',EE,' KA=',KA
      WRITE(6,*)' V0=',VZERO
      T=TZERO
      IF(X1.NE.1)X1=0
      IF(Y1.NE.1.AND.Y1.NE.2)Y1=0
      IF(INTERACTIVE) WRITE(6,*)'Option: X1=',X1,' Y1=',Y1
C X1=1: periode fixe
C Y1=1: ellipticite fixe
C Y1=2: orbite circulaire
      IF(Y1.EQ.2)THEN
       EE=0.
       OMEGA=0.
      ENDIF

C Nombre de mesures, nombre d'iterations:
      READ(LU_IN,*)NN,NITER
C1100  FORMAT(I3,6X,I2)
      IF(INTERACTIVE) 
     1     WRITE(6,*) 'Nombre de mesures:',NN,' Nombre d''iterations',NITER
      DO I=1,NN
C      READ(LU_IN,1200)TOBS(I),VIT(I),POI(I)
 1200  FORMAT(F9.3,F7.1,F6.2)
       READ(LU_IN,10)BUFFER
C When the weight is not present, assumes that weight=1
       POI(I)=-1.0
       READ(BUFFER,*)TOBS(I),VIT(I)
       READ(BUFFER,*,ERR=181,END=181)W1,W2,POI(I)
181    IF(POI(I).EQ.-1.0)POI(I)=1.0

      END DO
      CLOSE(LU_IN)

C Creation du fichier en sortie OBJECT.txt 
      LU_OUT_TXT=8
      OUT_RESULTS=OBJECT(1:INDEX(OBJECT,' ')-1)//'.txt'
      OPEN(LU_OUT_TXT,FILE=OUT_RESULTS,STATUS='UNKNOWN')
C Nombre de parametres
C 1=PP 2=T 3=OMEGA 4=EE 5=KA 6=VZERO
C X=1: periode fixe
C Y=1: ellipticite fixe 
C Y=2: orbite circulaire 
      NP=6-X1-Y1
      IF(INTERACTIVE) WRITE(LU_OUT_TXT,1500)PP,T,OMEGA,EE,KA,VZERO
 1500 FORMAT(13X,'ELEMENTS PROVISOIRES',/,16X,'P =',F15.7,/,16X,'T =',
     1F11.3,/,12X,'OMEGA =',F13.5,/,16X,'E =',F13.5,/,16X,'K =',F10.2,/,15X,
     2'VO =',F10.2,/,/)

C Minimisation par moindres carrés:
      CALL BS1_FIT(POI,TOBS,VIT,RES,ERR,EQM,CORR,
     1                   PP,T,TZERO,OMEGA,EE,KA,VZERO,
     1                   ASINI,FDEM,ERR_ASINI,ERR_FDEM,
     1                   LU_OUT_TXT,NP,X1,Y1,NN,NITER,INTERACTIVE)

C Calcul de l'ecart-type des residus
      CALL RMS_RESIDUS(RES,POI,NN,SIGMA)
      WRITE(6,3000)SIGMA
      IF(INTERACTIVE)THEN
         WRITE(LU_OUT_TXT,3000)SIGMA
      ELSE
         WRITE(LU_OUT_TXT,*)'SIGMA=',SIGMA
      ENDIF
3000  FORMAT('Ecart-type moyen des residus (Sigma O-C): ',
     1        F5.2,' km.s$^{-1}$',/)
      CLOSE(LU_OUT_TXT)

C********************************************************************
C BOF INTERACTIVE mode
      IF(INTERACTIVE)THEN
C Visualisation des courbes:
      WRITE(6,*) ' Courbe: 0=aucune  1=orbite  2=residus '
      WRITE(6,*) '         3=les deux  4=mode pour publi. ?'
      READ(5,*) IOPT
      PUBLI=(IOPT.EQ.4)
      IF(IOPT.EQ.1)THEN
C Vitesse/Phase
        CALL JLP_PLOT_BS1(TOBS,VIT,T,PP,NN,EE,OMEGA,KA,VZERO,OBJECT,
     1                    PUBLI) 
      ELSE IF(IOPT.EQ.2)THEN
C Residus:
        CALL JLP_PLOT_RES1(TOBS,RES,POI,NN,OBJECT,PUBLI)
      ELSE IF(IOPT.EQ.3.OR.IOPT.EQ.4)THEN
        CALL JLP_PLOT_BS1(TOBS,VIT,T,PP,NN,EE,OMEGA,KA,VZERO,OBJECT,
     1                    PUBLI) 
        CALL JLP_PLOT_RES1(TOBS,RES,POI,NN,OBJECT,PUBLI)
      ENDIF

C Latex format:
      OUT_LATEX=OBJECT(1:INDEX(OBJECT,' ')-1)//'.tex'
      CALL SB1_TO_LATEX(PP,T,OMEGA,EE,KA,VZERO,ASINI,FDEM,
     1                  SIGMA,ERR,ERR_ASINI,ERR_FDEM,
     1                  TOBS,VIT,RES,POI,NN,OUT_LATEX,OBJECT)
C Information sur l'ecran:
      WRITE(6,1300)CORR(1),PP,ERR(1),CORR(2),T,ERR(2),
     1CORR(3)*180/PI,OMEGA*180/PI,ERR(3)*180/PI,
     1CORR(3),OMEGA,ERR(3),
     1CORR(4),EE,ERR(4),CORR(5),KA,ERR(5),CORR(6),VZERO,ERR(6),
     1ASINI,ERR_ASINI,FDEM,ERR_FDEM
 1300 FORMAT(1X,'CORRECTIONS',15X,'NOUVEAUX ELEMENTS',
     115X,'ERREURS STD.',/,
     13X,E9.3,17X,'P =',F15.7,13X,E9.3,/,
     13X,E9.3,17X,'T =',F11.3,17X,E9.3,/,
     13X,E9.3,13X,'OMEGA =',F13.5,15X,E9.3,4X,'(deg.)',/,
     13X,E9.3,13X,'OMEGA =',F13.5,15X,E9.3,4X,'(rad.)',/,
     13X,E9.3,17X,'E =',F13.5,15X,E9.3,/,
     13X,E9.3,17X,'K =',F10.2,18X,E9.3,/,
     13X,E9.3,16X,'V0 =',F10.2,18X,E9.3,/,
     124X,'A.SINI =',6X,E10.4,12X,E9.3,/,
     126X,'F(M) =',6X,E10.4,12X,E9.3,/,/)

C Sauvegarde des residus:
      LU_OUT=8
      OPEN(LU_OUT,FILE='RESID.DAT',STATUS='UNKNOWN')
      WRITE(LU_OUT,*) NN
      WRITE(LU_OUT,1200)(TOBS(I),RES(I),POI(I),I=1,NN)
      CLOSE(LU_OUT)

      WRITE(6,*)' Fichiers de sortie: ',
     1         OUT_RESULTS(1:INDEX(OUT_RESULTS,' ')-1),', ',
     1         OUT_LATEX(1:INDEX(OUT_LATEX,' ')-1),', RESID.DAT,'
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
C EOF INTERACTIVE mode
      ENDIF
C********************************************************************

      STOP
      END
C*********************************************************************
C Pour visualiser les points et la courbe
C*********************************************************************
      SUBROUTINE JLP_PLOT_BS1(TOBS,VIT,T,PP,NN,EE,OMEGA,KA,VZERO,
     1                        OBJECT,PUBLI) 
      IMPLICIT NONE
      INTEGER IDIM,KCUR
      PARAMETER (IDIM=1024,KCUR=3)
      REAL TOBS(*),VIT(*)
      DOUBLE PRECISION T,PP,PHI
      REAL PI,EE,OMEGA,KA,VZERO,VIC,VV
      INTEGER I,NN,KK,NPTS(KCUR)
      REAL XPLOT(IDIM,KCUR),YPLOT(IDIM,KCUR)
      CHARACTER CHAR1*30,CHAR2*30,CHAR3*40,OBJECT*20
      CHARACTER PLOTDEV*40,IN_FILE*40,IN_COMMENTS*80
      CHARACTER*4 NCHAR(KCUR)
      LOGICAL PUBLI

      PI=3.1415926535
C On charge les mesures en K=1
      DO I=1,NN
        PHI=DMOD((TOBS(I)-T)/PP,1.D0)
        IF(PHI.LT.0)PHI=1.+PHI
        XPLOT(I,1)=PHI
        YPLOT(I,1)=VIT(I)
      END DO
C Cercles pleins (primaire)
      NCHAR(1)='913'
      NPTS(1)=NN

C On charge la courbe en K=2 
      DO I=1,1000
C Phase PHI entre 0 et 1:
        PHI = REAL(I)/REAL(1000)
        CALL VITESSE_CALCULEE(PHI,VIC,VV,EE,OMEGA,KA,VZERO)
        XPLOT(I,2)=PHI
        YPLOT(I,2)=VIC
      END DO
      NCHAR(2)='L0'
      NPTS(2)=1000
      KK=2
C For papers: '&LANDSCAPE' (no rotation)
C For work  : '&landscape' (rotation)
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
C KA: parametre K1
C VZERO: vitesse radiale du centre de gravite
C
C OUTPUT:
C VV: anomalie vraie
C VIC: vitesse radiale
C*************************************************************
      SUBROUTINE VITESSE_CALCULEE(PHI,VIC,VV,EE,OMEGA,KA,VZERO)
      IMPLICIT NONE
      DOUBLE PRECISION PHI
      REAL PI,M,EMME,ANOM,ANOM_OLD,VV,VIC
      REAL EE,OMEGA,KA,VZERO
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
12     VV=2.*ATAN(TAN(ANOM/2.)*SQRT((1+EE)/(1-EE)))

C Calcul de la vitesse radiale (formule fondamentale, p81):
      VIC=VZERO+KA*(EE*COS(OMEGA)+COS(VV+OMEGA))

      RETURN
      END
C*************************************************************
C Calcul de l'ecart-type des residus
C*************************************************************
      SUBROUTINE RMS_RESIDUS(RES,POI,NN,SIGMA)
      IMPLICIT NONE
      INTEGER I,NN
      REAL SUM,SUMSQ
      REAL RES(*),POI(*),SIGMA
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
C Methode du pivot (non optimisee)
C 
C INPUT:
C A(NP,NP): matrice a inverser
C
C OUTPUT:
C D(NP,NP): matrice inverse
C************************************************************************
      SUBROUTINE MAT_INV(A,D,NP)
      IMPLICIT NONE
      INTEGER NP
      REAL A(6,6),C(6,6),D(6,6),G,H
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
C Pour visualiser les residus 
C*********************************************************************
      SUBROUTINE JLP_PLOT_RES1(TOBS,RES1,POI1,NN,OBJECT,PUBLI)
      IMPLICIT NONE
      INTEGER DIM,KCUR
      PARAMETER (DIM=1024,KCUR=1)
      REAL TOBS(120),RES1(120),POI1(120)
      INTEGER KK,NPTS(KCUR),I,NN,I1
      REAL XPLOT(DIM,KCUR),YPLOT(DIM,KCUR)
      CHARACTER CHAR1*30,CHAR2*30,CHAR3*40,OBJECT*20
      CHARACTER PLOTDEV*40,IN_FILE*40,IN_COMMENTS*80
      CHARACTER*4 NCHAR(KCUR)
      LOGICAL PUBLI

C On charge les residus en K=1:
      I1 = 1
      DO I=1,NN
        IF(POI1(I).NE.0.)THEN
           XPLOT(I1,1)=TOBS(I)
           YPLOT(I1,1)=RES1(I)
           I1=I1+1
        ENDIF
      END DO
C Cercles pleins (primaire)
      NCHAR(1)='913'
      NPTS(1)=I1-1

      KK=1
C For papers: '&LANDSCAPE' (no rotation)
C For work  : '&landscape' (rotation)
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
      SUBROUTINE SB1_TO_LATEX(PP,T,OMEGA,EE,KA,VZERO,ASINI,FDEM,
     1                        SIGMA,ERR,ERR_ASINI,ERR_FDEM,
     1                        TOBS,VIT,RES,POI,NN,OUT_LATEX,OBJECT)
      IMPLICIT NONE
      DOUBLE PRECISION PP,T,PHI
      INTEGER NN
      REAL POI(NN),TOBS(NN),VIT(NN),RES(NN)
      REAL OMEGA,EE,KA,VZERO,ASINI,FDEM,SIGMA,ERR(6)
      REAL ERR_ASINI,ERR_FDEM
      CHARACTER OUT_LATEX*40,OBJECT*20
      INTEGER I,LU_OUT
      LU_OUT=8
      OPEN(LU_OUT,FILE=OUT_LATEX,STATUS='UNKNOWN')
      WRITE(LU_OUT,3000)OBJECT(1:INDEX(OBJECT,' ')-1)
3000  FORMAT('\\documentclass{article}',/,
     1 '\\usepackage{graphicx}',/,
     1 '\\voffset=-1cm',/,
     1 '\\hoffset=-3cm',/,
     1 '\\textwidth=15.6cm',/,
     1 '\\textheight=26cm',/,
     1 '\\newcommand{\\nodata}{\\ldots}',/,/,
     1 '\\begin{document}',/,/,
     1 '\\centerline{\\large \\bf ',A,'}',/,/,
     1 '\\vskip 1cm',/,
     1 '\\tabcolsep=1mm',/,
     1 '\\begin{tabular}{lcccccccccccc}',/,
     1 '\\hline',/,
     1 'Name & $P$ & $T_0$ (JD)& $\\omega$ & $e$ &',
     1 ' $K_1$ & $K_2$ & $V_0$ ',
     1 '& $a_1 \\sin i$ & $a_2 \\sin i$',/,
     1 '& $f(m_1)$ & $f(m_2)$ & $\\sigma_{(O-C)}$ \\\\',/,
     1 '& days & 2400000+ & deg. & & km.s$^{-1}$ ',
     1 '& km.s$^{-1}$ & km.s$^{-1}$',
     1 '& Gm & Gm & M$_\\odot$ & M$_\\odot$ & km.s$^{-1}$ \\\\',
     1 /,'\\hline')
      WRITE(LU_OUT,3001)OBJECT(1:INDEX(OBJECT,' ')-1),
     1                  PP,T,OMEGA*180/3.14159,EE,KA,VZERO,
     1                  ASINI/1.e6,FDEM,SIGMA
3001  FORMAT(A,' & ',F13.6,' & ',F13.3,' & ',F5.1,' & ',F8.3,
C K1:
     1       ' & ',F8.2, 
C K2:
     1       ' & \\nodata ',
C V0
     1       ' & ',F8.2,
C A1 SINI:
     1       ' & ',F8.2,
C A2 SINI:
     1       ' & \\nodata ',
C F(M1):
     1       ' & ',F10.6,
C F(M2):
     1       ' & \\nodata ',
     1       ' & ',F5.2,' \\\\')
      WRITE(LU_OUT,3002)ERR(1),ERR(2),ERR(3)*180/3.14159,ERR(4),
     1                  ERR(5),ERR(6),ERR_ASINI/1.e6,ERR_FDEM
3002  FORMAT(' & $\\pm',F10.6,'$ & $\\pm',F8.3,'$ & $\\pm',
     1       F5.1,'$ & $\\pm',F8.3,'$ & $\\pm',F8.2,
     1       '$ & \\nodata & $\\pm',F8.2,'$ & $\\pm',F8.2,
     1       '$ & \\nodata & $\\pm',F10.6,'$ & \\nodata & \\\\')
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
     1 '\\caption{',A,': radial velocity curve.}',/,
     1 '\\end{figure}',/)

      WRITE(LU_OUT,3005)OBJECT(1:INDEX(OBJECT,' ')-1),
     1 OBJECT(1:INDEX(OBJECT,' ')-1)
3005  FORMAT('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%',/,
     1 '\\begin{figure}[h]',/,
     1 '\\centering\\includegraphics*[width=12cm]{',A,'_resi.ps}',/,
     1 '\\caption{',A,': radial velocity residuals.}',/,
     1 '\\end{figure}',/)

C*****************************************
C O-C
C*****************************************
      WRITE(LU_OUT,3009) OBJECT(1:INDEX(OBJECT,' ')-1)
3009  FORMAT('\\twocolumn',/,/,
     1 '\\footnotesize',/,
     1 '\\begin{tabular}{rrrr}',/,
     1 '\\hline',/,
     1 '\\multicolumn{4}{c}{',A,'}\\\\',/,
     1 '\\hline',/,
     1 'Date (JD) & Cycle\\hfil & $RV$\\hfil & $(O-C)$\\hfil \\\\',/,
     1 ' 2400000+ & & km.s$^{-1}$ & km.s$^{-1}$ \\\\',/,
     1 '\\hline')
      DO I=1,NN
        IF(POI(I).NE.0.)THEN
C Cycle = phase PHI complète, Cf. Griffin):
C          PHI=DMOD((TOBS(I)-T)/PP,1.D0)
C          IF(PHI.LT.0)PHI=1.+PHI
          PHI=(TOBS(I)-T)/PP
          WRITE(LU_OUT,3010)TOBS(I),PHI,VIT(I),RES(I)
3010      FORMAT(F12.2,' & ',F12.2,' & ',F8.1,' & ',F8.1,' \\\\')
        ENDIF
      ENDDO
      WRITE(LU_OUT,3011)
3011  FORMAT('\\hline',/,'\\end{tabular}',/)

      WRITE(LU_OUT,3012)OBJECT(1:INDEX(OBJECT,' ')-1),
     1 OBJECT(1:INDEX(OBJECT,' ')-1)
3012  FORMAT('\\bigskip',/,
     1 '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%',/,
     1 '\\begin{figure}[h]',/,
     1 '\\centering\\includegraphics*[width=12cm]{',A,'_per.ps}',/,
     1 '\\caption{',A,': periodogram.}',/,
     1 '\\end{figure}',/)

      WRITE(LU_OUT,3013)
3013  FORMAT(/,'\\end{document}')
      CLOSE(LU_OUT)
      RETURN
      END
C*************************************************************************
C
C INPUT: 
C
C PP,T,TZERO,OMEGA,EE,KA,VZERO: elements provisoires
C
C OUTPUT:
C
C ASINI,ERR_ASINI,FDEM,ERR_FDEM
C*************************************************************************
      SUBROUTINE BS1_FIT(POI,TOBS,VIT,RES,ERR,EQM,CORR,
     1                   PP,T,TZERO,OMEGA,EE,KA,VZERO,
     1                   ASINI,FDEM,ERR_ASINI,ERR_FDEM,
     1                   LU_OUT_TXT,NP,X1,Y1,NN,NITER,INTERACTIVE)
      IMPLICIT NONE 
      INTEGER IDIM
      PARAMETER (IDIM=120)
      REAL POI(IDIM),TOBS(IDIM),VIT(IDIM),RES(IDIM)
      REAL AA(6,IDIM),B(6),DELT(6),ERR(6),EQM(6),CORR(6)
      REAL A(6,6),D(6,6),SOM_POI
      INTEGER LU_OUT_TXT,NP,X1,Y1,F,DELTA(6)
      INTEGER NN,NITER,I,J,K,ITER
      REAL EE,OMEGA,KA,VZERO,PI
      REAL TZERO,SOME,ASINI,DIST
      REAL VV,VIC,S,Q,SOMB,SOM,SOMD,FDEM,ERR_ASINI,ERR_FDEM
      DOUBLE PRECISION PP,T,PHI
      LOGICAL INTERACTIVE

      ITER=0
C
      PI=3.1415926535
      DIST=1.

C Calcul de la somme des poids: 
      SOM_POI=0.
      DO I=1,NN
        SOM_POI=SOM_POI+POI(I)
      END DO

C Boucle principale (ITER)
   89 DO J=1,NP
        B(J)=0
      END DO
      DO 9 J=1,NP
      DO 9 K=1,NP
        A(J,K)=0
    9 CONTINUE
C Sortie des RES(I) = O-C
C seulement pour les elements provisoires et pour la derniere iteration:
      IF(ITER.EQ.0.OR.ITER.EQ.NITER)THEN
        IF(INTERACTIVE)WRITE(LU_OUT_TXT,1600)
 1600   FORMAT('OBS.',2X,'DATE (JJ)',1X,'PHASE',3X,'VR OBS.',2X,'VR CALC.',
     1          3X,'O-C',4X,'POIDS')
      ENDIF
      DO 53 I=1,NN
C A voir plus tard:
C      PHI=FRAC(TOBS(I)-T)/PP)
C Calcul de la phase a partir de l'epoque d'observation et de la periode
C Phase PHI (entre 0 et 1):
      PHI=DMOD((TOBS(I)-T)/PP,1.D0)
      IF(PHI.LT.0)PHI=1.+PHI
C Calcul des coefficients des equations normales
C

C Calcul du residu:
      CALL VITESSE_CALCULEE(PHI,VIC,VV,EE,OMEGA,KA,VZERO)
      RES(I)=VIT(I)-VIC
C Sortie des RES(I) = O-C
C seulement pour les elements provisoires et pour la derniere iteration:
      IF(ITER.EQ.0.OR.ITER.EQ.NITER)THEN
        IF(INTERACTIVE) 
     1    WRITE(LU_OUT_TXT,1400)I,TOBS(I),PHI,VIT(I),VIC,RES(I),POI(I)
 1400 FORMAT(I3,3X,F9.3,2X,F4.3,3(2X,F7.2),3X,F5.2)
      ENDIF
C Critere de sortie de boucle:
      IF(ITER.EQ.NITER.OR.DIST.LT.1.E-7)GO TO 53
      S=SIN(VV+OMEGA)
      Q=S*(1+EE*COS(VV))**2/(1-EE**2)**1.5
C Correction de (-2.*PI/PP**2) par rapport a p83 (coeff "D"), 
C puisque qu'on utilise la periode PP au lieu de n (moyen mouvement)   
C X=1: periode fixe
      IF(X1.EQ.0)THEN
        AA(1,I)=(2.*PI/PP**2)*(TOBS(I)-T)*KA*Q
      ENDIF
C Coeff "F" p83
      AA(2-X1,I)=KA*Q*(2.*PI/PP)
C Y=2: orbite circulaire 
      IF(Y1.NE.2)THEN
C Coeff "B" p 83
        AA(3-X1,I)=-KA*(EE*SIN(OMEGA)+S)
C Y=0: ellipticite libre 
C Y=1: ellipticite fixe 
        IF(Y1.EQ.0)THEN
C Coeff "C" p 83
          AA(4-X1,I)=KA*(COS(OMEGA)-S*SIN(VV)*(2.+EE*COS(VV))/(1-EE**2))
        ENDIF
      ENDIF
C Coeff "A" p 83
      AA(5-X1-Y1,I)=EE*COS(OMEGA)+COS(VV+OMEGA)
      AA(6-X1-Y1,I)=1

C Calcul du 2eme membre de l'equation normale: POI(I) A^T RES
      DO J=1,NP
        B(J)=B(J)+POI(I)*AA(J,I)*RES(I)
      END DO
C Calcul du 1er membre de l'equation normale POI(I) A^T A
      DO 3 J=1,NP
      DO 3 K=J,NP
      A(J,K)=A(J,K)+POI(I)*AA(J,I)*AA(K,I)
    3 CONTINUE
C Fin de la boucle secondaire
   53 CONTINUE
C Test de convergence ou d'atteinte de NITER iterations
C (second test, en fin de boucle)
      IF(ITER.EQ.NITER.OR.DIST.LT.1.E-7)GO TO 81

C Symmetrise la matrice A 
C (puisque par construction elle est symmetrique, mais
C on n'avait pas encore affecte l'autre partie de cette matrice)
      DO 4 J=2,NP
        DO 4 K=1,J-1
         A(J,K)=A(K,J)
    4 CONTINUE

C Inversion de la matrice: 
C en sortie: D=matrice inverse
      CALL MAT_INV(A,D,NP)

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
C Calcul des erreurs avec les résidus precedents:
      DO I=1,NN
        SOMB=SOMB+POI(I)*RES(I)**2
      END DO
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
C Erreur estimee:
      SOME=SOMB-SOMD
      SOM=SQRT(SOME/(NN-NP))
      DO I=1,NP
        EQM(I)=SOM*SQRT(D(I,I))
      END DO

C Nouveaux elements:
      DIST=0
      F=0
      DO I=1,6
        DELTA(I)=0
      END DO
      DELTA(1)=X1
      DELTA(4)=Y1
        IF(Y1.EQ.2)THEN
          DELTA(3)=1
          DELTA(4)=1
        ENDIF

      DO 801 I=1,6
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
C JLP2002: j'ajoute une contrainte sur la positivite de EE:
       IF(EE+CORR(4).LT.0.)THEN
        EE=EE/2.
C Surtout ne pas rajouter OMEGA=OMEGA+PI, comme on pourrait le penser!
C (diverge avec HD100054B_fixed)
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
      ERR_ASINI=ASINI*(ERR(1)/PP + ERR(5)/KA + ERR(4)*EE/(1.-EE**2))
      ERR_FDEM=FDEM*(3.*EE*ERR(4)/(1.-EE**2) + ERR(1)/PP
     1         + 3.*ERR(5)/KA)
      ITER=ITER+1
      IF(INTERACTIVE)THEN
         WRITE(6,2200)ITER,DIST
         WRITE(LU_OUT_TXT,2200)ITER,DIST
 2200 FORMAT(3X,I2,'eme iteration',' dist=',E10.3)
         IF(X1.EQ.1) WRITE(LU_OUT_TXT,1700)
 1700 FORMAT(3X,'PERIODE FIXEE')
         IF(Y1.EQ.1) WRITE(LU_OUT_TXT,1800)
 1800 FORMAT(3X,'EXCENTRICITE FIXEE')
         IF(Y1.EQ.2) WRITE(LU_OUT_TXT,1801)
 1801 FORMAT(3X,'ORBITE CIRCULAIRE')
         WRITE(LU_OUT_TXT,1300)CORR(1),PP,ERR(1),CORR(2),T,ERR(2),
     1       CORR(3)*180./PI,OMEGA*180./PI,ERR(3)*180./PI,
     1       CORR(3),OMEGA,ERR(3),
     1       CORR(4),EE,ERR(4),CORR(5),KA,ERR(5),CORR(6),VZERO,ERR(6),
     1       ASINI,ERR_ASINI,FDEM,ERR_FDEM
 1300 FORMAT(1X,'CORRECTIONS',15X,'NOUVEAUX ELEMENTS',
     1   15X,'ERREURS STD.',/,
     1   3X,E9.3,17X,'P =',F15.7,13X,E9.3,/,
     1   3X,E9.3,17X,'T =',F11.3,17X,E9.3,/,
     1   3X,E9.3,13X,'OMEGA =',F13.5,15X,E9.3,4X,'(deg.)',/,
     1   3X,E9.3,13X,'OMEGA =',F13.5,15X,E9.3,4X,'(rad.)',/,
     1   3X,E9.3,17X,'E =',F13.5,15X,E9.3,/,
     1   3X,E9.3,17X,'K =',F10.2,18X,E9.3,/,
     1   3X,E9.3,16X,'V0 =',F10.2,18X,E9.3,/,
     1   24X,'A.SINI =',6X,E10.4,12X,E9.3,/,
     1   26X,'F(M) =',6X,E10.4,12X,E9.3,/,/)
      ENDIF
      IF(ABS(EE).GE.1.)THEN
C Erreur grave 
        WRITE(LU_OUT_TXT,2000)
        WRITE(6,2000)
 2000   FORMAT(5X,'ABS(EE).GE.1, POURSUITE DU CALCUL IMPOSSIBLE')
        CLOSE(LU_OUT_TXT)
        STOP
      ELSE
        GO TO 89
      ENDIF
C************************************************************************
C Sortie normale:
   81 RETURN
      END
