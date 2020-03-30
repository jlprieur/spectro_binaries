C************************************************************************
C Main program to call "BS1_PROCESS" routine which does the job
C
C JLP
C Version 30/05/2005
C************************************************************************
      PROGRAM BS1
      IMPLICIT NONE
      REAL PERIOD,SIGMA
      CHARACTER OBJECT*40

      WRITE(6,*)' BS1            JLP/Version 20/06/2005'
      WRITE(6,*)' Titre (sans blancs, car utilisé comme préfixe):'
      READ(5,10) OBJECT 
10    FORMAT(A)
      PERIOD=-1.
       CALL BS1_PROCESS(OBJECT,PERIOD,SIGMA)
      WRITE(6,*) 'PERIOD = ',PERIOD,' SIGMA = ',SIGMA

      STOP
      END
