      PROGRAMBS1
      IMPLICITNONE
      REALPERIOD,SIGMA
      CHARACTEROBJECT*40
      WRITE(6,*)' BS1            JLP/Version 20/06/2005'
      WRITE(6,*)' Titre (sans blancs, car utilisé comme préfixe):'
      READ(5,10)OBJECT
   10 FORMAT(A)
      PERIOD=-1.
      CALLBS1_PROCESS(OBJECT,PERIOD,SIGMA)
      WRITE(6,*)'PERIOD = ',PERIOD,' SIGMA = ',SIGMA
      STOP
      END
