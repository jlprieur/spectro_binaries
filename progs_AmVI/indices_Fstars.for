C****************************************************************
C Photométrie de Stromgren pour les étoiles F
C
C Calcul des corrections de rougissement et divers
C autres paramètres: delta m1, delta c1
C en vue de calculer [Fe/H] et Mv
C
C From E. Oblak, 27 janvier 2004
C
C JLP
C Version 27/01/2004
C****************************************************************
      program indice_oblak
      real c1zac,deltc1,deltm1,am1hyb,deltmb,deczac,deltb
      real am1hy,c1z,cro_m1,cro_c1,cro_ub
      real by,am1,c1,beta,cons1,cons2
      integer ihd 
      character filename*40, buffer*80

      write(6,122)
122   format(" Photométrie de Stromgren pour les étoiles F"/,
     1  " Format du fichier des mesures:",/,
     1  " HD, b-y, m1, c1 et beta mesurés",/,
     1  " un entier et 4 réels, séparés par des blancs",/,
     1  " (NB: la première ligne est réservée pour les commentaires") 
      write(6,*) ' Entrez le nom du fichier de mesures:'
      read(5,*) filename
      open(1,file=filename,status='unknown')

C Premiere ligne avec commentaires
      read(1,'(A)') buffer
       print *,'Commentaires: ',buffer 

C Boucle sur les lignes suivantes:
    1  read(1,*,end=99,err=1010)ihd,by,am1,c1,beta
       print *,'OK: ', ihd,by,am1,c1,beta
       cons1=-0.04
       cons2=-0.025
C
C calcul des indices crochets (non utilisés ensuite)
C
       cro_m1=am1+0.3*by
       cro_c1=c1-0.2*by
       cro_ub=cro_c1+2.*cro_m1
       print 109,cro_m1,cro_c1,cro_ub
  109  format(' Indices crochets: cro_m1= ',F7.3,1x,
     1        'cro_c1= ',F7.3, 1x,'cro_ub= '
     1,F7.3)
C
C
C calcul de c1Z(ZAMS) et de MVZAMS (Crawford 1975 Astron J 60, 955
C (étoiles F, p960)

       if(beta.gt.2.74) then
          c1z=1.93*beta-4.63
       else if(beta.le.2.74.and.beta.gt.2.72) then
          c1z=3.*beta-7.56
       else if(beta.le.2.72.and.beta.gt.2.69) then
          c1z=2.85*beta-7.168
       else if(beta.le.2.69.and.beta.gt.2.65) then
          c1z=2.6*beta-6.501
       else if(beta.le.2.65) then
          c1z=2.0*beta-4.91
       endif
C 
C m1hy(m1 pour les hyades Crawford 1970)
C c1zac(c1 pour la ZAMS Crawford 1970)
C

C
C calcul de m1 (Hyades) a partir des tables de Crawford et Barnes 1970
C et de Crawford 1975 pour les étoiles F (2.72 < beta < 2.50)
C
       if(beta.gt.2.855)then
         am1hy=0.207
       else if(beta.le.2.855.and.beta.gt.2.805)then
         am1hy=0.08*beta-0.2
       else if(beta.le.2.805.and.beta.gt.2.74)then
         am1hy=0.371*beta-0.837
       else if(beta.le.2.74.and.beta.gt.2.7)then
         am1hy=0.23*beta-0.449
       else if(beta.le.2.7.and.beta.gt.2.65)then
         am1hy=0.974*beta**(-1.758)
       else if(beta.le.2.65.and.beta.gt.2.40)then
         am1hy=4219.836*beta**(-10.349)
       else
         goto 20
       endif


       c1zac=0.81209*beta**2-1.99*beta
       deltc1=c1-c1z
       deltm1=am1hy-am1
       am1hyb=0.31788*beta**2-2.0548*beta+3.4
       deltmb=am1hyb-am1
       deczac=c1-c1zac
       deltb=2.72-beta
C
C corrections de Stromgren sur c1 et m1
C
       if(by.gt.0.23)then
         c1cs=c1+0.7*deltm1
         corc1=0.7*deltm1
         am1cs=am1
         corm1=9.99
       else if(by.le.0.23.and.by.gt.-0.05)then
         am1cs=am1+0.1*deltc1
         corm1=0.1*deltc1
         c1cs=c1
         corc1=9.99
       else
         goto 20
       endif

       deltc1s=c1cs-c1z
       deltm1s=am1hy-am1cs
       am1by=2.5*by**2-1.31514*by+0.34436
       dm1byc=am1by-am1
       do 423 i=1.5
         byc=by-0.5*dm1byc
         am1by=2.5*byc**2-1.31514*byc+0.34436
         dm1byc=am1by-am1cs
423    continue

       do 424 i=1,5

       if(byc.le.0.27)then
          c1zbyc=-2.4*byc+1.1165
       else if(byc.gt.0.27.and.byc.le.0.313)then
         c1zbyc=-1.851*byc+0.966
       else if(byc.gt.0.313)then
         c1zbyc=-1.211*byc+0.767
       else
         goto 20
       endif

       dc1byc=c1cs-c1zbyc
       byc=0.222+1.11*deltb+2.7*deltb**2-0.05*dc1byc-(0.1+3.6*deltb)*
     1     dm1byc
  424 continue

       byc0=byc
       Ebyc0=by-byc0
       am1byc=am1+0.3*Ebyc0
       c0byc=c1-0.25*Ebyc0
       avbyc=4.3*Ebyc0

       if(byc0.le.0.266)then
         amzbyc=10.86*byc0+0.424
       else if(byc0.gt.0.266.and.byc0.le.0.377)then
         amzbyc=9.877*byc0+0.799
       else if(byc0.gt.0.377)then
         amzbyc=6.671*byc0+1.657
       endif

       dcbyc0=dcbyc

   20 continue
C
C relation de Crawford 1975 astron J 80 955 (p 962, for F-type stars)
C
       by0cst=0.222+1.11*deltb+2.7*deltb**2-0.05*delc1s-(0.1+3.6*deltb)*
     1 deltm1s
       by0c=0.222+1.11*deltb+2.7*deltb**2-0.05*deltc1-(0.1+3.6*deltb)*
     1 deltm1
       c0=c1-0.25*(by-by0c)
       am0=am1+0.32*(by-by0c)
       Eby0=by-by0c
       avby0=4.3*Eby0
       avbycs=4.3*Eby0cs
       if(am0.gt.0.22) print 271
  271 format(1h ,'Etoile metallique')

      if(Eby0.le.cons1) print 274
  274 format(1h ,'Ap ou Fp')

      if((c0.gt.1.19).or.(deltc1.ge.0.28)) print 571
  571 format(1h ,'etoile a grande luminosite')

      if(deltm1.ge.0.01)then
        print 572
  572 format(1h ,'supergeante')
      else if(deltm1.ge.0.05.and.Eby0.ge.0.1)then
        print 572
      else if(by0c.le.0.08)then
        print 573
  573 format(1h ,'etoile du groupe intermediaire possible')
      endif

       by0s=2.946-0.15*deltc1-beta
       Eby0s=by-by0s
       Eby0cs=by-by0cst
C
C relation de Crawford pour les etoiles A (2.72<beta<2.89) pour (b-y)0
C
      by0a=2.943-beta-0.1*(deltc1+deltm1)
      Eby0a=by-by0a
      avbya=4.3*Eby0a
      c1crs=c1cs-0.2*by
      c1crz=c1z-0.2*by0a
      dc1cr=c1crs-c1crz
C
C relation de Crawford 1979, Astron. J, 84, 1858 pour les etoiles A
C
       if(deltm1.lt.0.)then
         bycr=2.946-beta-0.1*deltc1-0.25*deltm1
       else if(deltm1.ge.0.)then
         bycr=2.946-beta-0.1*deltc1
       endif
       Ebycr=by-bycr
       avbycr=4.3*Ebycr
       print 888
  888  format(1h ,
     1'(b-y)c',2x,'HD',9x,'(b-y)0',2x,'beta',4x,'am10',4x,'c0',5x,
     1'E(b-y)0',1x,
     2'Av(b-y)0',6x,'c1Z',3x,'deltam1',1x,'deltac1',1x,'deltamb',
     31x,'deczac')
       print 703,ihd,by0c,beta,am0,c0,Eby0,avby0,c1z,deltm1,deltc1,
     1 deltmb,deczac
  703 format(1h ,'by0c',2x,i8,2x,6f8.3,4x,6F8.3)
      print 704,ihd,byc0,beta,am1byc,c0byc,Ebyc0,avbyc,dm1byc,dc1byc
  704 format(1h ,'byc0',2x,i8,2x,6F8.3,20x,2F8.3)
      print 705,ihd,by0cst,beta,am1cs,c1cs,Eby0cs,avbycs,deltm1s,delc1s,
     1 corm1,corc1
  705 format(1h ,'by0cst',i8,2x,6F8.3,20x,2F8.3,4x,2F8.3)
      print 706,ihd,by0s,beta,Eby0s
  706 format(1h ,'by0s',2x,i8,2x,F8.3,1x,F8.3,16x,F8.3)
      print 707,ihd,by0a,beta,Eby0a
  707 format(1h ,'by0a',2x,i8,2x,F8.3,1x,F8.3,16X,F8.3)
      print 708,ihd,bycr,beta,Ebycr
  708 format(1h ,'bycr',2x,i8,2x,F8.3,1x,F8.3,16x,F8.3)
      GOTO 1
C Error opening file:
   99 continue
 1010 continue
C 999 format(1h ,'erreur')
      stop
      end
