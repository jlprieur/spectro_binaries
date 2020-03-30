C****************************************************************
C Photométrie de Stromgren pour des étoiles A
C
C Calcul des corrections de rougissement et divers
C autres paramètres: delta m1, delta c1
C en vue de calculer [Fe/H] et Mv
C
C From E. Oblak, 27 janvier 2004 (sa version est pour les étoiles F)
C
C JLP
C Version 27/01/2004
C****************************************************************
      program indice_oblak
      real deltc1,deltm1,deltmb,deltb
      real am1hy,c1z,cro_m1,cro_c1,cro_ub
      real by,am1,c1,beta,cons1,cons2
      integer ihd 
      character filename*40, buffer*80
      integer nn_a
      real beta_ZAMS(9), c1_ZAMS(9), m1_ZAMS(9)
C Calcul de c1Z (ZAMS) Crawford 1979, AJ, 84, 1858, table 1 de beta et c1 (page 1859)
      data beta_ZAMS/2.88, 2.86, 2.84, 2.82, 2.80, 2.78, 2.76, 2.74, 2.72/
      data c1_ZAMS/0.93, 0.89, 0.85, 0.82, 0.78, 0.74, 0.70, 0.66, 0.60/
      data m1_ZAMS/0.200, 0.205, 0.208, 0.206, 0.203, 0.196, 0.188, 
     1          0.182, 0.177/
       
      nn_a = 9

      write(6,122)
122   format(" Photométrie de Stromgren pour des étoiles A"/,
     1  " Format du fichier des mesures:",/,
     1  " HD, b-y, m1, c1 et beta mesurés",/,
     1  " un entier et 4 réels, séparés par des blancs",/,
     1  " (NB: la première ligne est réservée pour les commentaires") 
C      write(6,*) ' Entrez le nom du fichier de mesures:'
C      read(5,*) filename
      filename='photom.dat'
      open(1,file=filename,status='old')

C Premiere ligne avec commentaires
      read(1,'(A)') buffer
       print *,'Commentaires: ',buffer 

C Boucle sur les lignes suivantes:
    1  read(1,*,end=99,err=1010)ihd,by,am1,c1,beta
       print *,'***************************************************'
       print *,'New object: HD ', ihd,by,am1,c1,beta
       cons1=-0.04
       cons2=-0.025
C
C calcul des indices crochets (non utilisés ensuite)
C
       cro_m1=am1+0.3*by
       cro_c1=c1-0.2*by
       cro_ub=cro_c1+2.*cro_m1
C       print 109,cro_m1,cro_c1,cro_ub
  109  format(' Indices crochets: cro_m1= ',F7.3,1x,
     1        'cro_c1= ',F7.3, 1x,'cro_ub= ',F7.3)
C
C
C Calcul de c1Z (ZAMS) Crawford 1979, AJ, 84, 1858, table 1 p 1859)
       call interpol_decr(beta_ZAMS, c1_ZAMS, nn_a, beta, c1z)
       write(6,*) ' Interpolation: beta=',beta,' c1=',c1z

C 
C Calcul de am1hy: m1 (Hyades) a partir de Crawford 1979 pour les étoiles A
       call interpol_decr(beta_ZAMS, m1_ZAMS, nn_a, beta, am1hy)
       write(6,*) ' Interpolation: beta=',beta,' m1=',am1hy

C Delta c1 = c10 (observé, corrigé) - c1(standard)
       deltc1=c1-c1z
C Delta m1 = m1(standard) - m10 (observé, corrigé du rougissement)
       deltm1=am1hy-am1
C
C Corrections de Stromgren sur c1 et m1
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

C Delta c1 Stromgren (i.e., corrigé par Stromgren)
       deltc1s=c1cs-c1z
C Delta m1 Stromgren (i.e., corrigé par Stromgren)
       deltm1s=am1hy-am1cs

       am1by=2.5*by**2-1.31514*by+0.34436
       dm1byc=am1by-am1
       do 423 i=1.5
         byc=by-0.5*dm1byc
         am1by=2.5*byc**2-1.31514*byc+0.34436
         dm1byc=am1by-am1cs
423    continue

C Delta beta (Crawford, p 962, for F stars)
       deltb=2.72-beta

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

   20 continue
C
C relation de Crawford 1979 astron J 84 1858 (p 1860, for A-type stars)
C
C by0cst = (b-y)0 corrigé, Stromgren
       if(deltm1s.lt.0)then
         by0cst = 2.946 - beta - 0.1 * deltc1s - 0.25 * deltm1s
       else
         by0cst = 2.946 - beta - 0.1 * deltc1s 
       endif
C by0c = (b-y)0 corrigé:
       if(deltm1s.lt.0)then
         by0c = 2.946 - beta - 0.1 * deltc1 - 0.25 * deltm1
       else
         by0c = 2.946 - beta - 0.1 * deltc1
       endif
       c0=c1-0.25*(by-by0c)
       am0=am1+0.32*(by-by0c)

C E(b-y)0
       Eby0=by-by0c
       avby0=4.3*Eby0

c E(b-y)0 with Stromgren correction:
       Eby0cs=by-by0cst
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

C Sortie des resultats:
       print 888
  888  format(' (b-y)c',4x,'HD',4x,'(b-y)0',2x,'beta',3x,'am10',
     1        4x,'c0',2x,'E(b-y)0',1x,
     2'Av(b-y)0',1x,'c1Z',2x,'deltam1',1x,'deltac1',1x,'deltamb')
       print 703,ihd,by0c,beta,am0,c0,Eby0,avby0,c1z,deltm1,deltc1
  703 format(' by0c',2x,i8,1x,9(1x,f6.3))
      print 705,ihd,by0cst,beta,am1cs,c1z,Eby0cs,avbycs,c1zs,
     1 deltm1s,deltc1s
  705 format(' by0cst',i8,1x,9(1x,F6.3))
      GOTO 1
C Error opening file:
   99 continue
 1010 continue
C 999 format(1h ,'erreur')
      stop
      end
C************************************************************************
C
C INPUT:
C tab_x(nn): tableau en ordre décroissant!
C tab_y(nn)
C nn
C xx: valeur de x
C
C OUTPUT:
C  yy: valeur de y interpolée à partir de tab_x et tab_y
C
C************************************************************************
       subroutine interpol_decr(tab_x, tab_y, nn, xx, yy)
       implicit none
       integer nn,i,i0
       real tab_x(nn), tab_y(nn), xx, yy
       real x1, x2, x3, y1, y2, y3, aa, bb, cc

       i0 = 0
       do 12 i=nn,2,-1
          if(xx.gt.tab_x(i))then
            i0 = i
            goto 13
          endif
12     continue

13     if(i0.eq.0)then
        write(6,*)' Erreur fatale: xx=',xx,' en dehors du tableau'
        stop
       endif

C Cas du bord supérieur:
       if(i0.eq.nn) i0 = i0-1

       x1 = tab_x(i0-1)
       x2 = tab_x(i0)
       x3 = tab_x(i0+1)
       y1 = tab_y(i0-1)
       y2 = tab_y(i0)
       y3 = tab_y(i0+1)
       call parfit(x1, x2, x3, y1, y2, y3, aa, bb, cc)
       yy = aa * xx * xx + bb * xx + cc

       return
       end
C -------------------------------------------------------------------------
C Ajustement d'une parabole y=a x^2 + b x + c
C a partir de 3 couples de points (x1,y1), (x2,y2) et (x3,y3)
C (From Eric Aristidi)
C ------------------------------------------------------------------------- 
      subroutine parfit(x1, x2, x3, y1, y2, y3, aa, bb, cc)
      implicit none
      real x1, x2, x3, y1, y2, y3, aa, bb, cc
      real dd, da, db

       dd=(x3*x3-x1*x1)*(x2-x1)-(x2*x2-x1*x1)*(x3-x1)
       da=(y3-y1)*(x2-x1)-(y2-y1)*(x3-x1)
       db=(x3*x3-x1*x1)*(y2-y1)-(x2*x2-x1*x1)*(y3-y1)
      aa=da/dd
      bb=db/dd
      cc=y1-aa*x1*x1-bb*x1
      return
      end


