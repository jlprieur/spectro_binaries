C****************************************************************
C Photométrie de Stromgren pour des étoiles A
C
C Calcul des corrections de rougissement et divers
C autres paramètres: delta m1, delta c1
C en vue de calculer [Fe/H] et Mv
C
C From E. Oblak, 27 janvier 2004 (sa version est pour les étoiles F)
C 
C Input: photom.dat
C
C JLP
C Version 30/01/2004
C****************************************************************
      program indice_Astars
      implicit none
      real deltc1,deltm1,deltm1s,deltmb,deltb
      real am1_hy,c1_hy,cro_m1,cro_c1,cro_ub
      real by,am1,c1,beta,cons1,cons2,c1s,am1s
      real am1_0,c1_0,Eby0,Eby0cs,avby0,avbycs
      real am1_0s,c1_0s,delta1s,deltc1s
      real by0cs,by0c,Mv,Mvs,Mv_hy,FeH,FeHs,FeH_hy
      integer ihd,iter,i,k 
      character filename*40, buffer*80
      integer nn_a
      real beta_ZAMS(9), c1_ZAMS(9), m1_ZAMS(9), Mv_ZAMS(9)
C Table 1 de Crawford 1979, AJ, 84, 1858, beta, c1, m1 et Mv (page 1859)
C from observations of ZAMS stars of the Hyades
      data beta_ZAMS/2.88, 2.86, 2.84, 2.82, 2.80, 2.78, 
     1               2.76, 2.74, 2.72/
      data c1_ZAMS/0.93, 0.89, 0.85, 0.82, 0.78, 0.74, 
     1             0.70, 0.66, 0.60/
      data m1_ZAMS/0.200, 0.205, 0.208, 0.206, 0.203, 0.196, 0.188, 
     1             0.182, 0.177/
      data Mv_ZAMS/2.30, 2.50, 2.64, 2.70, 2.76, 2.82, 2.88, 2.96,
     1             3.10/
       
      nn_a = 9

      write(6,122)
122   format(" Photométrie de Stromgren pour des étoiles A"/,
     1  " Format du fichier des mesures:",/,
     1  " HD, b-y, m1, c1 et beta mesurés",/,
     1  " un entier et 4 réels, séparés par des blancs",/,
     1  " (NB: la première ligne est réservée pour les commentaires") 
C      write(6,*) ' Entrez le nom du fichier de mesures:'
C      read(5,*) filename
      write(6,*) ' Lit le fichier de mesures photom.dat'
      filename='photom.dat'
      open(1,file=filename,status='old')
C Logfile:
      open(2,file="indice.log",status='unknown')

C Premiere ligne avec commentaires
      read(1,'(A)') buffer
       print *,'Commentaires: ',buffer 

C Boucle sur les lignes suivantes:
      do 111 iter=1,1000
       read(1,*,end=99,err=1010)ihd,by,am1,c1,beta
       print *,'***************************************************'
C HD, (b-y), m_1, c_1, beta:
       print *,'Nouvel objet: HD ', ihd,by,am1,c1,beta
       cons1=-0.04
       cons2=-0.025
C
C calcul des indices crochets (non utilisés ensuite)
C
       cro_m1=am1+0.3*by
       cro_c1=c1-0.2*by
       cro_ub=cro_c1+2.*cro_m1
C       print 109,cro_m1,cro_c1,cro_ub
  109  format(' Indices crochets: cro_m1= ',F6.2,1x,
     1        'cro_c1= ',F6.2, 1x,'cro_ub= ',F6.2)

C Calcul de c1_hy (ZAMS), à partir du beta mesuré, en 
C utilisant la calibration de Crawford 1979, AJ, 84, 1858, table 1 p 1859)
C Interpolation des valeurs de la table (beta_ZAMS, c1_ZAMS)
C pour calculer c1_hy(beta) (Hyades):
       call interpol_decr(beta_ZAMS, c1_ZAMS, nn_a, beta, c1_hy)

C Calcul de am1_hy: m1 (Hyades) à partir du beta mesuré
C en utilisant la calibration de Crawford 1979 pour les étoiles A
C Interpolation des valeurs de la table (beta_ZAMS, m1_ZAMS)
C pour calculer am1_hy(beta) (Hyades):
       call interpol_decr(beta_ZAMS, m1_ZAMS, nn_a, beta, am1_hy)

C Calcul de Mv_hy: Mv (Hyades) a partir de Crawford 1979 pour les étoiles A
C Interpolation des valeurs de la table (beta_ZAMS, m1_ZAMS)
C pour calculer Mv_hy(beta) (Hyades):
       call interpol_decr(beta_ZAMS, Mv_ZAMS, nn_a, beta, Mv_hy)
       write(6,43) beta, c1_hy, am1_hy, Mv_hy 
43     format(' Interpolation: beta_ZAMS=',F5.2,' c1_ZAMS=',F5.2,
     1        ' m1_ZAMS=',F5.2,' Mv_ZAMS=',F5.2)

C Delta c1 = c10 (observé, corrigé) - c1(standard)
       deltc1=c1-c1_hy
C Delta m1 = m1(standard) - m10 (observé, corrigé du rougissement)
       deltm1=am1_hy-am1

C
C relation de Crawford 1979 astron J 84 1858 (p 1860, for A-type stars)
C
C by0c = (b-y)0 corrigé du rougissement, sans corrections de Stromgren
       if(deltm1.lt.0)then
         by0c = 2.946 - beta - 0.1 * deltc1 - 0.25 * deltm1
       else
         by0c = 2.946 - beta - 0.1 * deltc1
       endif
       write(6,*)'1ere iteration: E(b-y)=',by-by0c

C Crawford, 1975, AJ 80, 955  page 957 en haut à gauche
       c1_0 = c1 - 0.2 * (by-by0c)
       am1_0 = am1 + 0.3 * (by-by0c)

C Iterations for finding solution compatible with the two relations:
       do 35 k=1,8
          deltc1= c1_0 - c1_hy
          deltm1= am1_hy - am1_0
          if(deltm1.lt.0)then
            by0c = 2.946 - beta - 0.1 * deltc1 - 0.25 * deltm1
          else
            by0c = 2.946 - beta - 0.1 * deltc1
          endif

C Constraint: color excess E(b-y) cannot be negative
          if(by - by0c .lt. 0.) by0c = by
          c1_0 = c1 - 0.2 * (by-by0c)
          am1_0 = am1 + 0.3 * (by-by0c)
35     continue
       write(6,*) k-1,'eme iteration: E(b-y)=',by-by0c

C Color excess E(b-y) (for reddening correction)
       Eby0=by-by0c
       avby0=4.3*Eby0

C Stromgren (1963)'s corrections of c1 and m1 for absolute magnitudes 
       call stromgren_correc(c1_0,am1_0,deltm1,deltc1,by0c,
     1                       c1_0s,am1_0s)
C Change c1 and am1 with c1s and am1s to take into account this correction
C and find another E(b-y) (JLP2007):
       c1s = c1_0s + 0.2 * (by-by0c)
       am1s = am1_0s - 0.3 * (by-by0c)

C Iterations for finding a solution compatible with
C Stromgren's and Crawford's relations: 

      do 36 k=1,8

C Delta c1 Stromgren (i.e., corrigé par Stromgren)
       deltc1s = c1_0s - c1_hy
C Delta m1 Stromgren (i.e., corrigé par Stromgren)
       deltm1s = am1_hy - am1_0s

C by0cs = (b-y)0 corrigé du rougissement, avec corrections de Stromgren
       if(deltm1s.lt.0)then
         by0cs = 2.946 - beta - 0.1 * deltc1s - 0.25 * deltm1s
       else
         by0cs = 2.946 - beta - 0.1 * deltc1s 
       endif
C Constraint: color excess E(b-y) cannot be negative
       if(by - by0cs .lt. 0.) by0cs = by

C Color excess E(b-y) (for correction of reddening), 
C with Stromgren correction:
       Eby0cs = by - by0cs
       avbycs = 4.3 * Eby0cs
C       write(6,*) k-1,'eme iteration: E(b-y)stromgren=',by-by0c
C       write(6,*) 'am1_0s= ',am1_0s,' c1_0s= ',c1_0s

C Crawford, 1975, AJ 80, 955
       c1_0s = c1s - 0.2 * (by-by0cs)
       am1_0s = am1s + 0.3 * (by-by0cs)

36     enddo
       write(6,*) k-1,'eme iteration: E(b-y)stromgren=',by-by0c

C Fe/H (Crawford, 1975, AJ 80, 955) p969, from Cayrel 
       FeH = 0.20 - 10. * deltm1
       FeHs = 0.20 - 10. * deltm1s

C Mv (Crawford, 1979, AJ 84, 1858) for A stars, p1862
       Mv = Mv_hy - 9. * deltc1
       Mvs = Mv_hy - 9. * deltc1s

       if(am1_0.gt.0.22) print 271
  271 format(1h ,'Etoile metallique')

      if(Eby0.le.cons1) print 274
  274 format(1h ,'Ap ou Fp')

      if((c1_0.gt.1.19).or.(deltc1.ge.0.28)) print 571
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
       write(6,888)
       if(iter.eq.1) write(2,888)
  888  format(' L''indice 0 est ajoute pour les valeurs corriges du',
     1        ' rougissement.',/,
     1        ' Pour chaque objet il y a deux lignes: ',/,
     1        '   by0c: sans correction de Stromgren(1963)',/,
     1        '   by0cs: avec correction de am1 et c1 ',
     1        'selon Stromgren(1963)',/,
     1        ' (utile pour la détermination de Mv',
     1        ' avec les indices de Stromgren)',
     1        /,'(b-y)c',4x,'HD',2x,'beta',4x,'b-y',2x,'(b-y)0',
     1        3x,'m1',4x,'(m1)0',2x,'c1',2x,'(c1)0',1x,'E(b-y)',
     1        2x,'Av',4x,'dm1',3x,'dc1',3x,'Mv',3x,'Fe/H')
C Without Stromgren corrections (not well suited to Mv determination)
       write(6,703) "by0c ",ihd,beta,by,by0c,am1,am1_0,c1,c1_0,Eby0,
     1              avby0,deltm1,deltc1,Mv,FeH
       write(2,703) "by0c ",ihd,beta,by,by0c,am1,am1_0,c1,c1_0,Eby0,
     1              avby0,deltm1,deltc1,Mv,FeH
  703 format(a4,1x,i8,2(f6.3,1x),f5.2,1x,f6.3,1x,f5.2,
     1       3(1x,f6.3),4(1x,f5.2))
C With Stromgren corrections (well suited to Mv determination)
       write(6,703) "by0cs",ihd,beta,by,by0cs,am1s,am1_0s,c1s,c1_0s,
     1              Eby0cs,avbycs,deltm1s,deltc1s,Mvs,FeHs
       write(2,703) "by0cs",ihd,beta,by,by0cs,am1s,am1_0s,c1s,c1_0s,
     1              Eby0cs,avbycs,deltm1s,deltc1s,Mvs,FeHs

C End of loop on iter:
111   continue
C Error opening file:
   99 continue
 1010 continue
      write(6,*) " Sortie des résultats dans ""indice.log"""
      close(2)
      close(1)
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
          if(xx.lt.tab_x(i))then
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
C ------------------------------------------------------------------------- 
C Corrections de Stromgren (1963) de c1 et m1 
C Surtout utile pour les magnitudes absolues ? 
C
C INPUT:
C c1_0, am1_0: c1 et m1 corriges du rougissement
C deltm1_0, deltc1_0: delta m1 et delta c1, corriges du rougissement 
C by_0: b - y corrige du rougissement 
C
C OUPUT:
C c1_0s, am1_0s: c1 et m1 corriges selon Stromgren (1963)
C ------------------------------------------------------------------------- 
      subroutine stromgren_correc(c1_0,am1_0,deltm1_0,deltc1_0,by_0,
     1                            c1_0s,am1_0s)
      implicit none
      real c1_0,am1_0,deltm1_0,deltc1_0,c1_0s,am1_0s,by_0


C For Mv determination, calibration on the Hyades of late F to early G stars
C suggest a correction on c1:
C Stromgren, 1963, p23: if 0.25 < b - y < 0.4,  c1 = c1 + 0.75 deltm1
       if(by_0.gt.0.23.and.by_0.lt.0.42)then
         c1_0s = c1_0 + 0.75 * deltm1_0
         am1_0s = am1_0
C For Mv determination, the absorption line Hbeta introduce some
C effects that should be corrected.
C Stromgren, 1963, p20: if 0 < b - y < 0.2,  m1 = m1 + 0.1 deltc1
       else if(by_0.le.0.23.and.by_0.gt.-0.05)then
         am1_0s = am1_0 + 0.1 * deltc1_0
         c1_0s = c1_0
       else
         am1_0s = am1_0 
         c1_0s = c1_0 
       endif

       return
       end
