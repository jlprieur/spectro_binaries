C****************************************************************
C Calcul de l'influence du compagnon d'une étoile Am 
C sur la photométrie
C Attention, les etoiles Am n'ont pas forcement la meme correspondance
C T, B-V que dans Flower (1996), donc on procede par Delta T, Delta B-V
C
C JLP
C Version 05/02/2004
C****************************************************************
      program compagnon 
      implicit none
      real B_V,B_V_1, B_V_2,delta_B,delta_V,delta_B_V_1,delta_T1
      real B_V_1_Flo,T1_Am,T1_Flo_0,T1_Flo,T2_Flo,T2_T1,T1
      integer iter,i,k,n_Flo,ihd 
      real B_V_Flo(21), T_Flo(21)
      character filename*40, buffer*80

C Table 3 de Flower 1996, ApJ 469 355 
      data B_V_Flo/0.20, 0.22, 0.24, 0.26, 0.28, 0.30, 0.32, 0.34,
     1             0.36, 0.38, 0.40, 0.42, 0.44, 0.46, 0.48, 0.50,
     1             0.52, 0.54, 0.56, 0.58, 0.60/
      data T_Flo/7766, 7648, 7535, 7426, 7319, 7216, 7113, 7014,
     1           6916, 6820, 6725, 6632, 6541, 6453, 6366, 6282,
     1           6198, 6117, 6040, 5964, 5891/ 
       
      n_Flo = 21 

      write(6,122)
122   format(" Influence du compagnon pour une étoile A"/,
     1  " Format du fichier des mesures:",/,
     1  " HD, B-V global, Delta V, T1, T2/T1",/,
     1  " (NB: la première ligne est réservée pour les commentaires") 
C      write(6,*) ' Entrez le nom du fichier de mesures:'
C      read(5,*) filename
      filename='compagnon.dat'
      open(1,file=filename,status='old')
C Logfile:
      open(2,file="compagnon.log",status='unknown')

C Premiere ligne avec commentaires
      read(1,'(A)') buffer
       print *,'Commentaires: ',buffer 

C Boucle sur les lignes suivantes:
      do 111 iter=1,1000
       read(1,*,end=99,err=1010)ihd,B_V,delta_V,T1_Am,T2_T1
       print *,'***************************************************'
       print *,'Nouvel objet: HD ',ihd,B_V,delta_V,T1_Am,T2_T1,B_V_2

C Initialisation:
       B_V_1 = B_V
       T1 = T1_Am
C Calcul de T1_Flo_0
       call interpol_incr(B_V_Flo,T_Flo, n_Flo, B_V, T1_Flo_0)
       write(6,41) T1_Flo_0,B_V_1 
41     format(' Interpolation: T1_Flo_0=',F6.1,' B_V_1=',F5.3)

C Iterations:
C delta_V = V_2 - V_1
C delta_B = B_2 - B_1
       do 35 k=1,8
C Calcul de (B_V)2 a partir de T2_Flo et de Flower 1979 
         T2_Flo = T1 * T2_T1
         call interpol_decr(T_Flo, B_V_Flo, n_Flo, T2_Flo, B_V_2)
         write(6,43) T2_Flo,B_V_2 
43       format(' Interpolation: T2_Flo=',F6.1,' B_V_2=',F5.3)

C Calcul de (B-V)1
         delta_B = delta_V + B_V_2 - B_V_1
         B_V_1 = B_V - 2.5 * log10((1. + 10.**(-0.4 * delta_V))/
     1                             (1. + 10.**(-0.4 * delta_B)))
C Calcul de T1_Flo
         call interpol_incr(B_V_Flo,T_Flo, n_Flo, B_V_1, T1_Flo)
         write(6,44) T1_Flo,B_V_1 
44       format(' Interpolation: T1_Flo=',F6.1,' B_V_1=',F5.3)
         delta_B_V_1 = B_V - B_V_1
         delta_T1 = T1_Flo - T1_Flo_0 
         T1 = T1_Am + delta_T1
         write(6,*) k,'eme iteration: T1 =',T1
35     continue

C Sortie des resultats:
       write(6,888)
       if(iter.eq.1) write(2,888)
  888  format(4x,'HD',5x,'B-V',2x,'T2/T1',2x,'T1',5x,
     1          'T2',4x,'(B-V)1',1x,'(B-V)2')
       write(6,703) ihd,B_V,T2_T1,T1,T2_Flo,B_V_1,B_V_2
       write(2,703) ihd,B_V,T2_T1,T1,T2_Flo,B_V_1,B_V_2
  703 format(i8,1x,f6.3,1x,f4.2,1x,2(f6.1,1x),2(f6.3,1x))

C End of loop on iter:
111   continue
C Error opening file:
   99 continue
 1010 continue
      write(6,*) " Sortie des résultats dans ""compagnon.log"""
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
C************************************************************************
C
C INPUT:
C tab_x(nn): tableau en ordre croissant!
C tab_y(nn)
C nn
C xx: valeur de x
C
C OUTPUT:
C  yy: valeur de y interpolée à partir de tab_x et tab_y
C
C************************************************************************
       subroutine interpol_incr(tab_x, tab_y, nn, xx, yy)
       implicit none
       integer nn,i,i0
       real tab_x(nn), tab_y(nn), xx, yy
       real x1, x2, x3, y1, y2, y3, aa, bb, cc

       i0 = 0
       do 12 i=2,nn
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

C       write(6,*)'xx=',xx,' x1, x2, x3',x1,x2,x3
C       write(6,*)'yy=',yy,' y1, y2, y3',y1,y2,y3
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
