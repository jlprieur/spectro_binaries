/**************************************************************************
* To compute the characterictic times, t_sync and t_circ
* for synchronization and circularisation of the orbits
* in the case of radiative enveloppe
* using the formulae of Zahn (1975)
*
* cc tsync.c -lm
* a.out .....
*
* JLP
* Version 23/03/2005
**************************************************************************/
#include<stdio.h>
#include<math.h>

/* 
#define DEBUG
*/

/* Gravitation constant: */
static double G=6.67e-11;
/* Solar mass in kg */
static double Msol=1.989e30;
/* Solar radius in m */
static double Rsol=6.9599e8;

int compute_imr2(double mm, double *mr2);
int compute_e2(double mm, double *e2);
int compute_tsync(double mm, double qq, double e2, double rr, double aa,
                  double mr2, double *t_sync);
int compute_tcirc(double mm, double qq, double e2, double rr, double aa,
                  double mr2, double *t_circ);

int main(int argc, char *argv[])
{
double mm, qq, e2, rr, aa, mr2, t_sync, t_circ;

printf(" Program to compute t_sync and t_circ with Zahn's formulae (1975, A&A 41, 329)\n");
#ifdef DEBUG
 printf("argc=%d\n", argc);
 printf("argv[1]=%s\n", argv[1]);
#endif

if(argc != 5) {
  printf(" Usage: \n tsync M q R a \n");
  printf(" M = stellar mass (Msol) \n");
  printf(" q = mass ratio \n");
  printf(" R = radius (Rsol) \n");
  printf(" a = semi-major axis (Rsol) \n");
  exit(-1);
  }
else {
  sscanf(argv[1],"%lf",&mm);
  sscanf(argv[2],"%lf",&qq);
  sscanf(argv[3],"%lf",&rr);
  sscanf(argv[4],"%lf",&aa);
  }
printf("OK: M=%.3f Msol;  q=%.3f;  R=%.3f Rsol;  a=%.3f Rsol; \n",
       mm, qq, rr, aa);

/* Compute I/MR2 corresponding to mm Msol (stellar mass) */
compute_imr2(mm, &mr2);
printf(" I/MR2 = %e \n",mr2);

/* Compute E_2 corresponding to mm Msol (stellar mass) */
compute_e2(mm, &e2);
printf(" E_2 = %e \n",e2);

compute_tsync(mm, qq, e2, rr, aa, mr2, &t_sync);
printf(" t_sync = %e years \n log t_sync = %f \n", t_sync, log10(t_sync));

compute_tcirc(mm, qq, e2, rr, aa, mr2, &t_circ);
printf(" t_circ = %e years \n log t_circ = %f \n", t_circ, log10(t_circ));
return(0);
}
/***********************************************************************
*
***********************************************************************/
int compute_tsync(double mm, double qq, double e2, double rr, double aa,
                  double mr2, double *t_sync)
{
double tt;

tt = 5 * pow(2.,(double)(5./3.)) * sqrt(G * Msol / (Rsol*Rsol*Rsol)); 
tt *= sqrt(mm / (rr * rr * rr)) * qq * qq * pow((1. + qq),(double)(5./6.)); 
tt *= (e2 * pow((rr / aa),(double)8.5) / mr2);
*t_sync = 1./tt;

/* Conversion from seconds to years: */
*t_sync /= (3600. * 24. * 365.25);

return(0);
}
/**************************************************************************
* To compute t_circ  with Zahn (1975)'s formula
**************************************************************************/
int compute_tcirc(double mm, double qq, double e2, double rr, double aa,
                  double mr2, double *t_circ)
{
double tt;

tt = (21. / 2.) * sqrt(G * Msol / (Rsol*Rsol*Rsol)); 
tt *= sqrt(mm / (rr * rr * rr)) * qq * pow((1. + qq),(double)(11./6.)); 
tt *= e2 * pow((rr / aa),(double)10.5);
*t_circ = 1./tt;

/* Conversion from seconds to years: */
*t_circ /= (3600. * 24. * 365.25);
return(0);
}
/* -------------------------------------------------------------------------
 Ajustement d'une parabole y=a x^2 + b x + c
 a partir de 3 couples de points (x1,y1), (x2,y2) et (x3,y3)
(From Eric Aristidi)
 ------------------------------------------------------------------------- */
void parfit(double x1, double x2, double x3, double y1,
            double y2, double y3, double *a, double *b, double *c)

{
 double d,da,db;

  d=(x3*x3-x1*x1)*(x2-x1)-(x2*x2-x1*x1)*(x3-x1);
  da=(y3-y1)*(x2-x1)-(y2-y1)*(x3-x1);
  db=(x3*x3-x1*x1)*(y2-y1)-(x2*x2-x1*x1)*(y3-y1);
  *a=da/d ;
  *b=db/d;
  *c=y1-(*a)*x1*x1-(*b)*x1;
}

/* -------------------------------------------------------------------------
 Interpolation d'un tableau de points
 float *tab1,*tab2: tableaux de depart et d'arrivee
 dim1 est la dimension de depart, dim2 celle d'arrivee.
(From Eric Aristidi)
 ------------------------------------------------------------------------- */

void interpol(float *tab1, float *tab2, int dim1, int dim2)
{
 int i;
 float step;
 double x1,x2,x3,y1,y2,y3,a,b,c,x,y;

 step=((float) dim1-1)/(dim2-1);

/* Traitement des premiers points  */

 for (i=0;(int) (i*step)<(dim1-2);i++)
 {
  x=(double) i*step;
  x1=(double) ((int) x);
  x2=(double) ((int) x+1);
  x3=(double) ((int) x+2);
  y1=(double) *(tab1+(int) x1);
  y2=(double) *(tab1+(int) x2);
  y3=(double) *(tab1+(int) x3);
  parfit(x1,x2,x3,y1,y2,y3,&a,&b,&c);
  y=a*x*x+b*x+c;
  *(tab2+i)=y;
 }

/* Traitement des derniers points  */

 for (i=dim2-2;(int) (i*step)>=(dim1-2);i--)
 {
  x=(double) i*step;
  y=a*x*x+b*x+c;
  *(tab2+i)=y;
 }
 *(tab2+dim2-1)=*(tab2+dim1-1); // Dernier point

}
/*********************************************************
* Compute I/MR2 corresponding to mm Msol (stellar mass) 
*********************************************************/
int compute_imr2(double mm, double *mr2)
{
double imr2[7], imm[7], aa, bb, cc;
register int i;

imm[0] = 1.6; imm[1] = 2.; imm[2] = 3., imm[3] = 5.; imm[4] = 7.;
imm[5] = 10.; imm[6] = 15.;
imr2[0] = 7.57e-2; imr2[1] = 7.71e-2; imr2[2] = 7.95e-2; imr2[3] = 8.35e-2;
imr2[4] = 8.72e-2; imr2[5] = 9.13e-2; imr2[6] = 9.59e-2;

if(mm <= imm[0]) {
  *mr2 = imr2[0];
  printf("compute_imr2/Warning: mass smaller than minimum mass (= %f)!\n",
         imm[0]);
  return(-1);
  }
if(mm >= imm[6]) {
  *mr2 = imr2[6];
  printf("compute_imr2/Warning: mass larger than maximum mass (= %f)!\n",
         imm[6]);
  return(-1);
  }
/* Parabolic interpolation on 3 points: */
/* Special treatment of the end: */
if(mm >= imm[5])
  i = 4; 
else {
  for(i = 0; i <= 5; i++) if(mm < imm[i+1]) break; 
  }

#ifdef DEBUG
 printf(" OK: i=%d\n", i);
#endif

 parfit(imm[i], imm[i+1], imm[i+2], imr2[i], imr2[i+1], imr2[i+2], 
        &aa, &bb, &cc);
 *mr2 = aa * mm * mm + bb * mm + cc;

#ifdef DEBUG
 printf(" OK: mm=%f I/MR2=%e\n", mm, *mr2);
#endif

return(0);
}
/*********************************************************
* Compute E_2 corresponding to mm Msol (stellar mass) 
*********************************************************/
int compute_e2(double mm, double *e2)
{
double ie2[7], imm[7], aa, bb, cc;
register int i;

imm[0] = 1.6; imm[1] = 2.; imm[2] = 3., imm[3] = 5.; imm[4] = 7.;
imm[5] = 10.; imm[6] = 15.;
ie2[0] = 2.41e-9; 
ie2[1] = 1.45e-8; 
ie2[2] = 4.72e-8; 
ie2[3] = 1.53e-7; 
ie2[4] = 3.80e-7; 
ie2[5] = 1.02e-6; 
ie2[6] = 3.49e-6; 

if(mm <= imm[0]) {
  *e2 = ie2[0];
  printf("compute_imr2/Warning: mass smaller than minimum mass (= %f)!\n",
         imm[0]);
  return(-1);
  }
if(mm >= imm[6]) {
  *e2 = ie2[6];
  printf("compute_imr2/Warning: mass larger than maximum mass (= %f)!\n",
         imm[6]);
  return(-1);
  }
/* Parabolic interpolation on 3 points: */
/* Special treatment of the end: */
if(mm >= imm[5])
  i = 4; 
else {
  for(i = 0; i <= 5; i++) if(mm < imm[i+1]) break; 
  }

#ifdef DEBUG
 printf(" OK: i=%d\n", i);
#endif
 parfit(imm[i], imm[i+1], imm[i+2], ie2[i], ie2[i+1], ie2[i+2], 
        &aa, &bb, &cc);
 *e2 = aa * mm * mm + bb * mm + cc;
#ifdef DEBUG
 printf(" OK: mm=%f E_2=%e\n", mm, *e2);
#endif

return(0);
}
