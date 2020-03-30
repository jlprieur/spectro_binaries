/***************************************************************
* period_lomb.c
*
* To compute the period of a periodic signal
* in the case of non-uniform sampling
*
* JLP
* Version 27/05/2002
***************************************************************/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>

#define TWOPID 6.2831853071795865

#define DEBUG
#define STOP_ERROR { printf("Fatal error reading input file\n"); fclose(fp); exit(-1);}

/* Prototypes of functions defined here:
*/
void period_phase(float *x, float *y, int n, float period);

/* Prototypes of functions defined in "period_lomb1.c":
*/
void jlp_period_lomb0(float *x, float *y, int n, float ofac, float fmax,
                      float *x_ave, float *y_ave, float *y_var,
                      float *fstep, int *nout);
void jlp_period_lomb1(float *x, float *y, int n, float ofac,
                      float *px, float *py, int np, int nout,
                      float x_ave, float y_ave, float y_var, float fstep,
                      int *jmax, float *prob);

/* Prototypes of functions defined in "sb_read.c":
*/
int period_read_data(FILE *fp, float *x, float *y, int ic, int *n);

int main(int argc, char *argv[])
{
float *x, *y, ofac, *px, *py, prob;
float fstep, fmax, x_ave, y_ave, y_var, pmax;
int n, np, nout, jmax, ic; 
register int i, j;
char filename[60], buffer[80];
FILE *fp;

printf("period_lomb  Version 01/10/2004 \n");
if((argc != 1) && (argc != 4)) {
  printf("Syntax: period_lomb filename fmax curve_nber \n");
  printf("curve_nber: \n 0=primary (3rd column) \n");
  printf(" 1=primary (2nd column) \n 2=secondary (4th column with weights)\n");
  printf(" 3=prim. & sec. \n 4=3rd body \n");
  exit(-1);
  }
else if (argc == 4) {
  sscanf(argv[1],"%s",filename);
  sscanf(argv[2],"%f",&fmax);
  sscanf(argv[3],"%d",&ic);
  }
else {
/* Input file: */
  printf(" Filename: ");
  scanf("%s", filename);
/* maximum frequency: */
  printf(" Maximum frequency: ");
  scanf("%f", &fmax);
/* Curve number: */
  printf("curve nber: 1=primary 2=secondary 3=prim. & sec. 4=3rd body \n");
  printf(" Curve number : ");
  scanf("%d", &ic);
  }

#ifdef DEBUG
printf(" OK: fmax=%f ic=%d\n", fmax, ic);
#endif


/* Real data */ 
if((fp = fopen(filename,"r")) == NULL) {
  printf("Fatal error opening %s\n", filename);
  exit(-1);
 }

/* First lines is useless here: */
if(fgets(buffer,80,fp) == NULL) STOP_ERROR
if(fgets(buffer,80,fp) == NULL) STOP_ERROR
sscanf(buffer,"%d",&n);
printf("OK: nber of data points=%d\n", n);

x = (float *)malloc((n + 1) * sizeof(float));
y = (float *)malloc((n + 1) * sizeof(float));

period_read_data(fp, x, y, ic, &n);

/* Save input data: */
strcpy(filename,"period_in.dat");
if((fp = fopen(filename,"w")) == NULL) {
  printf("Fatal error opening %s\n", filename);
  exit(-1);
 }
else {
  for(i = 1; i <= n; i++) {
   fprintf(fp,"%f %f\n", x[i], y[i]);
   }
  fclose(fp);
  }

/* JLP2005: ofac = 4 => too many points */
ofac = 2.;
jlp_period_lomb0(x, y, n, ofac, fmax, &x_ave, &y_ave, &y_var, &fstep, &nout);

printf(" Number of points: nout = %d\n", nout);

np = nout;
px = (float *)malloc((np + 1) * sizeof(float));
py = (float *)malloc((np + 1) * sizeof(float));

/* Mean Nyquist frequency: nf = (1 / <Delta x>) = n / xdif
* fstep = 1.0 / (xdif * ofac);
*  hifac : fraction of (average) Nyquist frequency to be used
* fmax = fstep * nout = 0.5 * ofac * hifac * n / (xdif * ofac)
       = 0.5 * hifac * nf
*/

jlp_period_lomb1(x, y, n, ofac, px, py, np, nout, x_ave, y_ave, y_var,
                 fstep, &jmax, &prob);

printf("Output from period: max at py[%d]=%g prob=%e\n", jmax, py[jmax], prob);

/* Results: */
strcpy(filename,"period.dat");
if((fp = fopen(filename,"w")) == NULL) {
  printf("Fatal error opening %s\n", filename);
 }
else {
  for(i = 1; i <= nout; i++) {
   fprintf(fp,"%f %f\n", px[i], py[i]);
   }
  fclose(fp);
 }

/* Output phases of the 8 higher periods */
pmax = 1./px[jmax];
period_phase(x, y, n, pmax);
printf(" px[%d]=%f (max at %g) period=%f\n", jmax, px[jmax], py[jmax], pmax);

for(j = 0; j < 7; j++) {
   py[jmax] = 0.;
   pmax = 0.; 
   for(i = 1; i <= nout; i++) {
     if(pmax < py[i]) {
       pmax = py[i]; jmax = i;
       } 
     } 
   pmax = 1./px[jmax];
   printf(" px[%d]=%f (max at %g) period=%f\n", jmax, px[jmax], py[jmax], pmax);
   period_phase(x, y, n, pmax);
   }

free(x); 
free(y);
free(px); 
free(py);
return 0;
}
/*****************************************************************
*
* To compute the phase if the period is known.
*
*****************************************************************/
void period_phase(float *x, float *y, int n, float period)
{
float phi;
char filename[60];
FILE *fp;
register int i;

/* Output data */
sprintf(filename,"phase_%.2f.dat",period);
if((fp = fopen(filename,"w")) == NULL) {
  printf("Fatal error opening %s\n", filename);
  exit(-1);
 }

for(i = 1; i <= n; i++) {
  phi = x[i] - period * (float)((int)(x[i] / period));
  fprintf(fp,"%f %f\n", phi, y[i]);
  }

fclose(fp);
}
