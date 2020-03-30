/***************************************************************
* period_resid
*
* To compute the period of a periodic signal
* in the case of non-uniform sampling
* by minimizing the residuals 
*
* JLP
* Version 30/05/2005
***************************************************************/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>
#include <jlp_ftoc.h>

#define TWOPID 6.2831853071795865

/*
#define DEBUG
*/
#define STOP_ERROR { printf("Fatal error reading input file\n"); fclose(fp); exit(-1);}

/* Prototypes of functions defined here:
*/
static void output_phase(float *x, float *y, int n, float period);
static int compute_VR_residuals(float period_in1, float *resid1, 
                                float *period_out1, char *infile, int ic,
                                int n_order);

/* Prototypes of functions defined in "sb_read.c":
*/
int period_read_data(FILE *fp, float *x, float *y, int ic, int *n);
/* in "sb_initial.c" */
int sb_initial_process(char *infile, int n_order, int ic, int talk);


int main(int argc, char *argv[])
{
float *x, *y, *period_in, *period_out, *resi;
float pstep, pmin, pmax, resid_min, ww1, ww2;
int n, np, ic, jmin, n_order; 
register int i, j;
char infile[60], filename[60], buffer[80];
FILE *fp;

printf("period_resid  Version 01/10/2004 \n");
if(argc == 7){
  for(i = 0; i < 7 && argv[i][0]; i++);
   argc = i;
  }
if((argc != 1) && (argc != 4)) {
  printf("argc=%d\n", argc);
  printf("Syntax: period_resid infile pmin,pmax,pstep n_order \n");
  printf("n_order: polynomial order for decomposition of Fourier coeff.\n");
  printf("Should be larger than 1 (if equal to 1, eccentricity is set to zero)\n");
  exit(-1);
  }
else if (argc == 4) {
  sscanf(argv[1],"%s",infile);
  sscanf(argv[2],"%f,%f,%f",&pmin,&pmax,&pstep);
  sscanf(argv[3],"%d",&n_order);
  }
else {
/* Input file: */
  printf(" Input data file: ");
  scanf("%s", infile);
/* maximum frequency: */
  printf(" Minimum, maximum periods and step: (in days)");
  scanf("%f,%f,%f", &pmin, &pmax, &pstep);
  printf("n_order (polynomial order for decomposition of Fourier coeff.):\n");
  scanf("%d", &n_order);
  }

#ifdef DEBUG
printf(" OK: pmin=%f pmax=%f n_order=%d\n", pmin, pmax, n_order);
#endif


/* Real data */ 
if((fp = fopen(infile,"r")) == NULL) {
  printf("Fatal error opening %s\n", infile);
  exit(-1);
 }

/* First lines is useless here: */
if(fgets(buffer,80,fp) == NULL) {printf("infile=%s\n",infile); STOP_ERROR}
if(fgets(buffer,80,fp) == NULL) {printf("infile=%s\n",infile); STOP_ERROR}
sscanf(buffer,"%d",&n);
printf("OK: nber of data points=%d\n", n);

x = (float *)malloc((n + 1) * sizeof(float));
y = (float *)malloc((n + 1) * sizeof(float));

/* Work with first component (ic=1): */
ic = 1;
period_read_data(fp, x, y, ic, &n);
/* Input data was read in the [1, n] range
* Transfer to [0, n-1] range: */
for(i = 0; i < n; i++) {x[i] = x[i+1]; y[i] = y[i+1];}

/* Save input data: */
strcpy(filename,"period_input.dat");
if((fp = fopen(filename,"w")) == NULL) {
  printf("Fatal error opening %s\n", filename);
  exit(-1);
 }
else {
  for(i = 0; i < n; i++) {
   fprintf(fp,"%f %f\n", x[i], y[i]);
   }
  fclose(fp);
  }

/* Sampling of pstep days: */
np = 1 + (pmax - pmin) / pstep; 
period_in = (float *)malloc(np * sizeof(float));
period_out = (float *)malloc(np * sizeof(float));
resi = (float *)malloc(np * sizeof(float));

/* Compute period array: */
for(i = 0; i < np; i++) period_in[i] = pmin + i * pstep;

/* Initialisation: */
for(i = 0; i < np; i++) period_out[i] = 0.;
for(i = 0; i < np; i++) resi[i] = 10000.;

/* Compute residuals array: */
for(i = 0; i < np; i++) {
  compute_VR_residuals(period_in[i], &resi[i], &period_out[i], infile, ic,
                                             n_order); 
/* JLP2006, to avoid nonsensical values: */
  if(ABS(period_out[i] - period_in[i]) > pstep) period_out[i] = period_in[i];
  }

/* Results: */
strcpy(filename,"period_resid.dat");
if((fp = fopen(filename,"w")) == NULL) {
  printf("Fatal error opening %s\n", filename);
 }
else {
  for(i = 0; i < np; i++) {
   fprintf(fp,"%f %f %f\n", period_out[i], resi[i], period_in[i]);
   }
  fclose(fp);
 }

/* Output phases of the 8 lower periods */

jmin = 0;
for(j = 0; j < MINI(8,np); j++) {
   if(j) resi[jmin] = 1.e+6;
   resid_min = 1.e+6; 
   for(i = 0; i < np; i++) {
     if(resid_min > resi[i]) {
       resid_min = resi[i]; jmin = i;
       } 
     } 
   printf("Iter #%d/ Residuals are minimum at %g for period=%f\n", 
           j, resi[jmin], period_out[jmin]);
   output_phase(x, y, n, period_out[jmin]);
/* Generate VR_BS1.DAT for the minimum: */
   if(j == 0) {
               compute_VR_residuals(period_out[jmin], &ww1,
                                             &ww2, infile, ic,
                                             n_order); 
              }
   }

printf("Output in period_input.dat, period_resid.dat \n");
printf("and VR_BS1.DAT for the minimum \n");
free(x); 
free(y);
free(period_in); 
free(period_out); 
free(resi);
return(0);
}
/*****************************************************************
*
* To compute and output the phase if the period is known.
*
*****************************************************************/
static void output_phase(float *x, float *y, int n, float period)
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

return;
}
/*****************************************************************
* To compute the residuals between the fitted curve and the RV data 
* (in VR_BS1.DAT)
*
*****************************************************************/
static int compute_VR_residuals(float period_in1, float *resid1, 
                                float *period_out1, char *infile, int ic,
                                int n_order)
{
char buffer[80], outfile[40], object[40];
int talk;
FILE *fp_in, *fp_out;

/* Creating VR_BS1.TMP: */
#ifdef DEBUG
    printf("Opening %s\n", infile);
#endif

  if((fp_in = fopen(infile,"r")) == NULL) {
    printf("Fatal error opening %s\n", infile);
    exit(-1);
    }
  strcpy(outfile,"VR_BS1.TMP");
  if((fp_out = fopen(outfile,"w")) == NULL) {
    printf("Fatal error opening %s\n", infile);
    exit(-1);
    }
/* period is on the first line: */
/* Read dummy first line: */
  fgets(buffer,80,fp_in);

/* Write period only: */
  fprintf(fp_out,"%10.6f \n", period_in1);

/* Then transfer all lines: */
  while(fgets(buffer,80,fp_in) != NULL) fputs(buffer,fp_out);

  fclose(fp_in);
  fclose(fp_out);
/**********************************************************/

/* First step: compute the initial values for the fit: */
talk = 0;
sb_initial_process(outfile, n_order, ic, talk);

/* Second step: fit a model to the data and compute the mean residual: */
strcpy(object,"mode-automatique");
*period_out1 = -1;
bs1_process__(object, period_out1, resid1);

return(0);
}
