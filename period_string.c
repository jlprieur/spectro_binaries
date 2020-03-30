/***************************************************************
* period_string
*
* To compute the period of a periodic signal
* in the case of non-uniform sampling
* by minimising the length of the string in the phase diagram
* Cf. M.M. Dworetsky MNRAS 1983, 203, 917-924
* " A period-finding method for sparse randomly spaced observations
* or "How long is a piece of string ?" " 
*
* JLP
* Version 01/10/2004
***************************************************************/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>
#include <jlp_ftoc.h>

#define TWOPID 6.2831853071795865

#define DEBUG
#define STOP_ERROR { printf("Fatal error reading input file\n"); fclose(fp); exit(-1);}

/* Prototypes of functions defined here:
*/
static void output_phase(float *x, float *y, int n, float period);
static float phase_string_length(float *c_x, float *c_y, int n);

void jlp_period_string0(float *x, float *y, int n, float pmin, float pstep,
                      float *px, float *py, int np);

/* Prototypes of functions defined in "sb_read.c":
*/
int period_read_data(FILE *fp, float *x, float *y, int ic, int *n);

int main(int argc, char *argv[])
{
float *x, *y, *px, *py;
float pstep, pmin, pmax, phase_string_min;
int n, np, ic, jmin; 
register int i, j;
char filename[60], buffer[80];
FILE *fp;

printf("period_string  Version 01/10/2004 \n");
if(argc == 7){
  for(i = 0; i < 7 && argv[i][0]; i++);
   argc = i;
  }
if((argc != 1) && (argc != 4)) {
  printf("argc=%d\n", argc);
  printf("Syntax: period_string filename pmin,pmax,pstep curve_nber \n");
  printf("curve_nber: \n 0=primary (3rd column) \n");
  printf(" 1=primary (2nd column) \n 2=secondary (4th column with weights)\n");
  printf(" 3=prim. & sec. \n 4=3rd body \n");
  exit(-1);
  }
else if (argc == 4) {
  sscanf(argv[1],"%s",filename);
  sscanf(argv[2],"%f,%f,%f",&pmin,&pmax,&pstep);
  sscanf(argv[3],"%d",&ic);
  }
else {
/* Input file: */
  printf(" Filename: ");
  scanf("%s", filename);
/* maximum frequency: */
  printf(" Minimum, maximum periods and step: (in days)");
  scanf("%f,%f,%f", &pmin, &pmax, &pstep);
/* Curve number: */
  printf("curve nber: 1=primary 2=secondary 3=prim. & sec. 4=3rd body \n");
  printf(" Curve number : ");
  scanf("%d", &ic);
  }

#ifdef DEBUG
printf(" OK: pmin=%f pmax=%f ic=%d\n", pmin, pmax, ic);
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
/* Input data was read in the [1, n] range
* Transfer to [0, n-1] range: */
for(i = 0; i < n; i++) {x[i] = x[i+1]; y[i] = y[i+1];}

/* Save input data: */
strcpy(filename,"period_in.dat");
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
px = (float *)malloc(np * sizeof(float));
py = (float *)malloc(np * sizeof(float));

jlp_period_string0(x, y, n, pmin, pstep, px, py, np);

/* Results: */
strcpy(filename,"period_string.dat");
if((fp = fopen(filename,"w")) == NULL) {
  printf("Fatal error opening %s\n", filename);
 }
else {
  for(i = 0; i < np; i++) {
   fprintf(fp,"%f %f\n", px[i], py[i]);
   }
  fclose(fp);
 }

/* Output phases of the 8 lower periods */

jmin = 0;
for(j = 0; j < 7; j++) {
   if(j) py[jmin] = 1.e+6;
   phase_string_min = 1.e+6; 
   for(i = 0; i < np; i++) {
     if(phase_string_min > py[i]) {
       phase_string_min = py[i]; jmin = i;
       } 
     } 
   printf("Iter #%d/ Chord is minimum at %g for period=%f\n", 
           j, py[jmin], px[jmin]);
   output_phase(x, y, n, px[jmin]);
   }

printf("Output in period_string.dat \n");
free(x); 
free(y);
free(px); 
free(py);
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
*
* To compute and output the phase if the period is known.
*
* INPUT:
* x, y: input data (time and velocity)
* n: number of data points
* OUTPUT:
*
*****************************************************************/
void jlp_period_string0(float *x, float *y, int n, float pmin, float pstep,
                      float *px, float *py, int np)
{
float *c_x, *c_y, period;
int *index;
register int i, j;

if((c_x = (float *)malloc(n * sizeof(float))) == NULL) {
 printf("jlp_period_string0/Fatal error allocating memory (n=%d)\n", n);
 exit(-1);
 }
if((c_y = (float *)malloc(n * sizeof(float))) == NULL) {
 printf("jlp_period_string0/Fatal error allocating memory (n=%d)\n", n);
 exit(-1);
 }

if((index = (int *)malloc(n * sizeof(int))) == NULL) {
 printf("jlp_period_string0/Fatal error allocating memory (n=%d)\n", n);
 exit(-1);
 }

/* Compute period array: */
for(i = 0; i < np; i++) px[i] = pmin + i * pstep;

/* Compute phase string array: */
for(i = 0; i < np; i++) {

/* Compute cycle: */
  period = px[i]; 
  for(j = 0; j < n; j++) 
     c_x[j] = (x[j] - period * (float)((int)(x[j] / period))) / period;

/* Sort cycle members: */
  JLP_QSORT_INDX(c_x, index, &n);
  for(j = 0; j < n; j++) c_y[j] = y[index[j]];

/* Compute phase string length: */
  py[i] = phase_string_length(c_x, c_y, n);

/* EOF loop on i */
  }

free(c_x);
free(c_y);
free(index);
return;
}
/**********************************************************************
* Compute phase string length, i.e., length of cycle of a given period
*
* INPUT:
* c_x, c_y: arrays with cycle data (cycle, velocity)
* n: number of points
**********************************************************************/
static float phase_string_length(float *c_x, float *c_y, int n)
{
double sum;
register int i;

sum = 0;
for(i = 1; i < n; i++) {
 sum += sqrt(SQUARE(c_x[i] - c_x[i-1]) + SQUARE(c_y[i] - c_y[i-1]));
 }

return(sum);
}
