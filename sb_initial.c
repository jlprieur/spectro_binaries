/******************************************************
* sb_initial.c
* To compute Fourier coefficients with least-square minimization
*
* JLP
* Version 05-06-2002
******************************************************/
#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <jlp_ftoc.h>

/*
#define DEBUG 1
*/
int JLP_CGRAD(double *aa, double *psi, double *phi, INT4 *nx1, INT4 *ny1,
              INT4 *ifail);

/* Prototypes of functions included here: */
int sb_initial_process(char *infile, int n_order, int ic, int talk);
static void input_data_simu(double **psi, double **aa, double **x0,
                            double period, int *ny, int nx, int n_order);
static void input_data(double **psi, double **aa, double **x0,
                       double period, int *ny, int nx, int n_order, 
                       char *infile, int ic);
static int sb_initial_param(double *phi, int n_order, double period,
                            double t1, double t2, char *infile, int ic,
                            int talk);
static int compute_fterm(double x0, double *phi, int n_order, double *ww);
static int write_to_file(char *infile, char *first_line, int ic);

/* Prototypes of functions defined in sb_read.c:
*/
int period_read_data(FILE *fp, float *x, float *y, int ic, int *n);

#define TWOPI (2.*PI)

/* BOF MAIN_PROGRAM */

#ifdef MAIN_PROGRAM
int main(int argc, char *argv[])
{
int n_order, ic, talk;
register int i;
char infile[60];

if(argc == 7) {
for(i = 0; i <= 6 && argv[i][0]; i++) {}
argc = i;
}

if(argc != 4) {
  printf("Error: argc = %d (!= 4) \n",argc);
  printf("Error: syntax is \n sb_initial input_file n_order ic \n");
  printf(" n_order=order of the developement of Fourier Series ( >= 1)\n"); 
  printf(" ic=curve nber: 1=primary, 2=secondary, 3=prim. & sec., 4=3rd body)\n"); 
  printf(" if(n_order == 1) eccentricity is set to 0 \n");
  exit(-1);
  }
else {
  sscanf(argv[1], "%s", infile);
  sscanf(argv[2], "%d", &n_order);
  sscanf(argv[3], "%d", &ic);
  printf(" OK: infile=%s n_order=%d ic=%d\n", infile, n_order, ic);
  }

if(n_order < 1) {printf("Fatal error: n_order=%d < 1 \n", n_order);
  exit(-1);
  }
 
/* The job is done by "sb_initial_process": */
talk = 1;
sb_initial_process(infile, n_order, ic, talk);

return(0);
}
#endif
/* EOF MAIN_PROGRAM */
/*************************************************************
* Prepare aa matrix to fit Fourier coefficients by solving psi = aa phi
* phi: Fourier coefficients: a0, a1, b1, a2, b2
* psi = V(t)
* V(t) = a0/2 + a1 cos(2 PI t / P) + b1 sin(2 PI t / P)
*             + a2 cos(2 PI 2t / P) + b2 sin(2 PI 2t / P)
*
*************************************************************/
static void input_data(double **psi, double **aa, double **x0, 
                       double period, int *ny, int nx, int n_order, 
                       char *infile, int ic)
{
char buffer[80];
float *xx, *yy;
register int i, j;
FILE *fp;

/* Real data */
if((fp = fopen(infile,"r")) == NULL) {
  printf("Fatal error opening %s\n", infile);
  exit(-1);
 }
/* First lines is useless here: */
fgets(buffer,80,fp);
fgets(buffer,80,fp);
sscanf(buffer,"%d",ny);
#ifdef DEBUG
printf("Nber of data points should be: %d\n", *ny);
#endif

/* Indices start at 1 for xx and yy: */
xx = (float *)malloc(((*ny) + 1) * sizeof(float));
yy = (float *)malloc(((*ny) + 1) * sizeof(float));

/* Warning: fp is closed there! */
period_read_data(fp, xx, yy, ic, ny);

#ifdef DEBUG
printf("Nber of data points is actually: %d (nx=%d)\n", *ny, nx);
#endif

(*psi) = (double *)malloc((*ny) * sizeof(double));
(*x0) = (double *)malloc((*ny) * sizeof(double));
(*aa) = (double *)malloc((*ny) * nx * sizeof(double));

for(i = 0; i < *ny; i++) {
    (*psi)[i] = yy[i+1]; 
    (*x0)[i] = (xx[i+1] - period * (int)(xx[i+1] / period)) / period;
    (*aa)[0 + i * nx] = 1.0; 
    for(j = 1; j <= n_order; j++) {
/*
       (*aa)[(2*j) - 1 + i * nx] = cos(TWOPI * (double)j * xx[i+1] / period); 
       (*aa)[2*j + i * nx] = sin(TWOPI * (double)j * xx[i+1] / period); 
*/
       (*aa)[(2*j) - 1 + i * nx] = cos(TWOPI * (double)j * (*x0)[i]); 
       (*aa)[2*j + i * nx] = sin(TWOPI * (double)j * (*x0)[i]); 
       }
  }

free(xx);
free(yy);
}
/*************************************************************
* We examine the problem: psi = aa phi
* We simulate a problem whose solution is: phi = (1, 2, 2, 4, 5)
* V(t) = 1 + 2 cos(2 pi t / T) + 2 sin(2 pi t / T)
*         + 4 cos(2 pi 2 t / T) + 5 sin(2 pi 2 t / T)
*************************************************************/
static void input_data_simu(double **psi, double **aa, double **x0,
                            double period, int *ny, int nx, int n_order)
{
double xx;
register int i, j;

*psi = (double *)malloc((*ny) * sizeof(double));
*aa = (double *)malloc((*ny) * nx * sizeof(double));

for(i = 0; i < *ny; i++) {
    xx = (double)i; 
if(n_order >= 2) {
    (*psi)[i] = 1.0 + 2. * cos(TWOPI * xx / period) 
              + 2. * sin(TWOPI * xx / period)
              - 0.04 * cos(TWOPI * 2. * xx / period) 
              + 0.05 * sin(TWOPI * 2. * xx / period); 
  } else {
    (*psi)[i] = 1.0 + 2. * cos(TWOPI * xx / period) 
              + 2. * sin(TWOPI * xx / period);
  }
    (*aa)[0 + i * nx] = 1.0; 
    for(j = 1; j <= n_order; j++) {
       (*aa)[(2*j) - 1 + i * nx] = cos(TWOPI * (double)j * xx / period); 
       (*aa)[2*j + i * nx] = sin(TWOPI * (double)j * xx / period); 
       }
  }
printf(" Solving   aa phi = psi \n");
if(n_order >= 2) {
  printf(" phi solution should be: (1,2,2,-0.04,0.05) \n");
  } else {
  printf(" phi solution should be: (1,2,2) \n");
  }
}
/*********************************************************************
*
*********************************************************************/
static int sb_initial_param(double *phi, int n_order, double period,
                            double t1, double t2, char *infile, int ic,
                            int talk)
{
double beta1, beta2, r2, r1, ee, ee0, v0;
double gamma, alpha, mu_t, t0, tt;
double kk, k1, k2, tan_omega, omega;
int ee_is_zero;
register int i;
char first_line[80];

if(n_order <= 1) ee_is_zero = 1;
else ee_is_zero = 0;

if(!ee_is_zero && (phi[2] == 0 || phi[4] == 0))
    {
    printf("sb_initial_param/Error: sinus coefficient null for order 1 or 2\n");
    ee_is_zero = 1;
    }

if(!ee_is_zero) {
  beta1 = atan(phi[1] / phi[2]);
  r1 = phi[1] / sin(beta1);
  beta2 = atan(phi[3] / phi[4]);
  r2 = phi[3] / sin(beta2);
  if(r1 * r2 < 0.) {
    r2 *= -1;
    beta2 += PI;
    }
  gamma = 4. * beta1 - 2. * beta2;

/* Excentricity with successive approximations: */
/* Both methods leads to the same result but the second is more tolerant
* for large eccentricity */
  i = 0; ee = r2 / r1;
  if(talk) printf(" r1=%.2f r2=%.2f ee=%.2f\n", r1, r2, ee);
  if(ee < 1.0) {
#ifdef TTT
#define TTT
  ee0 = 0;
  while(ABS(ee - ee0) > 1.e-5 && i < 100) {
    ee0 = ee;
    ee = r2 / r1 + ee * ee * ee * ( (0.25 + cos(gamma) / 24.)
         + ee * ee * cos(gamma) / 96.);
/*
    printf("iter=%d/ ee = %.5f (Delta=%.4g)\n", i, ee, ABS(ee - ee0)); 
*/
    i++;
  }
#else
  ee0 = 1.; 
  while(ABS(ee0) > 1.e-5 && i < 100) {
    ee0 = - ee * ee * ee * ( (0.25 + cos(gamma) / 24.)
         - ee * ee * cos(gamma) / 96.) + ee - r2/r1;
    ee0 /= ee * ee * (5. * ee * ee * cos(gamma) / 96. 
           + 3. * (0.25 + cos(gamma) / 24.)) - 1.;
    ee += ee0;
/*
    printf("iter=%d/ ee = %.5f (Delta=%.4g)\n", i, ee, ee0); 
*/
    i++;
    }
#endif
/* End of ee < 1 */
  }
/* End of !ee_is_zero */
}

if(ee < 0. || ee >= 1.) {
 ee_is_zero = 1;
 if(talk) printf(" Error: e=%.2f !\n", ee);
 }

if(ee_is_zero) {
   ee = 0.;
   t0 = t2;
   omega = 0.;
   kk = sqrt(phi[1] * phi[1] + phi[2] * phi[2]);
 } else {
/* T0: */ 
  alpha = ee * ee * ee * (1./24. + ee * ee / 96.);
  alpha *= sin(gamma) * r1 / r2;
  mu_t = alpha + beta1 - beta2;

  tt = period * mu_t / TWOPI;

  t0 = tt + t1 + period * (int)((t2 - t1) / period);
  while( t0 < t2) t0 += period;

/* omega: */
/* k1 and k2 (cf equations 8) */
  k1 = r1 * sin(mu_t + beta1) / 
      ( 1. + ee * ee * (-1./ 8. + ee * ee / 192.)); 
  k2 = -r1 * cos(mu_t + beta1) / 
      ( 1. + ee * ee * (-3./ 8. + ee * ee * 5. / 192.)); 
  tan_omega = sqrt(1. - ee * ee) * k2 / k1; 
  if(k1 < 0) 
    omega = PI + atan(tan_omega); 
  else {
    omega = atan(tan_omega);
    if(omega < 0) omega += TWOPI;
    }

  kk = k2 / (sin(omega) * sqrt(1. - ee * ee));
}

v0 = phi[0];
if(talk) {
         printf(" Preliminary elements: \n");
         printf(" Period = %.6f days \n", period);
         printf(" V0 = %.2f km/s \n", v0);
         printf(" T0 = %.4f (JD) \n", t0);
         printf(" e = %.3f \n", ee);
         printf(" K = %.3f km/s \n", kk);
         printf(" omega = %.3f rad  (or %.3f deg)\n", omega, omega * 180. / PI);

         printf(" Now writing these elements in the first line of VR_BS1.DAT\n");
}
sprintf(first_line,"%10.4f%10.3f%10.5f%10.5f%9.2f %6.2f    ", 
        period, t0, omega, ee, kk, v0);
write_to_file(infile, first_line, ic); 

return(0);
}
/**************************************************************************
*
**************************************************************************/
static int compute_fterm(double x0, double *phi, int n_order, double *ww)
{
register int i;
*ww = phi[0]; 
for(i = 1; i <= n_order; i++) {
   *ww += phi[2*i - 1] * cos(TWOPI * (double)i * x0)
          + phi[2*i] * sin(TWOPI * (double)i * x0);
   }
return(0);
}
/******************************************************************
* READ input data from "infile" and write VR_BS1.DAT file
*
******************************************************************/
static int write_to_file(char *infile, char *first_line, int ic) 
{
float *xx, *yy;
char buffer[80];
int ny;
FILE *fp, *fp_out;
register int i;

if((fp = fopen(infile,"r")) == NULL) {
  printf("write_to_file/Fatal error opening %s \n", infile);
  exit(-1);
  }
/* First two lines: */
fgets(buffer,80,fp);
fgets(buffer,80,fp);
sscanf(buffer,"%d %*d", &ny);
#ifdef DEBUG
printf(" ny = %d\n", ny);
#endif

/* Indices start at 1 for xx and yy: */
xx = (float *)malloc((ny + 1) * sizeof(float));
yy = (float *)malloc((ny + 1) * sizeof(float));

/* Warning: fp is closed there! */
period_read_data(fp, xx, yy, ic, &ny);

if((fp_out = fopen("VR_BS1.DAT","w")) == NULL) {
  printf("write_to_file/Fatal error opening VR_BS1.DAT \n");
  }
else {
/* First line with initial values: */
   fprintf(fp_out, "%s\n", first_line);
   fprintf(fp_out, "%d  5\n", ny);
   for(i = 1; i <= ny; i++) 
     fprintf(fp_out, "%f %f 1.0\n", xx[i], yy[i]); 
  fclose(fp_out);
  }

free(xx);
free(yy);
return(0);
}
/**************************************************************************
* Main routine which computes the initial parameters 
* (directly called by period_resid)
*
**************************************************************************/
int sb_initial_process(char *infile, int n_order, int ic, int talk)
{
INT4 nx, ny, ifail;
int status;
register int i;
double *aa, *phi, *psi, *x0, period, t1, t2;
#ifdef DEBUG
double ww;
#endif
float cst;
char buffer[80];
FILE *fp;

nx = 1 + 2 * n_order;

#ifdef TEST
#define TEST
  ny = 30;
  period = (double)ny;
  input_data_simu(&psi, &aa, &x0, period, &ny, nx, n_order);
#else

  if((fp = fopen(infile,"r")) == NULL) {
    printf("Fatal error opening %s\n", infile);
    exit(-1);
    }
/* Format: same as VR_BS1.DAT or VR_BS2.DAT
**********************************************************/
/* period is on the first line: */
  fgets(buffer,80,fp);
  sscanf(buffer,"%lf", &period);
/* ny is on the second line: */
  fgets(buffer,80,fp);
  sscanf(buffer,"%d", &ny);
/* First date is taken as initial value for t2 */
  fgets(buffer,80,fp);
  sscanf(buffer,"%lf", &t2);
  fclose(fp);
  if(talk) printf(" Period: %f  ny = %d\n", period, ny);
  t1 = 0.;
/*
  printf(" Epoch of origin: T1=%.3f\n", t1);
*/
  if(talk) printf(" Epoch T0 that will be as close as possible to T2=%.3f\n", t2);
  input_data(&psi, &aa, &x0, period, &ny, nx, n_order, infile, ic);
#endif

/*
  for (i = 0; i < ny; i++) printf(" psi[%d] = %12.5e \n", i, psi[i]);
*/

/* First guess for phi */
phi = (double *)malloc(nx * sizeof(double));

/*
printf(" Enter initial guess for phi (one value only) := ");
scanf("%f",&cst);
*/

cst = 0.;
for (i = 0; i < nx; i++) phi[i] = cst;

/* DEBUG:
printf("Initial guess: phi=(");
for (i = 0; i < nx; i++) printf("%.2f ",phi[i]);
printf(")\n");
*/

/* Fit Fourier coefficients by solving psi = aa phi
* aa: matrix filled with sin/cos(2 PI t /P)
* phi: Fourier coefficients: a0, a1, b1, a2, b2
* psi = V(t)
* V(t) = a0/2 + a1 cos(2 PI t / P) + b1 sin(2 PI t / P)
*             + a2 cos(2 PI 2t / P) + b2 sin(2 PI 2t / P)
*/
status = JLP_CGRAD(aa,psi,phi,&nx,&ny,&ifail);

#ifdef DEBUG
printf("Solution: phi \n %.3f \n", phi[0]);
for (i = 1; i < nx; i+=2) printf("%.3f %.3f \n",phi[i],phi[i+1]);
#endif

sb_initial_param(phi, n_order, period, t1, t2, infile, ic, talk);

#ifdef DEBUG
fp = fopen("sb_initial.dat","w");
for(i = 0; i < ny; i++) {
  compute_fterm(x0[i], phi, n_order, &ww);
  fprintf(fp,"%f %f %f\n", x0[i], psi[i], ww);
  }
fclose(fp);
#endif

free(phi);
free(x0);
free(psi);
free(aa);

return(0);
}
