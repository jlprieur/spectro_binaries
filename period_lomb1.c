/***************************************************************
* period_lomb1.c
*
* To compute the period of a periodic signal
* in the case of a non-uniform sampling
* Called by main program contained in "period_lomb.c"
*
* Lomb-Scargle method
* From Numerical recipees p 582
* with a simpler version with JLP modifications
*
* WARNING: arrays indices start at 1!
*
* Version 01/10/2004
***************************************************************/
#include <stdio.h>
#include <stdlib.h> /* RAND_MAX is defined here! */
#include <math.h>
#include <malloc.h>

#define TWOPID 6.2831853071795865
#define MOD(a,b) while(a >= b) a -= b;
/* Square: */
#define SQR(a) ( (a) == 0.0 ? 0.0 : (a)*(a))

#define DEBUG

/* Prototypes of functions defined here:
*/
void jlp_period_lomb0(float *x, float *y, int n, float ofac, float fmax,
                      float *x_ave, float *y_ave, float *y_var, 
                      float *fstep, int *nout);
void jlp_period_lomb1(float *x, float *y, int n, float ofac,
                      float *px, float *py, int np, int nout, 
                      float x_ave, float y_ave, float y_var, float fstep,
                      int *jmax, float *prob);
void period(float *x, float *y, int n, float ofac, float hifac,
            float *px, float *py, int np, int *nout, int *jmax, float *prob);
static void fasper(float *x, float *y, unsigned long n, float ofac, 
                   float hifac, float *wk1, float *wk2, unsigned long nwk, 
                   unsigned long *nout, unsigned long *jmax, float *prob);
static void avevar(float *data, unsigned long n, float *ave, float *var);
/*
static void realft(float *data, unsigned long n, int isign);
static void spread(float y, float *yy, unsigned long n, float x, int m);
*/

/*
#define TEST_MAIN
*/
#ifdef TEST_MAIN
/*************************************************************
* Test version
*************************************************************/
int main(int argc, char *argv[])
{
float *x, *y, ofac, hifac, *px, *py, prob, xmax, fmax;
float x_ave, y_ave, y_var, fstep;
int n, np, nout, jmax; 
register int i;
char filename[60];
FILE *fp;

printf("period_lomb JLP/Test version 24/05/2002 \n");
xmax = 1.;
n = 100;
x = (float *)malloc((n + 1) * sizeof(float));
y = (float *)malloc((n + 1) * sizeof(float));

/* Simulation:  (cos 2 PI f t) */
for(i = 1; i <= n; i++) {
  x[i] = xmax * (float)rand()/(float)RAND_MAX;
  y[i] = cos(2.2 * TWOPID * x[i] / xmax)
   + 0.5 * sin(7. * TWOPID * x[i] / xmax);
 }

printf("frequencies: %f %f \n", 2.2 / xmax, 7. / xmax);

/* Results: */
strcpy(filename,"period_in.dat");
if((fp = fopen(filename,"w")) == NULL) {
  printf("Fatal error opening %s\n", filename);
 }
else {
  for(i = 1; i <= n; i++) {
   fprintf(fp,"%f %f\n", x[i], y[i]);
   }
  fclose(fp);
  }

/* maximum frequency at 10.: */
fmax = 10.;
ofac = 8.;
jlp_period_0(x, y, n, ofac, fmax, &x_ave, &y_ave, &y_var, &fstep, &nout);

np = nout;
px = (float *)malloc((np + 1) * sizeof(float));
py = (float *)malloc((np + 1) * sizeof(float));

/* Mean Nyquist frequency: nf = (1 / <Delta x>) = n / xdif
* fstep = 1.0 / (xdif * ofac);
*  hifac : fraction of (average) Nyquist frequency to be used
* fmax = fstep * nout = 0.5 * ofac * hifac * n / (xdif * ofac) 
       = 0.5 * hifac * nf 
*/


jlp_period_1(x, y, n, ofac, px, py, np, nout, x_ave, y_ave, y_var, 
             fstep, &jmax, &prob);

printf("Output from period: max at py[%d]=%g prob=%f\n", jmax, py[jmax], prob);

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

free(x); 
free(y);
free(px); 
free(py);
return 0;
}
#endif
/*********************************************************************
* jlp_period_lomb0
* To compute parameters necessary for jlp_period_lomb1
* INPUT:
***********************************************************************/
void jlp_period_lomb0(float *x, float *y, int n, float ofac, float fmax,
                      float *x_ave, float *y_ave, float *y_var, 
                      float *fstep, int *nout)
{
int j;
float xdif, xmax, xmin;

/* Get mean and variance of the input data: */
avevar(y, n , y_ave, y_var);

#ifdef DEBUG
printf("period/Input data (y) average=%f variance=%f\n", *y_ave, *y_var);
#endif

/* Go through data to get range of abscissa: */
xmax = xmin = x[1];
for (j = 1; j <= n; j++) {
    if(x[j] > xmax) xmax = x[j];
    if(x[j] < xmin) xmin = x[j];
  }
xdif = xmax - xmin;
*x_ave = 0.5 * (xmax + xmin);
#ifdef DEBUG
printf("period/Input data (x) average=%f range=%f\n",*x_ave, xdif);
#endif

/* Step frequency: */
*fstep = 1.0 / (xdif * ofac);

/* Number of frequency terms: */
*nout = fmax / (*fstep);

 printf(" First frequency: %e (= step in frequency domain)\n", 
         *fstep);
/* Mean Nyquist frequency: nf = (1 / <Delta x>) = n / xdif */

 printf(" Mean Nyquist frequency = %e \n", (float)n / xdif);

 printf(" Last frequency: %e \n", (*nout) * (*fstep));

}
/*********************************************************************
* jlp_period_lomb1
* Modified from Numrecipes in C, p579
*
* Given n data points with abscissa x[1...n] (which need not to be equally
* spaced) and ordinates y[1...n] and given a desired oversampling 
* factor ofac (a typical value being 4 or larger), this routine fills 
* array px[1...np] with an increasing sequence of frequencies (not
* angular frequencies) up to hifac times the "average" Nyquist frequency,
* and fills array ny[1...np] with the values of the Lomb normalized 
* periodogram at those frequencies. 
*
* The arrays x and y are not altered.
* The dimension np of nx and ny must be large enough to contain the output,
* or an error results.
* The routine also returns jmax such that py[jmax] is the maximum element
* in py, and prob, an estimate of the significance of that maximum against
* the hypothesis of random noise. A small value of prob indicates that 
* a significant periodic noise is present.
*
* WARNING: arrays indices start at 1!
***********************************************************************/
void jlp_period_lomb1(float *x, float *y, int n, float ofac,
                      float *px, float *py, int np, int nout, 
                      float x_ave, float y_ave, float y_var, float fstep,
                      int *jmax, float *prob)
{
int i, j;
float cc, effm, expy, pymax, yy;
float ss, sumc, sumcy, sums, sumsy;
double arg;

if(nout > np) {
  printf("period/fatal error: output arrays too short in period\n");
  exit(-1);
  }

pymax = 0.0;

/* Main loop over the frequencies to be evaluated.
* First, loop over the data to get tau and related quantities: */
for (i = 1; i <= nout; i++) {
   px[i] = (float)i * fstep;

/* Loop over the data again to get the periodogram value: */
   sums = sumc = sumsy = sumcy = 0.0;
   for (j = 1; j <= n; j++) {
       arg = TWOPID * ((x[j] - x_ave) * px[i]);
       cc = cos(arg);
       ss = sin(arg);
       sums += ss * ss;
       sumc += cc * cc;

       yy = y[j] - y_ave;
       sumsy += yy * ss;
       sumcy += yy * cc;
   }
/*
py[i] = (1/(2 sig^2)) * 
( [ sum_j (y_j - y_ave) * cos(2 PI f_i (x_j - x_ave))]^2
    / sum_j cos^2(2 PI f_i (x_j - x_ave))
*/
   py[i] = (sumcy * sumcy / sumc + sumsy * sumsy / sums) / ( 2. * y_var);
/*
#ifdef DEBUG
printf("period/px[%d]=%f py=%f \n", i, px[i], py[i]);
#endif
*/
   if(py[i] >= pymax) pymax = py[(*jmax=i)]; 

 }

/* Evaluate statistical significance of the maximum:  */
 expy = exp(-pymax);
/* effm = M parameter (cf. p 576-577) */
 effm = 2.0 * nout / ofac;
 *prob = effm * expy;
 if(*prob > 0.01) *prob = 1.0 - pow(1.0 - expy, effm); 

}
/*********************************************************************
* period
* Cf. Numrecipes in C, p579
*
* Given n data points with abscissa x[1...n] (which need not to be equally
* spaced) and ordinates y[1...n] and given a desired oversampling 
* factor ofac (a typical value being 4 or larger), this routine fills 
* array px[1...np] with an increasing sequence of frequencies (not
* angular frequencies) up to hifac times the "average" Nyquist frequency,
* and fills array ny[1...np] with the values of the Lomb normalized 
* periodogram at those frequencies. 
*
* The arrays x and y are not altered.
* The dimension np of nx and ny must be large enough to contain the output,
* or an error results.
* The routine also returns jmax such that py[jmax] is the maximum element
* in py, and prob, an estimate of the significance of that maximum against
* the hypothesis of random noise. A small value of prob indicates that 
* a significant periodic noise is present.
*
* WARNING: arrays indices start at 1!
***********************************************************************/
void period(float *x, float *y, int n, float ofac, float hifac,
            float *px, float *py, int np, int *nout, int *jmax, float *prob)
{
int i, j;
float y_ave, y_var, c, cc, cwtau, effm, expy, pymax;
float s , ss, sumc, sumcy, sums, sumsh, sumsy, swtau, wtau, x_ave;
float xdif, xmax, xmin, yy, fstep;
double arg, wtemp, *wi, *wpi, *wpr, *wr;

/* Double precision is needed for stability of trigonometric recurrences */
wi = (double *)malloc((n + 1) * sizeof(double));
wpi = (double *)malloc((n + 1) * sizeof(double));
wpr = (double *)malloc((n + 1) * sizeof(double));
wr = (double *)malloc((n + 1) * sizeof(double));

*nout = 0.5 * ofac * hifac * n;
if(*nout > np) {
  printf("period/fatal error: output arrays too short in period\n");
  exit(-1);
  }

/* Get mean and variance of the input data: */
avevar(y, n , &y_ave, &y_var);

#ifdef DEBUG
printf("period/Input data (y) average=%f variance=%f\n", y_ave, y_var);
#endif

/* Go through data to get range of abscissa: */
xmax = xmin = x[1];
for (j = 1; j <= n; j++) {
    if(x[j] > xmax) xmax = x[j];
    if(x[j] < xmin) xmin = x[j];
  }
xdif = xmax - xmin;
x_ave = 0.5 * (xmax + xmin);
#ifdef DEBUG
printf("period/Input data (x) average=%f range=%f\n",x_ave, xdif);
#endif

pymax = 0.0;
/* Starting frequency: */
fstep = 1.0 / (xdif * ofac);

/* Initialize values for the trigonometric recurrences (p 178) 
* at each data point. */
/* cos(u + step) = cos u + (alpha cos u - beta sin u)
*  sin(u + step) = sin u + (alpha sin u + beta cos u)
* with alpha = 2 sin^2(step/2) and  beta = sin(step) 
* here step[j] = 2 PI (x[j] - x_ave) * fstep 
*/
for (j = 1; j <= n; j++) {
    arg = TWOPID * ((x[j] - x_ave) * fstep);
/* wpr[j] = alpha and wpi[j] = beta for data point #j: */
    wpr[j] = - 2.0 * SQR( sin(0.5 * arg)); 
    wpi[j] = sin(arg);
/* First values for wr[j] = cos(j*step) and wi[j] = sin(j*step)
*/
    wr[j] = cos(arg);
    wi[j] = wpi[j];
 }

/* Main loop over the frequencies to be evaluated.
* First, loop over the data to get tau and related quantities: */
for (i = 1; i <= (*nout); i++) {
   px[i] = (float)i * fstep;
   sumsh = sumc = 0.0;
   for (j = 1; j <= n; j++) {
       c = wr[j];
       s = wi[j];
       sumsh += s * c;
       sumc += (c - s) * (c + s);
   }
   wtau = 0.5 * atan2(2.0 * sumsh, sumc);
   swtau = sin(wtau);
   cwtau = cos(wtau);

/* Loop over the data again to get the periodogram value: */
   sums = sumc = sumsy = sumcy = 0.0;
   for (j = 1; j <= n; j++) {
       c = wr[j];
       s = wi[j];
/* cc = cos(omega (t - tau)) */
/* cos(a-b) = cos a cos b + sin a sin b:
* Hence cos(omega (t - tau)) = cos (omega t) cos (omega tau)
*                              + sin (omega t) sin (omega tau) */
       cc = c * cwtau + s * swtau;
/* ss = sin(omega (t - tau)) */
/* sin(a-b) = sin a cos b - cos a sin b:
* Hence sin(omega (t - tau)) = sin (omega t) cos (omega tau)
*                              - cos (omega t) sin (omega tau) */
       ss = s * cwtau - c * swtau;
       sums += ss * ss;
       sumc += cc * cc;

       yy = y[j] - y_ave;
       sumsy += yy * ss;
       sumcy += yy * cc;
/* Update the trigonometric recurrences:
* cos(u + step) = cos u + (alpha cos u - beta sin u)
* sin(u + step) = sin u + (alpha sin u + beta cos u)
* wpr[j] = alpha and wpi[j] = beta for data point #j: 
* 
* i.e., wr[j] = cos(j*step) and wi[j] = sin(j*step)
*/
       wr[j] = ((wtemp = wr[j]) * wpr[j] - wi[j] * wpi[j]) + wr[j];
       wi[j] = (wi[j] * wpr[j] + wtemp * wpi[j]) + wi[j];
#ifdef DEBUG
if(i == -1)
{
arg = TWOPID * ((x[j] - x_ave) * fstep);
printf("period/wr[%d]=%f wi=%f cos(j*step)=%f sin(j*step)=%f\n", 
          j, wr[j], wi[j], cos((double)j * arg), sin((double)j * arg));
}
#endif
   }
/*
py[i] = (1/(2 sig^2)) * 
( [ sum_j (y_j - y_ave) * cos(2 PI f_i (x_j - x_ave))]^2
    / sum_j cos^2(2 PI f_i (x_j - x_ave))
*/
   py[i] = (sumcy * sumcy / sumc + sumsy * sumsy / sums) / ( 2. * y_var);
/*
#ifdef DEBUG
printf("period/px[%d]=%f py=%f \n", i, px[i], py[i]);
#endif
*/
   if(py[i] >= pymax) pymax = py[(*jmax=i)]; 

 }

/* Evaluate statistical significance of the maximum:  */
 expy = exp(-pymax);
/* effm = M parameter (cf. p 576-577) */
 effm = 2.0 * (*nout) / ofac;
 *prob = effm * expy;
 if(*prob > 0.01) *prob = 1.0 - pow(1.0 - expy, effm); 

 printf(" First frequency: %e (= step in frequency domain)\n", 
         1.0 / (ofac * xdif));
 printf(" Last frequency: %e \n", (*nout) / (ofac * xdif));

free(wi);
free(wpi);
free(wpr);
free(wr);
}
/*********************************************************************
* fasper
* Cf. Numrecipes in C, p582
* Same as "period" but faster (with use of FFT)
*
* Given n data points with abscissa x[1...n] (which need not to be equally
* spaced) and ordinates y[1...n] and given a desired oversampling 
* factor ofac (a typical value being 4 or larger), this routine fills 
* array wk1[1...nwk] with a sequence of nout increasing frequencies (not
* angular frequencies) up to hifac times the "average" Nyquist frequency,
* and fills array wk2[1...nwk] with the values of the Lomb normalized 
* periodogram at those frequencies. 
*
* The arrays x and y are not altered.
* The dimension of wk1 and wk2 must be large enough for intermediate 
* work space, or an error results.
* The routine also returns jmax such that wk2[jmax] is the maximum element
* in wk2, and prob, an estimate of the significance of that maximum against
* the hypothesis of random noise. A small value of prob indicates that 
* a significant periodic noise is present.
***********************************************************************/
void fasper(float *x, float *y, unsigned long n, float ofac, float hifac,
            float *wk1, float *wk2, unsigned long nwk, unsigned long *nout,
            unsigned long *jmax, float *prob)
{
/* Number of interpolation points per 1/4 cycle 
of highest frequency */
/*
#define MACC 4
unsigned long j, k, ndim, nfreq, nfreqt;
*/
printf("fasper/Fatal error: fasper not implemented yet \n");
exit(-1);
}

/******************************************************************
* avervar (p 617)
* Given array data[1...n], returns its mean and its variance var.
*
* WARNING: arrays indices start at 1!
******************************************************************/ 
static void avevar(float *data, unsigned long n, float *ave, float *var)
{
register int i;
float s, ep;

*ave = 0.0;
for(i = 1; i <= n; i++) *ave += data[i];
*ave /= n;
*var = ep = 0.0;
for(i = 1; i <= n; i++) {
   s= data[i] - (*ave);
   ep += s;
   *var += s * s;
 }

/* corrected two-pass formula (14.1.8, p 613): */
 *var = (*var - ep * ep / (float)n) / (float)(n - 1);

}
