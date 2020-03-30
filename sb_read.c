/***************************************************************
* To read spectroscopic binary data
*
* Version 27/05/2002
***************************************************************/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>

/*
#define DEBUG
*/
#define STOP_ERROR { printf("Fatal error reading input file\n"); fclose(fp); exit(-1);}

/* Prototypes of functions defined here:
*/
int period_read_data(FILE *fp, float *x, float *y, int ic, int *n);
void period_hrange(float *y, int n, float *half_range, float *mean);

/****************************************************************
* read data from Coravel VR arrays
*
* 1=primary 2=secondary 3=prim. & sec. 4=3rd body 
*
* Warning arrays start at one here! 
****************************************************************/
int period_read_data(FILE *fp, float *x, float *y, int ic, int *n)
{
float *y1, *y2, k1, k2, v01, v02, weight;
register int i, j;
char buffer[80];

#ifdef DEBUG
printf("period_read_data/ic = %d \n",ic);
#endif

/* Format x, w, y with zeros when no data */
if(ic == 0) {
/* Warning arrays start at one here! */
   j = 0;
   for(i = 1; i <= (*n); i++) {
     if(fgets(buffer,80,fp) == NULL) STOP_ERROR
     j++;
     sscanf(buffer,"%f %*f %f", &x[j], &y[j]);
/*
     printf("%d: %f %f \n", j, x[j], y[j]);
*/
     if(y[j] == 0.) j--;
     }
   (*n) = j;
   }
else if(ic == 1) {
/* Warning arrays start at one here! */
   j = 0;
   for(i = 1; i <= (*n); i++) {
     if(fgets(buffer,80,fp) == NULL) STOP_ERROR
     j++;
     sscanf(buffer,"%f %f", &x[j], &y[j]);
/*
     printf("%d: %f %f \n", j, x[j], y[j]);
*/
     }
   (*n) = j;
   }
else if(ic == 2) {
   j = 0;
   for(i = 1; i <= (*n); i++) {
     if(fgets(buffer,80,fp) == NULL) STOP_ERROR
/* Warning arrays start at one here! */
     j++;
     sscanf(buffer,"%f %*f %*f %f %f", &x[j], &y[j], &weight);
/*
     printf("%d: %f %f %f\n", j, x[j], y[j], weight);
*/
     if(weight == 0) j--;
     }
   (*n) = j;
   }
else {
   y1 = (float *)malloc(((*n) + 1) * sizeof(float));
   y2 = (float *)malloc(((*n) + 1) * sizeof(float));
   j = 0;
   for(i = 1; i <= (*n); i++) {
     if(fgets(buffer,80,fp) == NULL) STOP_ERROR
/* Warning arrays start at one here! */
     j++;
     sscanf(buffer,"%f %f %*f %f %f", &x[j], &y1[j], &y2[j], &weight);
/*
     printf("%d: x=%f y1=%f y2=%f w=%f\n", j, x[j], y1[j], y2[j], weight);
*/
     if(weight == 0) j--;
     }
   (*n) = j;
   period_hrange(y1, *n, &k1, &v01);
#ifdef DEBUG
   printf(" k1=%f v01=%f\n", k1, v01);
#endif
   period_hrange(y2, *n, &k2, &v02);
#ifdef DEBUG
   printf(" k2=%f v02=%f\n", k2, v02);
#endif
   if(ic == 3) {
     for(j = 1; j <= (*n); j++) 
          y[j] = (y1[j] - v01) / k1 - (y2[j] - v02) / k2;
     }
   else {
     for(j = 1; j <= (*n); j++) 
          y[j] = (y1[j] - v01) / k1 + (y2[j] - v02) / k2;
     }
   free(y1);
   free(y2);
  }

fclose(fp);

#ifdef DEBUG
printf("period_read_data/n=%d data points\n",(*n));
#endif

return(0);
}
/*****************************************************************
*
* To compute the half range and the mean of an array.
*
*****************************************************************/
void period_hrange(float *y, int n, float *half_range, float *mean)
{
double sum;
float max, min;
register int i;

if( n <= 0) {
   printf("period_range/Fatal error: n=%d\n", n);
   exit(-1);
   }

sum = 0.;
max = min = y[1];
for(i = 1; i <= n; i++)
  {
  if(min > y[i]) min = y[i];
  if(max < y[i]) max = y[i];
  sum += y[i];
  }
*mean = sum / (float)n;
*half_range = (max - min) / 2.;
}
