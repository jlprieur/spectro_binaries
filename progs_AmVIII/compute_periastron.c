/*************************************************************************
* Compute the periastron from an ASCII file containg 3 columns
*    xxx, a sin i, e 
* xxx is the HD number for instance
* JLP
* Version 16/03/2006
*************************************************************************/
#include <stdio.h>
#include <math.h>
#include <string.h>

int main(int argc, char *argv[])
{
float val, asini, ee, pp;
char infile[60], outfile[60], comments[100], buffer[80];
int nval;
register int i, k;
FILE *fp_in, *fp_out;

/* To handle "runs periastron" */
for(k = 7; k > 0; k--) if(argc == k && argv[k-1][0] == '\0') argc = k-1;
if(argc != 3) {
  printf("Error, argc=%d \n Syntax is:\n", argc);
  printf("runs compute_periastron input_list output_file\n");
  printf("Format of input file: 3 columns (xxx, a sin i, e)\n", argc);
  return(-1);
 }

strcpy(infile, argv[1]);
strcpy(outfile, argv[2]);
printf("OK: reading %s and writing %s\n", infile, outfile);
sprintf(comments,"From %s\n", infile);

if((fp_in = fopen(infile, "r")) == NULL) {
 fprintf(stderr, "Fatal error opening input file: %s\n", infile);
 return(-1);
 }

if((fp_out = fopen(outfile, "w")) == NULL) {
 fprintf(stderr, "Fatal error opening input file: %s\n", outfile);
 return(-1);
 }

while(!feof(fp_in)) {
  if(fgets(buffer, 80, fp_in)) {
#ifdef DEBUG
  printf("%s", buffer);
#endif
    if(buffer[0] == '%') {
      fprintf(fp_out,"%s\n", buffer); 
    } else {
    nval = sscanf(buffer,"%f %f %f\n", &val, &asini, &ee); 
#ifdef DEBUG
    printf("nval= %d val=%f asini=%f ee=%f\n", nval, val, asini, ee);
#endif
    if(nval == 3) {
      pp = asini * (1. - ee);
      fprintf(fp_out,"%f %f %f %f\n", val, asini, ee, pp); 
      }
    }
  }
}

fclose(fp_out);
fclose(fp_in);
return(0);
}
