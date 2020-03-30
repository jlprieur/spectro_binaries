/*********************************************************************
*
* JLP March 2014
*********************************************************************/
#include<stdio.h>
#include<math.h>

int main(int argc, char* argv[])
{
double mv, varpi, Dmv;
double Mv, L1, L2, ss, Mv1, Mv2;
double Cv = -0.1;

printf("Enter mV_global, parallax, Delta mV:\n");
scanf("%lf,%lf,%lf", &mv, &varpi, &Dmv);
printf("OK: mv=%.3f varpi=%f Dmv=%.3f \n", mv, varpi, Dmv);

Mv = mv + 5. * log10(varpi) + 5.;
printf("Mv = %.3f\n", Mv);

ss = 1. + pow(10., (-0.4 * Dmv));
L1 = pow(10., (-0.4 * (Mv -Cv - 4.75))) / ss;
L2 = pow(10., (-0.4 * (Mv + Dmv -Cv - 4.75))) / ss;

printf("L1 = %f    L2=%f  ss=%f \n", L1, L2, ss);

Mv1 = -2.5 * log10(L1) + 4.75 - Cv;
Mv2 = -2.5 * log10(L2) + 4.75 - Cv;

printf("Mv1=%f    Mv2=%f\n", Mv1, Mv2);

return(0);
}
