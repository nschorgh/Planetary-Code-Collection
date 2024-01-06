#include <stdio.h>
#include <stdlib.h>


void tridag(double a[], double b[], double c[], double r[],
	    double u[], unsigned long n)
/* Tridiagonal solver */
{
   unsigned long j;
   double bet, gam[n+1];

   if (b[1] == 0.0) {
      fprintf(stderr,"%s\n","tridag: rewrite equations");
      exit(1);
   }
   bet = b[1]; 
   u[1] = r[1] / bet;
   for (j=2; j<=n; j++) {
      gam[j] = c[j-1]/bet;
      bet = b[j]-a[j]*gam[j];
      if (bet == 0.0)	{
	 fprintf(stderr,"%s\n","tridag failed");
	 exit(1);
      }
      u[j] = (r[j]-a[j]*u[j-1]) / bet;
   }
   for (j=n-1; j>=1; j--)
      u[j] -= gam[j+1]*u[j+1];
}
/* based on Numerical Recipes' tridag.c, but heavily modified */
