#include <stdio.h>
#include <stdlib.h>


void tridag(double a[], double b[], double c[], double r[],
	    double u[], unsigned long N)
/* solves Ax = r, where A is a tridiagonal matrix consisting of vectors a, b, c
   N ... number of equations
   indices start with 1 (Fortran convention); index 0 is not used
   a ... subdiagonal, index 1 is not used
   b ... main diagonal, indexed from 1,...,N
   c ... superdiagonal, index N is not used
   Method: A=LU, first solve Ls=r and then Ux=s
   stable when |b_i| > |a_i|+|c_i|
*/
{
   unsigned long j;
   double bet, scratch[N+1];

   if (b[1] == 0.0) {
      fprintf(stderr,"%s\n","tridag: rewrite equations");
      exit(1);
   }
   bet = b[1];
   u[1] = r[1] / bet;
   for (j=2; j<=N; j++) {  // downward sweep
      scratch[j] = c[j-1] / bet;
      bet = b[j] - a[j]*scratch[j];
      if (bet == 0.0) {
	 fprintf(stderr,"%s\n","tridag failed");
	 exit(1);
      }
      u[j] = ( r[j] - a[j]*u[j-1]) / bet;
   }
   for (j=N-1; j>=1; j--)  // upward sweep
      u[j] -= scratch[j+1]*u[j+1];
}

