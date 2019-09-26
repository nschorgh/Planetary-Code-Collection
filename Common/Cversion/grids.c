#include <math.h>


void setgrid(int nz, double z[], double zmax, double zfac)
{
/*   construct regularly or geometrically spaced 1D grid 
     z[n]-z[1] = 2*z(1)*(zfac**(n-1)-1)/(zfac-1) 
     choice of z[1] and z[2] is compatible with conductionQ 
*/
   int i;
   double dz;

   dz = zmax / nz;
   for (i=1; i<=nz; i++)  z[i] = (i - 0.5) * dz;

   if (zfac > 1.) {
      dz = zmax / ( 3. + 2. * zfac * ( pow(zfac,nz-2) - 1.) / (zfac - 1.) );
      z[1] = dz;
      z[2] = 3*z[1];
      for (i=3; i<=nz; i++) {
	 z[i] = (1. + zfac) * z[i-1] - zfac * z[i-2];
	 // same as z[i+1] - z[i] = zfac * (z[i] - z[i-1]) 
      }
   }
}



