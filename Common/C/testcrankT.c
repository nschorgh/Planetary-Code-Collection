#include <math.h>
#include <stdio.h>
#include "thrmlLib.h"


int main(void)  {
/* *********************************************************************/
/*     test Crank-Nicolson subroutine                                  */
/* *********************************************************************/

  const double pi=3.1415926535897931, Period=88775.244*670;
  const double Fgeo = 0., Ta=30., Tm=190.;

  const int nz = 70;
  int n, i, STEPSPERSOL;
  double T[nz+1], time, dt;
  double thIn;  // thermalInertia
  double rhocv[nz+1]; // volumetric heat capacity, rho*c
  double Tsurf, Tsurfp1;
  double delta, ti[nz+1], Fsurf, z[nz+1];
  double zmax, zfac;
  FILE *fout0, *fout;
 
  STEPSPERSOL = 120;
  dt = Period/STEPSPERSOL;
  zmax = 2.5; zfac=1.02;
  thIn = 120.;

  for (i=1; i <= nz; i++) {
    rhocv[i] = 1200.*800.;
    ti[i] = thIn;
  }
  
  delta = thIn/rhocv[1]*sqrt(Period/pi);
  printf("Skin depth= %f\n", delta);
  printf("zmax=%f\n", zmax);
  printf("Time step=%f\n", dt);
  printf("Thermal inertia=%f  Period=%f\n", thIn, Period);
  
  /* Initialize */
  fout = fopen("Tprofile","wt"); // temperature profle

  for (i=0; i<=nz; i++) T[i] = 0.;

  setgrid(nz,z,zmax,zfac);
  fout0 = fopen("z","wt");
  for (i=1; i <= nz; i++) fprintf(fout0,"%f ",z[i]);
  fprintf(fout0,"\n");
  fclose(fout0);
	
  time = 0.;

  Tsurf = Tm + Ta*sin(2*pi*time/Period);

  for (n = 0; n <= 50000; n++) {
    // printf("%f %f\n",time/Period,Qn);
    time = (n+1)*dt;     // time at n+1;
    Tsurfp1 = Tm + Ta*sin(2*pi*time/Period);
    conductionT(nz, z, dt, T, Tsurf, Tsurfp1, ti, rhocv, Fgeo, &Fsurf);
    Tsurf = Tsurfp1;

    //printf("%f %f\n",time/Period,Tsurf);
    
    //write 12 profiles from the last sol
    if (n>50000-STEPSPERSOL) {
      if (n%10==9) {  // synched with test_Tprofile.m
	printf("%f %f\n",time/Period,Tsurf);
	fprintf(fout, "%7.2f ", Tsurf);
	for (i=1; i<=nz; i++)  fprintf(fout, "%7.2f ",T[i]);
	fprintf(fout, "\n");
	
      }
    }
  } // end of loop
  fclose(fout);
}
