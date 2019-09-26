#include <math.h>
#include <stdio.h>
#include "thrmlLib.h"


int main(void)  {
/* *********************************************************************/
/*     test Crank-Nicolson subroutine                                  */
/* *********************************************************************/
  
    const int nz = 60;
    const double pi=3.1415926535897931, zero=0., d2r=pi/180.;
    const double Period=88775.244*670, emiss=1., Fgeo = 0.2;
    int i, n, STEPSPERSOL;
    double latitude, albedo;
    double thIn; // thermal inertia
    double z[nz+1], zmax, zfac;
    double T[nz+1]; // T[0] is surface temperature
    double dt, time, Qn, Qnp1;
    double Rau, Decl, HA;
    double delta, rhocv[nz+1], ti[nz+1], Fsurf;
    FILE *fout0, *fout1, *fout2;

    STEPSPERSOL = 120;
    dt = Period / STEPSPERSOL;
    thIn = 120.;
    albedo = 0.2;
    latitude = 5.; // [degree]
    zmax = 2.5; zfac = 1.05;  /* set the spatial resolution */
    
    for (i=1; i <= nz; i++) {
       rhocv[i] = 1200.*800.;
       ti[i] = thIn;
    }
    delta = thIn / rhocv[1] * sqrt(Period/pi);
    printf("Skin depth= %f\n", delta);
    Rau = 1.52;
    Decl = 0.;
    
    printf("Time step=%f\n", dt);
    printf("zmax=%f\n", zmax);
    printf("Thermal inertia=%f  Period=%f\n", thIn, Period);
    printf("Heat Flux at bottom boundary=%f\n",-Fgeo);

    fout1 = fopen("Tsurface","wt"); // surface temperature
    fout2 = fopen("Tprofile","wt"); // temperature profle
    
    for (i=0; i <= nz; i++) T[i] = 210.;
    setgrid(nz, z, zmax, zfac);
    fout0 = fopen("z","wt");
    for (i=1; i <= nz; i++) fprintf(fout0,"%f ",z[i]);
    fprintf(fout0,"\n");
    fclose(fout0);
    
    // print *,z(1:nz) 
    latitude *= d2r;
    HA = 0.; // hour angle
    /*  solar insolation */
    Qn = (1-albedo)*flux_noatm(Rau, Decl, latitude, HA, zero, zero);
    
    for (n = 0; n <= 50000; n++) {
       time = (n + 1) * dt;  //  time at n+1
       //printf("%f %f\n",time,Qn);
       HA = 2*pi*fmod( time/Period, 1.); //  hour angle 
       Qnp1 = (1-albedo)*flux_noatm(Rau, Decl, latitude, HA, zero, zero);
       conductionQ(nz, z, dt, Qn, Qnp1, T, ti, rhocv, emiss, Fgeo, &Fsurf);
       Qn = Qnp1;
       if (n%3==0) {
	  fprintf(fout1, "%12.6f  %9.3f %9.3f\n", time/Period,T[0],T[nz]);
       }
       if (n > 50000 - STEPSPERSOL) {
	  if (n%10 == 0) {
	     fprintf(fout2, "%7.2f  ", T[0]);
	     for (i=1; i <= nz; i++)  fprintf(fout2, "%7.2f  ",T[i]);
	     fprintf(fout2, "\n");
	  }
       }
    } /* end the loop */
    
    fclose(fout1);
    fclose(fout2);
    
    return 0;
} /* testcrankQ */
