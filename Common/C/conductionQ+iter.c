#include <math.h>
#include "thrmlLib.h"


void conductionQ(int nz, double z[], double dt, double Qn, double Qnp1,
		double T[], double ti[], double rhoc[], double emiss,
		double Fgeotherm, double *Fsurf)  {
/************************************************************************ 
   conductionQ:  program to calculate the diffusion of temperature 
                 into the ground and thermal emission at the surface 
                 with variable thermal properties on irregular grid 
   Crank-Nicolson scheme, flux conservative 
                          uses Samar's radiation formula 
   Eqn: rhoc*T_t = (k*T_z)_z 
   BC (z=0): Q(t) + kT_z = em*sig*T^4 
   BC (z=L): heat flux = Fgeotherm 

   nz = number of grid points
   z[1..nz] = depths below surface
   dt = time step
   ti[1..nz] = thermal inertia [J m^-2 K^-1 s^-1/2]  
   rhoc[1..nz] = rho*c  where rho=density [kg m^-3] and 
                              c=specific heat [J K^-1 kg^-1] 
   ti and rhoc are not allowed to vary in the layers immediately adjacent 
               to the surface or the bottom 
   T[0..nz] = vertical temperature profile [K] (in- and output) 
   Qn,Qnp1 = net solar insolation at time steps n and n+1 [Watts/m^2] 
   emiss = emissivity 
   Fgeotherm = geothermal heat flux at bottom boundary [W/m^2] 
   Fsurf = heat flux at surface [W/m^2]  (output) 

   Grid: surface is at z=0 
         z[0]=0, z[2]=3*z[1], i.e., the top layer has half the width 
         T[0] is at z[0]=0; ...; T[i] is at z[i]
         rhoc[i], ti[i] are midway between z[i-1] and z[i]

   input arrays go from 1...nz (as in the Fortran version), 
         zeroth elements are not used, except that 
         T[0] is the surface temperature

   originally written by Samar Khatiwala, 2001 
   extended to variable thermal properties 
         and irregular grid by Norbert Schorghofer 
   converted from Fortran to C in 2019
************************************************************************/
  
  int i, iter;
  const double sigSB = 5.6704e-8;
  double a[nz+1], b[nz+1], c[nz+1], r[nz+1];
  double k[nz+1], k1dz, alpha[nz+1], gamma[nz+1], Tr;
  double arad, brad, ann, annp1, bn, buf, dz, beta;
  double Told[nz+1];
  
  /* set some constants */
  for (i=1; i<=nz; i++) {
     k[i] = ti[i] * ti[i] / rhoc[i];  // thermal conductivity
  }
  dz = 2.*z[1];
  beta = dt / rhoc[1] / (2.*dz*dz);  // assumes rhoc[0]=rhoc[1]
  alpha[1] = beta * k[2];
  gamma[1] = beta * k[1];
  for (i=2; i<nz; i++) {
     buf = 2. * dt / (rhoc[i] + rhoc[i+1]) / (z[i+1] - z[i-1]);
     alpha[i] = k[i+1] * buf / (z[i+1] - z[i]);
     gamma[i] = k[i] * buf / (z[i] - z[i-1]);
  }
  buf = dt / (z[nz] - z[nz-1]) / (z[nz] - z[nz-1]);
  gamma[nz] = k[nz] * buf / (2. * rhoc[nz]);  // assumes rhoc(nz+1)=rhoc(nz)
  
  k1dz = k[1] / dz;
  
  /* elements of tridiagonal matrix */
  for (i=1; i<nz; i++) {
     a[i] = -gamma[i];  //  a[1] is not used 
     b[i] = 1. + alpha[i] + gamma[i];  //  b[1] has to be reset at ever
     c[i] = -alpha[i];  //  c[nz] is not used 
  }
  a[nz] = -2.*gamma[nz];
  b[nz] = 1. + 2.*gamma[nz];
  
  Tr = T[0];    //   'reference' temperature
  iter = 0;
  for (i=1; i<=nz; i++) Told[i] = T[i];
 lbpredcorr: 

  /* Emission */
  arad = -3 * emiss * sigSB * Tr * Tr * Tr * Tr;
  brad = 2 * emiss * sigSB * Tr * Tr * Tr;
  ann = (Qn - arad) / (k1dz + brad);
  annp1 = (Qnp1 - arad) / (k1dz + brad);
  bn = (k1dz - brad) / (k1dz + brad);
  b[1] = 1. + alpha[1] + gamma[1] - gamma[1] * bn;
  
  /* Set RHS */
  r[1] = gamma[1]*(annp1 + ann) +
     (1.-alpha[1]-gamma[1]+gamma[1]*bn)*T[1] + alpha[1]*T[2];
  for (i=2; i<nz; i++) {
     r[i] = gamma[i]*T[i-1] + (1.-alpha[i]-gamma[i])*T[i] + alpha[i]*T[i+1];
  }
  r[nz] = 2.*gamma[nz]*T[nz-1] + (1.-2.*gamma[nz])*T[nz] +
     2.*dt/rhoc[nz]*Fgeotherm/(z[nz]-z[nz-1]);   // assumes rhoc[nz+1]=rhoc[nz]
  
  /*  Solve for T at n+1 */
  tridag(a, b, c, r, T, (unsigned long)nz);  // update by tridiagonal inversion 
  
  T[0] = 0.5 * (annp1 + bn * T[1] + T[1]);

  /* iterative predictor-corrector */
  if ((T[0] > 1.2*Tr || T[0] < 0.8*Tr) && iter<10) {  // linearization error expected
     /* redo until Tr is within 20% of new surface temperature
       (under most circumstances, the 20% threshold is never exceeded)
     */
     iter++;
     Tr = sqrt(Tr * T[0]);  // linearize around an intermediate temperature
     for (i=1; i<=nz; i++) T[i] = Told[i];
     goto lbpredcorr;
  }
  // if (iter>=5) printf("consider taking shorter time steps %d\n",iter);
  
  *Fsurf = -k[1] * (T[1] - T[0]) / z[1];  // heat flux into surface 

} /* conductionQ */
