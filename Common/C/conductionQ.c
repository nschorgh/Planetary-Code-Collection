#include <math.h>
#include "thrmlLib.h"


void cranknQ(int nz, double z[], double dt, double Qn, double Qnp1,
	     double T[], double ti[], double rhoc[], double emiss,
	     double Fgeotherm, double *Fsurf)  {
/************************************************************************ 
   cranknnQ:  program to calculate the diffusion of temperature into the
              ground and thermal emission at the surface with variable
              thermal properties on irregular grid 
   Crank-Nicolson scheme, flux conservative, uses Samar's radiation formula

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
   added Volterra predictor 1/2024 -norbert
************************************************************************/
  
  int i;
  const double sigSB = 5.6704e-8, pi=3.1415926535897931;
  double a[nz+1], b[nz+1], c[nz+1], r[nz+1];
  double k[nz+1], k1dz, alpha[nz+1], gamma[nz+1], Tr;
  double arad, brad, ann, annp1, bn, buf, dz, beta;
  double seb, Tpred;
  
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

  /* Volterra predictor (optional) */
  *Fsurf = - k[1] * ( T[1]-T[0] ) / z[1];  // heat flux;
  seb = -(*Fsurf) -emiss*sigSB*pow(T[0],4) + (2*Qnp1 + Qn)/3.;
  // Tpred = T[0] + sqrt(4*dt/pi) / ti(1) * seb;  ! 1st order
  Tpred = T[0] + seb / ( sqrt(pi/(4.*dt))*ti[1] + 8./3.*emiss*sigSB*pow(T[0],3) );
  Tr = (T[0]+Tpred)/2.;  // better reference temperature
  
  /* Emission */
  //Tr = T[0];    //   'reference' temperature
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

  *Fsurf = -k[1] * (T[1] - T[0]) / z[1];  // heat flux into surface 

} /* cranknQ */



void conductionQ(int nz, double z[], double dt, double Qn, double Qnp1,
		double T[], double ti[], double rhoc[], double emiss,
		double Fgeotherm, double *Fsurf)  {
  /*
    conductionQ:  wrapper for cranknQ, which improves stability
    Arguments and restrictions are the same as for function cranknQ above.
    created wrapper using flux smoothing 1/2024  
  */
  int i, j;
  const int Ni=5;  // for flux smoothing
  double Told[nz+1], Qartiold, Qarti, avFsurf;

  for (i=0; i<=nz; i++) Told[i] = T[i];
  
  cranknQ(nz, z, dt, Qn, Qnp1, T, ti, rhoc, emiss, Fgeotherm, Fsurf);

  /* artificial flux smoothing */
  if ( T[0]>1.2*Told[0] || T[0]<0.8*Told[0] ) { // linearization error
    for (i=0; i<=nz; i++) T[i] = Told[i];
    avFsurf = 0.;
    for (j=1; j<=Ni; j++) {
      Qartiold = ( (Ni-j+1)*Qn + (j-1)*Qnp1 ) / Ni;
      Qarti    = ( (Ni-j)*Qn + j*Qnp1 ) / Ni;
      cranknQ(nz,z,dt/Ni,Qartiold,Qarti,T,ti,rhoc,emiss,Fgeotherm,Fsurf);
      avFsurf += *Fsurf;
    }
    *Fsurf = avFsurf/Ni;
  }
  
} /* conductionQ */
