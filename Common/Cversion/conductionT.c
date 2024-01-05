#include <math.h>
#include "thrmlLib.h"

void conductionT(int nz, double z[], double dt, double T[], double Tsurf,
		 double Tsurfp1, double ti[], double rhoc[], double Fgeotherm,
		 double *Fsurf)  {
/***********************************************************************
    conductionT:  program to calculate the diffusion of temperature 
                  into the ground with prescribed surface temperature 
                  and variable thermal properties on irregular grid
    Crank-Nicholson scheme, flux conservative
 
    Eqn: rhoc*T_t = (k*T_z)_z 
    BC (z=0): T=T(t)
    BC (z=L): heat flux = Fgeotherm
 
    nz = number of grid points
    dt = time step
    T = vertical temperature profile [K]  (in- and output)
    Tsurf, Tsurfp1 = surface temperatures at times n and n+1  
    ti = thermal inertia [J m^-2 K^-1 s^-1/2]
    rhoc = rho*c  heat capacity per volume [J m^-3 K^-1]
    ti and rhoc are not allowed to vary in the layers immediately
                adjacent to the surface or the bottom
    Fgeotherm = geothermal heat flux at bottom boundary [W/m^2]
    Fsurf = heat flux at surface [W/m^2]  (output)
 
    Grid: surface is at z=0
          T[1] is at z(1); ...; T[i] is at z[i]
          k[i] is midway between z[i-1] and z[i]
          rhoc[i] is midway between z[i-1] and z[i]
 ***********************************************************************/

    int i;
    double alpha[nz+1], k[nz+1], gamma[nz+1], buf;
    double a[nz+1], b[nz+1], c[nz+1], r[nz+1];
  
    /* set some constants */
    for (i=1; i<=nz; i++) {
        k[i] = ti[i] * ti[i] / rhoc[i]; // thermal conductivity
    }
    alpha[1] = k[2]*dt / rhoc[1] / (z[2]-z[1]) / z[2];
    gamma[1] = k[1]*dt / rhoc[1] / z[1] / z[2];
    for (i=2; i<nz; i++) {
        buf = 2. * dt / (rhoc[i]+rhoc[i+1]) / (z[i+1]-z[i-1]);
	alpha[i] = k[i+1] * buf / (z[i+1]-z[i]);
	gamma[i] = k[i] * buf / (z[i]-z[i-1]);
    }
    buf = dt / (z[nz]-z[nz-1]) / (z[nz]-z[nz-1]);
    gamma[nz] = k[nz] * buf / (2. * rhoc[nz]);  // assumes rhoc(nz+1)=rhoc(nz)
    
    /* elements of tridiagonal matrix */
    for (i=1; i<nz; i++) {
        a[i] = -gamma[i];   //  a[1] is not used
	b[i] = 1. + alpha[i] + gamma[i];
	c[i] = -alpha[i];   //  c[nz] is not used
    }
    a[nz] = -2.*gamma[nz];
    b[nz] = 1. + 2.*gamma[nz];
  
    /* Set RHS */
    r[1] = alpha[1]*T[2] + (1.-alpha[1]-gamma[1])*T[1] + gamma[1]*(Tsurf+Tsurfp1);
    for (i=2; i<nz; i++) {
        r[i] = gamma[i]*T[i-1] + (1.-alpha[i]-gamma[i])*T[i] + alpha[i]*T[i+1];
    }
    r[nz] = 2.*gamma[nz]*T[nz-1] + (1.-2.*gamma[nz])*T[nz] + 
        2.*dt/rhoc[nz]*Fgeotherm/(z[nz]-z[nz-1]); // assumes rhoc[nz+1]=rhoc[nz]

    /* Solve for T at n+1 */
    tridag(a,b,c,r,T,nz); // update by tridiagonal inversion
    
    *Fsurf = -k[1] * (T[1]-Tsurfp1) / z[1]; // heat flux into surface
    
} /* conductionT */

