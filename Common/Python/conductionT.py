import numpy as np
import scipy.linalg


def conductionT(nz,z,dt,T,Tsurf,Tsurfp1,ti,rhoc,Fgeotherm,Fsurf):
#***********************************************************************
#   conductionT:  program to calculate the diffusion of temperature 
#                 into the ground with prescribed surface temperature 
#                 and variable thermal properties on irregular grid
#   Crank-Nicholson scheme, flux conservative
#
#   Eqn: rhoc*T_t = (k*T_z)_z 
#   BC (z=0): T=T(t)
#   BC (z=L): heat flux = Fgeotherm
#
#   nz = number of grid points
#   dt = time step
#   T = vertical temperature profile [K]  (in- and output)
#   Tsurf, Tsurfp1 = surface temperatures at times n and n+1  
#   ti = thermal inertia [J m^-2 K^-1 s^-1/2]  VECTOR
#   rhoc = rho*c  heat capacity per volume [J m^-3 K^-1]  VECTOR
#   ti and rhoc are not allowed to vary in the layers immediately
#               adjacent to the surface or the bottom
#   Fgeotherm = geothermal heat flux at bottom boundary [W/m^2]
#   Fsurf = heat flux at surface [W/m^2]  (output)
#
#   Grid: surface is at z=0
#         T[1] is at z[1]; ...; T[i] is at z[i]
#         k[i] is midway between z[i-1] and z[i]
#         rhoc[i] is midway between z[i-1] and z[i]
#
#   converted to Python 3/2021
#***********************************************************************

    alpha = np.zeros(nz+1)
    beta = np.zeros(nz+1)
    gamma = np.zeros(nz+1)
  
    # set some constants
    k = np.zeros(nz+1)
    k = ti[:]**2/rhoc[:] # thermal conductivity
    alpha[1] = k[2]*dt/rhoc[1]/(z[2]-z[1])/z[2] 
    gamma[1] = k[1]*dt/rhoc[1]/z[1]/z[2]
    for i in range(2,nz):
        buf = dt/(z[i+1]-z[i-1])
        alpha[i] = k[i+1]*buf*2./(rhoc[i]+rhoc[i+1])/(z[i+1]-z[i])
        gamma[i] = k[i]*buf*2./(rhoc[i]+rhoc[i+1])/(z[i]-z[i-1])
    buf = dt/(z[nz]-z[nz-1])**2
    gamma[nz] = k[nz]*buf/(rhoc[nz]+rhoc[nz]) # assumes rhoc[nz+1]=rhoc[nz]
  
    # elements of tridiagonal matrix
    a = -gamma[:]   #  a[1] is not used
    b = 1. + alpha[:] + gamma[:]
    c = -alpha[:]   #  c[nz] is not used
    b[nz] = 1. + gamma[nz]

    # Set RHS
    r = np.zeros(nz+1)
    r[1] = alpha[1]*T[2] + (1.-alpha[1]-gamma[1])*T[1] + gamma[1]*(Tsurf+Tsurfp1)
    for i in range(2,nz):
        r[i] = gamma[i]*T[i-1] + (1.-alpha[i]-gamma[i])*T[i] + alpha[i]*T[i+1]
    r[nz] = gamma[nz]*T[nz-1] + (1.-gamma[nz])*T[nz] + \
        dt/rhoc[nz]*Fgeotherm/(z[nz]-z[nz-1]) # assumes rhoc[nz+1]=rhoc[nz]

    # special matrix for solve_banded
    D = np.zeros((3,nz))
    D[0,1:] = c[1:-1]
    D[1,:] = b[1:]
    D[2,:-1] = a[2:]

    # Solve for T at n+1
    T[1:] = scipy.linalg.solve_banded((1,1), D, r[1:])
    T[0] = Tsurfp1
    
    Fsurf = -k[1]*(T[1]-Tsurfp1)/z[1] # heat flux into surface
  

