import numpy as np
import scipy.linalg


def conductionQ(nz, z, dt, Qn, Qnp1, T, ti, rhoc, emiss, Fgeotherm, Fsurf):
    """
    conductionQ: program to calculate the diffusion of temperature
                 into the ground and thermal emission at the surface
                 with variable thermal properties on irregular grid Crank-Nicolson scheme, flux conservative uses Samar's radiation formula

    Eqn: rhoc*T_t = (k*T_z)_z
    BC (z=0): Q(t) + kT_z = em*sig*T^4
    BC (z=L): heat flux = Fgeotherm

    nz = number of grid points (not counting the surface)
    dt = time step
    Qn,Qnp1 = net solar insolation at time steps n and n+1 [W/m^2]
    T = vertical temperature profile [K]  (in- and output)
    T[0] = surface temperature [K]  (in- and output)
    ti = thermal inertia [J m^-2 K^-1 s^-1/2]  VECTOR
    rhoc = rho*c  VECTOR where rho=density [kg/m^3] and
                              c=specific heat [J K^-1 kg^-1]
    ti and rhoc are not allowed to vary in the layers immediately
               adjacent to the surface or the bottom
    emiss = emissivity
    Fgeotherm = geothermal heat flux at bottom boundary [W/m^2]
    Fsurf = heat flux at surface [W/m^2]  (output)

    Grid: surface is at z=0
         z[0]=0, z[2]=3*z[1], i.e., the top layer has half the width
         T[1] is at z[1]; ...; T[i] is at z[i]
         k[i], rhoc[i], ti[i] are midway between z[i-1] and z[i]

    originally written by Samar Khatiwala, 2001
    extended to variable thermal properties
         and irregular grid by Norbert Schorghofer, 2004
    added predictor-corrector 9/2019
    converted to Python 3/2021
    Speed optimization by Cyril Mergny 5/2023
    """
    sigSB = 5.6704e-8

    # set some constants
    k = np.empty(nz+1)
    k = ti[:] ** 2 / rhoc[:]  # thermal conductivity
    alpha = np.empty(nz+1)
    gamma = np.empty(nz+1)

    dz = 2 * z[1]
    beta = dt / (rhoc[1] + rhoc[2]) / dz**2
    alpha[1] = beta * k[2]
    gamma[1] = beta * k[1]

    buf2 = dt * 2.0 / (rhoc[2:-1] + rhoc[3:]) / (z[3:] - z[1:-2])
    alpha[2:-1] = k[3:] * buf2 / (z[3:] - z[2:-1])
    gamma[2:-1] = k[2:-1] * buf2 / (z[2:-1] - z[1:-2])

    alpha[nz] = 0.0  
    gamma[nz] = dt * k[nz] / (2 * rhoc[nz]) / (z[nz]-z[nz-1])**2

    # for i in range(2,nz):  # 2 ... nz-1
    #    buf = dt / (z[i+1]-z[i-1])
    #    alpha[i] = 2*k[i+1]*buf/(rhoc[i]+rhoc[i+1])/(z[i+1]-z[i])
    #    gamma[i] = 2*k[i]*buf/(rhoc[i]+rhoc[i+1])/(z[i]-z[i-1])

    k1dz = k[1] / dz

    # elements of tridiagonal matrix
    # special matrix for solve_banded
    D = np.zeros((3, nz))
    D[0, 1:] = -alpha[1:-1]  # coefficient 'c'
    D[1, :] = 1.0 + alpha[1:] + gamma[1:]  # coefficient 'b'
    D[2, :-1] = -gamma[2:]  # coefficient 'a'
    # b[1] has to be reset at every timestep
    D[1,-1] = 1. + 2*gamma[nz]
    D[2,-2] = -2*gamma[nz]
    
    Tr = T[0]  # 'reference' temperature
    iter = 0
    Told = np.empty(nz+1)
    Told[:] = T[:]

    while iter < 10:
        if iter > 0:
            Tr = np.sqrt(Tr * T[0])  # linearize around an intermediate temperature
            T[1:] = Told[1:]

        # Emission
        arad = -3.0 * emiss * sigSB * Tr**4
        brad = 2.0 * emiss * sigSB * Tr**3
        ann = (Qn - arad) / (k1dz + brad)
        annp1 = (Qnp1 - arad) / (k1dz + brad)
        bn = (k1dz - brad) / (k1dz + brad)
        b1 = 1.0 + alpha[1] + gamma[1] - gamma[1] * bn  # b[1]

        # Set RHS
        r = np.empty(nz+1)
        r[1] = (
            gamma[1] * (annp1 + ann)
            + (1.0 - alpha[1] - gamma[1] + gamma[1] * bn) * T[1]
            + alpha[1] * T[2]
        )
        # for i in range(2,nz): # 2...nz-1
        #    r[i] = gamma[i]*T[i-1] + (1.-alpha[i]-gamma[i])*T[i]+ alpha[i]*T[i+1]
        r[2:-1] = (
            gamma[2:-1] * T[1:-2]
            + (1.0 - alpha[2:-1] - gamma[2:-1]) * T[2:-1]
            + alpha[2:-1] * T[3:]
        )
        r[nz] = (
            2*gamma[nz] * T[nz-1]
            + (1.0 - 2*gamma[nz]) * T[nz]
            + 2 * dt / rhoc[nz] * Fgeotherm / (z[nz] - z[nz-1])
        )

        D[1, 0] = b1  # coefficient b[1]

        # Solve for T at n+1
        T[1:] = scipy.linalg.solve_banded((1, 1), D, r[1:])
        T[0] = 0.5 * (annp1 + bn * T[1] + T[1])  # (T0+T1)/2

        # iterative predictor-corrector
        if T[0] < 1.2 * Tr and T[0] > 0.8 * Tr:
            break

        iter += 1
        # redo until Tr is within 20% of new surface temperature
        # (under most circumstances, the 20% threshold is never exceeded)

    Fsurf = -k[1] * (T[1] - T[0]) / z[1]  # heat flux into surface
