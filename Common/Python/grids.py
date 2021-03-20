import numpy as np


def setgrid(nz, zmax, zfac):
    # construct regularly or geometrically spaced 1D grid
    # z(n)-z(1) = 2*z(1)*(zfac**(n-1)-1)/(zfac-1)
    # choice of z(1) and z(2) is compatible with conductionQ
  
    z = np.zeros(nz+1)
  
    if zfac>1.:
        dz = zmax/(3.+2.*zfac*(zfac**(nz-2)-1.)/(zfac-1.))
        z[1] = dz
        z[2] = 3*z[1]
        for i in range(3,nz+1):
            z[i] = (1.+zfac)*z[i-1] - zfac*z[i-2]
            # z[i] = z[i-1] + zfac*(z[i-1]-z[i-2]) # equivalent
    else:
        dz = zmax/nz
        for i in range(1,nz+1):
            z[i] = (i-0.5)*dz  

    return z

