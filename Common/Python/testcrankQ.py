import numpy as np
import grids
from conductionQ import conductionQ
from flux_noatm import flux_noatm


if __name__ == '__main__':
#***********************************************************************
# test Crank-Nicolson subroutine
#***********************************************************************
    sigSB = 5.6704e-8
    Period = 88775.244*670 # [seconds]
    NSTEPS = 50000
    emiss = 1.
    Fgeo = 0.2  # [W/m^2]

    STEPSPERSOL = 120
    dt = Period / STEPSPERSOL
    thIn = 120.  # thermal inertia
    albedo = 0.2
    latitude = 5.  # [degree]

    nz = 60; zmax = 2.5; zfac = 1.05

    # rhoc = thIn * np.sqrt(Period/pi)  # skin depth = 1
    rhocv = np.zeros(nz+1)
    rhocv[:] = 1200.*800.  # (density) * (heat capacity)
    delta = thIn/rhocv[1] * np.sqrt(Period/np.pi)
    print('Skin depth= ',delta)
    ti = np.zeros(nz+1)
    ti[:]= thIn

    T = np.zeros(nz+1)
    Tmean = np.zeros(nz+1)
    
    Rau = 1.52
    Decl = 0.
  
    print('Time step=',dt)
    print('zmax=',zmax)
    print('Thermal inertia=',thIn,' Period=',Period)
    print('Heat Flux at bottom boundary=',-Fgeo)
    
    # Initialize
    fout1 = open('Tsurface',"w")  # surface temperature
    fout2 = open('Tprofile',"w")  # temperature profile
    
    T[:] = 210.
    z = grids.setgrid(nz,zmax,zfac)
    np.savetxt('z', z[1:], fmt='%g', newline=" ")
  
    latitude = np.deg2rad(latitude)

    time = 0.
    #fout1.write('%12.6f %9.3f %9.3f\n' % (0.,T[0],T[nz]) )
    #np.savetxt(fout2, np.column_stack(T[:]), fmt=" %7.2f"*(nz+1))
    
    Fmean = 0.
    HA = 0.
    Qn = (1-albedo) * flux_noatm(Rau,Decl,latitude,HA,0.,0.)
    Fsurf = 0.

    for n in range (0,NSTEPS+1):
        
        time = (n+1)*dt   #   time at n+1; 
        HA = 2 * np.pi * (time/Period % 1.) #  hour angle
        Qnp1 = (1-albedo) * flux_noatm(Rau,Decl,latitude,HA,0.,0.)
        conductionQ(nz,z,dt,Qn,Qnp1,T,ti,rhocv,emiss,Fgeo,Fsurf)
        Qn = Qnp1

        if n%3 == 0:
            fout1.write('%12.6f %9.3f %9.3f\n' % (time/Period,T[0],T[nz]) )

        if (n > NSTEPS-STEPSPERSOL):
            if n%10 == 0:
                np.savetxt(fout2,np.column_stack(T[:]),fmt="%7.2f"*(nz+1))

        Fmean += Fmean
        Tmean = Tmean[:] + T[:]
     
  
    Fmean = Fmean / STEPSPERSOL
    Tmean[:] = Tmean[:] / STEPSPERSOL

    fout1.close()
    fout2.close()

