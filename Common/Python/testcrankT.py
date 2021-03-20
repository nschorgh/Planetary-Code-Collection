import numpy as np
import grids
from conductionT import conductionT
from math import sin, pi, sqrt


if __name__ == '__main__':
#***********************************************************************
# test Crank-Nicolson subroutine
#***********************************************************************
      Period = 88775.244*670
      NSTEPS = 50000
      Fgeo = 0.
      Ta=30.; Tm=190.

      STEPSPERSOL = 120
      dt = Period/STEPSPERSOL
      nz = 70
      zmax = 2.5; zfac=1.02
      thIn = 120.

      rhocv = np.zeros(nz+1)
      rhocv[:] = 1200.*800.
      delta = thIn/rhocv[1]*sqrt(Period/pi)
      print('Skin depth= ',delta)
      print('zmax=',zmax)
      print('Time step=',dt)
      print('Thermal inertia=',thIn,' Period=',Period)

      #  Initialize
      fout = open('Tprofile','w') # temperature profile

      T = np.zeros(nz+1)
      ti = np.zeros(nz+1)
      ti[:] = thIn

      z = grids.setgrid(nz,zmax,zfac)
      np.savetxt('z', z[1:], fmt='%g', newline=" ")
      
      time = 0.
      Fsurf = 0.

      Tsurf = Tm + Ta*sin(2*pi*time/Period)

      for n in range(0,NSTEPS+1):
          time = (n+1)*dt     #   time at n+1;
          Tsurfp1 = Tm + Ta*sin(2*pi*time/Period) 
          conductionT(nz,z,dt,T,Tsurf,Tsurfp1,ti,rhocv,Fgeo,Fsurf)
          Tsurf = Tsurfp1

          # write 12 profiles from the last sol
          if n > NSTEPS-STEPSPERSOL:
              if n%10==9: # synched with test_Tprofile.m
                  print(time/Period, Tsurf)
                  np.savetxt(fout,np.column_stack(T[:]),fmt=' %7.2f'*(nz+1) )

      # end of time loop

      fout.close()
