# read and plot output of testcrankT.py

from matplotlib.pyplot import *
import numpy

wpth='./'

a = numpy.loadtxt(wpth + 'Tprofile')
T = a[:,1:]
z = numpy.loadtxt(wpth + 'z')


# plot profile
plot(T[0,:],z,'k-',label='Numerical')
for i in range(1,12):
    plot(T[i,:],z,'k-')
        
xlabel('Temperature (K)')
ylabel('z (m)')
gca().invert_yaxis()


# compare with analytical solution for sinusoidal surface temperature
#                                oscillation and semi-infinite domain
Ta=30.; Tm=190.; P=670.*88775.244
rhoc = 1200*800; thIn=120
delta = thIn/rhoc * numpy.sqrt(P/numpy.pi)   # skin depth
w = 2*numpy.pi/P
dt = P/12
TT = numpy.zeros(len(z))

for i in range(0,12):
    t = i*dt
    TT = Tm + Ta * numpy.exp(-z/delta) * numpy.sin(z/delta-w*t)
    if i==0:
        plot(TT,z,'r-',linewidth=1,label='Analytical')
    else:
        plot(TT,z,'r-',linewidth=1)
        


legend()

show()
# savefig('test_Tprofile.png',format='png')

