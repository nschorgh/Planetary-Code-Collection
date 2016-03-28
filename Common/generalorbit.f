C=============================================================
C Subroutine to return distance, longitude, and declination of 
C the sun in planetocentric coordinates from orbital elements
C
C INPUTS: 
C edays = time since perihelion (earth days)
C a = semimajor axis (AU)
C ecc = eccentricity
C omega = Ls of perihelion, relative to equinox (radians)
C eps = obliquity (radians)
C
C OUTPUTS:
C Ls = areocentric longitude (radians)
C dec = planetocentric solar declination (radians)
C r = heliocentric distance (AU)
C=============================================================

      SUBROUTINE generalorbit(edays,a,ecc,omega,eps,Ls,dec,r)
      implicit none
      real*8 edays,a,ecc,omega,eps  ! input
      real*8 Ls,dec,r  ! output
      real*8 pi,d2r
      parameter (pi=3.1415926535897932,d2r=pi/180.d0)
      integer j
      real*8 M,E,nu,T,Eold

c     T = orbital period (days)
      T = sqrt(4*pi**2/(6.674e-11*1.989e30)*(a*149.598e9)**3)/86400.

c     M = mean anomaly (radians)
      M = 2.*pi*edays/T  ! M=0 at perihelion

c     E = eccentric anomaly 
c     solve M = E - ecc*sin(E) by newton method
      E = M
      do j=1,10
         Eold = E
         E = E - (E - ecc*sin(E) - M)/(1.-ecc*cos(E))
         if (abs(E-Eold)<1.e-8) exit
      enddo

c     nu = true anomaly
      !nu = acos(cos(E)-ecc/(1.-ecc*cos(E)))
      !nu = sqrt(1-ecc^2)*sin(E)/(1.-ecc*cos(E))
      !nu = atan(sqrt(1-ecc^2)*sin(E)/(1-cos(E)))
      nu = 2.*atan(sqrt((1.+ecc)/(1.-ecc))*tan(E/2.))

      !r = a*(1.-ecc**2)/(1.+ecc*cos(nu))
      r = a*(1-ecc*cos(E))
      Ls = mod(nu + omega,2.*pi)   
      dec = asin(sin(eps)*sin(Ls))

      END

