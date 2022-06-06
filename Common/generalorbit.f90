!=============================================================
! Subroutine to return distance, longitude, and declination of 
! the sun in planetocentric coordinates from orbital elements
! 
! INPUTS: 
! edays = time since perihelion (earth days)
! a = semimajor axis (AU)
! ecc = eccentricity
! omega = Ls of perihelion, relative to equinox (radians)
! eps = obliquity (radians)
! 
! OUTPUTS:
! Ls = areocentric longitude (radians)
! dec = planetocentric solar declination (radians)
! r = heliocentric distance (AU)
!=============================================================

SUBROUTINE generalorbit(edays,a,ecc,omega,eps,Ls,dec,r)
  implicit none
  real(8), intent(IN) :: edays, a, ecc, omega, eps
  real(8), intent(OUT) :: Ls, dec, r 
  real(8), parameter :: pi=3.1415926535897932, d2r=pi/180.d0
  integer j
  real(8) M, E, nu, T, Eold

  ! T ... orbital period (days)
  T = sqrt(4*pi**2/(6.674e-11*1.989e30)*(a*149.598e9)**3)/86400.

  ! M ... mean anomaly (radians)
  M = 2.*pi*edays/T  ! M=0 at perihelion

  ! E ... eccentric anomaly 
  ! solve M = E - ecc*sin(E) by Newton method
  E = M
  do j=1,10
     Eold = E
     E = E - (E - ecc*sin(E) - M)/(1.-ecc*cos(E))
     if (abs(E-Eold)<1.e-8) exit
  enddo

  ! nu ... true anomaly
  !nu = acos(cos(E)-ecc/(1.-ecc*cos(E)))
  !nu = sqrt(1-ecc^2)*sin(E)/(1.-ecc*cos(E))
  !nu = atan(sqrt(1-ecc^2)*sin(E)/(1-cos(E)))
  nu = 2.*atan(sqrt((1.+ecc)/(1.-ecc))*tan(E/2.))

  !r = a*(1.-ecc**2)/(1.+ecc*cos(nu))
  r = a*(1-ecc*cos(E))
  Ls = mod(nu + omega,2.*pi)   
  dec = asin(sin(eps)*sin(Ls))

END SUBROUTINE generalorbit

