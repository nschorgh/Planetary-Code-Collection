function sublrate(T)
  ! sublimation rate of H2O in #molecules/m^2/s
  implicit none
  real(8),intent(IN) :: T
  real(8) sublrate
  real(8), parameter :: pi=3.141592653589793, kB = 1.38065e-23
  real(8), parameter :: mu = 18.015*1.6605e-27
  real(8), external :: psv

  ! crystalline ice
  sublrate = psv(T)/sqrt(2*pi*kB*T*mu);

  ! amorphous ice according to Sack and Baragiola, 1993
  !sublrate = 1.82e21*1e4*T**3.5*exp(-0.45*1.6022e-19/(kB*T))
end function sublrate


pure function desorptionrate(T)
  ! Hibbitts et al. (2011), Icarus
  implicit none
  real(8), intent(IN) :: T
  real(8) desorptionrate
  real(8), parameter :: nu = 1e13
  real(8), parameter :: k=1.38065e-23, m=18*1.66e-27
  real(8), parameter :: theta = 9.8947e+18
  desorptionrate = nu*theta*exp(-0.5*1.6022e-19/(k*T)) ! zeroth order
end function desorptionrate


pure function padsr(x)
  ! based on Cadenhead's isotherms
  implicit none
  real(8), intent(IN) :: x
  real(8) padsr, b, c
  real(8), parameter :: a=0.402, x0=2.45
  
  if (x<x0) then
     b = (exp(-a*x0)*(2+a*x0) - 2)/x0**3
     c = (3 - exp(-a*x0)*(3+a*x0))/x0**2
     padsr = b*x**3 + c*x**2
  else
     padsr = 1 - exp(-a*x)
  endif
  ! 0<= padsr <=1
end function padsr


real(8) function sublrate_argon(T)
  ! sublimation rate in #molecules/m^2/s
  ! The vapor pressure of He is, for all practical purposes, infinite.
  implicit none
  real(8),intent(IN) :: T
  real(8), parameter :: pi=3.1415926535897932
  real(8), parameter :: mu = 39.962*1.6605e-27, kB = 1.38065e-23
  ! Argon
  real(8), parameter :: A=-7814.5, B=+7.5741   ! Argon, Ict Vol 3
  real(8), parameter :: sigma0 = 8.42e18   ! (1623/(40*1.66e-27))**(2./3.)  
  real(8) psv, residence_time

  psv = 10**(0.05223*A/T + B)*133.32  ! mmHg -> Pa
  sublrate_argon = psv*sqrt(1./(2*pi*kB*T*mu));
  residence_time = sigma0/sublrate_argon
end function sublrate_argon
