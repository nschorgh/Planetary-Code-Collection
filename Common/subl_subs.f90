! Functions that parametrizes sublimation rates into vacuum

function sublrate(T)
  ! sublimation rate of H2O in #molecules/m^2/s
  implicit none
  real(8),intent(IN) :: T
  real(8) sublrate
  real(8), parameter :: pi=3.141592653589793, kB = 1.38065e-23
  real(8), parameter :: mu = 18.015*1.66054e-27
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
  real(8), parameter :: kB=1.38065e-23
  real(8), parameter :: theta = 9.8947e+18
  desorptionrate = nu*theta*exp(-0.5*1.6022e-19/(kB*T)) ! zeroth order
end function desorptionrate


pure function padsr(x)
  ! relative vapor pressure of adsorbed H2O
  ! based on adsorption isotherms by Cadenhead & Stetter (1974)
  ! see Schorghofer & Aharonson (2014)
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


pure real(8) function restime_species(T)
  ! sublimation rate in #molecules/m^2/s
  implicit none
  real(8),intent(IN) :: T
  real(8), parameter :: pi=3.1415926535897932, kB = 1.38065e-23
  real(8) psv, sublrate_species
  
  ! Argon
  !real(8), parameter :: mu = 39.962*1.66054e-27
  !real(8), parameter :: A=-7814.5, B=+7.5741   ! Ict Vol 3
  !real(8), parameter :: sigma0 = 8.42e18   ! (1623/(40*1.66e-27))**(2./3.)  
  !psv = 10**(0.05223*A/T + B)*133.32  ! mmHg -> Pa
  
  ! CO2
  real*8, parameter :: mu = 44.01*1.66054e-27
  real*8, parameter :: A=6.81228, B=1301.679, C=-3.494
  real*8, parameter :: sigma0 = 7.5e18   ! (1500/mu)**(2./3.)
  psv = 10**(A-(B/(T+C)))*1.e5   ! Antoine equation
  
  ! SO2
  !real*8, parameter :: mu = 64.06*1.66054e-27
  !real*8, parameter :: A=3.48586, B=668.225, C=-72.252  ! 178-263K
  !real*8, parameter :: sigma0 = 8.5e18   ! (2630/mu)**(2./3.)
  !psv = 10**(A-(B/(T+C)))*1.e5   ! Antoine equation

  sublrate_species = psv*sqrt(1./(2*pi*kB*T*mu));
  restime_species = sigma0/sublrate_species
end function restime_species
