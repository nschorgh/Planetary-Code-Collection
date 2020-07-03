! Functions that parametrize sublimation rates into vacuum


function sublrate(T)
  ! sublimation rate of H2O [#molecules/m^2/s]
  implicit none
  real(8), intent(IN) :: T
  real(8) sublrate
  real(8), parameter :: pi=3.1415926535897932, kB = 1.38065e-23
  real(8), parameter :: mu = 18.015*1.66054e-27
  real(8), external :: psv

  ! crystalline ice
  sublrate = psv(T)/sqrt(2*pi*kB*T*mu);

  ! amorphous ice according to Sack & Baragiola (1993)
  !sublrate = 1.82e21*1e4*T**3.5*exp(-0.45*1.6022e-19/(kB*T))
end function sublrate



pure real(8) function restime_species(T)
  ! sublimation rate [#molecules/m^2/s]
  implicit none
  real*8, intent(IN) :: T
  real*8, parameter :: pi=3.1415926535897932, kB = 1.38065e-23, u=1.66054e-27
  real*8 psv, sublrate_species
  
  ! Argon
  real*8, parameter :: mu = 39.962*u
  real*8, parameter :: A=-7814.5, B=+7.5741   ! Int. Crit. Tbls. Vol 3
  real*8, parameter :: sigma0 = 8.42e18   ! (1623/(40*1.66e-27))**(2./3.)  
  psv = 10**(0.05223*A/T + B)*133.32  ! mmHg -> Pa
  !real*8, parameter :: A=4.46903, B=481.012, C=22.156  ! 114-150K, NIST Webbook
  !psv = 10**(A-(B/(T+C)))*1.e5   ! Antoine equation
  
  ! CO2
  !real*8, parameter :: mu = 44.01*u
  !real*8, parameter :: A=6.81228, B=1301.679, C=-3.494  ! NIST Webbook
  !real*8, parameter :: sigma0 = 7.5e18   ! (1500/mu)**(2./3.)
  !psv = 10**(A-(B/(T+C)))*1.e5   ! Antoine equation
  
  ! CH4
  !real*8, parameter :: mu = 16.04*u
  !real*8, parameter :: A=3.9895, B=443.028, C=-0.49 ! 91-190K, NIST Webbook
  !real*8, parameter :: sigma0 = 7.1e19  ! (500/mu)**(2./3.) ! obscure source
  !psv = 10**(A-(B/(T+C)))*1.e5   ! Antoine equation
  
  sublrate_species = psv*sqrt(1./(2*pi*kB*T*mu));
  restime_species = sigma0/sublrate_species
end function restime_species
