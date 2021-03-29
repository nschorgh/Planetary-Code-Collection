elemental function evap_vacuum_species(T)
  ! sublimation rate of CO2 ice into vacuum
  implicit none
  real(8) evap_vacuum_species  ! [kg/m^2/s]
  real(8), intent(IN) :: T  ! [Kelvin]
  real(8), parameter :: pi=3.1415926535897932
  real(8), parameter :: Ru=8314.46  ! [J/kmol/K]
  real(8) psv

  !! CO2
  real(8), parameter :: mu = 44.01 
  
  !real(8), parameter :: A=6.81228, B=1301.679, C=-3.494  ! NIST Webbook
  !psv = 10**(A-(B/(T+C)))*1.e5   ! Antoine equation

  ! Fray & Schmitt (2009)
  !real(8), parameter :: A0=1.861e1, A1=-4.154e3, A2=1.041e5  ! 194.7<T<216.58
  !psv = exp(A0+A1/T+A2/T2**2)*1e5
  ! 40<T<194.7
  !real(8), parameter :: A0=1.476e1, A1=-2.571e3, A2=-7.781e4
  !real(8), parameter :: A3=4.325e6, A4=-1.207e8, A5=1.350e9
  !psv = exp(A0+A1/T+A2/T**2+A3/T**3+A4/T**4+A5/T**5)*1e5
  
  ! Bryson et al. (1974) after unit conversions
  real(8), parameter :: A=28.693, B=3270.9 ! 91-190K
  psv = exp(A-B/T)
  !! end of CO2

  evap_vacuum_species = psv*sqrt(mu/(2*pi*Ru*T))
  
end function evap_vacuum_species
