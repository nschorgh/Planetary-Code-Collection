real(8) function psv(T)
  ! saturation vapor pressure of H2O ice [Pascal]
  ! input is temperature [Kelvin]
  implicit none
  real(8), intent(IN) :: T

  !-----parametrization 1
  ! real(8), parameter :: DHmelt=6008., DHvap=45050.
  ! real(8), parameter :: DHsub=DHmelt+DHvap ! sublimation enthalpy [J/mol]
  ! real(8), parameter :: R=8.314, pt=6.11e2, Tt=273.16
  ! real(8) C
  ! C = (DHsub/R)*(1./T - 1./Tt)
  ! psv = pt*exp(-C)

  !-----parametrization 2
  ! eq. (2) in Murphy & Koop, Q. J. R. Meteor. Soc. 131, 1539 (2005)
  ! differs from parametrization 1 by only 0.1%
  real(8), parameter :: A=-6143.7, B=28.9074
  psv = exp(A/T+B)  ! Clapeyron

  !-----parametrization 3
  ! eq. (7) in Murphy & Koop, Q. J. R. Meteor. Soc. 131, 1539 (2005)
  ! psv = exp(9.550426 - 5723.265/T + 3.53068*log(T) - 0.00728332*T)
      
end function psv


      
real(8) function frostpoint(p)
  ! inverse of psv
  ! input is partial pressure [Pascal]
  ! output is temperature [Kelvin]
  implicit none
  real(8), intent(IN) :: p
      
  !-----inverse of parametrization 1
  ! real(8), parameter :: DHmelt=6008.,DHvap=45050.
  ! real(8), parameter :: DHsub=DHmelt+DHvap
  ! real(8), parameter :: R=8.314, pt=6.11e2, Tt=273.16
  ! frostpoint = 1./(1./Tt-R/DHsub*log(p/pt))
      
  !-----inverse of parametrization 2
  ! inverse of eq. (2) in Murphy & Koop (2005)
  real(8), parameter :: A=-6143.7, B=28.9074
  frostpoint = A / (log(p) - B)

  !-----approximate inverse of parametrization 3
  ! eq. (8) in Murphy & Koop (2005)
  ! frostpoint = (1.814625*log(p) + 6190.134)/(29.120 - log(p))
      
end function frostpoint
