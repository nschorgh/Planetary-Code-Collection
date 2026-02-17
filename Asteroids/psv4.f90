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
  !real(8), parameter :: A=-6143.7, B=28.9074
  !psv = exp(A/T+B)  ! Clapeyron

  !-----parametrization 3
  ! eq. (7) in Murphy & Koop, Q. J. R. Meteor. Soc. 131, 1539 (2005)
  psv = exp(9.550426 - 5723.265/T + 3.53068*log(T) - 0.00728332*T)

  !-----parametrization 4
  ! International Association for the Properties of Water and Steam R14-08(2011)
  ! http://www.iapws.org/relguide/MeltSub2011.pdf
  !real(8), parameter :: a(3)=(/ -0.212144006d2,  0.273203819d2, -0.610598130d1 /)
  !real(8), parameter :: b(3)=(/  0.333333333d-2, 0.120666667d1,  0.170333333d1 /)
  !real(8) theta
  !theta = T/273.16
  !psv = 611.657*exp(1./theta*(a(1)*theta**b(1)+a(2)*theta**b(2)+a(3)*theta**b(3)))
  
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
  !real(8), parameter :: A=-6143.7, B=28.9074
  !frostpoint = A / (log(p) - B)

  !-----approximate inverse of parametrization 3
  ! eq. (8) in Murphy & Koop (2005)
  frostpoint = (1.814625*log(p) + 6190.134)/(29.120 - log(p))
      
end function frostpoint
