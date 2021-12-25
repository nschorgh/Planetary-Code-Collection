function bet_isotherm(pp0)
  ! BET isotherm, amount adsorbed as a function of pressure
  ! parameter c was fitted to adsorption isotherms by Cadenhead & Stetter (1974)
  implicit none
  real(8) bet_isotherm, vvm
  real(8), intent(IN) :: pp0 ! p/p0 
  real(8), parameter :: c=8.231
  
  vvm = c * pp0 / (1-pp0) / (1 + (c-1)*pp0)
  bet_isotherm = vvm  ! [# monolayers]
end function bet_isotherm



function p_BET(v)
  ! relative vapor pressure of adsorbed H2O using BET isotherm
  ! parameter c was fitted to adsorption isotherms by Cadenhead & Stetter (1974)
  ! this is the inverse of function bet_isotherm
  implicit none
  real(8) p_BET  ! vapor pressure reletave to saturation pressure of ice (v=oo) 
  real(8), intent(IN) :: v  ! [# monolayers]
  real(8), parameter :: c=8.231

  p_BET = ( c*(v-1) - 2*v + sqrt( c**2*(v-1)**2 + 4*c*v ) ) / ( 2*(c-1)*v )
  if (v<1e-3) p_BET = v/c
  if (v>1e3) p_BET = 1. - 1./v
end function p_BET

function vfv_BET(v)
  ! v / f(v), where f(v) is p_BET(v) 
  ! well-behaved for v=0
  implicit none
  real(8) vfv_BET  ! = v/f(v)
  real(8), intent(IN) :: v  ! [# monolayers]
  real(8), parameter :: c=8.231  ! must match value in p_BET
  real(8), external :: p_BET

  if (v>1e-3) then
     vfv_BET = v / p_BET(v)
  else
     vfv_BET = c
  end if
end function vfv_BET
